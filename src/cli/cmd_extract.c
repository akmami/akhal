#include "akhal/gfa.h"
#include "akhal/call.h"
#include "akhal/fasta.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FASTA_WRAP 80

// wrap at wrap_len columns; col carries the column across calls
static void emit_wrapped(FILE *out, const char *seq, size_t len, int *col, int wrap_len) {
    size_t printed = 0;
    while (printed < len) {
        size_t space = (size_t)(wrap_len - *col);
        size_t take = (len - printed < space) ? (len - printed) : space;
        fwrite(seq + printed, 1, take, out);
        printed += take;
        *col += (int)take;
        if (*col >= wrap_len) {
            fputc('\n', out);
            *col = 0;
        }
    }
}

// print the extract usage lines
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal extract fa   <r/GFA> <out.fa|.fasta> [WRAP-LENGTH]");
    ak_log(AK_LOG_ERROR, NULL, "       akhal extract path <r/GFA> <out.fa|.fasta> <PATH-NAME> [PATH-NAME ...] [WRAP-LENGTH]");
    ak_log(AK_LOG_ERROR, NULL, "       akhal extract vcf  <r/GFA> <out.vcf> [--ref <NAME>] [--fasta <FILE>]");
}

// 1 when the extension fits, else 0 and the reason is logged
static int want_fasta(const char *fn) {
    if (ak_ends_with(fn, ".fa") || ak_ends_with(fn, ".fasta")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "output must be a .fa/.fasta file: %s", fn);
    return 0;
}

// 1 when the extension fits, else 0 and the reason is logged
static int want_gfa(const char *fn) {
    if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
    return 0;
}

// 1 when arg is a positive wrap length, else 0 and the reason is logged
static int want_wrap_len(const char *arg, int *out) {
    int v;
    if (!ak_str2int(arg, &v) || v <= 0) {
        ak_log(AK_LOG_ERROR, NULL, "wrap length must be a positive integer: %s", arg);
        return 0;
    }
    *out = v;
    return 1;
}

// one FASTA record: the header, then the bases of every segment in order
static void emit_record(FILE *out, const gfa_t *g, const char *name, const uint32_t *segs, int64_t n, int wrap_len) {
    fprintf(out, ">%s\n", name);

    int col = 0;
    for (int64_t i = 0; i < n; i++) {
        if (segs[i] == GFA_NIL) continue;
        const gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
        if (s->seq) {
            emit_wrapped(out, s->seq, s->len, &col, wrap_len);
        }
    }
    if (col) {
        fputc('\n', out);
    }
}

// A bare contig name selects a PanSN path too, so "chr22" finds
// "GRCh38#0#chr22" without the caller having to spell the whole thing
static int path_selected(const char *name, const char *key) {
    if (!key) return 1;
    size_t nl = strlen(name), kl = strlen(key);
    if (nl == kl && !strcmp(name, key)) return 1;
    for (size_t i = nl; i > 0; i--) {
        if (name[i - 1] != '#') continue;
        return nl - i == kl && !strcmp(name + i, key);
    }
    return 0;
}

// one record per P line. `key` picks the paths of that name, NULL takes them
// all. Returns the number of records written, or a negative AK_E* code
static int64_t emit_paths(const gfa_t *g, FILE *out, const char *key, int wrap_len) {
    int64_t n = 0;
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        const char *name = gfa_path_name(g, k);
        if (!path_selected(name, key)) continue;
        const uint32_t *segs;
        int ns = gfa_path_segs(g, k, &segs);
        emit_record(out, g, name, segs, ns, wrap_len);
        n++;
    }
    if (n == 0 && key) {
        ak_log(AK_LOG_ERROR, NULL, "no path named '%s' in the graph", key);
        return AK_EINVAL;
    }
    return n;
}

// `extract fa` - every P line as one FASTA record
static int extract_fa(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL;
    int wrap_len = FASTA_WRAP, seen_wrap = 0;

    for (int i = 3; i < argc; i++) {
        if (argv[i][0] == '-') {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[i]);
            usage();
            return 1;
        } else if (!in) {
            in = argv[i];
        } else if (!out_fn) {
            out_fn = argv[i];
        } else if (!seen_wrap) {
            if (!want_wrap_len(argv[i], &wrap_len)) return 1;
            seen_wrap = 1;
        } else {
            usage();
            return 1;
        }
    }
    if (!in || !out_fn) {
        usage();
        return 1;
    }
    if (!want_gfa(in) || !want_fasta(out_fn)) return 1;

    // writing the paths as they lie needs their steps and the bases, nothing else
    gfa_t *g = gfa_read(in, GFA_PATHS | GFA_SEQ);
    if (!g) return 1;

    if (gfa_n_path(g) == 0) {
        ak_log(AK_LOG_WARN, "extract", "graph has no P lines: there is nothing to write");
    }

    FILE *out = fopen(out_fn, "w");
    if (!out) {
        ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
        gfa_destroy(g);
        return 1;
    }

    int ret = 0;
    int64_t n = emit_paths(g, out, NULL, wrap_len);
    if (n < 0) {
        ret = 1;
    } else {
        ak_log(AK_LOG_INFO, NULL, "wrote %lld record(s) to %s, one per P line", (long long)n, out_fn);
    }

    fclose(out);
    gfa_destroy(g);
    return ret;
}

// `extract path` - the named paths, one FASTA record each. At least one name
// is required; `fa` is how to take them all
static int extract_path(int argc, char **argv) {
    if (argc < 6) {
        usage();
        return 1;
    }
    const char *in = argv[3], *out_fn = argv[4];

    // the names run to the end, except that a trailing number is the wrap
    // length - which needs a name in front of it to be told from one
    int first = 5, last = argc - 1;
    int wrap_len = FASTA_WRAP, v;
    if (last > first && ak_str2int(argv[last], &v)) {
        if (!want_wrap_len(argv[last], &wrap_len)) return 1;
        last--;
    }
    if (!want_gfa(in) || !want_fasta(out_fn)) return 1;

    gfa_t *g = gfa_read(in, GFA_PATHS | GFA_SEQ);
    if (!g) return 1;

    FILE *out = fopen(out_fn, "w");
    if (!out) {
        gfa_destroy(g);
        ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
        return 1;
    }

    // a name that matches nothing stops the whole thing rather than quietly
    // leaving a half-written file behind
    int64_t n_total = 0;
    int ret = 0;
    for (int i = first; i <= last; i++) {
        int64_t n = emit_paths(g, out, argv[i], wrap_len);
        if (n < 0) {
            ret = 1;
            break;
        }
        n_total += n;
    }

    if (!ret) {
        ak_log(AK_LOG_INFO, NULL, "wrote %lld record(s) to %s", (long long)n_total, out_fn);
    }

    fclose(out);
    gfa_destroy(g);
    return ret;
}

// `extract vcf` - every detour off the reference backbone becomes a VCF row
static int extract_vcf(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL, *ref_name = NULL, *fa_fn = NULL;

    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--ref") && i + 1 < argc) {
            ref_name = argv[++i];
        } else if (!strcmp(argv[i], "--fasta") && i + 1 < argc) {
            fa_fn = argv[++i];
        } else if (argv[i][0] == '-') {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[i]);
            usage();
            return 1;
        } else if (!in) {
            in = argv[i];
        } else if (!out_fn) {
            out_fn = argv[i];
        } else {
            usage();
            return 1;
        }
    }
    if (!in || !out_fn) {
        usage();
        return 1;
    }
    if (!want_gfa(in)) return 1;
    if (!ak_ends_with(out_fn, ".vcf")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .vcf file: %s", out_fn);
        return 1;
    }

    gfa_t *g = gfa_read(in, GFA_ALL);
    if (!g) return 1;

    // label the backbone; everything the graph adds around it is a variation
    call_ref_t *ref = NULL;
    if (fa_fn) {
        fasta_t *fa = fasta_read(fa_fn);
        if (!fa) {
            gfa_destroy(g);
            return 1;
        }
        ref = call_ref_fasta(g, fa, ref_name);
        fasta_destroy(fa);
    } else {
        ref = call_ref_path(g, ref_name);
    }
    if (!ref) {
        gfa_destroy(g);
        return 1;
    }

    int ret = 1;
    call_t *c = call_variants(g, ref);
    if (c) {
        if (call_write_vcf(c, ref, out_fn) == AK_OK) {
            ak_log(AK_LOG_INFO, NULL, "wrote %lld variant(s) to %s", (long long)call_n(c), out_fn);
            ret = 0;
        }
        call_destroy(c);
    }

    call_ref_destroy(ref);
    gfa_destroy(g);
    return ret;
}

// `extract` entry point; see cli.h
int cmd_extract(int argc, char **argv) {
    if (argc < 3) {
        usage();
        return 1;
    }

    if (!strcmp(argv[2], "fa"))   return extract_fa(argc, argv);
    if (!strcmp(argv[2], "path")) return extract_path(argc, argv);
    if (!strcmp(argv[2], "vcf"))  return extract_vcf(argc, argv);

    ak_log(AK_LOG_ERROR, NULL, "unknown extract target: %s", argv[2]);
    usage();
    return 1;
}
