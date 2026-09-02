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

// wrap at FASTA_WRAP columns; col carries the column across calls
static void emit_wrapped(FILE *out, const char *seq, size_t len, int *col) {
    size_t printed = 0;
    while (printed < len) {
        size_t space = (size_t)(FASTA_WRAP - *col);
        size_t take = (len - printed < space) ? (len - printed) : space;
        fwrite(seq + printed, 1, take, out);
        printed += take;
        *col += (int)take;
        if (*col >= FASTA_WRAP) {
            fputc('\n', out);
            *col = 0;
        }
    }
}

// print the extract usage lines
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal extract fa   <r/GFA> <out.fa|.fasta>");
    ak_log(AK_LOG_ERROR, NULL, "       akhal extract path <r/GFA> <out.fa|.fasta> [PATH-NAME]");
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

// `extract fa` - one FASTA record per P line, exactly as the graph stores them
static int extract_fa(int argc, char **argv) {
    if (argc != 5) {
        usage();
        return 1;
    }
    const char *in = argv[3], *out_fn = argv[4];
    if (!want_gfa(in) || !want_fasta(out_fn)) return 1;

    FILE *out = fopen(out_fn, "w");
    if (!out) {
        ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
        return 1;
    }

    gfa_t *g = gfa_read(in, GFA_PATHS);
    if (!g) {
        fclose(out);
        return 1;
    }

    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        fprintf(out, ">%s\n", gfa_path_name(g, k));

        const uint32_t *segs;
        int n = gfa_path_segs(g, k, &segs);
        int col = 0;
        for (int i = 0; i < n; i++) {
            if (segs[i] == GFA_NIL) continue;
            gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
            if (s->seq) {
                emit_wrapped(out, s->seq, s->len, &col);
            }
        }
        if (col) {
            fputc('\n', out);
        }
    }

    fclose(out);
    gfa_destroy(g);
    ak_log(AK_LOG_INFO, NULL, "extracted reference to %s", out_fn);
    return 0;
}

// `extract path` - fragments that the links join leave as one FASTA record
static int extract_path(int argc, char **argv) {
    if (argc < 5 || argc > 6) {
        usage();
        return 1;
    }
    const char *in = argv[3], *out_fn = argv[4];
    const char *key = (argc == 6) ? argv[5] : NULL;
    if (!want_gfa(in) || !want_fasta(out_fn)) return 1;

    gfa_t *g = gfa_read(in, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;

    gfa_merge_t *m = gfa_path_merge(g, key);
    if (!m) {
        gfa_destroy(g);
        return 1;
    }

    FILE *out = fopen(out_fn, "w");
    if (!out) {
        gfa_merge_destroy(m);
        gfa_destroy(g);
        ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
        return 1;
    }

    int ret = 0;
    int32_t n_merged = 0;
    for (int32_t k = 0; k < m->n; k++) {
        uint32_t *segs;
        int64_t ns = gfa_merge_segs(g, m, k, &segs);
        if (ns < 0) {
            ak_log(AK_LOG_ERROR, NULL, "cannot assemble path %s: %s", m->name[k], ak_strerror((int)ns));
            ret = 1;
            break;
        }
        if (m->off[k + 1] - m->off[k] > 1) {
            n_merged++;
        }

        fprintf(out, ">%s\n", m->name[k]);
        int col = 0;
        for (int64_t i = 0; i < ns; i++) {
            gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
            if (s->seq) {
                emit_wrapped(out, s->seq, s->len, &col);
            }
        }
        if (col) {
            fputc('\n', out);
        }
        free(segs);
    }

    fclose(out);
    if (!ret) {
        ak_log(AK_LOG_INFO, NULL, "wrote %d path(s) to %s: %d stitched from fragments, %d left as they were", m->n, out_fn, n_merged, m->n - n_merged);
    }

    gfa_merge_destroy(m);
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

    gfa_t *g = gfa_read(in, GFA_LINKS | GFA_PATHS);
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
