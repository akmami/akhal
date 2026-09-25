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
    ak_log(AK_LOG_ERROR, NULL, "       akhal extract vcf  <r/GFA> <out.vcf> [--ref <NAME>[,<NAME>...] | --ref all] [--fasta <FILE>]");
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

// one record per P line. `key` picks the paths of that name, NULL takes them all. 
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

// `extract path` - the named paths, one FASTA record each. At least one name is required; `fa` is how to take them all
static int extract_path(int argc, char **argv) {
    if (argc < 6) {
        usage();
        return 1;
    }
    const char *in = argv[3], *out_fn = argv[4];

    // the names run to the end, except that a trailing number is the wrap length
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

    // a name that matches nothing stops the whole thing rather than quietly leaving a half-written file behind
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

// Split a --ref value on commas, in place: "chr1,chr2, chr3" is three names
static int split_refs(char *list, char ***out, int *n_out) {
    int n = 0, m = 8;
    char **names = (char **)malloc((size_t)m * sizeof(*names));
    if (!names) return AK_ENOMEM;

    for (char *p = list;;) {
        char *comma = strchr(p, ',');
        if (comma) *comma = '\0';
        while (*p == ' ' || *p == '\t') p++;
        char *e = p + strlen(p);
        while (e > p && (e[-1] == ' ' || e[-1] == '\t')) *--e = '\0';

        if (!*p) {
            ak_log(AK_LOG_ERROR, NULL, "--ref has an empty name in its list");
            free(names);
            return AK_EINVAL;
        }
        for (int i = 0; i < n; i++) {
            if (!strcmp(names[i], p)) {
                ak_log(AK_LOG_ERROR, NULL, "--ref names '%s' more than once", p);
                free(names);
                return AK_EINVAL;
            }
        }
        if (n == m) {
            m <<= 1;
            char **r = (char **)realloc(names, (size_t)m * sizeof(*names));
            if (!r) {
                free(names);
                return AK_ENOMEM;
            }
            names = r;
        }
        names[n++] = p;

        if (!comma) break;
        p = comma + 1;
    }
    *out = names;
    *n_out = n;
    return AK_OK;
}

// this is needed for qsort
static inline int name_cmp(const void *a, const void *b) {
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

// Every name --ref all stands for: each P line, or with --fasta each record, in file order
static int all_names(const gfa_t *g, const fasta_t *fa, const char *src, char ***out, int *n_out) {
    int n = fa ? (int)fasta_n(fa) : (int)gfa_n_path(g);
    if (n == 0) {
        ak_log(AK_LOG_ERROR, NULL, fa ? "--ref all: %s holds no sequence" : "--ref all: %s has no P lines", src);
        return AK_EINVAL;
    }

    char **names  = (char **)malloc((size_t)n * sizeof(char *));
    char **sorted = (char **)malloc((size_t)n * sizeof(char *));
    if (!names || !sorted) {
        free(names);
        free(sorted);
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return AK_ENOMEM;
    }
    for (int i = 0; i < n; i++) {
        names[i] = fa ? fa->rec[i].name : (char *)gfa_path_name(g, i);
    }

    // check if there is any name duplicates in fasta or P
    memcpy(sorted, names, (size_t)n * sizeof(char *));
    qsort(sorted, (size_t)n, sizeof(char *), name_cmp);
    int rc = AK_OK;
    for (int i = 1; i < n && rc == AK_OK; i++) {
        if (strcmp(sorted[i - 1], sorted[i]) != 0) continue;
        int run = 2;
        while (i + run - 1 < n && !strcmp(sorted[i - 1], sorted[i + run - 1])) run++;
        ak_log(AK_LOG_ERROR, NULL, "--ref all: '%s' names %d %s in %s; each needs a name of its own - a path written as fragments has to be joined first", sorted[i], run, fa ? "records" : "P lines", src);
        rc = AK_EINVAL;
    }
    free(sorted);
    if (rc != AK_OK) {
        free(names);
        return rc;
    }

    ak_log(AK_LOG_INFO, NULL, "--ref all: %d backbone(s)", n);
    *out = names;
    *n_out = n;
    return AK_OK;
}

// With no --ref, the one backbone is the library's own default: the graph's first P line, or the FASTA's first record
static int default_name(const gfa_t *g, const fasta_t *fa, char ***out, int *n_out) {
    if (fa ? fasta_n(fa) == 0 : gfa_n_path(g) == 0) {
        ak_log(AK_LOG_ERROR, NULL, fa ? "the FASTA holds no sequence" : "graph has no P lines to use as a backbone");
        return AK_EINVAL;
    }
    char **names = (char **)malloc(sizeof(char *));
    if (!names) {
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return AK_ENOMEM;
    }
    names[0] = fa ? fa->rec[0].name : (char *)gfa_path_name(g, 0);
    *out = names;
    *n_out = 1;
    return AK_OK;
}

// the first P line carrying exactly this name, or -1 - the same one call_ref_path() takes when the name is there
static int32_t path_named(const gfa_t *g, const char *name) {
    for (int32_t k = 0; k < gfa_n_path(g); k++)
        if (!strcmp(gfa_path_name(g, k), name)) return k;
    return -1;
}

// Check every name and find each contig's length, before any backbone is built
static int ref_lengths(const gfa_t *g, const fasta_t *fa, const char *src, char *const *names, int n, int64_t **out) {
    int64_t *lens = (int64_t *)malloc((size_t)n * sizeof(int64_t));
    if (!lens) {
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return AK_ENOMEM;
    }
    for (int i = 0; i < n; i++) {
        if (fa) {
            const fasta_rec_t *fr = fasta_get(fa, names[i]);
            if (!fr) {
                ak_log(AK_LOG_ERROR, NULL, "no sequence named '%s' in %s", names[i], src);
                free(lens);
                return AK_EINVAL;
            }
            lens[i] = fr->len;
        } else {
            int32_t k = path_named(g, names[i]);
            if (k < 0) {
                ak_log(AK_LOG_ERROR, NULL, "no P line named '%s' in %s", names[i], src);
                free(lens);
                return AK_EINVAL;
            }
            lens[i] = (int64_t)gfa_path_len(g, k);
        }
    }
    *out = lens;
    return AK_OK;
}

// one contig: its backbone built, its variants called and appended, and both released before the next
static int write_contig(const gfa_t *g, const fasta_t *fa, const char *name, FILE *fp, int64_t *total) {
    call_ref_t *ref = fa ? call_ref_fasta(g, fa, name) : call_ref_path(g, name);
    if (!ref) return AK_EINVAL;

    call_t *c = call_variants(g, ref);
    if (!c) {
        call_ref_destroy(ref);
        return AK_ENOMEM;
    }
    int rc = call_vcf_records(fp, c, ref);
    *total += call_n(c);

    call_destroy(c);
    call_ref_destroy(ref);
    return rc;
}

// the whole file: the header naming every contig, then each contig's records in the same order. 
static int write_vcf(const gfa_t *g, const fasta_t *fa, char *const *names, const int64_t *lens, int n, const char *out_fn) {
    FILE *fp = fopen(out_fn, "w");
    if (!fp) {
        ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
        return AK_EOPEN;
    }

    int64_t total = 0;
    int rc = call_vcf_header(fp, (const char *const *)names, lens, n);
    for (int i = 0; i < n && rc == AK_OK; i++) {
        rc = write_contig(g, fa, names[i], fp, &total);
    }
    if (fclose(fp) != 0 && rc == AK_OK) {
        rc = AK_EIO;
    }

    if (rc != AK_OK) {
        if (rc == AK_EIO) ak_log(AK_LOG_ERROR, NULL, "write failed on %s", out_fn);
        remove(out_fn);
        return rc;
    }
    ak_log(AK_LOG_INFO, NULL, "wrote %lld variant(s) over %d contig(s) to %s", (long long)total, n, out_fn);
    return AK_OK;
}

// `extract vcf` - every detour off the reference backbone becomes a VCF row.
// --ref takes one name, several comma-separated, or "all"; each becomes a contig, in the order given
static int extract_vcf(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL, *fa_fn = NULL;
    char *ref_arg = NULL;

    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--ref") && i + 1 < argc) {
            ref_arg = argv[++i];
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

    // A list is checked for empty and repeated names now, before the graph is read; "all" stands for names only the input can supply
    int want_all = ref_arg && !strcmp(ref_arg, "all");
    char **names = NULL;
    int n_names = 0;
    if (ref_arg && !(want_all) && split_refs(ref_arg, &names, &n_names) != AK_OK) return 1;

    gfa_t *g = gfa_read(in, GFA_ALL);
    fasta_t *fa = NULL;
    int64_t *lens = NULL;

    int rc = g ? AK_OK : AK_EOPEN;
    if (rc == AK_OK && fa_fn) {
        fa = fasta_read(fa_fn);
        if (!fa) rc = AK_EOPEN;
    }
    if (rc == AK_OK && !names) {
        rc = want_all ? all_names(g, fa, fa ? fa_fn : in, &names, &n_names) : default_name(g, fa, &names, &n_names);
    }
    if (rc == AK_OK) rc = ref_lengths(g, fa, fa ? fa_fn : in, names, n_names, &lens);
    if (rc == AK_OK) rc = write_vcf(g, fa, names, lens, n_names, out_fn);

    free(lens);
    free(names);
    fasta_destroy(fa);
    gfa_destroy(g);
    return rc == AK_OK ? 0 : 1;
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
