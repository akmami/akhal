#include "akhal/gfa.h"
#include "akhal/diff.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <string.h>

// print the compare usage lines
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal compare gfa <A.gfa> <B.gfa> [--verbose]");
    ak_log(AK_LOG_ERROR, NULL, "       akhal compare gaf <A.gaf> <B.gaf> [--verbose]");
}

// 1 when the extension fits, else 0 and the reason is logged
static int want_gfa(const char *fn) {
    if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
    return 0;
}

// 1 when the extension fits, else 0 and the reason is logged
static int want_gaf(const char *fn) {
    if (ak_ends_with(fn, ".gaf")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "expected a .gaf file: %s", fn);
    return 0;
}

// the two inputs and an optional --verbose, which is all either target takes.
// Returns 1 on success, else 0 and the usage is printed
static int parse_args(int argc, char **argv, const char **a, const char **b, int *verbose) {
    *a = *b = NULL;
    *verbose = 0;

    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--verbose")) {
            *verbose = 1;
        } else if (argv[i][0] == '-') {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[i]);
            usage();
            return 0;
        } else if (!*a) {
            *a = argv[i];
        } else if (!*b) {
            *b = argv[i];
        } else {
            usage();
            return 0;
        }
    }
    if (!*a || !*b) {
        usage();
        return 0;
    }
    return 1;
}

// graphs

// the counts both graphs are measured by
static void print_counts(const gfa_t *a, const gfa_t *b, const diff_t *d) {
    printf("Segments A: %ld\n", (long)gfa_n_seg(a));
    printf("Segments B: %ld\n", (long)gfa_n_seg(b));
    printf("Segments shared: %ld\n", (long)d->n_seg_shared);
    printf("Segments only in A: %ld\n", (long)d->a.n_seg);
    printf("Segments only in B: %ld\n", (long)d->b.n_seg);

    printf("Links A: %ld\n", (long)gfa_n_link(a));
    printf("Links B: %ld\n", (long)gfa_n_link(b));
    printf("Links shared: %ld\n", (long)d->n_link_shared);
    printf("Links only in A: %ld\n", (long)d->a.n_link);
    printf("Links only in B: %ld\n", (long)d->b.n_link);

    printf("Paths A: %ld\n", (long)(d->n_path_same + d->n_path_differ + d->n_path_a_only));
    printf("Paths B: %ld\n", (long)(d->n_path_same + d->n_path_differ + d->n_path_b_only));
    printf("Paths identical: %ld\n", (long)d->n_path_same);
    printf("Paths differing: %ld\n", (long)d->n_path_differ);
    printf("Paths only in A: %ld\n", (long)d->n_path_a_only);
    printf("Paths only in B: %ld\n", (long)d->n_path_b_only);
}

// one line per path name, since that is the whole verdict for a path
static void print_paths(const diff_t *d) {
    for (int32_t i = 0; i < d->n_path; i++) {
        const diff_path_t *p = &d->path[i];
        switch (p->state) {
            case DIFF_SAME:
                printf("Path %s: identical (%llu bp)\n", p->name, (unsigned long long)p->len_a);
                break;
            case DIFF_DIFFER:
                printf("Path %s: differs (%llu bp in A, %llu bp in B)\n", p->name, (unsigned long long)p->len_a, (unsigned long long)p->len_b);
                break;
            case DIFF_A_ONLY:
                printf("Path %s: only in A (%llu bp)\n", p->name, (unsigned long long)p->len_a);
                break;
            default:
                printf("Path %s: only in B (%llu bp)\n", p->name, (unsigned long long)p->len_b);
                break;
        }
    }
}

// what one graph alone carries, id by id
static void print_side(const diff_side_t *s, const char *tag) {
    for (int32_t i = 0; i < s->n_seg; i++) {
        printf("Segment only in %s: %llu\n", tag, (unsigned long long)s->seg[i]);
    }
    for (int32_t i = 0; i < s->n_link; i++) {
        const diff_link_t *l = &s->link[i];
        printf("Link only in %s: %llu%c -> %llu%c (%luM)\n", tag,
               (unsigned long long)l->from, l->from_orient,
               (unsigned long long)l->to, l->to_orient,
               (unsigned long)l->overlap);
    }
}

// `compare gfa` - two graphs that need not agree on segment ids
static int compare_gfa(int argc, char **argv) {
    const char *fn_a, *fn_b;
    int verbose;
    if (!parse_args(argc, argv, &fn_a, &fn_b, &verbose)) return 2;
    if (!want_gfa(fn_a) || !want_gfa(fn_b)) return 2;

    gfa_t *a = gfa_read(fn_a, GFA_LINKS | GFA_PATHS);
    if (!a) return 2;
    gfa_t *b = gfa_read(fn_b, GFA_LINKS | GFA_PATHS);
    if (!b) {
        gfa_destroy(a);
        return 2;
    }

    // 0 the graphs match, 1 they differ, 2 the comparison could not be made
    int ret = 2;
    diff_t *d = diff_graphs(a, b);
    if (d) {
        print_counts(a, b, d);
        print_paths(d);
        if (verbose) {
            print_side(&d->a, "A");
            print_side(&d->b, "B");
        }

        ret = diff_identical(d) ? 0 : 1;
        if (ret) {
            ak_log(AK_LOG_INFO, NULL, "the graphs differ");
        } else {
            ak_log(AK_LOG_INFO, NULL, "the graphs are identical");
        }
        diff_destroy(d);
    }

    gfa_destroy(a);
    gfa_destroy(b);
    return ret;
}

// alignments

// the counts both alignment sets are measured by
static void print_gaf_counts(const diff_gaf_t *d) {
    printf("Alignments A: %lld\n", (long long)d->n_aln_a);
    printf("Alignments B: %lld\n", (long long)d->n_aln_b);
    printf("Alignments shared: %lld\n", (long long)d->n_aln_shared);
    printf("Alignments only in A: %lld\n", (long long)d->n_aln_a_only);
    printf("Alignments only in B: %lld\n", (long long)d->n_aln_b_only);

    printf("Reads A: %lld\n", (long long)d->n_read_a);
    printf("Reads B: %lld\n", (long long)d->n_read_b);
    printf("Reads shared: %lld\n", (long long)d->n_read_shared);
    printf("Reads only in A: %lld\n", (long long)d->n_read_a_only);
    printf("Reads only in B: %lld\n", (long long)d->n_read_b_only);
    printf("Reads on the same path(s): %lld\n", (long long)d->n_read_all_same);
    printf("Reads partly on the same path(s): %lld\n", (long long)d->n_read_partial);
    printf("Reads on no shared path: %lld\n", (long long)d->n_read_none);
}

// every alignment one file has and the other does not, walk and all
static void print_gaf_alns(const diff_gaf_t *d) {
    for (int64_t i = 0; i < d->n_aln; i++) {
        const diff_aln_t *a = &d->aln[i];
        if (a->state == DIFF_ALN_SHARED) continue;

        // a read the other file never aligned is a different thing from a read
        // it aligned elsewhere, and the second is the interesting one
        printf("Alignment only in %s: %s %s%s\n",
               a->state == DIFF_ALN_A_ONLY ? "A" : "B",
               a->qname, a->path,
               a->read_both ? "" : " (read absent from the other file)");
    }
}

// `compare gaf` - two alignment sets over the same graph
static int compare_gaf(int argc, char **argv) {
    const char *fn_a, *fn_b;
    int verbose;
    if (!parse_args(argc, argv, &fn_a, &fn_b, &verbose)) return 2;
    if (!want_gaf(fn_a) || !want_gaf(fn_b)) return 2;

    diff_gaf_t *d = diff_gaf(fn_a, fn_b);
    if (!d) return 2;

    print_gaf_counts(d);
    if (verbose) {
        print_gaf_alns(d);
    }

    // 0 the files agree, 1 they differ, 2 the comparison could not be made
    int ret = diff_gaf_identical(d) ? 0 : 1;
    if (ret) {
        ak_log(AK_LOG_INFO, NULL, "the alignments differ");
    } else {
        ak_log(AK_LOG_INFO, NULL, "the alignments are identical");
    }

    diff_gaf_destroy(d);
    return ret;
}

// `compare` entry point; see cli.h
int cmd_compare(int argc, char **argv) {
    if (argc < 3) {
        usage();
        return 2;
    }

    if (!strcmp(argv[2], "gfa")) return compare_gfa(argc, argv);
    if (!strcmp(argv[2], "gaf")) return compare_gaf(argc, argv);

    ak_log(AK_LOG_ERROR, NULL, "unknown compare target: %s", argv[2]);
    usage();
    return 2;
}
