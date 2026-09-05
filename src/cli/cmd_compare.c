#include "akhal/gfa.h"
#include "akhal/diff.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <string.h>

// print the compare usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal compare <A.gfa> <B.gfa> [--verbose]");
}

// 1 when the extension fits, else 0 and the reason is logged
static int want_gfa(const char *fn) {
    if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
    return 0;
}

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

// `compare` entry point; see cli.h
int cmd_compare(int argc, char **argv) {
    if (argc < 4 || argc > 5) {
        usage();
        return 2;
    }
    const char *fn_a = argv[2], *fn_b = argv[3];

    int verbose = 0;
    if (argc == 5) {
        if (strcmp(argv[4], "--verbose") != 0) {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[4]);
            usage();
            return 2;
        }
        verbose = 1;
    }
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
