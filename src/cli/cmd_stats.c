#include "akhal/gfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/**
 * Min and max degree over segments with non-zero degree
 * @param g Graph to scan
 * @param in Non-zero to inspect in-degree, zero for out-degree
 * @param min_out Set to the minimum degree (-1 if the graph is empty)
 * @param max_out Set to the maximum degree (-1 if the graph is empty)
 */
static void degree_range(const gfa_t *g, int in, int *min_out, int *max_out) {
    if (gfa_n_seg(g) == 0) { *min_out = -1; *max_out = -1; return; }

    int mn = INT_MAX, mx = 0;
    for (int32_t i = 0; i < gfa_n_seg(g); i++) {
        int d = in ? g->seg[i].in_degree : g->seg[i].out_degree;
        if (d) {
            if (d < mn) mn = d;
            if (d > mx) mx = d;
        }
    }
    *min_out = mn;
    *max_out = mx;
}

// `stats` entry point; see cli.h
int cmd_stats(int argc, char **argv) {
    if (argc < 3) {
        ak_log(AK_LOG_ERROR, NULL, "usage: akhal stats <r/GFA>");
        return 1;
    }
    const char *fn = argv[2];
    if (!ak_ends_with(fn, ".gfa") && !ak_ends_with(fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
        return 1;
    }

    gfa_t *g = gfa_read(fn, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;

    int32_t n_seg  = gfa_n_seg(g);
    int32_t n_link = gfa_n_link(g);

    // segment length distribution
    double smean = 0.0, sstd = 0.0;
    size_t smin = 0, smax = 0;
    if (n_seg > 0) {
        size_t *slen = (size_t *)malloc((size_t)n_seg * sizeof(size_t));
        if (!slen) { gfa_destroy(g); ak_log(AK_LOG_ERROR, NULL, "out of memory"); return 1; }
        smin = smax = g->seg[0].len;
        for (int32_t i = 0; i < n_seg; i++) {
            size_t l = g->seg[i].len;
            slen[i] = l;
            if (l < smin) smin = l;
            if (l > smax) smax = l;
        }
        smean = ak_mean(slen, (size_t)n_seg);
        sstd  = ak_stddev(ak_variance(slen, (size_t)n_seg, smean));
        free(slen);
    }

    // link overlap distribution
    double omean = 0.0, ostd = 0.0;
    if (n_link > 0) {
        size_t *ov = (size_t *)malloc((size_t)n_link * sizeof(size_t));
        if (!ov) { gfa_destroy(g); ak_log(AK_LOG_ERROR, NULL, "out of memory"); return 1; }
        for (int32_t i = 0; i < n_link; i++) ov[i] = g->link[i].overlap;
        omean = ak_mean(ov, (size_t)n_link);
        ostd  = ak_stddev(ak_variance(ov, (size_t)n_link, omean));
        free(ov);
    }

    int min_in, max_in, min_out, max_out;
    degree_range(g, 1, &min_in, &max_in);
    degree_range(g, 0, &min_out, &max_out);

    printf("Segment count: %ld\n", (long)n_seg);
    printf("Rank 0 segment count: %ld\n", (long)g->n_path_seg);
    printf("Rank 0< segment count: %ld\n", (long)((long long)n_seg - (long long)g->n_path_seg));
    printf("Segment avg length: %f\n", smean);
    printf("Segment std length: %f\n", sstd);
    printf("Segment min. length %lu\n", (unsigned long)smin);
    printf("Segment max. length %lu\n", (unsigned long)smax);
    printf("Link count: %ld\n", (long)n_link);
    printf("Link overlapping avg length: %f\n", omean);
    printf("Link overlapping std length: %f\n", ostd);
    printf("Minimum in degree: %d\n", min_in);
    printf("Maximum in degree: %d\n", max_in);
    printf("Minimum out degree: %d\n", min_out);
    printf("Maximum out degree: %d\n", max_out);

    gfa_destroy(g);
    return 0;
}
