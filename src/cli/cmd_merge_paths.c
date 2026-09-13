#include "akhal/gfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// print the merge-paths usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal merge-paths <in.gfa> [out.gfa]");
}

// one chain flattened out of the merge set, ready to become a single P line
typedef struct {
    char     *name;   // owned
    uint32_t *seg;    // owned: segment indices, GFA_NIL entries dropped
    char     *ori;    // owned: the orientation each step carries
    int64_t   n;
} chain_t;

static void chains_free(chain_t *c, int32_t n) {
    if (!c) return;
    for (int32_t i = 0; i < n; i++) {
        free(c[i].name);
        free(c[i].seg);
        free(c[i].ori);
    }
    free(c);
}

// Flatten one chain's fragments into a single ordered walk.
//
// gfa_merge_segs() hands back the segment indices but not the orientation each
// step was visited with, and gfa_add_path() reads a NULL orientation array as
// all-'+' - so a path traversing a node in reverse would come back out
// forward. The orientations sit beside the segments in the path CSR, so take
// both.
static int chain_take(const gfa_t *g, const gfa_merge_t *m, int32_t k, chain_t *out) {
    int64_t n = 0;
    for (int32_t f = m->off[k]; f < m->off[k + 1]; f++) {
        const uint32_t *segs;
        int ns = gfa_path_segs(g, m->frag[f], &segs);
        for (int t = 0; t < ns; t++) {
            if (segs[t] != GFA_NIL) n++;
        }
    }

    out->name = strdup(m->name[k]);
    out->seg  = (uint32_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(uint32_t));
    out->ori  = (char *)malloc((size_t)(n > 0 ? n : 1));
    if (!out->name || !out->seg || !out->ori) return AK_ENOMEM;

    int64_t i = 0;
    for (int32_t f = m->off[k]; f < m->off[k + 1]; f++) {
        int32_t pi = m->frag[f];
        const uint32_t *segs;
        int ns = gfa_path_segs(g, pi, &segs);
        const char *ori = g->path_ori + g->path_off[pi];
        for (int t = 0; t < ns; t++) {
            if (segs[t] == GFA_NIL) continue;
            out->seg[i] = segs[t];
            out->ori[i] = ori[t];
            i++;
        }
    }
    out->n = n;
    return AK_OK;
}

// `merge-paths` entry point; see cli.h
int cmd_merge_paths(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL;

    for (int i = 2; i < argc; i++) {
        if (argv[i][0] == '-') {
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
    if (!in) {
        usage();
        return 1;
    }
    if (!ak_ends_with(in, ".gfa") && !ak_ends_with(in, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", in);
        return 1;
    }
    if (out_fn && !ak_ends_with(out_fn, ".gfa") && !ak_ends_with(out_fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .gfa/.rgfa file: %s", out_fn);
        return 1;
    }

    // chaining walks the L lines, so the links are needed as well as the paths
    gfa_t *g = gfa_read(in, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;

    if (gfa_n_path(g) == 0) {
        ak_log(AK_LOG_ERROR, NULL, "%s has no P lines to merge", in);
        gfa_destroy(g);
        return 1;
    }

    // NULL takes every path and groups each base name separately, which is
    // exactly "merge the fragments that share a name"
    gfa_merge_t *m = gfa_path_merge(g, NULL);
    if (!m) {
        ak_log(AK_LOG_ERROR, NULL, "cannot group the paths of %s", in);
        gfa_destroy(g);
        return 1;
    }

    int32_t n_before = gfa_n_path(g), n_after = m->n;

    // flatten every chain before touching the path block, since flattening
    // reads the very paths that clearing it would free
    chain_t *c = (chain_t *)calloc((size_t)(m->n > 0 ? m->n : 1), sizeof(chain_t));
    int rc = c ? AK_OK : AK_ENOMEM;
    for (int32_t k = 0; rc == AK_OK && k < m->n; k++) {
        rc = chain_take(g, m, k, &c[k]);
    }

    if (rc == AK_OK) {
        gfa_clear_paths(g);
        for (int32_t k = 0; rc == AK_OK && k < m->n; k++) {
            rc = gfa_add_path(g, c[k].name, c[k].seg, c[k].ori, c[k].n);
        }
    }

    chains_free(c, m->n);
    gfa_merge_destroy(m);

    if (rc != AK_OK) {
        ak_log(AK_LOG_ERROR, NULL, "cannot merge the paths of %s: %s", in, ak_strerror(rc));
        gfa_destroy(g);
        return 1;
    }

    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) {
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            gfa_destroy(g);
            return 1;
        }
    }

    rc = gfa_write(g, out);
    if (out_fn && fclose(out) != 0) {
        rc = AK_EIO;
    }

    if (rc == AK_OK) {
        ak_log(AK_LOG_INFO, NULL, "%d P line(s) merged into %d path(s)", n_before, n_after);
    } else {
        ak_log(AK_LOG_ERROR, NULL, "write failed on %s: %s", out_fn ? out_fn : "stdout", ak_strerror(rc));
    }

    gfa_destroy(g);
    return rc == AK_OK ? 0 : 1;
}
