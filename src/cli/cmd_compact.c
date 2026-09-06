#include "akhal/gfa.h"
#include "akhal/compact.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>

// 1 when the extension fits, else 0 and the reason is logged
static int want_gfa(const char *fn) {
    if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
    return 0;
}

// `compact` entry point; see cli.h
int cmd_compact(int argc, char **argv) {
    if (argc < 3 || argc > 4) {
        ak_log(AK_LOG_ERROR, NULL, "usage: akhal compact <in.gfa> [out.gfa]");
        return 1;
    }
    const char *in = argv[2];
    const char *out_fn = (argc >= 4) ? argv[3] : NULL;

    if (!want_gfa(in)) return 1;
    if (out_fn && !want_gfa(out_fn)) return 1;

    gfa_t *g = gfa_read(in, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;

    compact_t *c = compact_runs(g);
    if (!c) {
        gfa_destroy(g);
        return 1;
    }

    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) {
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            compact_destroy(c);
            gfa_destroy(g);
            return 1;
        }
    }

    int ret = compact_write(g, c, out) == AK_OK ? 0 : 1;
    if (out_fn) {
        fclose(out);
    }

    if (!ret) {
        ak_log(AK_LOG_INFO, NULL, "%d segment(s) in, %d out: %d folded into a neighbour", c->n_seg, c->n_run, c->n_merged);
    }

    compact_destroy(c);
    gfa_destroy(g);
    return ret;
}
