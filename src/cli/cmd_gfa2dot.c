#include "akhal/gfa.h"
#include "akhal/dot.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <string.h>

// print the gfa2dot usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal gfa2dot <in.gfa> [out.dot] [--ids]");
}

// `gfa2dot` entry point; see cli.h
int cmd_gfa2dot(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL;
    int flags = 0;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--ids")) {
            flags |= DOT_IDS;
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
    if (!in) {
        usage();
        return 1;
    }
    if (!ak_ends_with(in, ".gfa") && !ak_ends_with(in, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", in);
        return 1;
    }
    if (out_fn && !ak_ends_with(out_fn, ".dot") && !ak_ends_with(out_fn, ".gv")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .dot/.gv file: %s", out_fn);
        return 1;
    }

    gfa_t *g = gfa_read(in, GFA_LINKS);
    if (!g) return 1;

    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) {
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            gfa_destroy(g);
            return 1;
        }
    }

    int ret = dot_write(g, out, flags) == AK_OK ? 0 : 1;
    if (out_fn) {
        fclose(out);
    }

    if (!ret && out_fn) {
        ak_log(AK_LOG_INFO, NULL, "wrote %d node(s) and %d edge(s) to %s; render with: dot -Tpng %s -o graph.png", gfa_n_seg(g), gfa_n_link(g), out_fn, out_fn);
    }

    gfa_destroy(g);
    return ret;
}
