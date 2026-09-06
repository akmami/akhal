#include "akhal/gfa.h"
#include "akhal/rgfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <string.h>

// print the gfa2rgfa usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal gfa2rgfa <in.gfa> [out.rgfa] [--ref <NAME>]");
}

// 1 when the extension fits, else 0 and the reason is logged
static int want_gfa(const char *fn) {
    if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
    return 0;
}

// `gfa2rgfa` entry point; see cli.h
int cmd_gfa2rgfa(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL, *ref_name = NULL;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--ref") && i + 1 < argc) {
            ref_name = argv[++i];
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
    if (!want_gfa(in)) return 1;
    if (out_fn && !want_gfa(out_fn)) return 1;

    gfa_t *g = gfa_read(in, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;

    rgfa_stat_t st;
    if (rgfa_build(g, ref_name, &st) != AK_OK) {
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

    int ret = gfa_write_rgfa(g, out) == AK_OK ? 0 : 1;
    if (out_fn) {
        fclose(out);
    }

    if (!ret) {
        ak_log(AK_LOG_INFO, NULL, "%d path(s) after merging; %d node(s) on the backbone, %d placed in all, up to rank %d", st.n_path, st.n_rank0, st.n_labelled, st.max_rank);
        if (st.n_ambiguous) {
            ak_log(AK_LOG_WARN, "gfa2rgfa", "%d node(s) left without SN/SO: the paths reaching them disagree on where they sit", st.n_ambiguous);
        }
        if (st.n_unreached) {
            ak_log(AK_LOG_WARN, "gfa2rgfa", "%d node(s) left untagged: no path visits them", st.n_unreached);
        }
    }

    gfa_destroy(g);
    return ret;
}
