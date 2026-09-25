#include "akhal/gfa.h"
#include "akhal/rgfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// print the gfa2rgfa usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal gfa2rgfa <in.gfa> [out.rgfa] [--ref <NAME>[,<NAME>...] | --ref all] [--verbose] [--footprint]");
}

// 1 when the extension fits, else 0 and the reason is logged
static int want_gfa(const char *fn) {
    if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa")) return 1;
    ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
    return 0;
}

// The backbones --ref selects, as indices into the graph's paths, picked by ak_select()
static int pick_backbones(const gfa_t *g, const char *spec, const char *src, int32_t **out, int32_t *n_out) {
    int32_t n_path = gfa_n_path(g);
    int32_t *bb = (int32_t *)malloc((size_t)(n_path > 0 ? n_path : 1) * sizeof(*bb));
    if (!bb) {
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return AK_ENOMEM;
    }

    const char *what = "";
    int what_len = 0;
    int32_t n = 0;
    int verdict = ak_select((const char *const *)g->path, n_path, spec, bb, &n, &what, &what_len);
    if (verdict != AK_SELECT_OK) {
        char msg[512];
        ak_log(AK_LOG_ERROR, NULL, "--ref: %s in %s%s", ak_select_msg(msg, sizeof(msg), verdict, what, what_len, "P line"), src,
               verdict == AK_SELECT_SHARED ? " - a path written as fragments has to be joined first" : "");
        free(bb);
        return AK_EINVAL;
    }

    *out = bb;
    *n_out = n;
    return AK_OK;
}

// the labelled graph, to out_fn or to stdout when there is none
static int write_rgfa(const gfa_t *g, const char *out_fn) {
    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) {
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            return AK_EOPEN;
        }
    }
    int rc = gfa_write_rgfa(g, out);
    if (out_fn && fclose(out) != 0 && rc == AK_OK) {
        rc = AK_EIO;
    }
    return rc;
}

// `gfa2rgfa` entry point; see cli.h
int cmd_gfa2rgfa(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL, *ref_arg = NULL;
    int report = cli_take_report(&argc, argv, 2);

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--ref") && i + 1 < argc) {
            ref_arg = argv[++i];
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

    // A list's own shape - no empty or repeated name - is checked before the graph is read
    const char *what = "";
    int what_len = 0;
    int verdict = ak_select_check(ref_arg, &what, &what_len);
    if (verdict != AK_SELECT_OK) {
        char msg[512];
        ak_log(AK_LOG_ERROR, NULL, "--ref: %s", ak_select_msg(msg, sizeof(msg), verdict, what, what_len, "name"));
        return 1;
    }

    gfa_t *g = gfa_read(in, GFA_ALL | report);
    if (!g) return 1;

    int32_t *bb = NULL, n_bb = 0;
    rgfa_stat_t st;
    int rc = pick_backbones(g, ref_arg, in, &bb, &n_bb);
    if (rc == AK_OK) rc = rgfa_build(g, bb, n_bb, &st);
    if (rc == AK_OK) rc = write_rgfa(g, out_fn);

    if (rc == AK_OK) {
        ak_log(AK_LOG_INFO, NULL, "%d path(s), %d backbone(s); %d node(s) at rank 0, %d placed in all, up to rank %d",
               st.n_path, n_bb, st.n_rank0, st.n_labelled, st.max_rank);
        if (st.n_ambiguous) {
            ak_log(AK_LOG_WARN, "gfa2rgfa", "%d node(s) left without SN/SO: the paths reaching them disagree on where they sit", st.n_ambiguous);
        }
        if (st.n_unreached) {
            ak_log(AK_LOG_WARN, "gfa2rgfa", "%d node(s) left untagged: no path visits them", st.n_unreached);
        }
    }

    free(bb);
    gfa_destroy(g);
    return rc == AK_OK ? 0 : 1;
}
