/**
 * `akhal vg2gfa <in.vg> [out.gfa]` - convert vg's native format to GFA.
 *
 * The whole graph is read from the .vg (see the vg module: no protobuf runtime
 * required, only zlib) and emitted as GFA1: an H line, one S per node, one L
 * per edge, and one P per path. This mirrors `vg view -g`. If no output path is
 * given, the GFA is written to stdout.
 */

#include "akhal/vg.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>

/** 
 * @brief Emit a graph as GFA1 (H, S, L, P). 
 */
static void write_gfa(FILE *out, const vg_graph_t *g) {
    fprintf(out, "H\tVN:Z:1.0\n");

    for (int32_t i = 0; i < g->n_node; i++)
        fprintf(out, "S\t%lld\t%s\n", (long long)g->node[i].id, g->node[i].seq ? g->node[i].seq : "*");

    for (int32_t i = 0; i < g->n_edge; i++) {
        const vg_edge_t *e = &g->edge[i];
        fprintf(out, "L\t%lld\t%c\t%lld\t%c\t%dM\n",
                (long long)e->from, e->from_start ? '-' : '+',
                (long long)e->to,   e->to_end     ? '-' : '+',
                e->overlap);
    }

    for (int32_t i = 0; i < g->n_path; i++) {
        const vg_path_t *p = &g->path[i];
        fprintf(out, "P\t%s\t", p->name ? p->name : "");
        for (int32_t s = 0; s < p->n_step; s++)
            fprintf(out, "%s%lld%c", s ? "," : "", (long long)p->step[s].node_id, p->step[s].is_reverse ? '-' : '+');
        fprintf(out, "\t*\n");
    }
}

/** `vg2gfa` entry point; see cli.h. */
int cmd_vg2gfa(int argc, char **argv) {
    if (argc < 3) {
        ak_log(AK_LOG_ERROR, NULL, "usage: akhal vg2gfa <in.vg> [out.gfa]");
        return 1;
    }
    const char *in = argv[2];
    const char *out_fn = (argc >= 4) ? argv[3] : NULL;

    if (!ak_ends_with(in, ".vg")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .vg file: %s", in);
        return 1;
    }
    if (out_fn && !ak_ends_with(out_fn, ".gfa") && !ak_ends_with(out_fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .gfa/.rgfa file: %s", out_fn);
        return 1;
    }

    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) { ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn); return 1; }
    }

    vg_graph_t *g = vg_read(in);
    if (!g) { if (out_fn) fclose(out); return 1; }

    write_gfa(out, g);

    ak_log(AK_LOG_INFO, NULL, "converted %s (%d nodes, %d edges, %d paths)", in, g->n_node, g->n_edge, g->n_path);

    vg_graph_destroy(g);
    if (out_fn) fclose(out);
    return 0;
}
