#include "akhal/vg.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>

// emit a graph as GFA1 (H, S, L, P)
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

// `vg2gfa` entry point; see cli.h
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
        if (!out) {
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            return 1;
        }
    }

    vg_graph_t *g = vg_read(in);
    if (!g) {
        if (out_fn) {
            fclose(out);
        }
        return 1;
    }

    write_gfa(out, g);

    // Flush and close before tearing the graph down, not after: the teardown
    // walks the whole graph, and a process that dies in there would otherwise
    // leave the last buffered lines unwritten and the file never closed.
    int rc = 0;
    if (out_fn) {
        if (fclose(out) != 0) {
            ak_log(AK_LOG_ERROR, NULL, "write failed on %s", out_fn);
            rc = 1;
        }
    } else if (fflush(stdout) != 0) {
        ak_log(AK_LOG_ERROR, NULL, "write failed on stdout");
        rc = 1;
    }

    if (rc == 0) {
        ak_log(AK_LOG_INFO, NULL, "converted %s (%d nodes, %d edges, %d paths)", in, g->n_node, g->n_edge, g->n_path);
    }

    vg_graph_destroy(g);
    return rc;
}
