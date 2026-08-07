/**
 * `akhal sort <in.gfa> [out.gfa]` - topologically sort and renumber a graph.
 *
 * The nodes are ordered with gfa_toposort (Kahn's algorithm; ties broken
 * alphabetically by node sequence content), then renumbered 1..N in that order.
 * Ordering by content keeps the result independent of the input's node
 * numbering. The graph is
 * re-emitted as GFA1 with the new ids: S lines in sorted order, L lines with
 * remapped endpoints (orientation + overlap preserved), and P lines with
 * remapped, oriented steps. If no output path is given, stdout is used.
 *
 * Note: only the SR tag is carried over on S lines; SN/SO are not preserved
 * (the in-memory model does not keep them verbatim).
 */

#include "akhal/gfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>

// One remapped L line, used to emit links in a deterministic order.
typedef struct {
    long     from, to;
    char     from_orient, to_orient;
    unsigned overlap;
} lrec_t;

/** Order links by (from, to) so the output is deterministic. */
static int lrec_cmp(const void *A, const void *B) {
    const lrec_t *a = (const lrec_t *)A, *b = (const lrec_t *)B;
    if (a->from != b->from) return a->from < b->from ? -1 : 1;
    if (a->to   != b->to)   return a->to   < b->to   ? -1 : 1;
    return 0;
}

/** `sort` entry point; see cli.h. */
int cmd_sort(int argc, char **argv) {
    if (argc < 3) {
        ak_log(AK_LOG_ERROR, NULL, "usage: akhal sort <in.gfa> [out.gfa]");
        return 1;
    }
    const char *in = argv[2];
    const char *out_fn = (argc >= 4) ? argv[3] : NULL;

    if (!ak_ends_with(in, ".gfa") && !ak_ends_with(in, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", in);
        return 1;
    }
    if (out_fn && !ak_ends_with(out_fn, ".gfa") && !ak_ends_with(out_fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .gfa/.rgfa file: %s", out_fn);
        return 1;
    }

    gfa_t *g = gfa_read(in, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;

    int32_t n = gfa_n_seg(g);

    int32_t *order = (int32_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int32_t));
    int32_t *newid = (int32_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int32_t));
    if (!order || !newid) {
        free(order); free(newid); gfa_destroy(g);
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return 1;
    }

    int32_t placed = gfa_toposort(g, order);
    if (placed < 0) {
        free(order); free(newid); gfa_destroy(g);
        return 1;
    }
    if (placed < n)
        ak_log(AK_LOG_WARN, "sort", "graph is cyclic; %d node(s) in cycles appended after the acyclic prefix", n - placed);

    // new id (1..N) for each old segment index
    for (int32_t pos = 0; pos < n; pos++) newid[order[pos]] = pos + 1;

    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) {
            free(order); free(newid); gfa_destroy(g);
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            return 1;
        }
    }

    fprintf(out, "H\tVN:Z:1.0\n");

    // S lines in sorted order.
    for (int32_t pos = 0; pos < n; pos++) {
        gfa_seg_t *s = gfa_seg_at(g, order[pos]);
        fprintf(out, "S\t%d\t%s", pos + 1, s->seq ? s->seq : "*");
        if (s->rank >= 0) fprintf(out, "\tSR:i:%d", s->rank);
        fputc('\n', out);
    }

    // L lines, remapped and sorted by (from, to).
    int32_t n_link = gfa_n_link(g);
    if (n_link > 0) {
        lrec_t *ls = (lrec_t *)malloc((size_t)n_link * sizeof(lrec_t));
        if (!ls) {
            if (out_fn) fclose(out);
            free(order); free(newid); gfa_destroy(g);
            ak_log(AK_LOG_ERROR, NULL, "out of memory");
            return 1;
        }
        for (int32_t i = 0; i < n_link; i++) {
            gfa_link_t *e = gfa_link_at(g, i);
            ls[i].from = newid[e->v];
            ls[i].to   = newid[e->w];
            ls[i].from_orient = e->from_orient;
            ls[i].to_orient   = e->to_orient;
            ls[i].overlap = e->overlap;
        }
        qsort(ls, (size_t)n_link, sizeof(lrec_t), lrec_cmp);
        for (int32_t i = 0; i < n_link; i++)
            fprintf(out, "L\t%ld\t%c\t%ld\t%c\t%uM\n", ls[i].from, ls[i].from_orient, ls[i].to, ls[i].to_orient, ls[i].overlap);
        free(ls);
    }

    // P lines, remapped and oriented.
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        fprintf(out, "P\t%s\t", gfa_path_name(g, k));
        const uint32_t *segs;
        int ns = gfa_path_segs(g, k, &segs);
        const char *ori = g->path_ori + g->path_off[k];
        int written = 0;
        for (int i = 0; i < ns; i++) {
            if (segs[i] == GFA_NIL) continue;
            fprintf(out, "%s%d%c", written ? "," : "", newid[segs[i]], ori[i]);
            written = 1;
        }
        fprintf(out, "\t*\n");
    }

    if (out_fn) fclose(out);
    ak_log(AK_LOG_INFO, NULL, "sorted %s (%d nodes, %d links, %d paths)", in, n, n_link, gfa_n_path(g));

    free(order);
    free(newid);
    gfa_destroy(g);
    return 0;
}
