#include "akhal/gfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include "ksort.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Links are emitted in the order their endpoints came out of the sort, so the L block follows the S block
#define LINK_KEY(e) (((uint64_t)(e).v << 33) | ((uint64_t)(e).w << 2) | ((e).from_orient == '-' ? 2u : 0u) | ((e).to_orient == '-' ? 1u : 0u))

KRADIX_SORT_INIT(link, gfa_link_t, LINK_KEY, 8)

// sort the links in place, starting from the highest byte the keys reach
// rather than the top of a 64-bit word, which on any real graph is empty
static void sort_links(gfa_link_t *link, int32_t n_link, int32_t n_seg) {
    if (n_link <= RS_MIN_SIZE) {
        if (n_link > 1) rs_insertsort_link(link, link + n_link);
        return;
    }
    uint64_t max_key = ((uint64_t)n_seg << 33) | ((uint64_t)n_seg << 2) | 3u;
    int top = 0;
    while (top < 7 && (max_key >> (8 * (top + 1))) != 0) top++;
    rs_sort_link(link, link + n_link, RS_MAX_BITS, top * RS_MAX_BITS);
}

static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal sort <in.gfa> [out.gfa] [--tie seq|id] [--no-renumber] [--verbose] [--footprint]");
}

// `sort` entry point; see cli.h
int cmd_sort(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL;
    int renumber = 1, report = 0, tie = GFA_TIE_SEQ;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--no-renumber")) {
            renumber = 0;
        } else if (!strcmp(argv[i], "--verbose")) {
            report |= GFA_VERBOSE;
        } else if (!strcmp(argv[i], "--tie")) {
            if (++i >= argc) {
                usage();
                return 1;
            }
            if (!strcmp(argv[i], "seq")) {
                tie = GFA_TIE_SEQ;
            } else if (!strcmp(argv[i], "id")) {
                tie = GFA_TIE_ID;
            } else {
                ak_log(AK_LOG_ERROR, NULL, "--tie takes seq or id, not %s", argv[i]);
                return 1;
            }
        } else if (!strcmp(argv[i], "--footprint")) {
            report |= GFA_FOOTPRINT;
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
    if (out_fn && !ak_ends_with(out_fn, ".gfa") && !ak_ends_with(out_fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .gfa/.rgfa file: %s", out_fn);
        return 1;
    }

    // No GFA_DEGREES: the sort counts the in-degrees into the one array it
    // consumes anyway, so the graph's own pair would only be a second copy.
    // GFA_SEQ is here for the S lines, which carry the bases whatever breaks
    // a tie - so --tie id spends less time, not less memory
    gfa_t *g = gfa_read(in, GFA_SEGS | GFA_SEQ | GFA_LINKS | GFA_ARCS | GFA_PATHS | report);
    if (!g) return 1;

    int32_t n = gfa_n_seg(g), n_link = gfa_n_link(g), n_path = gfa_n_path(g);

    // order: the new sequence of segments. pos: where each landed in it,
    // 1-based, which is the id it takes when the graph is renumbered
    int32_t *order = (int32_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int32_t));
    if (!order) {
        gfa_destroy(g);
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return 1;
    }

    int32_t placed = gfa_toposort(g, order, tie);
    if (placed < 0) {
        free(order);
        gfa_destroy(g);
        return 1;
    }
    if (placed < n) {
        ak_log(AK_LOG_WARN, "sort", "graph is cyclic; %d node(s) in cycles appended after the acyclic prefix", n - placed);
    }

    int32_t *pos = (int32_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int32_t));
    if (!pos) {
        free(order);
        gfa_destroy(g);
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return 1;
    }

    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) {
            free(order);
            free(pos);
            gfa_destroy(g);
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            return 1;
        }
    }

    fprintf(out, "H\tVN:Z:1.0\n");

    // S lines in the new order, filling pos as they go
    for (int32_t p = 0; p < n; p++) {
        int32_t i = order[p];
        pos[i] = p + 1;
        const gfa_seg_t *s = gfa_seg_at(g, i);
        fprintf(out, "S\t%llu\t%s", renumber ? (unsigned long long)(p + 1) : (unsigned long long)s->id, s->seq ? s->seq : "*");
        // only round-trip ranks the input actually carried; ranks the reader
        // derived for a plain GFA are not the file's own and are not emitted
        if (g->has_sr && s->rank >= 0) {
            fprintf(out, "\tSR:i:%d", s->rank);
        }
        fputc('\n', out);
    }

    // the bases are written, and on a whole-genome graph they are the largest thing held; the ids and ranks that the L and P lines still need stay
    gfa_drop(g, GFA_SEQ);

    // L lines. The endpoints are replaced by their positions in place and the
    // links sorted where they lie, so the block costs no second copy of them
    if (n_link > 0) {
        gfa_link_t *link = g->link;
        for (int32_t k = 0; k < n_link; k++) {
            link[k].v = (uint32_t)pos[link[k].v];
            link[k].w = (uint32_t)pos[link[k].w];
        }
        sort_links(link, n_link, n);
        for (int32_t k = 0; k < n_link; k++) {
            const gfa_link_t *e = &link[k];
            // the positions are back to ids here: the new one is the position
            // itself, the old one is on the segment that landed there
            unsigned long long from = renumber ? (unsigned long long)e->v : (unsigned long long)gfa_seg_at(g, order[e->v - 1])->id;
            unsigned long long to   = renumber ? (unsigned long long)e->w : (unsigned long long)gfa_seg_at(g, order[e->w - 1])->id;
            fprintf(out, "L\t%llu\t%c\t%llu\t%c\t%uM\n", from, e->from_orient, to, e->to_orient, e->overlap);
        }
        gfa_drop(g, GFA_LINKS);
    }
    free(order);

    // P lines, in the order the file gave them, with their steps remapped
    for (int32_t k = 0; k < n_path; k++) {
        fprintf(out, "P\t%s\t", gfa_path_name(g, k));
        const uint32_t *segs;
        int ns = gfa_path_segs(g, k, &segs);
        const char *ori = g->path_ori + g->path_off[k];
        int written = 0;
        for (int i = 0; i < ns; i++) {
            if (segs[i] == GFA_NIL) continue;
            unsigned long long id = renumber ? (unsigned long long)pos[segs[i]] : (unsigned long long)gfa_seg_at(g, (int32_t)segs[i])->id;
            fprintf(out, "%s%llu%c", written ? "," : "", id, ori[i]);
            written = 1;
        }
        fprintf(out, "\t*\n");
    }

    free(pos);
    if (out_fn) {
        fclose(out);
    }
    ak_log(AK_LOG_INFO, NULL, "sorted %s (%d nodes, %d links, %d paths)", in, n, n_link, n_path);

    gfa_destroy(g);
    return 0;
}
