#include "akhal/rgfa.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

// the path to label rank 0: the one asked for by name, else the first
static int32_t backbone_path(const gfa_t *g, const char *ref_name) {
    if (!ref_name) return gfa_n_path(g) > 0 ? 0 : -1;
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        if (!strcmp(gfa_path_name(g, k), ref_name)) return k;
    }
    return -1;
}

// how many segments a path actually resolves; a P line whose ids the file
// never defined spells nothing and cannot be a backbone
static int64_t path_resolved(const gfa_t *g, int32_t k) {
    const uint32_t *segs;
    int ns = gfa_path_segs(g, k, &segs);
    int64_t n = 0;
    for (int i = 0; i < ns; i++) if (segs[i] != GFA_NIL) n++;
    return n;
}

// labelling

// drop whatever the reader and gfa_add_path() left behind: the walks below are
// the only thing allowed to place a segment
static void unlabel(gfa_t *g) {
    for (int32_t i = 0; i < gfa_n_seg(g); i++) {
        gfa_seg_t *s = gfa_seg_at(g, i);
        s->rank     = -1;
        s->ref_path = -1;
        s->start    = -1;
    }
}

// rank r, on the path at index `pi` at offset `off`
static void place(gfa_seg_t *s, int r, int32_t pi, int64_t off) {
    // start is int32_t, so an offset past its reach is no offset at all
    if (pi < 0 || off < 0 || off + s->len > INT32_MAX) {
        s->rank     = r;
        s->ref_path = -1;
        s->start    = -1;
        return;
    }
    s->rank     = r;
    s->ref_path = pi;
    s->start    = (int32_t)off;
}

// rank r and nothing else: reached, but with no offset anyone can stand behind
static void place_ranked(gfa_seg_t *s, int r) {
    place(s, r, -1, -1);
}

// the backbone: rank 0, offsets running the length of the walk. A segment the
// walk comes back to keeps the offset of its first visit
static void label_backbone(gfa_t *g, int32_t bb) {
    const int32_t pi = bb;
    const uint32_t *segs;
    int ns = gfa_path_segs(g, bb, &segs);

    int64_t off = 0;
    for (int i = 0; i < ns; i++) {
        if (segs[i] == GFA_NIL) continue;
        gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
        if (s->rank < 0) place(s, 0, pi, off);
        off += s->len;
    }
}

// One path, labelling only what it alone explains. `cur` is the offset the next
// new segment takes and `floor` the rank of the ground the walk last stood on.
// `anchored` says `cur` came from a segment somebody had placed rather than
// from this path's own start, and `amb` says the walk has lost its place: it
// may still hand out ranks, but no offsets, until it reaches somewhere another
// path has already pinned down
static void label_path(gfa_t *g, int32_t k) {
    const int32_t pi = k;
    const uint32_t *segs;
    int ns = gfa_path_segs(g, k, &segs);

    int64_t cur = 0;   // a path starting off the backbone counts from its own 0
    int floor = 0, amb = 0, anchored = 0;

    for (int i = 0; i < ns; i++) {
        if (segs[i] == GFA_NIL) continue;
        gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);

        if (s->rank >= 0) {
            // ground an earlier walk already covered
            if (s->start < 0) {
                amb = 1;                    // stepped through something unplaced
            } else if (amb || !anchored || s->rank == 0 || s->ref_path == pi) {
                // the backbone settles it, a segment this same path placed is
                // its own business, and before the first anchor there is
                // nothing to disagree with
                amb = 0;
            } else if (s->start != cur) {
                // another path put this somewhere else, and nothing
                // authoritative settles which offset is the real one
                place_ranked(s, s->rank);
                amb = 1;
            }
            floor = s->rank;
            if (s->start >= 0) {
                cur = (int64_t)s->start + s->len;
                anchored = 1;
            }
        } else if (amb) {
            place_ranked(s, floor + 1);     // a detour off ground we cannot place
        } else {
            place(s, floor + 1, pi, cur); // one rank deeper, carrying on from
            cur += s->len;                  // where the path left that ground
        }
    }
}

static void tally(const gfa_t *g, rgfa_stat_t *st) {
    memset(st, 0, sizeof(*st));
    st->n_path = gfa_n_path(g);

    for (int32_t i = 0; i < gfa_n_seg(g); i++) {
        const gfa_seg_t *s = gfa_seg_at(g, i);
        if (s->rank < 0) {
            st->n_unreached++;
            continue;
        }
        if (s->rank == 0) {
            st->n_rank0++;
        }
        if (s->rank > st->max_rank) {
            st->max_rank = s->rank;
        }
        if (s->start >= 0) {
            st->n_labelled++;
        } else {
            st->n_ambiguous++;
        }
    }
}

// label a graph as rGFA; see akhal/rgfa.h
int rgfa_build(gfa_t *g, const char *ref_name, rgfa_stat_t *st) {
    if (!(g->flags & GFA_PATHS)) {
        ak_log(AK_LOG_ERROR, "rgfa", "labelling requires the graph to be read with GFA_PATHS");
        return AK_EINVAL;
    }
    if (gfa_n_path(g) == 0) {
        ak_log(AK_LOG_ERROR, "rgfa", "graph has no P lines, so there is no backbone to label against");
        return AK_EINVAL;
    }

    int32_t bb = backbone_path(g, ref_name);
    if (bb < 0) {
        ak_log(AK_LOG_ERROR, "rgfa", "no path named '%s' in the graph", ref_name ? ref_name : "");
        return AK_EINVAL;
    }
    if (path_resolved(g, bb) == 0) {
        ak_log(AK_LOG_ERROR, "rgfa", "the backbone path resolves to no segments");
        return AK_EINVAL;
    }
    if (g->has_sr) {
        ak_log(AK_LOG_WARN, "rgfa", "the file carries its own SR tags; they are replaced by what the paths say");
    }

    unlabel(g);
    label_backbone(g, bb);
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        if (k != bb) label_path(g, k);
    }

    if (st) tally(g, st);
    return AK_OK;
}
