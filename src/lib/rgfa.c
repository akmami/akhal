#include "akhal/rgfa.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

// path consolidation

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

// flatten one chain's fragments into a single ordered walk
static int chain_take(const gfa_t *g, const gfa_merge_t *m, int32_t k, chain_t *out) {
    int64_t n = 0;
    for (int32_t f = m->off[k]; f < m->off[k + 1]; f++) {
        const uint32_t *segs;
        int ns = gfa_path_segs(g, m->frag[f], &segs);
        for (int t = 0; t < ns; t++)
            if (segs[t] != GFA_NIL) n++;
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

// the chain to label rank 0: the one asked for by name - as a whole chain or
// as one of the fragments it swallowed - or the one holding the first P line
static int32_t backbone_chain(const gfa_t *g, const gfa_merge_t *m, const char *ref_name) {
    for (int32_t k = 0; k < m->n; k++) {
        if (ref_name && !strcmp(m->name[k], ref_name)) return k;
        for (int32_t f = m->off[k]; f < m->off[k + 1]; f++) {
            if (!ref_name) {
                if (m->frag[f] == 0) return k;
            } else if (!strcmp(g->path[m->frag[f]], ref_name)) {
                return k;
            }
        }
    }
    return -1;
}

// swap the graph's P lines for one line per chain. A chain that spells nothing
// - every entry an id the file never defined - is dropped rather than written
// out as an empty P line, so `bb` is shifted down past any that went missing
static int install_chains(gfa_t *g, const chain_t *c, int32_t n, int32_t *bb) {
    if (c[*bb].n == 0) {
        ak_log(AK_LOG_ERROR, "rgfa", "the backbone path resolves to no segments");
        return AK_EINVAL;
    }
    int32_t skipped = 0;
    for (int32_t k = 0; k < *bb; k++) {
        if (c[k].n == 0) skipped++;
    }
    *bb -= skipped;

    gfa_clear_paths(g);
    for (int32_t k = 0; k < n; k++) {
        if (c[k].n == 0) continue;
        int rc = gfa_add_path(g, c[k].name, c[k].seg, c[k].ori, c[k].n);
        if (rc != AK_OK) return rc;
    }
    return AK_OK;
}

// labelling

// drop whatever the reader and gfa_add_path() left behind: the walks below are
// the only thing allowed to place a segment
static void unlabel(gfa_t *g) {
    for (int32_t i = 0; i < gfa_n_seg(g); i++) {
        gfa_seg_t *s = gfa_seg_at(g, i);
        s->rank     = -1;
        s->ref_name = NULL;
        s->start    = -1;
        s->end      = -1;
    }
}

// rank r, on sequence `name` at offset `off`
static void place(gfa_seg_t *s, int r, const char *name, int64_t off) {
    // start/end are int32_t, so an offset past their reach is no offset at all
    if (off < 0 || off + s->len > INT32_MAX) {
        s->rank     = r;
        s->ref_name = NULL;
        s->start    = -1;
        s->end      = -1;
        return;
    }
    s->rank     = r;
    s->ref_name = name;
    s->start    = (int32_t)off;
    s->end      = (int32_t)(off + s->len);
}

// rank r and nothing else: reached, but with no offset anyone can stand behind
static void place_ranked(gfa_seg_t *s, int r) {
    place(s, r, NULL, -1);
}

// the backbone: rank 0, offsets running the length of the walk. A segment the
// walk comes back to keeps the offset of its first visit
static void label_backbone(gfa_t *g, int32_t bb) {
    const char *name = gfa_path_name(g, bb);
    const uint32_t *segs;
    int ns = gfa_path_segs(g, bb, &segs);

    int64_t off = 0;
    for (int i = 0; i < ns; i++) {
        if (segs[i] == GFA_NIL) continue;
        gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
        if (s->rank < 0) place(s, 0, name, off);
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
    const char *name = gfa_path_name(g, k);
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
            } else if (amb || !anchored || s->rank == 0 || s->ref_name == name) {
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
            place(s, floor + 1, name, cur); // one rank deeper, carrying on from
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
    int need = GFA_LINKS | GFA_PATHS;
    if ((g->flags & need) != need) {
        ak_log(AK_LOG_ERROR, "rgfa", "labelling requires the graph to be read with GFA_LINKS | GFA_PATHS");
        return AK_EINVAL;
    }
    if (gfa_n_path(g) == 0) {
        ak_log(AK_LOG_ERROR, "rgfa", "graph has no P lines, so there is no backbone to label against");
        return AK_EINVAL;
    }

    gfa_merge_t *m = gfa_path_merge(g, NULL);
    if (!m) return AK_EINVAL;

    int32_t bb = backbone_chain(g, m, ref_name);
    if (bb < 0) {
        ak_log(AK_LOG_ERROR, "rgfa", "no path named '%s' in the graph", ref_name ? ref_name : "");
        gfa_merge_destroy(m);
        return AK_EINVAL;
    }
    if (g->has_sr) {
        ak_log(AK_LOG_WARN, "rgfa", "the file carries its own SR tags; they are replaced by what the paths say");
    }

    // take every chain before touching the path block: flattening one reads the
    // very P lines that clearing it would free
    int32_t nc = m->n;
    chain_t *c = (chain_t *)calloc((size_t)nc, sizeof(chain_t));
    int rc = c ? AK_OK : AK_ENOMEM;
    for (int32_t k = 0; rc == AK_OK && k < nc; k++) rc = chain_take(g, m, k, &c[k]);
    gfa_merge_destroy(m);

    if (rc == AK_OK) rc = install_chains(g, c, nc, &bb);
    chains_free(c, nc);
    if (rc != AK_OK) {
        ak_log(AK_LOG_ERROR, "rgfa", "cannot consolidate the paths: %s", ak_strerror(rc));
        return rc;
    }

    // the chains went in in order, so the backbone kept its index
    unlabel(g);
    label_backbone(g, bb);
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        if (k != bb) label_path(g, k);
    }

    if (st) tally(g, st);
    return AK_OK;
}
