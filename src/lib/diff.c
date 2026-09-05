#include "akhal/diff.h"
#include "akhal/kstr.h"
#include "akhal/util.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

// keeps the buffer NUL-terminated
static int ks_add(kstring_t *ks, const char *s, size_t len) {
    if (ks_resize(ks, ks->l + len + 1) != AK_OK) return AK_ENOMEM;
    memcpy(ks->s + ks->l, s, len);
    ks->l += len;
    ks->s[ks->l] = '\0';
    return AK_OK;
}

// segment labelling

// one segment, ordered by what it spells rather than by the id it was given
typedef struct {
    const char *seq;   // borrowed; "" for a segment carrying none
    uint32_t    idx;   // segment index in its graph
} sseg_t;

// content order, with the index as a stable tie-break between equal sequences
static int sseg_cmp(const void *A, const void *B) {
    const sseg_t *a = (const sseg_t *)A, *b = (const sseg_t *)B;
    int c = strcmp(a->seq, b->seq);
    if (c) return c;
    return a->idx < b->idx ? -1 : (a->idx > b->idx);
}

// every segment of a graph, sorted by sequence; NULL on allocation failure
static sseg_t *sorted_segs(const gfa_t *g) {
    int32_t n = gfa_n_seg(g);
    sseg_t *s = (sseg_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(sseg_t));
    if (!s) return NULL;

    for (int32_t i = 0; i < n; i++) {
        const gfa_seg_t *seg = gfa_seg_at(g, i);
        s[i].seq = seg->seq ? seg->seq : "";
        s[i].idx = (uint32_t)i;
    }
    qsort(s, (size_t)n, sizeof(sseg_t), sseg_cmp);
    return s;
}

// label both graphs' segments by content; see akhal/diff.h
diff_map_t *diff_map(const gfa_t *a, const gfa_t *b) {
    int32_t na = gfa_n_seg(a), nb = gfa_n_seg(b);

    diff_map_t *m = (diff_map_t *)calloc(1, sizeof(diff_map_t));
    sseg_t *sa = sorted_segs(a), *sb = sorted_segs(b);
    if (m) {
        m->a = (uint32_t *)malloc((size_t)(na > 0 ? na : 1) * sizeof(uint32_t));
        m->b = (uint32_t *)malloc((size_t)(nb > 0 ? nb : 1) * sizeof(uint32_t));
    }
    if (!m || !sa || !sb || !m->a || !m->b) {
        ak_log(AK_LOG_ERROR, "diff", "out of memory");
        free(sa);
        free(sb);
        diff_map_destroy(m);
        return NULL;
    }
    m->n_a = na;
    m->n_b = nb;

    // one pass over both sorted arrays, a run of equal sequences at a time. A
    // whole run takes one class, rather than pairing its members off one by
    // one: which copy of "A" pairs with which is not something the files
    // agree on, and every link touching them would inherit that arbitrary
    // choice. How many copies each side holds is still recorded, through
    // n_shared
    int32_t i = 0, j = 0;
    uint32_t cls = 0;
    while (i < na || j < nb) {
        int c;
        if (i >= na)      c =  1;
        else if (j >= nb) c = -1;
        else              c = strcmp(sa[i].seq, sb[j].seq);

        if (c < 0) {
            const char *s = sa[i].seq;
            while (i < na && strcmp(sa[i].seq, s) == 0) m->a[sa[i++].idx] = cls;
        } else if (c > 0) {
            const char *s = sb[j].seq;
            while (j < nb && strcmp(sb[j].seq, s) == 0) m->b[sb[j++].idx] = cls;
        } else {
            const char *s = sa[i].seq;
            int32_t ka = 0, kb = 0;
            while (i < na && strcmp(sa[i].seq, s) == 0) {
                m->a[sa[i++].idx] = cls;
                ka++;
            }
            while (j < nb && strcmp(sb[j].seq, s) == 0) {
                m->b[sb[j++].idx] = cls;
                kb++;
            }
            m->n_shared += ka < kb ? ka : kb;   // the surplus stays unmatched
        }
        cls++;
    }
    m->n_class = (int32_t)cls;

    free(sa);
    free(sb);
    return m;
}

// free a labelling; see akhal/diff.h
void diff_map_destroy(diff_map_t *m) {
    if (!m) return;
    free(m->a);
    free(m->b);
    free(m);
}

// segments

// the ids each graph alone carries, and how many matched
static int compare_segs(const gfa_t *ga, const gfa_t *gb, const diff_map_t *m, diff_t *d) {
    size_t nc = (size_t)(m->n_class > 0 ? m->n_class : 1);
    int32_t *cnt_a = (int32_t *)calloc(nc, sizeof(int32_t));
    int32_t *cnt_b = (int32_t *)calloc(nc, sizeof(int32_t));
    d->a.seg = (uint64_t *)malloc((size_t)(m->n_a > 0 ? m->n_a : 1) * sizeof(uint64_t));
    d->b.seg = (uint64_t *)malloc((size_t)(m->n_b > 0 ? m->n_b : 1) * sizeof(uint64_t));
    if (!cnt_a || !cnt_b || !d->a.seg || !d->b.seg) {
        free(cnt_a);
        free(cnt_b);
        return AK_ENOMEM;
    }

    // how many segments each graph puts in each class
    for (int32_t i = 0; i < m->n_a; i++) cnt_a[m->a[i]]++;
    for (int32_t j = 0; j < m->n_b; j++) cnt_b[m->b[j]]++;

    // spend the other graph's copies as we walk: the first min(ka, kb) of a
    // class are matched and the surplus is what only this graph carries. Which
    // copies end up as the surplus is decided by file order and means nothing;
    // how many there are is the answer
    for (int32_t i = 0; i < m->n_a; i++) {
        uint32_t k = m->a[i];
        if (cnt_b[k] > 0) {
            cnt_b[k]--;
        } else {
            d->a.seg[d->a.n_seg++] = gfa_seg_at(ga, i)->id;
        }
    }
    for (int32_t j = 0; j < m->n_b; j++) {
        uint32_t k = m->b[j];
        if (cnt_a[k] > 0) {
            cnt_a[k]--;
        } else {
            d->b.seg[d->b.n_seg++] = gfa_seg_at(gb, j)->id;
        }
    }

    d->n_seg_shared = m->n_shared;
    free(cnt_a);
    free(cnt_b);
    return AK_OK;
}

// links

// one link on the shared labelling, in the canonical of its two spellings
typedef struct {
    uint32_t v, w;      // labels of the two ends
    char     vo, wo;    // orientation each end carries
    uint32_t overlap;
    int32_t  idx;       // link index in its graph, so it can be reported
} klink_t;

static char flip(char orient) {
    return orient == '-' ? '+' : '-';
}

// order on everything that makes two links equal; the index is not part of it
static int klink_ord(const klink_t *a, const klink_t *b) {
    if (a->v  != b->v)  return a->v  < b->v  ? -1 : 1;
    if (a->vo != b->vo) return a->vo < b->vo ? -1 : 1;
    if (a->w  != b->w)  return a->w  < b->w  ? -1 : 1;
    if (a->wo != b->wo) return a->wo < b->wo ? -1 : 1;
    if (a->overlap != b->overlap) return a->overlap < b->overlap ? -1 : 1;
    return 0;
}

static int klink_cmp(const void *A, const void *B) {
    const klink_t *a = (const klink_t *)A, *b = (const klink_t *)B;
    int c = klink_ord(a, b);
    if (c) return c;
    return a->idx < b->idx ? -1 : (a->idx > b->idx);
}

// an edge has two spellings - `L a + b +` is `L b - a -` read from the other
// end - so both graphs are put on the smaller of the two before comparing
static klink_t link_canon(const klink_t *k) {
    klink_t r;
    r.v  = k->w;
    r.vo = flip(k->wo);
    r.w  = k->v;
    r.wo = flip(k->vo);
    r.overlap = k->overlap;
    r.idx = k->idx;
    return klink_ord(&r, k) < 0 ? r : *k;
}

// every link of a graph, relabelled and sorted; NULL on allocation failure
static klink_t *sorted_links(const gfa_t *g, const uint32_t *label) {
    int32_t n = gfa_n_link(g);
    klink_t *k = (klink_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(klink_t));
    if (!k) return NULL;

    for (int32_t i = 0; i < n; i++) {
        const gfa_link_t *e = gfa_link_at(g, i);
        klink_t t;
        t.v  = label[e->v];
        t.vo = e->from_orient;
        t.w  = label[e->w];
        t.wo = e->to_orient;
        t.overlap = e->overlap;
        t.idx = i;
        k[i] = link_canon(&t);
    }
    qsort(k, (size_t)n, sizeof(klink_t), klink_cmp);
    return k;
}

// report an unmatched link the way its own file spells it
static void record_link(diff_side_t *s, const gfa_t *g, int32_t li) {
    const gfa_link_t *e = gfa_link_at(g, li);
    diff_link_t *o = &s->link[s->n_link++];
    o->from = gfa_seg_at(g, (int32_t)e->v)->id;
    o->to   = gfa_seg_at(g, (int32_t)e->w)->id;
    o->from_orient = e->from_orient;
    o->to_orient   = e->to_orient;
    o->overlap     = e->overlap;
}

// the same merge pass as the segments, over the relabelled links
static int compare_links(const gfa_t *ga, const gfa_t *gb, const diff_map_t *m, diff_t *d) {
    int32_t na = gfa_n_link(ga), nb = gfa_n_link(gb);

    klink_t *ka = sorted_links(ga, m->a), *kb = sorted_links(gb, m->b);
    d->a.link = (diff_link_t *)malloc((size_t)(na > 0 ? na : 1) * sizeof(diff_link_t));
    d->b.link = (diff_link_t *)malloc((size_t)(nb > 0 ? nb : 1) * sizeof(diff_link_t));
    if (!ka || !kb || !d->a.link || !d->b.link) {
        free(ka);
        free(kb);
        return AK_ENOMEM;
    }

    int32_t i = 0, j = 0;
    while (i < na || j < nb) {
        int c;
        if (i >= na)      c =  1;
        else if (j >= nb) c = -1;
        else              c = klink_ord(&ka[i], &kb[j]);

        if (c == 0) {
            d->n_link_shared++;
            i++;
            j++;
        } else if (c < 0) {
            record_link(&d->a, ga, ka[i++].idx);
        } else {
            record_link(&d->b, gb, kb[j++].idx);
        }
    }

    free(ka);
    free(kb);
    return AK_OK;
}

// paths

// one merged chain, ordered by name so the two graphs' chains can be paired
typedef struct {
    const char *name;   // borrowed from the chain set
    int32_t     c;      // chain index
} chain_t;

static int chain_cmp(const void *A, const void *B) {
    const chain_t *a = (const chain_t *)A, *b = (const chain_t *)B;
    int c = strcmp(a->name, b->name);
    if (c) return c;
    return a->c < b->c ? -1 : (a->c > b->c);
}

// every chain of a merge set, sorted by name; NULL on allocation failure
static chain_t *sorted_chains(const gfa_merge_t *m) {
    chain_t *c = (chain_t *)malloc((size_t)(m->n > 0 ? m->n : 1) * sizeof(chain_t));
    if (!c) return NULL;

    for (int32_t i = 0; i < m->n; i++) {
        c[i].name = m->name[i];
        c[i].c = i;
    }
    qsort(c, (size_t)m->n, sizeof(chain_t), chain_cmp);
    return c;
}

// bases a chain spells, without materializing them
static uint64_t chain_len(const gfa_t *g, const gfa_merge_t *m, int32_t c) {
    uint64_t len = 0;
    for (int32_t f = m->off[c]; f < m->off[c + 1]; f++) {
        const uint32_t *segs;
        int ns = gfa_path_segs(g, m->frag[f], &segs);
        for (int t = 0; t < ns; t++)
            if (segs[t] != GFA_NIL) len += gfa_seg_at(g, (int32_t)segs[t])->len;
    }
    return len;
}

// the bases themselves, a '-' step contributing its reverse complement. Link
// overlaps are not trimmed off, exactly as `extract path` leaves them
static int chain_seq(const gfa_t *g, const gfa_merge_t *m, int32_t c, kstring_t *out) {
    ks_clear(out);
    for (int32_t f = m->off[c]; f < m->off[c + 1]; f++) {
        int32_t pi = m->frag[f];
        const uint32_t *segs;
        int ns = gfa_path_segs(g, pi, &segs);
        const char *ori = g->path_ori + g->path_off[pi];

        for (int t = 0; t < ns; t++) {
            if (segs[t] == GFA_NIL) continue;
            const gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[t]);
            if (!s->seq || s->len == 0) continue;

            size_t at = out->l;
            if (ks_add(out, s->seq, s->len) != AK_OK) return AK_ENOMEM;
            if (ori[t] == '-') ak_revcomp(out->s + at, s->len);
        }
    }
    return AK_OK;
}

// chain both graphs' P-line fragments, pair the chains by name, and compare
// what each pair spells
static int compare_paths(const gfa_t *ga, const gfa_t *gb, diff_t *d) {
    // a graph with no P lines has no chains, which is nothing to fail over -
    // its segments and links still compare
    gfa_merge_t *ma = NULL, *mb = NULL;
    chain_t *ca = NULL, *cb = NULL;
    int rc = AK_OK;

    if (gfa_n_path(ga) > 0 && !(ma = gfa_path_merge(ga, NULL))) rc = AK_EINVAL;
    if (rc == AK_OK && gfa_n_path(gb) > 0 && !(mb = gfa_path_merge(gb, NULL))) rc = AK_EINVAL;

    int32_t na = ma ? ma->n : 0, nb = mb ? mb->n : 0;
    if (rc == AK_OK && na > 0 && !(ca = sorted_chains(ma))) rc = AK_ENOMEM;
    if (rc == AK_OK && nb > 0 && !(cb = sorted_chains(mb))) rc = AK_ENOMEM;

    if (rc == AK_OK) {
        d->path = (diff_path_t *)calloc((size_t)(na + nb > 0 ? na + nb : 1), sizeof(diff_path_t));
        if (!d->path) rc = AK_ENOMEM;

        kstring_t sa = KS_INIT, sb = KS_INIT;
        int32_t i = 0, j = 0;
        while (rc == AK_OK && (i < na || j < nb)) {
            int c;
            if (i >= na)      c =  1;
            else if (j >= nb) c = -1;
            else              c = strcmp(ca[i].name, cb[j].name);

            diff_path_t *p = &d->path[d->n_path];
            p->name = strdup(c <= 0 ? ca[i].name : cb[j].name);
            if (!p->name) {
                rc = AK_ENOMEM;
                break;
            }
            d->n_path++;

            if (c < 0) {
                p->state = DIFF_A_ONLY;
                p->len_a = chain_len(ga, ma, ca[i++].c);
                d->n_path_a_only++;
            } else if (c > 0) {
                p->state = DIFF_B_ONLY;
                p->len_b = chain_len(gb, mb, cb[j++].c);
                d->n_path_b_only++;
            } else {
                p->len_a = chain_len(ga, ma, ca[i].c);
                p->len_b = chain_len(gb, mb, cb[j].c);

                // different lengths settle it; equal ones need the bases
                int same = 0;
                if (p->len_a == p->len_b) {
                    rc = chain_seq(ga, ma, ca[i].c, &sa);
                    if (rc == AK_OK) rc = chain_seq(gb, mb, cb[j].c, &sb);
                    if (rc != AK_OK) break;
                    same = sa.l == 0 || memcmp(sa.s, sb.s, sa.l) == 0;
                }
                p->state = same ? DIFF_SAME : DIFF_DIFFER;
                if (same) {
                    d->n_path_same++;
                } else {
                    d->n_path_differ++;
                }
                i++;
                j++;
            }
        }
        ks_free(&sa);
        ks_free(&sb);
    }

    free(ca);
    free(cb);
    gfa_merge_destroy(ma);
    gfa_merge_destroy(mb);
    return rc;
}

// compare two graphs; see akhal/diff.h
diff_t *diff_graphs(const gfa_t *a, const gfa_t *b) {
    int need = GFA_LINKS | GFA_PATHS;
    if ((a->flags & need) != need || (b->flags & need) != need) {
        ak_log(AK_LOG_ERROR, "diff", "comparison requires both graphs to be read with GFA_LINKS | GFA_PATHS");
        return NULL;
    }

    diff_t *d = (diff_t *)calloc(1, sizeof(diff_t));
    diff_map_t *m = diff_map(a, b);
    if (!d || !m) {
        diff_map_destroy(m);
        diff_destroy(d);
        return NULL;
    }

    int rc = compare_segs(a, b, m, d);
    if (rc == AK_OK) rc = compare_links(a, b, m, d);
    diff_map_destroy(m);
    if (rc == AK_OK) rc = compare_paths(a, b, d);

    if (rc != AK_OK) {
        ak_log(AK_LOG_ERROR, "diff", "cannot compare the graphs: %s", ak_strerror(rc));
        diff_destroy(d);
        return NULL;
    }
    return d;
}

// free a comparison; see akhal/diff.h
void diff_destroy(diff_t *d) {
    if (!d) return;
    free(d->a.seg);
    free(d->a.link);
    free(d->b.seg);
    free(d->b.link);
    if (d->path) {
        for (int32_t i = 0; i < d->n_path; i++) free(d->path[i].name);
        free(d->path);
    }
    free(d);
}
