#include "akhal/diff.h"
#include "akhal/gaf.h"
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

// GAF alignment comparison

// one alignment, reduced to what the comparison asks about
typedef struct {
    char    *qname;   // owned: read name
    char    *path;    // owned: canonical spelling of the walk
    uint64_t first;   // first node id of that spelling; 0 for a named path
} galn_t;

// one oriented step of a walk
typedef struct {
    uint64_t id;
    char     orient;   // '>' or '<'
} pstep_t;

// '>' before '<', so the two spellings of a walk order deterministically
static inline int orient_rank(char c) {
    return c == '>' ? 0 : 1;
}

// digits in an id, so the canonical spelling can be sized before it is built
static inline size_t id_width(uint64_t v) {
    size_t n = 1;
    while (v >= 10) {
        v /= 10;
        n++;
    }
    return n;
}

// A walk and the same walk read from its other end are one alignment: an
// aligner that hit the read's reverse complement writes `>3>2<1` where another
// writes `>1<2<3`, and neither spelling is more correct. So both are built and
// the smaller one is kept, which makes "same walk" a strcmp
static int steps_reversed(const pstep_t *v, int64_t n) {
    for (int64_t k = 0; k < n; k++) {
        const pstep_t *f = &v[k];
        uint64_t rid = v[n - 1 - k].id;
        int rrank = orient_rank(v[n - 1 - k].orient == '>' ? '<' : '>');

        if (f->id != rid) return f->id > rid;
        int frank = orient_rank(f->orient);
        if (frank != rrank) return frank > rrank;
    }
    return 0;   // a palindrome spells the same either way
}

// Canonicalize a GAF path field. A path naming a stable sequence rather than
// walking nodes ("chr1", as minigraph writes for an unplaced alignment) has no
// orientation to flip, so it is kept verbatim and sorts under id 0
static int path_canon(const char *path, char **out, uint64_t *first) {
    *out = NULL;
    *first = 0;

    if (path[0] != '>' && path[0] != '<') {
        *out = strdup(path);
        return *out ? AK_OK : AK_ENOMEM;
    }

    pstep_t stack[64], *v = stack;
    int64_t n = 0, cap = 64;

    const char *p = path;
    int used;
    uint64_t id;
    char orient;
    while ((used = gaf_path_next(p, &id, &orient)) > 0) {
        if (n == cap) {
            int64_t ncap = cap << 1;
            pstep_t *nv = (pstep_t *)malloc((size_t)ncap * sizeof(pstep_t));
            if (!nv) {
                if (v != stack) free(v);
                return AK_ENOMEM;
            }
            memcpy(nv, v, (size_t)n * sizeof(pstep_t));
            if (v != stack) free(v);
            v = nv;
            cap = ncap;
        }
        v[n].id = id;
        v[n].orient = orient;
        n++;
        p += used;
    }

    int rev = steps_reversed(v, n);

    size_t len = 0;
    for (int64_t k = 0; k < n; k++) len += 1 + id_width(v[k].id);

    char *s = (char *)malloc(len + 1);
    if (!s) {
        if (v != stack) free(v);
        return AK_ENOMEM;
    }

    size_t o = 0;
    for (int64_t k = 0; k < n; k++) {
        const pstep_t *st = &v[rev ? n - 1 - k : k];
        char c = st->orient;
        if (rev) c = (c == '>') ? '<' : '>';
        s[o++] = c;

        size_t w = id_width(st->id);
        uint64_t val = st->id;
        for (size_t d = w; d > 0; d--) {
            s[o + d - 1] = (char)('0' + (val % 10));
            val /= 10;
        }
        o += w;
    }
    s[o] = '\0';

    *first = n ? v[rev ? n - 1 : 0].id : 0;
    *out = s;
    if (v != stack) free(v);
    return AK_OK;
}

// read name first, then the walk: the first node id as the user reads it, and
// the whole spelling to break ties. Both files are ordered this way, so the
// comparison is one walk down the two arrays side by side
static int galn_cmp(const void *A, const void *B) {
    const galn_t *a = (const galn_t *)A, *b = (const galn_t *)B;
    int c = strcmp(a->qname, b->qname);
    if (c) return c;
    if (a->first != b->first) return a->first < b->first ? -1 : 1;
    return strcmp(a->path, b->path);
}

// the walk alone, for pairing within one read's block
static int galn_cmp_path(const galn_t *a, const galn_t *b) {
    if (a->first != b->first) return a->first < b->first ? -1 : 1;
    return strcmp(a->path, b->path);
}

static void galn_free(galn_t *v, int64_t n) {
    if (!v) return;
    for (int64_t i = 0; i < n; i++) {
        free(v[i].qname);
        free(v[i].path);
    }
    free(v);
}

// every alignment of one file, sorted. Streamed rather than slurped: only the
// read name and the canonical walk are kept, so a file of long CIGARs costs
// nothing beyond the line it is on. NULL on failure (logged), never for an
// empty file
static galn_t *galn_load(const char *fn, int64_t *n_out) {
    *n_out = 0;

    gaf_reader_t *r = gaf_open(fn);
    if (!r) return NULL;

    galn_t *v = NULL;
    int64_t n = 0, cap = 0;
    gaf_rec_t rec;
    gaf_rec_init(&rec);

    int rc;
    while ((rc = gaf_read1(r, &rec)) == 1) {
        if (n == cap) {
            int64_t ncap = cap ? cap << 1 : 4096;
            galn_t *nv = (galn_t *)realloc(v, (size_t)ncap * sizeof(galn_t));
            if (!nv) {
                rc = AK_ENOMEM;
                break;
            }
            v = nv;
            cap = ncap;
        }

        // the record's name is handed over rather than copied, and detached so
        // the next read does not free it
        v[n].qname = rec.qname;
        rec.qname = NULL;
        rc = path_canon(rec.path ? rec.path : "", &v[n].path, &v[n].first);
        if (rc != AK_OK) {
            free(v[n].qname);
            break;
        }
        n++;
    }

    gaf_rec_clear(&rec);
    gaf_close(r);

    if (rc < 0) {
        galn_free(v, n);
        ak_log(AK_LOG_ERROR, "diff", "cannot read %s: %s", fn, ak_strerror(rc));
        return NULL;
    }
    if (!v) {
        // an empty file is not a failure, but the caller tells NULL from it
        v = (galn_t *)malloc(sizeof(galn_t));
        if (!v) {
            ak_log(AK_LOG_ERROR, "diff", "out of memory");
            return NULL;
        }
    }

    qsort(v, (size_t)n, sizeof(galn_t), galn_cmp);
    *n_out = n;
    return v;
}

// end of the run of alignments sharing v[i]'s read name
static int64_t block_end(const galn_t *v, int64_t i, int64_t n) {
    int64_t j = i + 1;
    while (j < n && strcmp(v[j].qname, v[i].qname) == 0) j++;
    return j;
}

// record one verdict, taking the strings off the alignment that produced it
static void emit_aln(diff_gaf_t *d, galn_t *g, int state, int read_both) {
    diff_aln_t *a = &d->aln[d->n_aln++];
    a->qname = g->qname;
    a->path = g->path;
    a->first = g->first;
    a->state = state;
    a->read_both = read_both;
    g->qname = NULL;
    g->path = NULL;
}

// compare two GAF files; see akhal/diff.h
diff_gaf_t *diff_gaf(const char *fn_a, const char *fn_b) {
    int64_t na = 0, nb = 0;
    galn_t *va = galn_load(fn_a, &na);
    if (!va) return NULL;

    galn_t *vb = galn_load(fn_b, &nb);
    if (!vb) {
        galn_free(va, na);
        return NULL;
    }

    diff_gaf_t *d = (diff_gaf_t *)calloc(1, sizeof(diff_gaf_t));
    if (d) {
        // a pair takes one entry and a leftover one, so the two files together
        // are the bound; one allocation covers the whole walk
        d->aln = (diff_aln_t *)malloc((size_t)(na + nb > 0 ? na + nb : 1) * sizeof(diff_aln_t));
    }
    if (!d || !d->aln) {
        ak_log(AK_LOG_ERROR, "diff", "out of memory");
        diff_gaf_destroy(d);
        galn_free(va, na);
        galn_free(vb, nb);
        return NULL;
    }
    d->n_aln_a = na;
    d->n_aln_b = nb;

    // a read name at a time down both files. Names are in the same order on
    // both sides, so a name that is not next in the other file is not in it
    int64_t i = 0, j = 0;
    while (i < na || j < nb) {
        int c;
        if (i >= na)      c =  1;
        else if (j >= nb) c = -1;
        else              c = strcmp(va[i].qname, vb[j].qname);

        if (c < 0) {
            int64_t i2 = block_end(va, i, na);
            d->n_read_a++;
            d->n_read_a_only++;
            for (; i < i2; i++) {
                emit_aln(d, &va[i], DIFF_ALN_A_ONLY, 0);
                d->n_aln_a_only++;
            }
        } else if (c > 0) {
            int64_t j2 = block_end(vb, j, nb);
            d->n_read_b++;
            d->n_read_b_only++;
            for (; j < j2; j++) {
                emit_aln(d, &vb[j], DIFF_ALN_B_ONLY, 0);
                d->n_aln_b_only++;
            }
        } else {
            int64_t i2 = block_end(va, i, na), j2 = block_end(vb, j, nb);
            d->n_read_a++;
            d->n_read_b++;
            d->n_read_shared++;

            // within one name, the alignments pair off one-to-one on their
            // walk: a read aligned three ways here and twice there is two
            // pairs and a leftover, not a match
            int64_t matched = 0, left = 0;
            while (i < i2 || j < j2) {
                int k;
                if (i >= i2)      k =  1;
                else if (j >= j2) k = -1;
                else              k = galn_cmp_path(&va[i], &vb[j]);

                if (k < 0) {
                    emit_aln(d, &va[i++], DIFF_ALN_A_ONLY, 1);
                    d->n_aln_a_only++;
                    left++;
                } else if (k > 0) {
                    emit_aln(d, &vb[j++], DIFF_ALN_B_ONLY, 1);
                    d->n_aln_b_only++;
                    left++;
                } else {
                    emit_aln(d, &va[i++], DIFF_ALN_SHARED, 1);
                    j++;
                    d->n_aln_shared++;
                    matched++;
                }
            }

            if (matched && !left)  d->n_read_all_same++;
            else if (matched)      d->n_read_partial++;
            else                   d->n_read_none++;
        }
    }

    galn_free(va, na);
    galn_free(vb, nb);
    return d;
}

// free a GAF comparison; see akhal/diff.h
void diff_gaf_destroy(diff_gaf_t *d) {
    if (!d) return;
    if (d->aln) {
        for (int64_t i = 0; i < d->n_aln; i++) {
            free(d->aln[i].qname);
            free(d->aln[i].path);
        }
        free(d->aln);
    }
    free(d);
}
