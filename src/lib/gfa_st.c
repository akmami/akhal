// Summary statistics for a GFA file: a single streaming pass that never builds
// a graph (gfa_read_stats), plus the graph-backed fallback it uses when the
// segment ids are too scattered to index by value. Kept apart from gfa.c,
// which owns the in-memory model and its reader; this file only shares the
// public gfa_stat_t with it.

#include "akhal/gfa.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/util.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

// What the summary keeps per segment id: its two degrees and three bits
// saying where the id was seen. Indexed by id - lo, so a chunk numbered from
// forty million costs no more than one numbered from one
#define SEEN_SEG  0x1
#define SEEN_PATH 0x2
#define SEEN_LINK 0x4

typedef struct {
    uint32_t *in, *out;
    uint8_t  *seen;
    uint64_t  lo, n;     // ids lo .. lo + n - 1 are addressable
} sidx_t;

static void sidx_free(sidx_t *x) {
    free(x->in);
    free(x->out);
    free(x->seen);
    x->in = NULL;
    x->out = NULL;
    x->seen = NULL;
    x->lo = x->n = 0;
}

// Make room for `id`, growing in whichever direction it lies. Returns 0 when
// the range asked for is too sparse to be worth an array - `budget` is the
// caller's idea of how many entries the file has earned so far
static int sidx_fit(sidx_t *x, uint64_t id, int64_t budget) {
    if (x->n && id >= x->lo && id < x->lo + x->n) return 1;

    uint64_t lo = x->n ? (id < x->lo ? id : x->lo) : id;
    uint64_t hi = x->n ? (id >= x->lo + x->n ? id : x->lo + x->n - 1) : id;
    uint64_t need = hi - lo + 1;

    uint64_t cap = (uint64_t)budget * 4 + (1u << 20);
    if (need > cap) return 0;

    // grow geometrically so a file numbered 1..N reallocates log N times
    uint64_t n = x->n ? x->n : (1u << 16);
    while (n < need) {
        uint64_t next = n << 1;
        if (next < n) return 0;
        n = next;
    }
    if (n > cap) n = need;
    // extend downwards only as far as asked; the common case grows upwards
    if (lo < x->lo || !x->n) {
        // keep the new low end, and let the array run up from there
        if (x->n && lo + n < x->lo + x->n) n = x->lo + x->n - lo;
    } else {
        lo = x->lo;
    }
    if (n > cap) return 0;

    uint32_t *ni = (uint32_t *)calloc((size_t)n, sizeof(*ni));
    uint32_t *no = (uint32_t *)calloc((size_t)n, sizeof(*no));
    uint8_t  *ns = (uint8_t *)calloc((size_t)n, 1);
    if (!ni || !no || !ns) {
        free(ni);
        free(no);
        free(ns);
        return -1;
    }
    if (x->n) {
        size_t at = (size_t)(x->lo - lo);
        memcpy(ni + at, x->in,   (size_t)x->n * sizeof(*ni));
        memcpy(no + at, x->out,  (size_t)x->n * sizeof(*no));
        memcpy(ns + at, x->seen, (size_t)x->n);
        free(x->in);
        free(x->out);
        free(x->seen);
    }
    x->in = ni;
    x->out = no;
    x->seen = ns;
    x->lo = lo;
    x->n = n;
    return 1;
}

// the last field of an S line that looks like SR:i:<n>, or -1 when there is none
static int32_t sr_tag(const char *line) {
    for (const char *p = line; (p = strstr(p, "SR:i:")) != NULL; p += 5) {
        if (p != line && p[-1] != '\t') continue;   // a tag starts a field
        return (int32_t)strtol(p + 5, NULL, 10);
    }
    return -1;
}

// fall back to the graph when the ids defeat the index; the answer is the same
static int stats_via_graph(const char *fn, gfa_stat_t *st) {
    gfa_t *g = gfa_read(fn, GFA_LINKS | GFA_PATHS | GFA_DEGREES);
    if (!g) return AK_EOPEN;

    ak_dist_t sl = {0}, ov = {0};
    for (int32_t i = 0; i < g->n_seg; i++) ak_dist_add(&sl, (double)g->seg[i].len);
    for (int32_t i = 0; i < g->n_link; i++) ak_dist_add(&ov, (double)g->link[i].overlap);

    st->n_seg  = g->n_seg;
    st->n_link = g->n_link;
    st->n_path = g->n_path;
    st->has_sr = g->has_sr;
    st->seg_mean = sl.mean;
    st->seg_sd   = ak_dist_sd(&sl);
    st->seg_min  = g->n_seg ? (uint64_t)sl.min : 0;
    st->seg_max  = g->n_seg ? (uint64_t)sl.max : 0;
    st->ov_mean  = ov.mean;
    st->ov_sd    = ak_dist_sd(&ov);

    st->min_in = st->max_in = st->min_out = st->max_out = -1;
    st->n_rank0 = 0;
    for (int32_t i = 0; i < g->n_seg; i++) {
        const gfa_seg_t *s = &g->seg[i];
        if (s->rank == 0) st->n_rank0++;
        if (g->in_degree[i]) {
            if (st->min_in < 0 || g->in_degree[i] < st->min_in) st->min_in = g->in_degree[i];
            if (g->in_degree[i] > st->max_in) st->max_in = g->in_degree[i];
        }
        if (g->out_degree[i]) {
            if (st->min_out < 0 || g->out_degree[i] < st->min_out) st->min_out = g->out_degree[i];
            if (g->out_degree[i] > st->max_out) st->max_out = g->out_degree[i];
        }
    }
    st->n_undefined = 0;   // gfa_read() drops those lines rather than counting them
    gfa_destroy(g);
    return AK_OK;
}

// summarize a file in one pass; see akhal/gfa.h
int gfa_read_stats(const char *fn, gfa_stat_t *st) {
    ak_file *f = ak_open(fn);
    if (!f) return AK_EOPEN;

    ak_dist_t sl = {0}, ov = {0};
    sidx_t x = {0};
    int64_t n_seg = 0, n_link = 0, n_path = 0, sr0 = 0;
    int has_sr = 0, sparse = 0, oom = 0;

    kstring_t ks = KS_INIT;
    long len;
    while ((len = ak_getline(f, &ks)) >= 0 && !sparse && !oom) {
        if (len == 0) continue;
        char *p = ks.s;

        if (p[0] == 'S' && p[1] == '\t') {
            p += 2;
            uint64_t id = strtoull(p, &p, 10);
            if (*p == '\t') p++;
            const char *seq = p;
            while (*p && *p != '\t') p++;
            n_seg++;
            ak_dist_add(&sl, (double)(p - seq));

            int32_t r = sr_tag(ks.s);
            if (r >= 0) {
                has_sr = 1;
                if (r == 0) sr0++;
            }
            int ok = sidx_fit(&x, id, n_seg);
            if (ok < 0) { oom = 1; break; }
            if (!ok)    { sparse = 1; break; }
            x.seen[id - x.lo] |= SEEN_SEG;

        } else if (p[0] == 'L' && p[1] == '\t') {
            p += 2;
            char *q;
            uint64_t a = strtoull(p, &q, 10);
            if (q == p || *q != '\t') continue;          // malformed; gfa_read warns and skips
            p = q + 1;
            if (!*p || p[1] != '\t') continue;
            p += 2;
            uint64_t b = strtoull(p, &q, 10);
            // an L line needs both ids and both orientations; the overlap may
            // be absent or a "*", which counts as 0 - the same four fields
            // gfa_read()'s sscanf insists on
            if (q == p || *q != '\t' || !q[1]) continue;
            p = q + 1;                                   // the to-orientation
            uint64_t o = 0;
            if (p[1] == '\t') {
                p += 2;
                o = strtoull(p, &q, 10);
                if (q == p) o = 0;
            }
            n_link++;
            ak_dist_add(&ov, (double)o);

            int ok = sidx_fit(&x, a, n_seg);
            if (ok > 0) ok = sidx_fit(&x, b, n_seg);
            if (ok < 0) { oom = 1; break; }
            if (!ok)    { sparse = 1; break; }
            if (x.out[a - x.lo] != UINT32_MAX) x.out[a - x.lo]++;
            if (x.in[b - x.lo]  != UINT32_MAX) x.in[b - x.lo]++;
            x.seen[a - x.lo] |= SEEN_LINK;
            x.seen[b - x.lo] |= SEEN_LINK;

        } else if (p[0] == 'P' && p[1] == '\t') {
            p += 2;
            while (*p && *p != '\t') p++;                // past the name
            if (!*p) continue;
            p++;
            n_path++;
            while (*p && *p != '\t') {
                char *q;
                uint64_t id = strtoull(p, &q, 10);
                if (q == p) break;
                int ok = sidx_fit(&x, id, n_seg);
                if (ok < 0) { oom = 1; break; }
                if (!ok)    { sparse = 1; break; }
                x.seen[id - x.lo] |= SEEN_PATH;
                p = q;
                while (*p && *p != ',' && *p != '\t') p++;
                if (*p == ',') p++;
            }
        }
    }
    ks_free(&ks);
    ak_close(f);

    if (oom) {
        sidx_free(&x);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory summarizing %s", fn);
        return AK_ENOMEM;
    }
    if (sparse) {
        sidx_free(&x);
        ak_log(AK_LOG_INFO, "gfa", "segment ids are too scattered to index; summarizing through the graph instead");
        return stats_via_graph(fn, st);
    }

    memset(st, 0, sizeof(*st));
    st->n_seg  = n_seg;
    st->n_link = n_link;
    st->n_path = n_path;
    st->has_sr = has_sr;
    st->seg_mean = sl.mean;
    st->seg_sd   = ak_dist_sd(&sl);
    st->seg_min  = n_seg ? (uint64_t)sl.min : 0;
    st->seg_max  = n_seg ? (uint64_t)sl.max : 0;
    st->ov_mean  = ov.mean;
    st->ov_sd    = ak_dist_sd(&ov);
    st->min_in = st->max_in = st->min_out = st->max_out = -1;

    for (uint64_t i = 0; i < x.n; i++) {
        uint8_t s = x.seen[i];
        if (!s) continue;
        if (!(s & SEEN_SEG)) {
            if (s & (SEEN_LINK | SEEN_PATH)) st->n_undefined++;
            continue;   // only a real segment carries a degree or a rank
        }
        uint32_t di = x.in[i], dou = x.out[i];
        if (di) {
            if (st->min_in < 0 || (int32_t)di < st->min_in) st->min_in = (int32_t)di;
            if ((int32_t)di > st->max_in) st->max_in = (int32_t)di;
        }
        if (dou) {
            if (st->min_out < 0 || (int32_t)dou < st->min_out) st->min_out = (int32_t)dou;
            if ((int32_t)dou > st->max_out) st->max_out = (int32_t)dou;
        }
        if (!has_sr && (s & SEEN_PATH)) st->n_rank0++;
    }
    if (has_sr) st->n_rank0 = sr0;

    sidx_free(&x);
    return AK_OK;
}

