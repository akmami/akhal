// Summary statistics for a GFA file, without building a graph.
//
// Everything is a running total or a Welford accumulator except three things,
// which are properties of the file as a whole rather than of any one line: a
// segment's degree is the number of L lines that name it, "undefined" means
// named by an L or P line but by no S line, and (in a file without SR tags)
// "rank 0" means named by some P line. Those are answered by id arrays rather
// than per-id tables: the ids are appended to flat arrays as they stream
// past, radix-sorted in place once the file is read, and read off as run
// lengths and merges. An array costs 4 bytes per entry (8 if an id ever exceeds
// 32 bits) whatever the numbering looks like, so a graph whose ids are
// scattered costs exactly what a densely numbered one does, and there is no
// fallback path.
//
// Kept apart from gfa.c, which owns the in-memory model and its reader; this
// file only shares the public gfa_stat_t with it.

#include "akhal/gfa.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/util.h"
#include "akhal/error.h"

#include "ksort.h"

#include <stdlib.h>
#include <string.h>

// id arrays

// A flat array of segment ids, one per line that named one. It starts 32-bit
// wide and widens to 64 the first time an id needs it, so the common case
// pays 4 bytes per entry and the general case still works
typedef struct {
    void   *id;
    size_t  n_id, m_id;
    int     wide;        // 0: uint32_t entries, 1: uint64_t entries
} ids_t;

#define IDS_INIT { NULL, 0, 0, 0 }

static inline uint64_t ids_at(const ids_t *ids, size_t i) {
    return ids->wide ? ((const uint64_t *)ids->id)[i] : ((const uint32_t *)ids->id)[i];
}

// widen in place, spreading the entries from the back so nothing is overwritten
static int ids_widen(ids_t *ids) {
    if (ids->m_id == 0) {           // nothing allocated yet: just start out wide
        ids->wide = 1;
        return AK_OK;
    }
    void *id = realloc(ids->id, ids->m_id * sizeof(uint64_t));
    if (!id) return AK_ENOMEM;
    const uint32_t *s = (const uint32_t *)id;
    uint64_t *d = (uint64_t *)id;
    for (size_t i = ids->n_id; i-- > 0;) d[i] = s[i];
    ids->id = id;
    ids->wide = 1;
    return AK_OK;
}

// grow by half rather than doubling: an id array is the largest thing this pass
// holds, and the slack past the last growth is what the peak is charged for
static int ids_push(ids_t *ids, uint64_t id) {
    if (!ids->wide && id > UINT32_MAX) {
        if (ids_widen(ids) != AK_OK) return AK_ENOMEM;
    }
    if (ids->n_id == ids->m_id) {
        size_t m_id = ids->m_id ? ids->m_id + ids->m_id / 2 : (1u << 16);
        size_t esz = ids->wide ? sizeof(uint64_t) : sizeof(uint32_t);
        void *id = realloc(ids->id, m_id * esz);
        if (!id) return AK_ENOMEM;
        ids->id = id;
        ids->m_id = m_id;
    }
    if (ids->wide) ((uint64_t *)ids->id)[ids->n_id++] = id;
    else           ((uint32_t *)ids->id)[ids->n_id++] = (uint32_t)id;
    return AK_OK;
}

static void ids_free(ids_t *ids) {
    free(ids->id);
    ids->id = NULL;
    ids->n_id = ids->m_id = 0;
    ids->wide = 0;
}

KRADIX_SORT_INIT(u32, uint32_t, , 4)
KRADIX_SORT_INIT(u64, uint64_t, , 8)

// in-place radix sort. The klib sort starts from the key's top byte; on ids
// that never reach that high the early passes just move everything into one
// bucket, so start from the highest byte the largest id actually uses
static void ids_sort(ids_t *ids, uint64_t max_id) {
    if (ids->n_id <= RS_MIN_SIZE) {
        if (ids->wide) rs_insertsort_u64((uint64_t *)ids->id, (uint64_t *)ids->id + ids->n_id);
        else         rs_insertsort_u32((uint32_t *)ids->id, (uint32_t *)ids->id + ids->n_id);
        return;
    }
    int top = 0;                                   // highest byte in use
    while (top < 7 && (max_id >> (8 * (top + 1))) != 0) top++;
    if (ids->wide) rs_sort_u64((uint64_t *)ids->id, (uint64_t *)ids->id + ids->n_id, RS_MAX_BITS, top * RS_MAX_BITS);
    else         rs_sort_u32((uint32_t *)ids->id, (uint32_t *)ids->id + ids->n_id, RS_MAX_BITS, top * RS_MAX_BITS);
}

// A cursor over a sorted id array that yields each distinct id once, with how
// many times it occurred. ids_run_init() leaves it on the first run, so a walk
// is `for (; r.live; ids_run_next(&r))`; `live` is 0 once the array is spent
typedef struct {
    const ids_t *ids;
    size_t   i;
    uint64_t id;
    size_t   run;
    int      live;
} ids_run_t;

static int ids_run_next(ids_run_t *r) {
    if (r->i >= r->ids->n_id) return r->live = 0;
    r->id = ids_at(r->ids, r->i);
    size_t j = r->i + 1;
    while (j < r->ids->n_id && ids_at(r->ids, j) == r->id) j++;
    r->run = j - r->i;
    r->i = j;
    return r->live = 1;
}

static void ids_run_init(ids_run_t *r, const ids_t *ids) {
    r->ids = ids;
    r->i = 0;
    r->id = r->run = 0;
    r->live = ids ? ids_run_next(r) : 0;
}

// advance until the cursor is at or past `id`; returns 1 when it is on it
static int ids_run_seek(ids_run_t *r, uint64_t id) {
    while (r->live && r->id < id) ids_run_next(r);
    return r->live && r->id == id;
}

// degree distribution over one end of the links: the run length of each id is
// its degree. Ids no S line defines are skipped - only a real segment carries
// a degree - and so the extremes and moments cover segments with a non-zero
// degree, as they always have
static void count_degrees(const ids_t *ends, const ids_t *segs, ak_dist_t *d) {
    ids_run_t r, s;
    ids_run_init(&r, ends);
    ids_run_init(&s, segs);
    for (; r.live; ids_run_next(&r)) {
        if (ids_run_seek(&s, r.id)) ak_dist_add(d, (double)r.run);
    }
}

// distinct ids across up to three sorted arrays that no S line defines
static int64_t count_undefined(const ids_t *named[3], const ids_t *segs) {
    ids_run_t r[3], s;
    for (int k = 0; k < 3; k++) ids_run_init(&r[k], named[k]);
    ids_run_init(&s, segs);

    int64_t n = 0;
    for (;;) {
        // the smallest id any live cursor is on
        int pick = -1;
        for (int k = 0; k < 3; k++)
            if (r[k].live && (pick < 0 || r[k].id < r[pick].id)) pick = k;
        if (pick < 0) break;
        uint64_t id = r[pick].id;
        if (!ids_run_seek(&s, id)) n++;
        for (int k = 0; k < 3; k++)
            if (r[k].live && r[k].id == id) ids_run_next(&r[k]);
    }
    return n;
}

// distinct ids in a sorted array that an S line defines
static int64_t count_defined(const ids_t *ids, const ids_t *segs) {
    ids_run_t r, s;
    ids_run_init(&r, ids);
    ids_run_init(&s, segs);
    int64_t n = 0;
    for (; r.live; ids_run_next(&r))
        if (ids_run_seek(&s, r.id)) n++;
    return n;
}

// the last field of an S line that looks like SR:i:<n>, or -1 when there is none
static int32_t sr_tag(const char *line) {
    for (const char *p = line; (p = strstr(p, "SR:i:")) != NULL; p += 5) {
        if (p != line && p[-1] != '\t') continue;   // a tag starts a field
        return (int32_t)strtol(p + 5, NULL, 10);
    }
    return -1;
}

// summarize a file in one pass; see akhal/gfa.h
int gfa_read_stats(const char *fn, gfa_stat_t *st, int flags) {
    ak_file *f = ak_open(fn);
    if (!f) return AK_EOPEN;

    int want_deg = (flags & GFA_STAT_DEGREES) != 0;
    int want_rank = (flags & GFA_STAT_RANKS) != 0;

    ak_dist_t sl = {0}, ov = {0};
    ids_t segs = IDS_INIT, src = IDS_INIT, dst = IDS_INIT, steps = IDS_INIT;
    int64_t n_seg = 0, n_link = 0, n_path = 0, sr0 = 0;
    int has_sr = 0, rc = AK_OK;
    uint64_t max_id = 0;

    kstring_t ks = KS_INIT;
    long len;
    while ((len = ak_getline(f, &ks)) >= 0 && rc == AK_OK) {
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
            if (want_deg || want_rank) {
                if (id > max_id) max_id = id;
                rc = ids_push(&segs, id);
            }

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

            if (want_deg) {
                if (a > max_id) max_id = a;
                if (b > max_id) max_id = b;
                rc = ids_push(&src, a);
                if (rc == AK_OK) rc = ids_push(&dst, b);
            }

        } else if (p[0] == 'P' && p[1] == '\t') {
            n_path++;
            // the steps are only walked when rank 0 has to be derived from
            // them; a file that ranks itself needs nothing but the line count
            if (!want_rank || has_sr) continue;
            p += 2;
            while (*p && *p != '\t') p++;                // past the name
            if (!*p) continue;
            p++;
            while (*p && *p != '\t' && rc == AK_OK) {
                char *q;
                uint64_t id = strtoull(p, &q, 10);
                if (q == p) break;
                if (id > max_id) max_id = id;
                rc = ids_push(&steps, id);
                p = q;
                while (*p && *p != ',' && *p != '\t') p++;
                if (*p == ',') p++;
            }
        }
    }
    ks_free(&ks);
    ak_close(f);

    if (rc != AK_OK) {
        ids_free(&segs);
        ids_free(&src);
        ids_free(&dst);
        ids_free(&steps);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory summarizing %s", fn);
        return rc;
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
    st->n_rank0  = has_sr ? sr0 : -1;               // -1: not derived
    st->n_undefined = -1;                           // -1: not checked
    st->in_mean = st->out_mean = st->in_sd = st->out_sd = 0.0;
    st->min_in = st->max_in = st->min_out = st->max_out = -1;

    if (!want_deg && !want_rank) return AK_OK;

    // A file whose SR tags were read has already answered the rank question,
    // and its steps were never collected; that array is empty in that case
    ids_sort(&segs, max_id);

    if (want_deg) {
        ak_dist_t din = {0}, dout = {0};
        ids_sort(&src, max_id);
        count_degrees(&src, &segs, &dout);
        ids_sort(&dst, max_id);
        count_degrees(&dst, &segs, &din);
        if (din.n) {
            st->in_mean = din.mean;
            st->in_sd   = ak_dist_sd(&din);
            st->min_in  = (int32_t)din.min;
            st->max_in  = (int32_t)din.max;
        }
        if (dout.n) {
            st->out_mean = dout.mean;
            st->out_sd   = ak_dist_sd(&dout);
            st->min_out  = (int32_t)dout.min;
            st->max_out  = (int32_t)dout.max;
        }
    }
    if (want_rank && !has_sr) {
        ids_sort(&steps, max_id);
        // every segment some P line names sits on the backbone - the same
        // count gfa_rank_paths() makes, which can only mark defined segments
        st->n_rank0 = count_defined(&steps, &segs);
    }

    // undefined: named by an L line (when the link ends were kept) or a P
    // line (when the steps were), but by no S line
    const ids_t *named[3] = { want_deg ? &src : NULL, want_deg ? &dst : NULL, steps.n_id ? &steps : NULL };
    st->n_undefined = count_undefined(named, &segs);

    ids_free(&segs);
    ids_free(&src);
    ids_free(&dst);
    ids_free(&steps);
    return AK_OK;
}
