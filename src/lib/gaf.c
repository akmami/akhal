#include "akhal/gaf.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

struct gaf_reader {
    ak_file  *f;
    kstring_t line;
};

/** Zero-initialize a record. */
void gaf_rec_init(gaf_rec_t *r) {
    memset(r, 0, sizeof(*r));
}

/** Free owned strings and reset a record for reuse. */
void gaf_rec_clear(gaf_rec_t *r) {
    free(r->qname);
    free(r->path);
    free(r->cigar);
    memset(r, 0, sizeof(*r));
}

// line parsing

#if defined(__GNUC__) || defined(__clang__)
#  define GAF_LIKELY(x)   __builtin_expect(!!(x), 1)
#  define GAF_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#  define GAF_LIKELY(x)   (x)
#  define GAF_UNLIKELY(x) (x)
#endif

/** 
 * @brief End of the field starting at p: the next tab, or the NUL terminator. 
 */
static inline const char *gaf_field_end(const char *p) {
    const char *e = strchr(p, '\t');
    return e ? e : p + strlen(p);
}

/** 
 * @brief Copy n bytes into a fresh NUL-terminated string; NULL on failure. 
 */
static inline char *gaf_dup(const char *s, size_t n) {
    char *d = (char *)malloc(n + 1);
    if (GAF_UNLIKELY(!d)) return NULL;
    memcpy(d, s, n);
    d[n] = '\0';
    return d;
}

/** 
 * @brief Scan a decimal integer, leaving *pp on the first non-digit. 
 */
static inline int64_t gaf_scan_i64(const char **pp) {
    const char *p = *pp;
    unsigned c = (unsigned char)*p;
    int neg = (c == '-');
    p += (neg | (c == '+'));

    uint64_t v = 0;
    unsigned d;
    while ((d = (unsigned)((unsigned char)*p - '0')) <= 9u) {
        v = v * 10u + d;
        p++;
    }
    *pp = p;
    return neg ? -(int64_t)v : (int64_t)v;
}

// Negative powers of ten, indexed by fractional digit count.
static const double gaf_p10n[19] = {
    1e-0,  1e-1,  1e-2,  1e-3,  1e-4,  1e-5,  1e-6,  1e-7,  1e-8, 1e-9,
    1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18
};

/** 
 * @brief Scan a decimal float. Defers to strtod for exponents or >18 digits. 
 */
static inline double gaf_scan_f64(const char *p) {
    const char *orig = p;
    unsigned c = (unsigned char)*p;
    int neg = (c == '-');
    p += (neg | (c == '+'));

    uint64_t m = 0;
    unsigned d;
    int nd = 0, frac = 0;

    while ((d = (unsigned)((unsigned char)*p - '0')) <= 9u) {
        m = m * 10u + d; p++; nd++;
    }
    if (*p == '.') {
        p++;
        while ((d = (unsigned)((unsigned char)*p - '0')) <= 9u) {
            m = m * 10u + d; p++; nd++; frac++;
        }
    }
    if (GAF_UNLIKELY(nd > 18 || *p == 'e' || *p == 'E')) return strtod(orig, NULL);

    double r = (double)m * gaf_p10n[frac];
    return neg ? -r : r;
}

// Two-letter tag key, so tag dispatch is one switch instead of a strncmp chain.
#define GAF_TAG2(a, b) ((unsigned)(((unsigned char)(a) << 8) | (unsigned char)(b)))

// Step over the tab ending the current field; bail out if the line ran short.
#define GAF_EAT_TAB()                                   \
    do {                                                \
        if (GAF_UNLIKELY(*p != '\t')) {                 \
            while (*p && *p != '\t') p++;               \
            if (GAF_UNLIKELY(!*p)) goto short_line;     \
        }                                               \
        p++;                                            \
    } while (0)

/**
 * Parse one GAF line into a cleared record. The line is read but not modified.
 * @param line The line, NUL-terminated (a trailing newline is tolerated).
 * @param rec Destination record (cleared first).
 * @return AK_OK, AK_EFORMAT (too few fields), or AK_ENOMEM.
 */
static int gaf_parse_line(char *line, gaf_rec_t *rec) {
    gaf_rec_clear(rec);

    const char *p = line;
    const char *s, *e;
    size_t len;

    // qname
    e = gaf_field_end(p);
    if (GAF_UNLIKELY(*e != '\t')) goto short_line;
    rec->qname = gaf_dup(p, (size_t)(e - p));
    if (GAF_UNLIKELY(!rec->qname)) return AK_ENOMEM;
    p = e + 1;

    // qlen, qstart, qend
    rec->qlen = gaf_scan_i64(&p); GAF_EAT_TAB();
    rec->qstart = gaf_scan_i64(&p); GAF_EAT_TAB();
    rec->qend = gaf_scan_i64(&p); GAF_EAT_TAB();

    // strand
    if (GAF_LIKELY(*p != '\t' && *p != '\0')) { rec->strand = *p; p++; }
    GAF_EAT_TAB();

    // path
    e = gaf_field_end(p);
    if (GAF_UNLIKELY(*e != '\t')) goto short_line;
    rec->path = gaf_dup(p, (size_t)(e - p));
    if (GAF_UNLIKELY(!rec->path)) return AK_ENOMEM;
    p = e + 1;

    // plen, pstart, pend, matches, block_len
    rec->plen = gaf_scan_i64(&p); GAF_EAT_TAB();
    rec->pstart = gaf_scan_i64(&p); GAF_EAT_TAB();
    rec->pend = gaf_scan_i64(&p); GAF_EAT_TAB();
    rec->matches = gaf_scan_i64(&p); GAF_EAT_TAB();
    rec->block_len = gaf_scan_i64(&p); GAF_EAT_TAB();

    // mapq, the last mandatory field
    rec->mapq = (int)gaf_scan_i64(&p);
    while (*p && *p != '\t') p++;
    if (!*p) return AK_OK;      // exactly 12 fields, no tags
    p++;

    // optional tags
    for (;;) {
        s = p;
        p = gaf_field_end(p);
        len = (size_t)(p - s);
        while (len && (s[len - 1] == '\n' || s[len - 1] == '\r')) len--;

        if (GAF_LIKELY(len >= 5 && s[2] == ':' && s[4] == ':')) {
            switch (GAF_TAG2(s[0], s[1])) {
                case GAF_TAG2('N', 'M'):
                    if (s[3] == 'i') {
                        const char *q = s + 5;
                        rec->nm = (int)gaf_scan_i64(&q);
                        rec->has_nm = 1;
                    }
                    break;
                case GAF_TAG2('A', 'S'):
                    // minimap2 and GraphAligner both emit AS:i as well as AS:f
                    if (s[3] == 'f' || s[3] == 'i') {
                        rec->as = gaf_scan_f64(s + 5);
                        rec->has_as = 1;
                    }
                    break;
                case GAF_TAG2('d', 'v'):
                    if (s[3] == 'f') { rec->dv = gaf_scan_f64(s + 5); rec->has_dv = 1; }
                    break;
                case GAF_TAG2('i', 'd'):
                    if (s[3] == 'f') { rec->id = gaf_scan_f64(s + 5); rec->has_id = 1; }
                    break;
                case GAF_TAG2('c', 'g'):
                    if (s[3] == 'Z') {
                        rec->cigar = gaf_dup(s + 5, len - 5);
                        if (GAF_UNLIKELY(!rec->cigar)) return AK_ENOMEM;
                    }
                    break;
                default:
                    break;
            }
        }

        if (!*p) break;
        p++;
    }
    return AK_OK;

short_line: {
        int n = (*line != '\0');
        for (const char *q = line; *q; q++) n += (*q == '\t');
        ak_log(AK_LOG_WARN, "gaf", "line has only %d fields (need >= 12)", n);
        return AK_EFORMAT;
    }
}

#undef GAF_EAT_TAB

// streaming

/** Open a GAF file for streaming; see akhal/gaf.h. */
gaf_reader_t *gaf_open(const char *fn) {
    ak_file *f = ak_open(fn);
    if (!f) return NULL;

    gaf_reader_t *r = (gaf_reader_t *)calloc(1, sizeof(gaf_reader_t));
    if (!r) {
        ak_close(f);
        ak_log(AK_LOG_ERROR, "gaf", "out of memory");
        return NULL;
    }
    r->f = f;
    return r;
}

/**
 * Read the next well-formed record, skipping malformed lines.
 * @param r Open reader.
 * @param rec Destination record (reused each call).
 * @return 1 on success, 0 at EOF, AK_EINVAL/AK_ENOMEM on error.
 */
int gaf_read1(gaf_reader_t *r, gaf_rec_t *rec) {
    if (!r || !rec) return AK_EINVAL;

    long len;
    while ((len = ak_getline(r->f, &r->line)) >= 0) {
        if (len == 0) continue;
        int rc = gaf_parse_line(r->line.s, rec);
        if (rc == AK_ENOMEM) return AK_ENOMEM;
        if (rc == AK_EFORMAT) continue;   // skip malformed lines
        return 1;
    }
    return 0;   // EOF
}

/** Close a streaming reader. Safe with NULL. */
void gaf_close(gaf_reader_t *r) {
    if (!r) return;
    ak_close(r->f);
    ks_free(&r->line);
    free(r);
}

// batch

/** Load an entire GAF file into an array; see akhal/gaf.h. */
gaf_rec_t *gaf_slurp(const char *fn, int64_t *n) {
    if (n) *n = 0;

    gaf_reader_t *r = gaf_open(fn);
    if (!r) return NULL;

    gaf_rec_t *arr = NULL;
    int64_t cnt = 0, cap = 0;
    gaf_rec_t rec;
    gaf_rec_init(&rec);

    int rc;
    while ((rc = gaf_read1(r, &rec)) == 1) {
        if (cnt == cap) {
            int64_t ncap = cap ? cap << 1 : 1024;
            gaf_rec_t *p = (gaf_rec_t *)realloc(arr, (size_t)ncap * sizeof(*p));
            if (!p) { rc = AK_ENOMEM; break; }
            arr = p;
            cap = ncap;
        }
        // Move ownership of the record's strings into the array, then detach
        // them from `rec` so the next read does not free them.
        arr[cnt++] = rec;
        gaf_rec_init(&rec);
    }

    gaf_rec_clear(&rec);
    gaf_close(r);

    if (rc == AK_ENOMEM) {
        gaf_free(arr, cnt);
        ak_log(AK_LOG_ERROR, "gaf", "out of memory");
        return NULL;
    }

    if (n) *n = cnt;
    return arr;
}

/** Free an array from gaf_slurp(). */
void gaf_free(gaf_rec_t *recs, int64_t n) {
    if (!recs) return;
    for (int64_t i = 0; i < n; i++) gaf_rec_clear(&recs[i]);
    free(recs);
}

/** Parse the next oriented node from a GAF path string; see akhal/gaf.h. */
int gaf_path_next(const char *p, uint64_t *id, char *strand) {
    if (*p == '\0') return 0;
    *strand = *p;
    p++;
    int consumed = 1;
    uint64_t v = 0;
    while (*p >= '0' && *p <= '9') {
        v = v * 10 + (uint64_t)(*p - '0');
        p++;
        consumed++;
    }
    *id = v;
    return consumed;
}
