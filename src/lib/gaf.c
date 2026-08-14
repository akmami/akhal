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

// zero-initialize a record
void gaf_rec_init(gaf_rec_t *r) {
    memset(r, 0, sizeof(*r));
}

// free owned strings and reset a record for reuse
void gaf_rec_clear(gaf_rec_t *r) {
    free(r->qname);
    free(r->path);
    free(r->cigar);
    memset(r, 0, sizeof(*r));
}

/**
 * Parse one GAF line into a cleared record
 * @param line The line (modified in place by tokenizing)
 * @param rec Destination record (cleared first)
 * @return AK_OK, AK_EFORMAT (too few fields), or AK_ENOMEM
 */
static int gaf_parse_line(char *line, gaf_rec_t *rec) {
    gaf_rec_clear(rec);

    char *save;
    char *tok = strtok_r(line, "\t", &save);
    int field = 0;

    while (tok) {
        switch (field) {
            case 0: rec->qname = strdup(tok); if (!rec->qname) return AK_ENOMEM; break;
            case 1: rec->qlen      = strtoll(tok, NULL, 10); break;
            case 2: rec->qstart    = strtoll(tok, NULL, 10); break;
            case 3: rec->qend      = strtoll(tok, NULL, 10); break;
            case 4: rec->strand    = tok[0]; break;
            case 5: rec->path = strdup(tok); if (!rec->path) return AK_ENOMEM; break;
            case 6: rec->plen      = strtoll(tok, NULL, 10); break;
            case 7: rec->pstart    = strtoll(tok, NULL, 10); break;
            case 8: rec->pend      = strtoll(tok, NULL, 10); break;
            case 9: rec->matches   = strtoll(tok, NULL, 10); break;
            case 10: rec->block_len = strtoll(tok, NULL, 10); break;
            case 11: rec->mapq     = atoi(tok); break;
            default:
                if      (!strncmp(tok, "NM:i:", 5)) { rec->nm = atoi(tok + 5);       rec->has_nm = 1; }
                else if (!strncmp(tok, "AS:f:", 5)) { rec->as = atof(tok + 5);       rec->has_as = 1; }
                else if (!strncmp(tok, "dv:f:", 5)) { rec->dv = atof(tok + 5);       rec->has_dv = 1; }
                else if (!strncmp(tok, "id:f:", 5)) { rec->id = atof(tok + 5);       rec->has_id = 1; }
                else if (!strncmp(tok, "cg:Z:", 5)) {
                    rec->cigar = strdup(tok + 5);
                    if (!rec->cigar) return AK_ENOMEM;
                }
                break;
        }
        field++;
        tok = strtok_r(NULL, "\t", &save);
    }

    if (field < 12) {
        ak_log(AK_LOG_WARN, "gaf", "line has only %d fields (need >= 12)", field);
        return AK_EFORMAT;
    }
    return AK_OK;
}

// streaming

// open a GAF file for streaming; see akhal/gaf.h
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
 * Read the next well-formed record, skipping malformed lines
 * @param r Open reader
 * @param rec Destination record (reused each call)
 * @return 1 on success, 0 at EOF, AK_EINVAL/AK_ENOMEM on error
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

// close a streaming reader. Safe with NULL
void gaf_close(gaf_reader_t *r) {
    if (!r) return;
    ak_close(r->f);
    ks_free(&r->line);
    free(r);
}

// batch

// load an entire GAF file into an array; see akhal/gaf.h
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
        // move ownership of the record's strings into the array, then detach them from `rec` so the next read does not free them
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

// free an array from gaf_slurp()
void gaf_free(gaf_rec_t *recs, int64_t n) {
    if (!recs) return;
    for (int64_t i = 0; i < n; i++) gaf_rec_clear(&recs[i]);
    free(recs);
}

// parse the next oriented node from a GAF path string; see akhal/gaf.h
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
