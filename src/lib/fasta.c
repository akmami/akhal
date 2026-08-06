#include "akhal/fasta.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/error.h"

#include "khashl.h"

#include <stdlib.h>
#include <string.h>

// name -> record index
KHASHL_MAP_INIT(KH_LOCAL, fmap_t, fmap, const char *, int64_t, kh_hash_str, kh_eq_str)

/** Ensure room for one more record. @return AK_OK or AK_ENOMEM. */
static int reserve_rec(fasta_t *fa) {
    if (fa->n < fa->m) return AK_OK;
    int64_t m = fa->m ? fa->m << 1 : 256;
    fasta_rec_t *p = (fasta_rec_t *)realloc(fa->rec, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    fa->rec = p;
    fa->m = m;
    return AK_OK;
}

/**
 * Append bytes to a kstring, keeping it NUL-terminated.
 * @param ks Destination buffer.
 * @param s Bytes to append.
 * @param n Number of bytes.
 * @return AK_OK or AK_ENOMEM.
 */
static int str_append(kstring_t *ks, const char *s, size_t n) {
    if (ks_resize(ks, ks->l + n + 1) != AK_OK) return AK_ENOMEM;
    memcpy(ks->s + ks->l, s, n);
    ks->l += n;
    ks->s[ks->l] = '\0';
    return AK_OK;
}

/**
 * Copy a header (without '>') up to the first whitespace.
 * @param hdr Header text following '>'.
 * @return Newly allocated name (caller frees), or NULL on OOM.
 */
static char *first_token(const char *hdr) {
    size_t i = 0;
    while (hdr[i] && hdr[i] != ' ' && hdr[i] != '\t') i++;
    char *s = (char *)malloc(i + 1);
    if (!s) return NULL;
    memcpy(s, hdr, i);
    s[i] = '\0';
    return s;
}

/**
 * Store the accumulated (name, seq) as a record and index it.
 * @param fa Set being built.
 * @param h The name -> index map.
 * @param name Record name; ownership is taken (freed here on failure).
 * @param seq Accumulated sequence, copied into the record.
 * @return AK_OK (including when name is NULL), or AK_ENOMEM.
 */
static int flush_record(fasta_t *fa, fmap_t *h, char *name, kstring_t *seq) {
    if (!name) return AK_OK;
    if (reserve_rec(fa) != AK_OK) { free(name); return AK_ENOMEM; }

    fasta_rec_t *r = &fa->rec[fa->n];
    r->name = name;
    r->len  = (int64_t)seq->l;
    r->seq  = (char *)malloc(seq->l + 1);
    if (!r->seq) { free(name); return AK_ENOMEM; }
    memcpy(r->seq, seq->s ? seq->s : "", seq->l);
    r->seq[seq->l] = '\0';

    int absent;
    khint_t k = fmap_put(h, r->name, &absent);
    if (absent) kh_val(h, k) = fa->n;
    else ak_log(AK_LOG_WARN, "fasta", "duplicate sequence name %s", r->name);

    fa->n++;
    return AK_OK;
}

/** Load every record from a FASTA file; see akhal/fasta.h. */
fasta_t *fasta_read(const char *fn) {
    ak_file *f = ak_open(fn);
    if (!f) return NULL;

    fasta_t *fa = (fasta_t *)calloc(1, sizeof(fasta_t));
    if (!fa) { ak_close(f); ak_log(AK_LOG_ERROR, "fasta", "out of memory"); return NULL; }

    fmap_t *h = fmap_init();
    if (!h) { free(fa); ak_close(f); ak_log(AK_LOG_ERROR, "fasta", "out of memory"); return NULL; }
    fa->idx = h;

    kstring_t line = KS_INIT, seq = KS_INIT;
    char *name = NULL;
    int rc = AK_OK;
    long len;

    while ((len = ak_getline(f, &line)) >= 0) {
        if (len == 0) continue;
        if (line.s[0] == '>') {
            rc = flush_record(fa, h, name, &seq);
            name = NULL;
            if (rc != AK_OK) break;
            ks_clear(&seq);
            name = first_token(line.s + 1);
            if (!name) { rc = AK_ENOMEM; break; }
        } else if (name) {
            rc = str_append(&seq, line.s, (size_t)len);
            if (rc != AK_OK) break;
        }
    }
    if (rc == AK_OK) rc = flush_record(fa, h, name, &seq);
    else free(name);

    ks_free(&line);
    ks_free(&seq);
    ak_close(f);

    if (rc == AK_ENOMEM) {
        fasta_destroy(fa); 
        ak_log(AK_LOG_ERROR, "fasta", "out of memory");
        return NULL;
    }
    return fa;
}

/** Look up a record by name, or NULL if absent. */
const fasta_rec_t *fasta_get(const fasta_t *fa, const char *name) {
    fmap_t *h = (fmap_t *)fa->idx;
    khint_t k = fmap_get(h, name);
    return (k < kh_end(h)) ? &fa->rec[kh_val(h, k)] : NULL;
}

/** Free a loaded FASTA set. Safe with NULL. */
void fasta_destroy(fasta_t *fa) {
    if (!fa) return;
    for (int64_t i = 0; i < fa->n; i++) {
        free(fa->rec[i].name);
        free(fa->rec[i].seq);
    }
    free(fa->rec);
    if (fa->idx) fmap_destroy((fmap_t *)fa->idx);
    free(fa);
}
