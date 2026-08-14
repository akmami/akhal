#ifndef AKHAL_FASTA_H
#define AKHAL_FASTA_H

#include <stdint.h>

/**
 * A whole FASTA file loaded into memory, indexed by sequence name.
 *
 * Records live in a contiguous array and a hash table maps the sequence name
 * (the first whitespace-delimited token of the header) to its array index,
 * following the same "array + dict" pattern as the graph. This backs both the
 * read store used by gaf2sam and the reference used by sampoke.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char   *name;   // owned: sequence name (header up to first whitespace)
    char   *seq;    // owned: concatenated sequence, case preserved
    int64_t len;    // sequence length
} fasta_rec_t;

typedef struct {
    fasta_rec_t *rec;
    int64_t      n, m;
    void        *idx;   // opaque name -> index hash table
} fasta_t;

/**
 * Load every record from a FASTA file
 * @param fn Path to the FASTA file
 * @return The loaded set (release with fasta_destroy), or NULL on failure
 */
fasta_t *fasta_read(const char *fn);

/**
 * Look up a record by sequence name
 * @param fa Loaded FASTA set
 * @param name Sequence name to find
 * @return Pointer to the record, or NULL if absent
 */
const fasta_rec_t *fasta_get(const fasta_t *fa, const char *name);

/**
 * Release everything. Safe with NULL
 * @param fa Set to destroy
 */
void fasta_destroy(fasta_t *fa);

/** 
 * @return Number of records loaded
 */
static inline int64_t fasta_n(const fasta_t *fa) { return fa->n; }

#ifdef __cplusplus
}
#endif

#endif
