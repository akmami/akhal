#ifndef AKHAL_GAF_H
#define AKHAL_GAF_H

#include <stdint.h>

/**
 * GAF (Graph Alignment Format) records and readers.
 *
 * Two consumption styles are offered, and callers pick whichever suits them:
 *
 *   streaming  gaf_open() + gaf_read1() pull one record at a time into a
 *              caller-owned, reusable gaf_rec_t (cf. htslib sam_read1), which
 *              keeps memory flat for arbitrarily large alignment files.
 *
 *   batch      gaf_slurp() loads an entire file into an array in one call,
 *              for when random access or multiple passes are convenient.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char   *qname;       // query/read name (owned)
    int64_t qlen;        // query length
    int64_t qstart;      // query start (0-based)
    int64_t qend;        // query end
    char    strand;      // '+' or '-'
    char   *path;        // path string, e.g. ">1>2<3" (owned)
    int64_t plen;        // path length in bp
    int64_t pstart;      // path start
    int64_t pend;        // path end
    int64_t matches;     // number of matching bases
    int64_t block_len;   // alignment block length
    int     mapq;        // mapping quality (255 = unavailable)

    // optional tags; the has_* flags say whether the tag was present
    int     nm;   int has_nm;   // NM:i  mismatches + indels
    double  as;   int has_as;   // AS:f  alignment score
    double  dv;   int has_dv;   // dv:f  divergence
    double  id;   int has_id;   // id:f  identity
    char   *cigar;              // cg:Z difference CIGAR (owned, NULL if absent)
} gaf_rec_t;

/**
 * Zero-initialize a record. Does not allocate
 * @param r Record to initialize
 */
void gaf_rec_init(gaf_rec_t *r);

/**
 * Free any owned strings and reset the record for reuse
 * @param r Record to clear
 */
void gaf_rec_clear(gaf_rec_t *r);

// Streaming reader

typedef struct gaf_reader gaf_reader_t;

/**
 * Open a GAF file for streaming
 * @param fn Path to the .gaf file
 * @return A reader, or NULL on failure (logged)
 */
gaf_reader_t *gaf_open(const char *fn);

/**
 * Parse the next record. The record is cleared first, so the same one can be
 * reused across the whole file
 * @param r Open reader
 * @param rec Destination record (reused each call)
 * @return 1 on success, 0 at end of file, or a negative AK_E* code on error
 */
int gaf_read1(gaf_reader_t *r, gaf_rec_t *rec);

/**
 * Close a reader. Safe with NULL. Does not touch caller-owned records
 * @param r Reader to close
 */
void gaf_close(gaf_reader_t *r);

// Batch reader

/**
 * Load an entire GAF file into a newly allocated array
 * @param fn Path to the .gaf file
 * @param n Set to the number of records read
 * @return The record array (release with gaf_free), or NULL on failure
 */
gaf_rec_t *gaf_slurp(const char *fn, int64_t *n);

/**
 * Free an array returned by gaf_slurp(); clears each record, then the array
 * @param recs Array to free
 * @param n Number of records in the array
 */
void gaf_free(gaf_rec_t *recs, int64_t n);

// Path string traversal

/**
 * Parse the next oriented node from a GAF path string such as ">12<7>3"
 * Advance with:  p += gaf_path_next(p, &id, &strand);
 * @param p Pointer into a path string
 * @param id Set to the node id
 * @param strand Set to the orientation char ('>' or '<')
 * @return Bytes consumed, or 0 at the end of the string
 */
int gaf_path_next(const char *p, uint64_t *id, char *strand);

#ifdef __cplusplus
}
#endif

#endif
