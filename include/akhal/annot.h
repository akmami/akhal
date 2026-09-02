#ifndef AKHAL_ANNOT_H
#define AKHAL_ANNOT_H

#include <stdint.h>
#include "akhal/gfa.h"
#include "akhal/fasta.h"

/**
 * Node annotation store: where does each graph node come from?
 *
 * Every annotated node maps to one info string ("SNP chr1:12345 A>G rs99",
 * "REF chr1 0-1200", "SEQ sampleA 4501", ...). All info strings live packed
 * in a single NUL-separated buffer; a record holds only the node id, the
 * byte offset of its string in that buffer, its length, and a kind. That
 * keeps the store one allocation per component ("array + dict", like the
 * graph itself), and a lookup is a hash probe plus a pointer into the buffer
 * - no per-node string allocations.
 *
 * Nodes that were never annotated simply have no record: annot_get() reports
 * them as ANNOT_UNKNOWN. Backbone (reference path) nodes are ANNOT_BACKBONE;
 * everything traced to a variant or an external sequence is ANNOT_INFO.
 *
 * Builders fill the store from a graph plus a VCF (bubble/alt-path matching
 * against the backbone) or a FASTA (greedy walk of each sequence through the
 * graph). Both are optional: with neither input, nodes are left unknown.
 *
 * The store round-trips through a compact binary file (see annot_write),
 * so annotation is computed once and queried later without the graph.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Annotation kinds returned by annot_get()
enum {
    ANNOT_UNKNOWN  = 0,   // no record for this node
    ANNOT_BACKBONE = 1,   // node lies on the backbone reference path
    ANNOT_INFO     = 2    // node carries variant / origin information
};

// One annotated node. `off`/`len` locate the info string inside annot_t.buf:
// the string is buf[off .. off+len), always NUL-terminated at off+len.
typedef struct {
    uint64_t id;     // node (segment) id
    uint64_t off;    // byte offset of the info string in the shared buffer
    uint32_t len;    // info string length, excluding the NUL terminator
    uint8_t  kind;   // ANNOT_BACKBONE or ANNOT_INFO
} annot_rec_t;

typedef struct {
    annot_rec_t *rec;      // records, one per annotated node
    int64_t      n, m;     // used / allocated record count
    char        *buf;      // packed NUL-separated info strings
    uint64_t     buf_l;    // used bytes in buf (including NUL terminators)
    uint64_t     buf_m;    // allocated bytes in buf
    void        *idx;      // opaque id -> record-index hash table
} annot_t;

/**
 * Allocate an empty annotation store
 * @return The store (release with annot_destroy), or NULL on failure
 */
annot_t *annot_init(void);

/**
 * Release a store and everything it owns. Safe to call with NULL
 * @param a Store to destroy.
 */
void annot_destroy(annot_t *a);

// Writing annotations

/**
 * Set (or replace) a node's annotation. The info string is copied into the
 * shared buffer; replacing leaves the old bytes orphaned, which is the price
 * of keeping the buffer append-only.
 * @param a Store to modify
 * @param id Node id
 * @param kind ANNOT_BACKBONE or ANNOT_INFO
 * @param info Info string to copy; NULL or "" stores an empty annotation
 * @return AK_OK or AK_ENOMEM
 */
int annot_set(annot_t *a, uint64_t id, int kind, const char *info);

/**
 * Add to a node's annotation: like annot_set(), but if the node already has
 * a non-empty info string the new one is appended after "; " so multiple
 * origins (e.g. overlapping variants) accumulate instead of overwriting.
 * @param a Store to modify
 * @param id Node id
 * @param kind Kind to record (an existing kind is kept)
 * @param info Info string to add
 * @return AK_OK or AK_ENOMEM
 */
int annot_add(annot_t *a, uint64_t id, int kind, const char *info);

// Querying

/**
 * Look up a node's annotation.
 * Pointers returned through `info` are into the shared buffer and are
 * invalidated by the next annot_set()/annot_add() call.
 * @param a Store to query
 * @param id Node id
 * @param info If non-NULL: set to the info string, or NULL when the node is
 *             unknown or has an empty annotation
 * @return ANNOT_BACKBONE or ANNOT_INFO if a record exists, else ANNOT_UNKNOWN
 */
int annot_get(const annot_t *a, uint64_t id, const char **info);

/**
 * Human-readable name of an annotation kind.
 * @param kind An ANNOT_* value.
 * @return "backbone", "annot", or "unknown"; never NULL.
 */
const char *annot_kind_name(int kind);

/** 
 * @return Number of annotated nodes (records) in the store
 */
static inline int64_t annot_n(const annot_t *a) {
    return a->n;
}

/** 
 * @return Record at array index i (for iterating over all annotations)
 */
static inline const annot_rec_t *annot_at(const annot_t *a, int64_t i) {
    return &a->rec[i];
}

/** 
 * @return Info string of record i, or NULL if it is empty
 */
static inline const char *annot_info_at(const annot_t *a, int64_t i) {
    return a->rec[i].len ? a->buf + a->rec[i].off : NULL;
}

// Builders

/**
 * Identify the backbone reference path of a graph and mark every segment on
 * it as ANNOT_BACKBONE with info "REF <path> <start>-<end>".
 * @param a Store to fill
 * @param g Graph read with GFA_PATHS
 * @param ref_path Path name to use as backbone, or NULL to take the graph's
 *                 first path (logged when several are present)
 * @return The backbone path index (>= 0), or a negative AK_E* code when the
 *         graph has no paths / no path of that name (nothing is marked)
 */
int32_t annot_backbone(annot_t *a, const gfa_t *g, const char *ref_path);

/**
 * Annotate alternative-path nodes from a VCF.
 *
 * For each variant whose CHROM matches the backbone path name, the shared
 * REF/ALT prefix is stripped and the branch point located on the backbone
 * (the backbone segment ending at the variant's 0-based position). Each
 * non-backbone successor is then walked - node by node, following the
 * single non-backbone chain - and if the concatenated sequence equals the
 * alternate allele, every node on that walk is annotated
 * "<TYPE> <chrom>:<pos> <ref>><alt> [<id>]" (TYPE one of SNP/MNP/INS/COMPLEX).
 * Pure deletions produce only an edge, so they have no node to annotate and
 * are skipped. Variants that cannot be matched are logged at debug level and
 * left alone - their nodes simply stay unknown.
 *
 * @param a Store to fill (backbone must already be marked, see annot_backbone)
 * @param g Graph read with GFA_LINKS | GFA_PATHS
 * @param ref_path Backbone path index returned by annot_backbone()
 * @param vcf_fn Path to the VCF file
 * @return Number of alleles annotated (>= 0), or a negative AK_E* code
 */
int64_t annot_build_vcf(annot_t *a, const gfa_t *g, int32_t ref_path, const char *vcf_fn);

/**
 * Annotate nodes by walking external sequences through the graph.
 *
 * For each FASTA record, a start node is located (a source node whose
 * sequence prefixes the record, falling back to any matching node) and the
 * walk greedily follows the out-link whose sequence continues the record
 * (link overlaps are honoured). Every visited non-backbone node is annotated
 * "SEQ <name> <offset>", where offset is the node's 0-based position in the
 * record. A walk that stalls before consuming the whole sequence is logged
 * and keeps the annotations gathered so far.
 *
 * @param a Store to fill
 * @param g Graph read with GFA_LINKS (GFA_PATHS if a backbone is marked)
 * @param fa Loaded FASTA set to trace
 * @return Number of sequences fully traversed (>= 0), or a negative AK_E* code
 */
int64_t annot_build_fasta(annot_t *a, const gfa_t *g, const fasta_t *fa);

// File I/O

/**
 * Write a store to a binary .annot file.
 *
 * Layout (native endianness):
 *   magic "AKANNOT1" (8 bytes)
 *   uint64 n_rec, uint64 buf_len
 *   n_rec x { uint64 id, uint64 off, uint32 len, uint8 kind }
 *   buf_len bytes of NUL-separated info strings
 *
 * @param a Store to write
 * @param fn Output path
 * @return AK_OK, or a negative AK_E* code
 */
int annot_write(const annot_t *a, const char *fn);

/**
 * Load a store previously written by annot_write()
 * @param fn Path to the .annot file
 * @return The store (release with annot_destroy), or NULL on failure (logged)
 */
annot_t *annot_read(const char *fn);

#ifdef __cplusplus
}
#endif

#endif
