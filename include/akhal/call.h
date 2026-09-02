#ifndef AKHAL_CALL_H
#define AKHAL_CALL_H

#include <stdint.h>
#include "akhal/gfa.h"
#include "akhal/fasta.h"

/**
 * Variant discovery: read a graph as a reference plus its differences.
 *
 * One walk through the graph is declared the backbone - a P line, or an
 * external sequence traced through the nodes - and every node on it gets a
 * coordinate on that reference. Everything the backbone does not cover is, by
 * definition, a variation: each detour that leaves the backbone at one node
 * and rejoins it at a later one spells an alternate allele over the reference
 * span between them, and a link that skips forward along the backbone itself
 * spells a deletion of the bases it jumps over.
 *
 * Detours that leave from the same node and rejoin at the same node share a
 * REF span, so they are grouped into one multi-allelic record. The result is
 * an ordered array of records that maps directly onto VCF rows, which
 * call_write_vcf() emits.
 *
 * Coordinates assume blunt joins (0 bp link overlaps), which is what graph
 * builders like vg produce; a backbone carrying real overlaps is reported and
 * its coordinates shift by the overlap length.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Caps that keep a tangled region from exploding; exceeding one is logged
#define CALL_MAX_ALT   16        // alternate alleles collected per branch point
#define CALL_MAX_WALK  1024      // nodes visited along one alternate walk
#define CALL_MAX_LEN   1000000   // bases spelled by one alternate allele

/**
 * A labelled backbone: the reference sequence plus, for every segment, whether
 * it lies on that reference and where.
 */
typedef struct {
    char    *name;    // owned: reference name, used as the VCF CHROM column
    char    *seq;     // owned: backbone sequence, 5' -> 3'
    int64_t  len;     // length of seq
    uint8_t *on;      // length n_seg: 1 when the segment lies on the backbone
    int64_t *pos;     // length n_seg: 0-based start offset, -1 when off it
    int32_t  n_seg;   // segment count of the graph this was built from
    uint32_t *walk;   // owned: the backbone's segments in the order it visits
    int64_t  n_walk;  // length of walk; a repeated segment appears each time
} call_ref_t;

/**
 * One variation, already grouped over its alternate alleles. Offsets are
 * 0-based and half-open, so the REF allele is backbone[pos .. end)
 */
typedef struct {
    int64_t  pos;     // 0-based offset of the first REF base
    int64_t  end;     // 0-based, exclusive: one past the last REF base
    char    *ref;     // owned: REF allele
    char    *alt;     // owned: ALT alleles, comma-separated
    char    *type;    // owned: SNP/MNP/INS/DEL/COMPLEX tags, one per ALT
} call_var_t;

typedef struct {
    call_var_t *var;
    int64_t     n, m;   // used / allocated record count
} call_t;

// Backbone

/**
 * Label the backbone from the graph's P lines.
 *
 * Because builders such as vg split one reference over several consecutive P
 * lines, the candidates are the chains gfa_path_merge() stitches out of them
 * rather than the raw paths - so `path_name` may name either a whole reference
 * ("chr22") or one of its fragments. Their segments are concatenated in order
 * to give the reference sequence, and each is stamped with its offset in it.
 * A segment the backbone visits more than once keeps its first occurrence
 * @param g Graph read with GFA_LINKS | GFA_PATHS
 * @param path_name Reference to use, or NULL to take the chain holding the
 *                  graph's first path (logged when several are present)
 * @return The backbone (release with call_ref_destroy), or NULL on failure
 */
call_ref_t *call_ref_path(const gfa_t *g, const char *path_name);

/**
 * Label the backbone by tracing an external sequence through the graph.
 *
 * A start node is located (a source node whose sequence prefixes the record,
 * falling back to any matching node) and the walk greedily follows the
 * out-link whose sequence continues the record, honouring link overlaps.
 * A walk that stalls is reported and the backbone is truncated there, so only
 * the covered prefix yields variants
 * @param g Graph read with GFA_LINKS
 * @param fa Loaded FASTA set holding the reference
 * @param seq_name Record to trace, or NULL to take the first record (logged
 *                 when several are present)
 * @return The backbone (release with call_ref_destroy), or NULL on failure
 */
call_ref_t *call_ref_fasta(const gfa_t *g, const fasta_t *fa, const char *seq_name);

/**
 * Release a backbone and everything it owns. Safe to call with NULL
 * @param r Backbone to destroy
 */
void call_ref_destroy(call_ref_t *r);

// Variants

/**
 * Collect every variation the graph carries over a labelled backbone.
 *
 * Each backbone node is examined in reference order. A link into a node off
 * the backbone starts a bounded walk that ends where it rejoins, and the bases
 * it spells become an alternate allele; a link that lands further along the
 * backbone than the node it left is a deletion of the span it skips. Detours
 * sharing a REF span are merged into one multi-allelic record, duplicates
 * dropped, and the result is ordered by (pos, end).
 *
 * Records where either side of the difference is empty carry the preceding
 * reference base on both alleles, as VCF requires
 * @param g Graph read with GFA_LINKS
 * @param ref Backbone from call_ref_path() or call_ref_fasta()
 * @return The variants (release with call_destroy), or NULL on failure
 */
call_t *call_variants(const gfa_t *g, const call_ref_t *ref);

/**
 * Release a variant set and everything it owns. Safe to call with NULL
 * @param c Variant set to destroy
 */
void call_destroy(call_t *c);

/**
 * @return Number of variant records
 */
static inline int64_t call_n(const call_t *c) {
    return c->n;
}

/**
 * @return Record at array index i
 */
static inline const call_var_t *call_at(const call_t *c, int64_t i) {
    return &c->var[i];
}

/**
 * Write a variant set as VCF 4.2.
 *
 * One row per record, with POS converted to 1-based, ID/QUAL/FILTER left as
 * '.', and INFO carrying END (1-based inclusive end of the REF allele) and
 * TYPE (one tag per ALT). The backbone supplies the CHROM name and the
 * ##contig header line
 * @param c Variant set to write
 * @param ref Backbone the variants were called against
 * @param fn Output path
 * @return AK_OK, or a negative AK_E* code
 */
int call_write_vcf(const call_t *c, const call_ref_t *ref, const char *fn);

#ifdef __cplusplus
}
#endif

#endif
