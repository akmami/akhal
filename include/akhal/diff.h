#ifndef AKHAL_DIFF_H
#define AKHAL_DIFF_H

#include <stdint.h>
#include "akhal/gfa.h"

/**
 * Structural comparison of two graphs that do not agree on segment ids.
 *
 * Two files can spell the same graph and still share not a single id: builders
 * number nodes as they emit them, and `akhal sort` renumbers them again. So the
 * comparison never looks at an id. Segments are matched on their sequence
 * content, both graphs are relabelled onto one shared numbering, and everything
 * else is compared through those labels.
 *
 * The matching is a merge pass over both segment arrays sorted by sequence, so
 * it costs one sort per graph and a single linear walk. A label is a class of
 * equal sequences, not one segment: a graph full of single-base SNP nodes
 * carries "A" many times over, and which copy in one file "is" which copy in
 * the other is not something the two files agree on. Counting is still
 * multiset-style - with three copies here and two there, two match and the
 * third is reported as unmatched - but every copy shares one label, so no link
 * inherits an arbitrary pairing.
 *
 * Links are then relabelled and compared as (from, orientation, to,
 * orientation, overlap) tuples, canonicalized so that the two spellings of one
 * edge - `L a + b +` and `L b - a -` - are recognized as the same link. Where a
 * sequence repeats, this asks whether the other graph has that edge between
 * those sequences rather than between those particular nodes.
 *
 * Paths are compared by what they spell rather than by what they walk over:
 * P lines are paired by name across the two graphs and each pair's sequences
 * are compared. Two graphs that chop one reference into different nodes
 * therefore still report that path as identical.
 *
 * The same header also carries the alignment-side comparison, diff_gaf(),
 * which asks the corresponding question of two GAF files: how many of their
 * reads are put along the same walk. That one does look at ids, because two
 * GAF files are only comparable when they were aligned against the same graph.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Shared labelling

/**
 * A shared numbering of two graphs' segments, one label per distinct sequence.
 *
 * `a[i]` is the label of segment index i in the first graph and `b[j]` the
 * label of segment index j in the second. Both graphs draw from one counter, so
 * equal labels mean equal sequences and unequal labels mean unequal ones - in
 * either graph, or across the two. Several segments of one graph share a label
 * exactly when they spell the same bases.
 */
typedef struct {
    uint32_t *a;         // length n_a: segment index in A -> shared label
    uint32_t *b;         // length n_b: segment index in B -> shared label
    int32_t   n_a, n_b;  // segment counts of the two graphs
    int32_t   n_class;   // distinct sequences; labels run 0 .. n_class - 1
    int32_t   n_shared;  // segments that pair off, summed over the classes
} diff_map_t;

/**
 * Match two graphs' segments by sequence content and label them together.
 *
 * A segment with no sequence is treated as carrying the empty one, so two such
 * segments match each other. Needs no read flags: only the S lines are used
 * @param a First graph
 * @param b Second graph
 * @return The labelling (release with diff_map_destroy), or NULL on failure
 */
diff_map_t *diff_map(const gfa_t *a, const gfa_t *b);

/**
 * Release a labelling and everything it owns. Safe to call with NULL
 * @param m Labelling to destroy
 */
void diff_map_destroy(diff_map_t *m);

// Comparison

/**
 * One link, reported with the segment ids of the graph it came from
 */
typedef struct {
    uint64_t from, to;       // segment ids as that graph spells them
    char     from_orient;    // orientation of `from` on the L line ('+' or '-')
    char     to_orient;      // orientation of `to` on the L line ('+' or '-')
    uint32_t overlap;        // overlap length in bp
} diff_link_t;

// How a path name fared. A name present in both graphs is one entry, not two
enum {
    DIFF_SAME   = 0,   // both graphs spell it, identically
    DIFF_DIFFER = 1,   // both graphs spell it, differently
    DIFF_A_ONLY = 2,   // only the first graph has it
    DIFF_B_ONLY = 3    // only the second graph has it
};

/**
 * One path name's verdict
 */
typedef struct {
    char    *name;     // owned: the path name, as the P line spells it
    int      state;    // DIFF_SAME / DIFF_DIFFER / DIFF_A_ONLY / DIFF_B_ONLY
    uint64_t len_a;    // bases the first graph spells for it, 0 when absent
    uint64_t len_b;    // bases the second graph spells for it, 0 when absent
} diff_path_t;

/**
 * What one graph carries that the other does not
 */
typedef struct {
    uint64_t    *seg;      // owned: ids of segments only this graph has
    int32_t      n_seg;
    diff_link_t *link;     // owned: links only this graph has
    int32_t      n_link;
} diff_side_t;

/**
 * The result of comparing two graphs: what they share, and what each holds
 * alone. Shared counts are pair counts - `n_seg_shared` segments matched means
 * that many segments on each side.
 *
 * Where a sequence repeats, the ids listed as unmatched are whichever copies
 * file order left over; it is their number that is meaningful, not which ones
 */
typedef struct {
    diff_side_t  a, b;             // what each graph alone carries
    int32_t      n_seg_shared;     // segments matched one-to-one on sequence
    int32_t      n_link_shared;    // links matched one-to-one after relabelling

    diff_path_t *path;             // owned: one entry per name, ordered by name
    int32_t      n_path;
    int32_t      n_path_same;      // paired by name, spelling the same bases
    int32_t      n_path_differ;    // paired by name, spelling different bases
    int32_t      n_path_a_only;
    int32_t      n_path_b_only;
} diff_t;

/**
 * Compare two graphs: segments by sequence, links through the resulting
 * labelling, paths by the bases their chains spell.
 *
 * Both graphs must be read with GFA_LINKS | GFA_PATHS, though a graph carrying
 * no P lines at all is fine: it simply contributes no chains, and its segments
 * and links still compare. Path sequences are spelled a pair at a time rather
 * than all at once, so the peak cost is the two longest chains that share a
 * name, and link overlaps are not trimmed off the bases - the same blunt-join
 * assumption `extract path` makes
 * @param a First graph
 * @param b Second graph
 * @return The comparison (release with diff_destroy), or NULL on failure
 */
diff_t *diff_graphs(const gfa_t *a, const gfa_t *b);

/**
 * Release a comparison and everything it owns. Safe to call with NULL
 * @param d Comparison to destroy
 */
void diff_destroy(diff_t *d);

/**
 * Whether the two graphs came out equal: nothing unmatched on either side and
 * every path name spelling the same bases. Ids, node numbering, line order and
 * SR ranks are not part of it
 * @param d Comparison to test; must not be NULL
 * @return 1 when the graphs match, otherwise 0
 */
static inline int diff_identical(const diff_t *d) {
    return d->a.n_seg == 0 && d->b.n_seg == 0 &&
           d->a.n_link == 0 && d->b.n_link == 0 &&
           d->n_path_differ == 0 && d->n_path_a_only == 0 && d->n_path_b_only == 0;
}

// GAF alignment comparison

/**
 * Comparison of two GAF files: how much of one file's alignment set the other
 * one also carries.
 *
 * The question is whether the two aligners put the same read along the same
 * walk, so an alignment is reduced to the pair (read name, path) and nothing
 * else - coordinates, scores, mapping quality and CIGARs are not part of it.
 * The walk is compared in a canonical spelling, since `>1>2<3` and `>3>2<1`
 * are the same walk read from its two ends and no aligner is obliged to pick
 * one of them.
 *
 * Both files are loaded and sorted by (read name, first node id of the
 * canonical path, whole path), which lets the comparison be a single merge
 * walk over the two sorted arrays rather than a lookup per alignment. Reads
 * are handled a name-block at a time: within a block that both files carry,
 * alignments pair off one-to-one on their path, so a read aligned three ways
 * here and twice there reports two pairs and one leftover instead of simply
 * "matching".
 */

// How one alignment fared against the other file
enum {
    DIFF_ALN_SHARED = 0,   // the other file walks this read along the same path
    DIFF_ALN_A_ONLY = 1,   // only the first file has it
    DIFF_ALN_B_ONLY = 2    // only the second file has it
};

/**
 * One alignment's verdict. A pair that matched is one entry, not two
 */
typedef struct {
    char    *qname;      // owned: read name
    char    *path;       // owned: the walk, in its canonical spelling
    uint64_t first;      // first node id of that spelling; 0 for a named path
    int      state;      // DIFF_ALN_SHARED / DIFF_ALN_A_ONLY / DIFF_ALN_B_ONLY
    int      read_both;  // 1 when the read name occurs in both files
} diff_aln_t;

/**
 * The result of comparing two GAF files, counted per alignment and per read.
 *
 * Alignment counts are pair counts on the shared side: `n_aln_shared` pairs
 * means that many alignments in each file. Read counts partition the names -
 * `n_read_shared` splits into `n_read_all_same`, `n_read_partial` and
 * `n_read_none` - so the three read categories plus the two one-sided ones
 * account for every name in either file
 */
typedef struct {
    diff_aln_t *aln;             // owned: one entry per pair and per leftover
    int64_t     n_aln;           // entries in `aln`

    int64_t     n_aln_a, n_aln_b;        // alignments read from each file
    int64_t     n_aln_shared;            // pairs matched on (read, path)
    int64_t     n_aln_a_only;            // alignments the second file lacks
    int64_t     n_aln_b_only;            // alignments the first file lacks

    int64_t     n_read_a, n_read_b;      // distinct read names in each file
    int64_t     n_read_shared;           // names both files align
    int64_t     n_read_a_only;           // names only the first file aligns
    int64_t     n_read_b_only;           // names only the second file aligns
    int64_t     n_read_all_same;         // shared names whose alignments all pair off
    int64_t     n_read_partial;          // shared names with a pair and a leftover
    int64_t     n_read_none;             // shared names with no pairing at all
} diff_gaf_t;

/**
 * Compare two GAF files by the walks they put their reads along.
 *
 * Malformed lines are skipped with a warning, as the reader does everywhere
 * else; an empty file is not an error and simply contributes nothing. Both
 * files are held in memory for the duration, though only each alignment's read
 * name and canonical path are kept, not its whole record
 * @param fn_a First GAF file
 * @param fn_b Second GAF file
 * @return The comparison (release with diff_gaf_destroy), or NULL on failure
 */
diff_gaf_t *diff_gaf(const char *fn_a, const char *fn_b);

/**
 * Release a GAF comparison and everything it owns. Safe to call with NULL
 * @param d Comparison to destroy
 */
void diff_gaf_destroy(diff_gaf_t *d);

/**
 * Whether the two files align the same reads along the same walks, leaving
 * nothing over on either side. Coordinates, scores and record order are not
 * part of it
 * @param d Comparison to test; must not be NULL
 * @return 1 when the two files agree, otherwise 0
 */
static inline int diff_gaf_identical(const diff_gaf_t *d) {
    return d->n_aln_a_only == 0 && d->n_aln_b_only == 0;
}

#ifdef __cplusplus
}
#endif

#endif
