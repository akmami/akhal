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
 * fragmented P lines are chained with gfa_path_merge() first, chains are paired
 * by name across the two graphs, and each pair's sequences are compared. Two
 * graphs that chop one reference into different nodes therefore still report
 * that path as identical.
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
 * One path name's verdict, after fragmented P lines were chained together
 */
typedef struct {
    char    *name;     // owned: chain name, as gfa_path_merge() named it
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

#ifdef __cplusplus
}
#endif

#endif
