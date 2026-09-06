#ifndef AKHAL_COMPACT_H
#define AKHAL_COMPACT_H

#include <stdint.h>
#include <stdio.h>
#include "akhal/gfa.h"

/**
 * Compaction: a run of nodes with nowhere else to go becomes one node.
 *
 * Builders leave graphs full of chains - node after node joined by a single
 * link, with no branch anywhere along the way. Nothing is expressed by keeping
 * them apart, so a run like that folds into one segment carrying the same
 * bases, which is what "unitig" means in a de Bruijn graph and what `vg mod -u`
 * and `odgi unchop` do elsewhere.
 *
 * A join u -> v is taken when all of this holds:
 *
 *   - the link leaves u's right end and enters v's left end on the same strand,
 *     written `L u + v +` or, the same edge read from the other side,
 *     `L v - u -`. A join that flips strand is left alone
 *   - it is the only link on either of those two ends, so nothing branches
 *   - the link is a blunt join (0 bp overlap): bases the two share would have
 *     to be dropped, and this never rewrites sequence it was not given
 *   - u is not v, so a self-loop stays a self-loop
 *   - no path starts or ends at either of those ends. A path stopping in the
 *     middle of what would become one node cannot be written down afterwards,
 *     so that run is cut there instead
 *   - on an rGFA, the two agree on the stable sequence they sit on and on their
 *     rank, and their offsets are contiguous, so the merged node inherits a
 *     coordinate that still means something. On a plain GFA there is nothing to
 *     check, and nothing is tagged on the way out
 *
 * Every link and every path survives unchanged: links between runs are
 * repointed at the nodes that swallowed their endpoints, and a path walking a
 * run step by step comes out walking the one node instead.
 *
 * A caution on those tags. The reader does not keep a file's SN, and a graph
 * read with GFA_PATHS has its SO overwritten by the layout of the path that
 * owns each segment (see gfa_read). So SN and SO here are the owning path and
 * the offset along it, not whatever the input file wrote - the output is a
 * consistent rGFA, but for anything but the backbone its offsets are the path's
 * own rather than the ones rgfa_build() works out. Compacting first and
 * labelling afterwards is the order that leaves both saying what they mean.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Which run each segment ended up in.
 *
 * A run is a maximal chain of joins; a segment nothing could be merged with is
 * a run of one. Runs are numbered in the order their first segment appears, so
 * `n_run` is the segment count after compaction
 */
typedef struct {
    uint32_t *run;       // length n_seg: the run a segment belongs to
    int32_t  *pos;       // length n_seg: how far along that run it sits
    uint32_t *next;      // length n_seg: the segment after it, GFA_NIL at the end
    uint32_t *first;     // length n_run: the segment each run starts at
    int32_t  *count;     // length n_run: how many segments it gathered
    int32_t   n_seg;     // segments before compaction
    int32_t   n_run;     // segments after it
    int32_t   n_merged;  // segments folded into an earlier one: n_seg - n_run
} compact_t;

/**
 * Work out the runs. The graph is not modified.
 *
 * Requires GFA_LINKS; GFA_PATHS as well if the graph has any paths worth
 * protecting, since without it a path stopping mid-run cannot be seen
 * @param g Graph to examine
 * @return The runs (release with compact_destroy), or NULL on failure
 */
compact_t *compact_runs(const gfa_t *g);

/**
 * Release a run set. Safe to call with NULL
 * @param c Run set to destroy
 */
void compact_destroy(compact_t *c);

/**
 * Write the compacted graph: one S line per run, carrying the concatenated
 * sequence and the first segment's id, then the links and paths that survive.
 *
 * Ids are the first segment of each run, so a node that merged nothing keeps
 * the id it had. rGFA tags are emitted when the file carried its own SR tags,
 * the same rule `sort` follows - a plain GFA comes back as a plain GFA rather
 * than gaining ranks it never had
 * @param g Graph the runs came from
 * @param c Run set from compact_runs()
 * @param out Destination stream
 * @return AK_OK, or a negative AK_E* code
 */
int compact_write(const gfa_t *g, const compact_t *c, FILE *out);

#ifdef __cplusplus
}
#endif

#endif
