#ifndef AKHAL_RGFA_H
#define AKHAL_RGFA_H

#include <stdint.h>
#include "akhal/gfa.h"

/**
 * Turning a plain GFA into an rGFA: stable-sequence tags for every segment.
 *
 * rGFA says where each segment sits on a real sequence, through three tags -
 * SN:Z: the name of that sequence, SO:i: the offset on it, SR:i: how far the
 * segment is from the linear reference. A GFA carries none of them, but its P
 * lines say everything needed to work them out.
 *
 * One path is declared the backbone. Its segments are rank 0 and take their
 * offsets from walking it end to end. Every other path is then walked in turn:
 * while it runs over ground that is already labelled it only follows along, and
 * where it leaves that ground the segments it visits alone are one rank deeper,
 * named after that path, and offset onward from the point it left - so a bubble
 * carries the coordinate of the reference stretch it detours around, and the
 * counting stops as soon as the path merges back down to a lower rank.
 *
 * This reads a variation graph, where paths are haplotypes over a shared
 * backbone. An assembly graph, whose paths are unrelated contigs, will label
 * badly - there is no meaningful backbone for the rest to be an offset from.
 *
 * Where two paths disagree, nothing is invented. A segment that one path places
 * at one offset and another path at a different one has no single answer, so it
 * keeps its rank and loses SN and SO, as does everything the walk reaches after
 * it - until the walk arrives somewhere an authoritative path, the backbone
 * above all, has already pinned down.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * What the labelling did, for the caller to report
 */
typedef struct {
    int32_t n_path;        // paths after fragments were consolidated
    int32_t n_rank0;       // segments on the backbone
    int32_t n_labelled;    // segments that came out with SN and SO
    int32_t n_ambiguous;   // segments a walk reached but could not place
    int32_t n_unreached;   // segments no path visits at all
    int32_t max_rank;      // deepest rank handed out
} rgfa_stat_t;

/**
 * Label a graph as rGFA, in place.
 *
 * Fragmented P lines are consolidated first, exactly as gfa_path_merge()
 * chains them, so a reference that arrived as "chr22:0-1000",
 * "chr22:1000-2000", ... leaves as one path and one stable sequence name. The backbone is then
 * labelled rank 0 and every other path walked as described above.
 *
 * Segments are labelled through the fields gfa_write_rgfa() emits: `rank` is
 * SR, `ref_name` is SN, and `start`/`end` are the SO offset and its end. A
 * segment with no SN/SO carries -1 in `start` and NULL in `ref_name`, which is
 * how an ambiguous or unvisited one is told from a placed one.
 *
 * Requires the graph to be read with GFA_LINKS | GFA_PATHS
 * @param g Graph to label, modified in place
 * @param ref_name Path to use as the backbone, or NULL for the graph's first
 * @param st Filled in with what was labelled; may be NULL
 * @return AK_OK, or a negative AK_E* code
 */
int rgfa_build(gfa_t *g, const char *ref_name, rgfa_stat_t *st);

#ifdef __cplusplus
}
#endif

#endif
