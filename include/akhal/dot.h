#ifndef AKHAL_DOT_H
#define AKHAL_DOT_H

#include <stdio.h>
#include "akhal/gfa.h"

/**
 * Graphviz output, for looking at a graph rather than computing on it.
 *
 * One DOT node per segment and one arrow per link. GFA is bidirected and DOT is
 * not, so the strands a link joins are written on the arrow - `+/-` where it
 * flips - and left off where both ends are forward, which is nearly always.
 * Nothing about the graph is lost, but the picture is a directed reading of it.
 *
 * Meant for graphs small enough to look at: a bubble, a test case, a region
 * pulled out of something bigger. Graphviz will take a graph of any size and
 * spend a very long time on it, so dot_write() says so when handed one.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Sequences up to this many bases are printed on the node; longer segments
// show their length instead
#define DOT_SEQ_MAX 20

// Node count past which laying the graph out stops being worth waiting for
#define DOT_BIG 10000

// Flags for dot_write()
#define DOT_IDS 0x1     // label nodes with ids alone, however short the bases

/**
 * Write a graph as a Graphviz digraph.
 *
 * Nodes carry their id and, when the segment is at most DOT_SEQ_MAX bases, the
 * bases themselves; anything longer shows a length. On a file that ranks itself
 * the backbone is shaded apart from what hangs off it, so an rGFA's structure
 * reads at a glance. A link carries a label when it flips strand or when it has
 * an overlap to declare, and none at all otherwise
 * @param g Graph to draw; no read flags beyond the segments are required,
 *          though without GFA_LINKS there is nothing to draw between them
 * @param out Destination stream
 * @param flags DOT_IDS, or 0
 * @return AK_OK, or a negative AK_E* code
 */
int dot_write(const gfa_t *g, FILE *out, int flags);

#ifdef __cplusplus
}
#endif

#endif
