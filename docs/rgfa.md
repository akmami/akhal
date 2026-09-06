# `rgfa` - stable-sequence labelling, a GFA turned into an rGFA

Source: [`src/lib/rgfa.c`](../src/lib/rgfa.c) &middot; Header: [`include/akhal/rgfa.h`](../include/akhal/rgfa.h)

rGFA says where each segment sits on a real sequence, through three tags:

| Tag  | Type | Says                                                        |
| ---- | ---- | ----------------------------------------------------------- |
| `SN` | `Z`  | name of the stable sequence the segment came from           |
| `SO` | `i`  | offset on that sequence                                     |
| `SR` | `i`  | rank: 0 on the linear reference, higher the further from it |

A plain GFA carries none of them. Its `P` lines, though, say everything needed
to work them out, and that is the whole of this module: one call that reads the
paths and writes the three tags onto every segment they explain.

This module backs [`akhal gfa2rgfa`](../README.md#8-gfa2rgfa).

```c
#include "akhal/rgfa.h"    // rgfa_stat_t, rgfa_build()
#include "akhal/gfa.h"     // gfa_read(), gfa_write_rgfa() (rgfa.h pulls this in)
#include "akhal/error.h"   // AK_OK / AK_E* return codes, ak_log()
```

## Contents

- [What the walk does](#what-the-walk-does) - [`rgfa_build`](#rgfa_build)
- [Where the tags live](#where-the-tags-live)

## What the walk does

One path is declared the backbone. Its segments are rank 0, named after it,
with offsets running the length of the walk - a segment the backbone comes back
to keeps the offset of its first visit.

Every other path is then walked in turn. While it runs over ground some earlier
walk already labelled it only follows along, keeping track of where it is. Where
it leaves that ground, the segments it visits alone are **one rank deeper**,
named after that path, and offset **onward from the point it left** - so a
bubble carries the coordinate of the reference stretch it detours around, and
the counting stops as soon as the path merges back down to a lower rank.

```text
ref   1(AAAA) --- 2(C) --- 4(TTTT)                 rank 0, SN=ref, SO 0, 4, 5
samp  1(AAAA) --- 3(G) --- 4(TTTT) --- 5(GG)
                   `-- rank 1, SN=samp, SO=4        the offset node 2 sits at
                                          `-- rank 1, SN=samp, SO=9   past the end
```

Rank is depth, not path order: a detour hanging off a rank-1 stretch is rank 2.
Ranks and offsets therefore describe a position in the graph rather than which
path happened to be walked first - though *which* path lends its name to a
segment several could explain is decided by walk order, which is the order
[`gfa_path_merge`](gfa.md#gfa_path_merge) returns chains in.

A path that begins away from the backbone has no split point to count from, so
it starts at 0 on its own name, which is a coordinate on that path and nothing
else.

**Where two paths disagree, nothing is invented.** A segment one path places at
one offset and another at a different one has no single answer: it keeps its
rank and loses `SN` and `SO`, and so does everything that walk reaches
afterwards - until it arrives somewhere an authoritative path has already
pinned down. The backbone is authoritative always; a rank-0 segment is never
unplaced by a later disagreement, which is what keeps one tangled sample from
poisoning the reference coordinates.

Segments no path visits at all are left with no tags whatsoever, rank included.
Both they and the ambiguous ones are counted, so the caller can report them.

```c
typedef struct {
    int32_t n_path;        // paths after fragments were consolidated
    int32_t n_rank0;       // segments on the backbone
    int32_t n_labelled;    // segments that came out with SN and SO
    int32_t n_ambiguous;   // segments a walk reached but could not place
    int32_t n_unreached;   // segments no path visits at all
    int32_t max_rank;      // deepest rank handed out
} rgfa_stat_t;
```

This all assumes a **variation graph**, whose paths are haplotypes over a shared
backbone. An assembly graph, whose paths are unrelated contigs, will label
badly - there is no meaningful backbone for the rest to be an offset from, and
almost everything comes out either rank 1 from an arbitrary starting path or
ambiguous.

### `rgfa_build`

```c
int rgfa_build(gfa_t *g, const char *ref_name, rgfa_stat_t *st);
```

Labels the graph in place and returns `AK_OK`, or a negative `AK_E*` code with
the reason logged. Requires `GFA_LINKS | GFA_PATHS`, and a graph with no `P`
lines is refused - there is nothing to label against.

Fragmented `P` lines are consolidated first, exactly as
[`gfa_path_merge`](gfa.md#gfa_path_merge) chains them, so the graph comes back
with one `P` line per chain and one stable sequence name per reference. That
rewrite happens whether or not the labelling finds anything to say.

`ref_name` picks the backbone by chain name, or by the name of any fragment
that went into a chain, so either `chr22` or `chr22[0]` selects the same one.
`NULL` takes the chain holding the graph's first `P` line.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

rgfa_stat_t st;
if (rgfa_build(g, "chr22", &st) != AK_OK) {
    gfa_destroy(g);
    return 1;
}

printf("%d on the backbone, %d placed, %d ambiguous, up to rank %d\n",
       st.n_rank0, st.n_labelled, st.n_ambiguous, st.max_rank);

// the tags only reach a file through the rGFA writer
gfa_write_rgfa(g, stdout);
gfa_destroy(g);
```

## Where the tags live

There is no new storage for any of this: the three tags are three fields
[`gfa_seg_t`](gfa.md#segments) already had, which is what lets the ordinary
graph writer emit them.

| Tag | Field | Unset means |
| --- | --- | --- |
| `SN` | `ref_name` | `NULL` - borrowed from the path block, so it dies with the paths |
| `SO` | `start` (and `end`) | `-1` |
| `SR` | `rank` | `-1` |

[`gfa_write_rgfa`](gfa.md#gfa_write_rgfa) emits each tag only where the segment
carries it, so an ambiguous segment comes out with `SR` alone and an unvisited
one with a bare `S` line. [`gfa_write`](gfa.md#gfa_write) ignores `SN` and `SO`
entirely and emits plain GFA, which is why `sort` and `rank` are unaffected by
any of this.

Because `ref_name` borrows from the path block, the labels last exactly as long
as the paths do: [`gfa_clear_paths`](gfa.md#gfa_clear_paths) drops every `SN`,
and [`gfa_add_path`](gfa.md#gfa_add_path) relabels the segments it lays out.
Write the graph before rewriting its paths, not after.

---

[Back to the library index](README.md)
