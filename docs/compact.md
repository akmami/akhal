# `compact` - folding non-branching runs into single nodes

Source: [`src/lib/compact.c`](../src/lib/compact.c) &middot; Header: [`include/akhal/compact.h`](../include/akhal/compact.h)

Builders leave graphs full of chains: node after node joined by a single link,
with nothing branching anywhere along the way. Nothing is expressed by keeping
them apart, so a run like that folds into one segment carrying the same bases -
"unitig" in de Bruijn terms, and what `vg mod -u` and `odgi unchop` do
elsewhere.

The module finds the runs and writes the folded graph; it never modifies the
graph it was given, which is what lets the caller decide whether to keep the
original around.

This module backs [`akhal compact`](../README.md#5-compact).

```c
#include "akhal/compact.h"   // compact_t, compact_runs(), compact_write()
#include "akhal/gfa.h"       // gfa_read() and the GFA_* flags (compact.h pulls this in)
#include "akhal/error.h"     // AK_OK / AK_E* return codes, ak_log()
```

## Contents

- [When a join is taken](#when-a-join-is-taken)
- [The runs](#the-runs) - [`compact_runs`](#compact_runs), [`compact_destroy`](#compact_destroy)
- [Writing it out](#writing-it-out) - [`compact_write`](#compact_write)

## When a join is taken

A join `u -> v` is taken when all of this holds.

| | |
| --- | --- |
| **Same strand** | the link leaves `u`'s right end and enters `v`'s left end - `L u + v +`, or the same edge read from the other side, `L v - u -`. A join that flips strand is left alone |
| **Nothing branches** | it is the only link on either of those two ends |
| **Blunt** | the link carries a 0 bp overlap. Bases the two share would have to be dropped, and this never rewrites sequence it was not given |
| **Not itself** | `u` is not `v`, so a self-loop stays a self-loop |
| **No path stops there** | no path starts or ends at either end. A `P` line cannot stop halfway through a node, so a run is cut wherever one would have to |
| **Tags still mean something** | on an rGFA, `u` and `v` sit on the same stable sequence at the same rank with contiguous offsets. A plain GFA has nothing to check |

One caution on that last rule. The reader does not keep a file's `SN` at all, and
a graph read with `GFA_PATHS` has its `SO` overwritten by the layout of the path
that owns each segment. So the `SN`/`SO` this module compares - and the ones it
writes - are the owning path and the offset along it, not what the input file
said. The output is a consistent rGFA, but away from the backbone its offsets
are the path's own rather than the split-point coordinates
[`rgfa_build`](rgfa.md#rgfa_build) works out. **Compact first, label afterwards**
is the order that leaves both saying what they mean.

Ends rather than nodes are what the first two rules count, which is what makes
the two spellings of one edge the same thing. A link records the end it uses at
each of its segments: leaving `v` it uses `v`'s right end on `+` and its left on
`-`; arriving at `w` it uses `w`'s left end on `+` and its right on `-`.

Because an interior node of a run has exactly one link on each end - the two
joins that put it there - no other link can touch it. Every link that survives
therefore lands on a run's first or last segment, and lands on the same side of
it as before, which is why the orientations on the `L` lines come through
unchanged.

## The runs

```c
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
```

A run is a maximal chain of joins, and a segment nothing could be merged with is
a run of one - so every segment is in exactly one run and `n_run` is the segment
count after compaction. Runs are numbered in the order their first segment
appears, and `first` + `next` walk one from end to end.

A cycle of joins has no first segment at all. It is cut at its lowest segment
index, so a circular contig comes out as one node carrying a self-loop rather
than being left alone.

### `compact_runs`

```c
compact_t *compact_runs(const gfa_t *g);
```

Works out the runs and returns them, or `NULL` with the reason logged.
`GFA_LINKS` is required. `GFA_PATHS` is not, but without it the paths are
invisible, and a run will happily swallow a segment some `P` line stops at -
read with both whenever the graph has paths.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

compact_t *c = compact_runs(g);
if (!c) {
    gfa_destroy(g);
    return 1;
}

printf("%d segment(s) fold into %d\n", c->n_seg, c->n_run);

// what each run gathered, in order
for (int32_t r = 0; r < c->n_run; r++) {
    if (c->count[r] < 2) continue;
    printf("run %d:", r);
    for (uint32_t s = c->first[r]; s != GFA_NIL; s = c->next[s])
        printf(" %llu", (unsigned long long)gfa_seg_at(g, (int32_t)s)->id);
    printf("\n");
}

compact_destroy(c);
gfa_destroy(g);
```

### `compact_destroy`

```c
void compact_destroy(compact_t *c);
```

Releases the run set and its five arrays. Safe to call with `NULL`. The graph is
untouched throughout - nothing in a `compact_t` points into it.

## Writing it out

### `compact_write`

```c
int compact_write(const gfa_t *g, const compact_t *c, FILE *out);
```

Writes the folded graph: one `S` line per run carrying the concatenated
sequence, then the links and paths that survive. Returns `AK_OK`, or a negative
`AK_E*` code.

Each run keeps the id of its first segment, so a node that merged nothing comes
out exactly as it went in. rGFA tags are emitted only when the input carried its
own `SR` tags - the same rule [`sort`](../README.md#6-sort) follows, so a plain
GFA is not handed ranks the reader merely derived for it.

A `P` line that steps into the middle of a run from somewhere else did not
follow the links, and collapsing it would rewrite what it spells. That is
reported at `WARN` and the step is written out on its own rather than folded,
which is the least the file can be turned into; [`akhal parse`](../README.md#1-parse)
is what reports the problem properly. A run walked twice round - a cycle
returning to where it began - is not that case, and comes out as two steps.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

compact_t *c = compact_runs(g);
int rc = c ? compact_write(g, c, stdout) : AK_ENOMEM;

compact_destroy(c);
gfa_destroy(g);
if (rc != AK_OK) return 1;
```

---

[Back to the library index](README.md)
