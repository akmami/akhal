# `diff` - comparing two graphs that disagree on ids

Source: [`src/lib/diff.c`](../src/lib/diff.c) &middot; Header: [`include/akhal/diff.h`](../include/akhal/diff.h)

Two files can spell the same graph and share not a single segment id: builders
number nodes as they emit them, and [`akhal sort`](../README.md#5-sort)
renumbers them again. So nothing here looks at an id. Segments are matched on
their **sequence content**, both graphs are relabelled onto one shared
numbering, and links are compared through those labels - which turns "is this
edge in the other graph?" into an integer comparison.

Paths are the exception: they are compared by what they spell rather than by
what they walk over, so two graphs that chop one reference into different nodes
still agree on it.

This module backs [`akhal compare`](../README.md#4-compare), which is the whole
of it: two graphs in, one verdict out.

```c
#include "akhal/diff.h"    // diff_map_t, diff_t and the diff_* API
#include "akhal/gfa.h"     // gfa_read() and the GFA_* flags (diff.h pulls this in)
#include "akhal/error.h"   // AK_OK / AK_E* return codes, ak_log()
```

## Contents

- [Matching by content](#matching-by-content) - [`diff_map`](#diff_map), [`diff_map_destroy`](#diff_map_destroy)
- [The comparison](#the-comparison) - [`diff_graphs`](#diff_graphs), [`diff_destroy`](#diff_destroy), [`diff_identical`](#diff_identical)

## Matching by content

The labelling is a shared numbering over both graphs' segments.

```c
typedef struct {
    uint32_t *a;         // length n_a: segment index in A -> shared label
    uint32_t *b;         // length n_b: segment index in B -> shared label
    int32_t   n_a, n_b;  // segment counts of the two graphs
    int32_t   n_class;   // distinct sequences; labels run 0 .. n_class - 1
    int32_t   n_shared;  // segments that pair off, summed over the classes
} diff_map_t;
```

Both graphs draw labels from one counter, so **equal labels mean equal
sequences** - within one graph as much as across the two. There is no sentinel:
every segment gets a label, matched or not.

Matching is a merge pass over both segment arrays sorted by sequence, so it
costs one sort per graph and a single linear walk, not a hash lookup per node.
Sorting is by content, the same thing [`gfa_toposort`](gfa.md#gfa_toposort)
breaks its ties with.

**A label is a class of equal sequences, not one segment.** Variation graphs are
full of single-base nodes, so `A` is carried dozens of times over, and which
copy in one file "is" which copy in the other is not something the two files
agree on. Pairing them off one by one would be a coin toss, and every link
touching them would inherit it - which is how a comparison ends up reporting
hundreds of spurious link differences between two graphs that are the same. So
a whole run of equal sequences takes one label instead.

Counting is still multiset-style, through `n_shared`: with three copies of a
sequence in one graph and two in the other, two match and the third is
unmatched. What is given up is only *which* copy is the odd one out, which was
never meaningful.

A segment with no sequence at all is treated as carrying the empty one, so two
such segments match each other - and, being equal sequences, they all share one
label as well.

### `diff_map`

```c
diff_map_t *diff_map(const gfa_t *a, const gfa_t *b);
```

Labels both graphs' segments and returns the mapping, or `NULL` on allocation
failure. Needs no read flags: only the `S` lines are used, so a graph read with
flags of 0 is enough.

```c
// Ids are never compared, so the flags can be 0 - only sequences matter here.
gfa_t *a = gfa_read("one.gfa", 0);
gfa_t *b = gfa_read("two.gfa", 0);
if (!a || !b) return 1;

diff_map_t *m = diff_map(a, b);
if (!m) return 1;

printf("%d matched, %d only in A, %d only in B\n",
       m->n_shared, m->n_a - m->n_shared, m->n_b - m->n_shared);

// Counting B's segments per class turns the labelling into a lookup: does the
// other graph spell this, and how many times?
int32_t *in_b = (int32_t *)calloc((size_t)m->n_class, sizeof(int32_t));
if (!in_b) return 1;
for (int32_t j = 0; j < m->n_b; j++) in_b[m->b[j]]++;

for (int32_t i = 0; i < m->n_a; i++)
    if (in_b[m->a[i]] == 0)
        printf("B spells nothing like segment %llu\n",
               (unsigned long long)gfa_seg_at(a, i)->id);

free(in_b);
diff_map_destroy(m);
gfa_destroy(a);
gfa_destroy(b);
```

### `diff_map_destroy`

```c
void diff_map_destroy(diff_map_t *m);
```

Releases the labelling and both arrays. Safe to call with `NULL`. The graphs it
was built from are untouched - nothing in a `diff_map_t` points into them.

## The comparison

[`diff_graphs`](#diff_graphs) does the whole job and reports it as one object:
what matched, and what each graph holds alone.

```c
typedef struct {
    uint64_t    *seg;      // owned: ids of segments only this graph has
    int32_t      n_seg;
    diff_link_t *link;     // owned: links only this graph has
    int32_t      n_link;
} diff_side_t;

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
```

Shared counts are **pair** counts: `n_seg_shared` segments matched means that
many on each side, so `gfa_n_seg(a) == d->n_seg_shared + d->a.n_seg` and the
same holds for B.

Unmatched links are reported the way their own file spells them, with that
graph's ids rather than the internal labels.

```c
typedef struct {
    uint64_t from, to;       // segment ids as that graph spells them
    char     from_orient;    // orientation of `from` on the L line ('+' or '-')
    char     to_orient;      // orientation of `to` on the L line ('+' or '-')
    uint32_t overlap;        // overlap length in bp
} diff_link_t;
```

Two links are the same when their relabelled endpoints, both orientations and
the overlap all agree. The comparison is canonical over the two spellings of
one edge, so `L a + b +` in one file and `L b - a -` in the other - the same
edge read from the other end - count as one link, not two differences.

Because a label is a class of equal sequences, a link between repeated
sequences asks whether the other graph has that edge *between those sequences*,
not between those particular nodes. That is the deliberate trade: a graph whose
`A` nodes are numbered differently compares equal, at the cost of not
distinguishing which `A` an edge landed on.

Paths get one entry per name, whether one graph has it or both.

```c
typedef struct {
    char    *name;     // owned: chain name, as gfa_path_merge() named it
    int      state;    // DIFF_SAME / DIFF_DIFFER / DIFF_A_ONLY / DIFF_B_ONLY
    uint64_t len_a;    // bases the first graph spells for it, 0 when absent
    uint64_t len_b;    // bases the second graph spells for it, 0 when absent
} diff_path_t;
```

Fragmented `P` lines are chained with [`gfa_path_merge`](gfa.md#gfa_path_merge)
before anything is compared, so a reference that arrived as `chr22:0-1000`,
`chr22:1000-2000`, ... is one entry named `chr22` on both sides. Chains are then
paired by that name and each pair's **sequences** compared, a `-` step
contributing its reverse complement. Lengths are compared first and the bases
only spelled when they match, so a path that differs in length costs nothing to
reject.

Two caveats come with pairing by name. The name is the one
[`gfa_path_merge`](gfa.md#gfa_path_merge) settled on, which depends on how many
chains a base name yielded - so a reference that chains into one piece here
(`chr1`) but two there (`chr1_1`, `chr1_2`, because a joining `L` line is
missing) pairs with nothing, and reports as three unmatched names rather than
one difference. And the comparison is byte-for-byte, so it is case-sensitive;
`ak_revcomp` uppercases as it complements, which means a soft-masked graph
compared against a case-normalized one can differ on case alone.

### `diff_graphs`

```c
diff_t *diff_graphs(const gfa_t *a, const gfa_t *b);
```

Compares two graphs, both of which must have been read with `GFA_LINKS |
GFA_PATHS`; anything less is refused through `ak_log()` with a `NULL` return,
as is an allocation failure along the way. A graph that carries no `P` lines is
not a failure - it contributes no chains, and its segments and links compare as
usual - so a plain GFA is a perfectly good input.

Path sequences are spelled a pair at a time rather than all at once, so peak
memory is the two longest chains that share a name, not every path in both
files. Link overlaps are not trimmed off those bases - the same blunt-join
assumption [`extract path`](../README.md#path) makes.

```c
gfa_t *a = gfa_read("one.gfa", GFA_LINKS | GFA_PATHS);
gfa_t *b = gfa_read("two.gfa", GFA_LINKS | GFA_PATHS);
if (!a || !b) return 1;

diff_t *d = diff_graphs(a, b);
if (!d) return 1;

printf("%d segments and %d links in common\n", d->n_seg_shared, d->n_link_shared);

// Each side lists what only it carries, in its own ids.
for (int32_t i = 0; i < d->a.n_seg; i++)
    printf("only in A: segment %llu\n", (unsigned long long)d->a.seg[i]);
for (int32_t i = 0; i < d->b.n_link; i++)
    printf("only in B: link %llu%c -> %llu%c\n",
           (unsigned long long)d->b.link[i].from, d->b.link[i].from_orient,
           (unsigned long long)d->b.link[i].to,   d->b.link[i].to_orient);

// A path entry covers both graphs, so one pass reports every verdict.
for (int32_t i = 0; i < d->n_path; i++)
    if (d->path[i].state == DIFF_DIFFER)
        printf("%s: %llu bp here, %llu bp there\n", d->path[i].name,
               (unsigned long long)d->path[i].len_a,
               (unsigned long long)d->path[i].len_b);

diff_destroy(d);
gfa_destroy(a);
gfa_destroy(b);
```

### `diff_destroy`

```c
void diff_destroy(diff_t *d);
```

Releases the comparison, both sides' arrays and every path name. Safe to call
with `NULL`.

### `diff_identical`

```c
static inline int diff_identical(const diff_t *d);
```

Returns 1 when nothing is unmatched on either side and no path name differs.
This is the whole verdict in one call, and what `akhal compare` turns into its
exit status. `d` must not be NULL - unlike the destroys, this one does not
check.

It is a statement about structure and content, not about the files: two graphs
that pass this can still have different ids, different node numbering, a
different line order, and different `SR` ranks - none of which is compared.

```c
gfa_t *a = gfa_read("one.gfa", GFA_LINKS | GFA_PATHS);
gfa_t *b = gfa_read("two.gfa", GFA_LINKS | GFA_PATHS);
if (!a || !b) return 1;

diff_t *d = diff_graphs(a, b);
if (!d) return 1;

// True for a graph and its own `akhal sort` output, which renumbers every node.
if (diff_identical(d)) printf("same graph, different numbering\n");

diff_destroy(d);
gfa_destroy(a);
gfa_destroy(b);
```

---

[Back to the library index](README.md)
