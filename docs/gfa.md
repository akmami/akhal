# `gfa` - (r)GFA graph model, reader and traversal

Source: [`src/lib/gfa.c`](../src/lib/gfa.c), [`src/lib/gfa_st.c`](../src/lib/gfa_st.c) &middot; Header: [`include/akhal/gfa.h`](../include/akhal/gfa.h)

The in-memory model of an assembly graph, and the single reader every command
uses. Storage follows the "array + dict" design: segments and links live in
contiguous arrays, and a hash table maps the external segment id to its array
index. Everything that cross-references a segment does so by index.

Edges and paths are both stored in CSR form (a flat array plus per-owner
offsets), so a node's out-edges and a path's ordered segments are each a
contiguous slice.

```c
#include "akhal/gfa.h"
#include "akhal/error.h"   // AK_OK / AK_E* return codes, ak_log()
```

## Contents

- [Reading and releasing](#reading-and-releasing) - [`gfa_read`](#gfa_read), [`gfa_write`](#gfa_write), [`gfa_seg_set_seq`](#gfa_seg_set_seq), [`gfa_destroy`](#gfa_destroy)
- [Lookup and accessors](#lookup-and-accessors) - [`gfa_idx`](#gfa_idx), [`gfa_get`](#gfa_get), [counts and element accessors](#counts-and-element-accessors)
- [Traversal](#traversal) - [`gfa_arcs`](#gfa_arcs), [`gfa_has_arc`](#gfa_has_arc), [`gfa_path_segs`](#gfa_path_segs)
- [Ranks](#ranks) - [`gfa_rank_paths`](#gfa_rank_paths), [`gfa_rank_mark`](#gfa_rank_mark)
- [Rewriting the path block](#rewriting-the-path-block) - [`gfa_clear_paths`](#gfa_clear_paths), [`gfa_add_path`](#gfa_add_path)
- [Ordering](#ordering) - [`gfa_toposort`](#gfa_toposort)
- [Statistics without a graph](#statistics-without-a-graph) - [`gfa_read_stats`](#gfa_read_stats)

## Read flags

`gfa_read()` only does the work you ask for. Pass the bitwise OR of:

| Flag | Effect | Fills |
| --- | --- | --- |
| `GFA_LINKS` | record edges, degrees and out-adjacency | `link`, `arc`, `arc_off`, `in_degree`, `out_degree` |
| `GFA_PATHS` | build path membership and reference layout | `path`, `path_len`, `path_off`, `path_seg`, `path_ori` |
| `GFA_VALIDATE` | check overlap consistency and integrity | nothing; reports through `ak_log()` |

Passing `0` reads segments only. A function that needs a flag says so, and
returns an error rather than misbehaving when it is missing.

`GFA_NIL` (`UINT32_MAX`) marks a path entry whose segment id was not found in
the file. Every loop over a path's segments must skip it.

## Reading and releasing

### `gfa_read`

```c
gfa_t *gfa_read(const char *fn, int flags);
```

Reads an `.gfa` / `.rgfa` file into a freshly allocated graph. Soft problems
found under `GFA_VALIDATE` are logged but do not fail the read; only an
unreadable file or an allocation failure returns `NULL`.

It also settles the ranks. A file's own `SR:i:` tags are authoritative and are
never touched; `g->has_sr` records whether it saw any. Only when it saw none,
*and* `GFA_PATHS` was requested, does it derive them the way
[`gfa_rank_paths`](#gfa_rank_paths) would - so a plain GFA with `P` lines comes
back with a rank-0 backbone, one without them comes back entirely rank 1, and a
real minigraph rGFA comes back exactly as it shipped. Read without `GFA_PATHS`
there is nothing to derive from and ranks stay absent (`-1`).

```c
// Ask for links and paths: the reader skips work you do not request, so a
// command that only needs sequences should pass 0 instead.
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;   // the reason is already logged

printf("%d segments, %d links, %d paths\n",
       gfa_n_seg(g), gfa_n_link(g), gfa_n_path(g));

// has_sr distinguishes ranks the file supplied from ranks the reader derived,
// which is what lets `sort` round-trip an rGFA without inventing SR tags
printf("ranks are %s\n", g->has_sr ? "the file's own" : "derived from the paths");

gfa_destroy(g);
```

### `gfa_write`

```c
int gfa_write(const gfa_t *g, FILE *out);
```

Emits the graph as GFA: an `H` line, one `S` per segment (carrying `SR:i:`
where a rank is set), one `L` per link and one `P` per path. Segment ids are
the graph's own, so this round-trips a graph you have modified in memory.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

// modify in memory, then write the whole graph back out
gfa_rank_paths(g);

FILE *out = fopen("ranked.gfa", "w");
if (!out) {
    gfa_destroy(g);
    return 1;
}
int rc = gfa_write(g, out);   // AK_OK, or AK_EIO if the stream went bad
fclose(out);

gfa_destroy(g);
return rc == AK_OK ? 0 : 1;
```

### `gfa_seg_set_seq`

```c
int gfa_seg_set_seq(gfa_t *g, gfa_seg_t *s, const char *seq, size_t len);
```

Copies `seq` into the graph's arena and points the segment at it, setting
`len` to match. Returns `AK_OK`, `AK_EINVAL` on a `NULL` graph or segment, or
`AK_ENOMEM`. An empty sequence leaves the segment with `seq == NULL` and
`len == 0`, which is what the reader records for an `S` line carrying none.

This is the only supported way to set `gfa_seg_t.seq`. Assigning a `malloc`'d
buffer to the field would be freed twice - once by you and once when the arena
goes - and `free()`ing the field would corrupt the arena. The field's `const`
makes both a compile-time error rather than a crash at teardown.

```c
gfa_t *g = gfa_read("graph.gfa", 0);
if (!g) return 1;

// Hard-mask a segment: same length, all Ns. The old bytes stay in the arena
// until the graph goes, which is a few bytes and not worth reclaiming.
gfa_seg_t *s = gfa_seg_at(g, 0);
char *masked = (char *)malloc(s->len + 1);
if (masked) {
    memset(masked, 'N', s->len);
    masked[s->len] = '\0';
    gfa_seg_set_seq(g, s, masked, s->len);
    free(masked);            // the arena has its own copy now
}

gfa_destroy(g);
```

### `gfa_destroy`

```c
void gfa_destroy(gfa_t *g);
```

Releases a graph and everything it owns - segment sequences, links, the CSR
arrays, path names and the id index. Safe to call with `NULL`, which makes it
usable on every error path.

Segment sequences are not released one at a time: they live in a single
[arena](arena.md) the graph holds, so a graph with 17.9 M segments costs 73
`free()`s rather than 17.9 M. That is also why `gfa_seg_t.seq` is a
`const char *` - it points into storage the graph owns, so it must never be
`free()`d or assigned directly. [`gfa_seg_set_seq`](#gfa_seg_set_seq) is the
only supported way to give a segment its sequence.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_PATHS);
if (!g) return 1;

FILE *out = fopen("out.fa", "w");
if (!out) {
    gfa_destroy(g);   // release before bailing out
    return 1;
}

fclose(out);
gfa_destroy(g);
```

## Lookup and accessors

### `gfa_idx`

```c
int32_t gfa_idx(const gfa_t *g, uint64_t id);
```

Turns a segment id as written in the file into its array index, in O(1)
through the hash table. Returns `-1` when the id is absent. Indices, not ids,
are what the rest of the API takes.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS);
if (!g) return 1;

// Segment ids come from the file and need not be contiguous; array indices
// are 0..n_seg-1 and are what gfa_seg_at() / gfa_arcs() expect.
int32_t i = gfa_idx(g, 42);
if (i < 0) {
    ak_log(AK_LOG_ERROR, NULL, "no segment 42 in the graph");
} else {
    printf("segment 42 is at index %d\n", i);
}

gfa_destroy(g);
```

### `gfa_get`

```c
gfa_seg_t *gfa_get(const gfa_t *g, uint64_t id);
```

The same lookup, but returning the segment itself, or `NULL` when the id is
absent. Use this when you want the node's fields and not its position.

```c
gfa_t *g = gfa_read("graph.gfa", 0);
if (!g) return 1;

gfa_seg_t *s = gfa_get(g, 42);
if (s) {
    // seq is NUL-terminated but may be NULL when the S line had no sequence;
    // len is the cached length, so there is no need to call strlen().
    printf("id %llu, %u bp, rank %d, ref span %d-%d\n",
           (unsigned long long)s->id, s->len, s->rank, s->start, gfa_seg_end(s));
}

gfa_destroy(g);
```

### Counts and element accessors

```c
int32_t     gfa_n_seg(const gfa_t *g);
int32_t     gfa_n_link(const gfa_t *g);
int32_t     gfa_n_path(const gfa_t *g);
gfa_seg_t  *gfa_seg_at(const gfa_t *g, int32_t i);
gfa_link_t *gfa_link_at(const gfa_t *g, int32_t i);
const char *gfa_path_name(const gfa_t *g, int32_t k);
uint64_t    gfa_path_len(const gfa_t *g, int32_t k);
```

Inline, unchecked accessors over the three arrays. They do no bounds checking,
so keep indices inside the matching count. `gfa_path_name()` borrows - the
graph still owns the string.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

// Walk every node. These accessors are inline and unchecked, so the loop
// bound must come from the matching gfa_n_*() call.
uint64_t total = 0;
for (int32_t i = 0; i < gfa_n_seg(g); i++)
    total += gfa_seg_at(g, i)->len;
printf("%llu bp of sequence\n", (unsigned long long)total);

for (int32_t k = 0; k < gfa_n_path(g); k++)
    printf("%s\t%llu bp\n", gfa_path_name(g, k),
           (unsigned long long)gfa_path_len(g, k));

gfa_destroy(g);
```

## Traversal

### `gfa_arcs`

```c
int gfa_arcs(const gfa_t *g, int32_t v, const uint32_t **arcs);
```

Hands back the out-edges of segment index `v` as a slice of the CSR adjacency:
`arcs` points at an array of *link indices*, which you feed to `gfa_link_at()`.
Returns 0 when `v` has no out-edges or the graph was read without `GFA_LINKS`.

```c
// Adjacency only exists under GFA_LINKS; without it every node looks isolated.
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS);
if (!g) return 1;

const uint32_t *arcs;
int na = gfa_arcs(g, 0, &arcs);          // out-edges of segment index 0
for (int i = 0; i < na; i++) {
    // arcs[i] is a link index, not a segment index
    const gfa_link_t *e = gfa_link_at(g, (int32_t)arcs[i]);
    printf("%llu%c -> %llu%c (%u bp overlap)\n",
           (unsigned long long)gfa_seg_at(g, (int32_t)e->v)->id, e->from_orient,
           (unsigned long long)gfa_seg_at(g, (int32_t)e->w)->id, e->to_orient,
           e->overlap);
}

gfa_destroy(g);
```

### `gfa_has_arc`

```c
int gfa_has_arc(const gfa_t *g, int32_t v, int32_t w);
```

Tests for a directed link `v -> w`, in O(out-degree of `v`). Requires
`GFA_LINKS`. Use it for the occasional membership question; scanning
`gfa_arcs()` yourself is better when you need the link's fields.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS);
if (!g) return 1;

// Both arguments are array indices, so translate the file's ids first.
int32_t v = gfa_idx(g, 3), w = gfa_idx(g, 5);
if (v >= 0 && w >= 0 && gfa_has_arc(g, v, w))
    printf("3 -> 5 exists\n");

gfa_destroy(g);
```

### `gfa_path_segs`

```c
int gfa_path_segs(const gfa_t *g, int32_t k, const uint32_t **segs);
```

The ordered segment indices of path `k`. The matching orientation characters
sit at the same offsets in `g->path_ori + g->path_off[k]`. Entries may be
`GFA_NIL` for ids the file referenced but never defined.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_PATHS);
if (!g) return 1;

const uint32_t *segs;
int ns = gfa_path_segs(g, 0, &segs);            // segments of the first path
const char *ori = g->path_ori + g->path_off[0]; // parallel '+'/'-' array

for (int i = 0; i < ns; i++) {
    if (segs[i] == GFA_NIL) continue;   // an id the P line named but no S line defined
    const gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
    printf("%llu%c ", (unsigned long long)s->id, ori[i]);
}
putchar('\n');

gfa_destroy(g);
```

## Ranks

rGFA's `SR:i:` tag says how far a segment sits from the reference: 0 is the
backbone, anything higher came from a sample. A file's own tags are
authoritative - `gfa_read()` records whether it saw any in `g->has_sr` and only
derives ranks when it saw none, so a real minigraph rGFA keeps exactly the
ranks it shipped with while a plain GFA with `P` lines comes back with a rank-0
backbone. These two functions are how you overwrite them deliberately.

### `gfa_rank_paths`

```c
int64_t gfa_rank_paths(gfa_t *g);
```

Ranks against the graph's own paths: every segment any `P` line visits becomes
rank 0, everything else rank 1. A graph with no paths comes back entirely rank
1. Requires `GFA_PATHS`; returns the rank-0 count or a negative `AK_E*` code.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

// re-rank even when the file already carried SR tags: gfa_read() leaves those
// alone, so this is the explicit way to replace them
int64_t n0 = gfa_rank_paths(g);
if (n0 < 0) {
    gfa_destroy(g);
    return 1;
}
printf("%lld backbone node(s), %lld off it\n",
       (long long)n0, (long long)((int64_t)gfa_n_seg(g) - n0));

gfa_destroy(g);
```

### `gfa_rank_mark`

```c
int64_t gfa_rank_mark(gfa_t *g, const uint8_t *on);
```

The general form: rank 0 wherever `on` is set, rank 1 everywhere else. Use it
for a backbone that did not come from the `P` lines - a reference sequence
traced through the graph, for instance, where
[`call_ref_fasta`](call.md#call_ref_fasta) hands you exactly that flag array.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

// mark only the source nodes as backbone, as a stand-in for any labelling
uint8_t *on = (uint8_t *)calloc((size_t)gfa_n_seg(g), 1);
if (!on) {
    gfa_destroy(g);
    return 1;
}
for (int32_t i = 0; i < gfa_n_seg(g); i++) {
    on[i] = gfa_seg_at(g, i)->in_degree == 0;
}

int64_t n0 = gfa_rank_mark(g, on);   // flags array must be n_seg long
printf("%lld node(s) at rank 0\n", (long long)n0);

free(on);
gfa_destroy(g);
```

## Rewriting the path block

### `gfa_clear_paths`

```c
void gfa_clear_paths(gfa_t *g);
```

Drops every path. Path names are owned by the graph but *borrowed* by each
segment's `ref_path` index, so this resets every `ref_path` before freeing them - no
segment is left pointing at a freed name. Segments and ranks are untouched.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

gfa_clear_paths(g);   // also resets every seg->ref_path to -1

printf("%d path(s) left\n", gfa_n_path(g));   // 0
gfa_destroy(g);
```

### `gfa_add_path`

```c
int gfa_add_path(gfa_t *g, const char *name, const uint32_t *segs, const char *ori, int64_t n);
```

Appends one path, laying it out as the reader would: each segment's reference
`start` is recomputed along it and its `ref_path` repointed at the new
name. `ori` may be `NULL` to treat every step as `'+'`, and `GFA_NIL` entries
are skipped. Pair it with `gfa_clear_paths()` to replace the path block
outright - which, with `gfa_rank_mark()`, is how an external reference becomes
the graph's backbone.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_PATHS);
if (!g) return 1;

// keep only the path named "chr22", walked exactly as it was read
for (int32_t k = 0; k < gfa_n_path(g); k++) {
    if (strcmp(gfa_path_name(g, k), "chr22") != 0) continue;

    const uint32_t *segs;
    int ns = gfa_path_segs(g, k, &segs);
    const char *ori = g->path_ori + g->path_off[k];

    // copy both arrays out first: clearing the block frees what they point at
    uint32_t *keep = (uint32_t *)malloc((size_t)ns * sizeof(*keep));
    char     *dir  = (char *)malloc((size_t)ns);
    if (keep && dir) {
        memcpy(keep, segs, (size_t)ns * sizeof(*keep));
        memcpy(dir,  ori,  (size_t)ns);
        gfa_clear_paths(g);
        gfa_add_path(g, "chr22", keep, dir, ns);
    }
    free(keep);
    free(dir);
    break;
}

gfa_destroy(g);
```

## Ordering

### `gfa_toposort`

```c
int gfa_toposort(const gfa_t *g, int32_t *order);
```

Kahn's algorithm over the links. Ties in the ready set are broken by node
*sequence content*, alphabetically, so the result does not depend on the input's
node numbering. Nodes left inside cycles are appended after the acyclic prefix,
so `order` is always a full permutation of `0..n_seg-1`. Requires `GFA_LINKS`.

Returns the number of nodes placed before any cycle (`n_seg` when the graph is
acyclic), or a negative `AK_E*` code.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS);   // required for the in-degrees
if (!g) return 1;

int32_t n = gfa_n_seg(g);
int32_t *order = (int32_t *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int32_t));
if (!order) { gfa_destroy(g); return 1; }

int32_t placed = gfa_toposort(g, order);       // order[] is caller-allocated
if (placed < 0)      ak_log(AK_LOG_ERROR, NULL, "%s", ak_strerror(placed));
else if (placed < n) ak_log(AK_LOG_WARN, NULL, "%d node(s) sit in cycles", n - placed);
else                 printf("acyclic; order[0] is index %d\n", order[0]);

free(order);
gfa_destroy(g);
```

---

[Back to the library index](README.md)

## Statistics without a graph

Summarizing a file does not need one. 
The counts are running totals and the distributions are Welford's online mean and variance ([`ak_dist_t`](util.md#ak_dist_t)), so the streaming part of `gfa_read_stats()` holds nothing but the line buffer. 
Three figures are properties of the file as a whole rather than of any line - a segment's degree is how many `L` lines name it, "undefined" is named but never defined, and without `SR` tags "rank 0" is named by some `P` line - and those are answered from flat **id arrays** rather than tables: the ids are appended to flat arrays as the lines stream past, radix-sorted in place once the file is read (klib's `KRADIX_SORT_INIT`, which needs no scratch buffer), and read off as run lengths and merges. 
An array costs 4 bytes per entry, or 8 once an id exceeds 32 bits, whatever the numbering looks like, so scattered ids cost exactly what dense ones do and there is no fallback path.

### `gfa_read_stats`

```c
#define GFA_STAT_DEGREES 0x1   // degree distributions and n_undefined
#define GFA_STAT_RANKS   0x2   // n_rank0 from the P lines when the file has no SR tags

int gfa_read_stats(const char *fn, gfa_stat_t *st, int flags);
```

Fills `st` from one pass over `fn`. 
Returns `AK_OK`, or a negative `AK_E*` code with the reason logged.

| Field | Type | Notes |
| --- | --- | --- |
| `n_seg`, `n_link`, `n_path` | `int64_t` | lines of each kind |
| `n_bp` | `uint64_t` | total sequence length, the sum of every segment's length |
| `n_rank0` | `int64_t` | segments at rank 0; `-1` when not derived |
| `has_sr` | `int` | whether the file carried its own `SR` tags |
| `seg_mean`, `seg_sd`, `seg_min`, `seg_max` | `double`, `uint64_t` | segment length |
| `ov_mean`, `ov_sd` | `double` | link overlap |
| `in_mean`, `in_sd`, `out_mean`, `out_sd` | `double` | degrees, over segments that have any |
| `min_in`, `max_in`, `min_out`, `max_out` | `int32_t` | degree extremes, `-1` when no segment has one |
| `n_undefined` | `int64_t` | distinct ids an `L` or `P` line names that no `S` line defines; `-1` when not checked |

What each flag costs. 
`GFA_STAT_DEGREES` keeps three arrays - the `S` ids, the `L` sources and the `L` targets - about `4 * (n_seg + 2 * n_link)` bytes, so 8 GB on a graph of 700 million segments and links, and it is what fills the degree fields and `n_undefined`. 
`GFA_STAT_RANKS` keeps the `S` ids and every `P` step, which on a file carrying all its haplotypes as `P` lines is far larger than the graph itself; a file that ranks itself with `SR` tags needs nothing and answers `n_rank0` with either setting, since those tags are authoritative and the `P` lines are then not consulted. 
With `flags == 0` the pass holds nothing per line, and the fields the flags would fill are `-1`.

Three things are worth knowing about the numbers. 
The standard deviations are population figures, dividing by n. 
The degree figures cover only segments with a non-zero degree, so a graph of unlinked segments reports `-1` for the four extremes rather than `0`; an `L` line naming an undefined id still counts toward its defined endpoint's degree, as it does toward `n_link`. 
And `n_rank0` without `SR` tags counts defined segments only, the same set [`gfa_rank_paths`](#gfa_rank_paths) can mark.

```c
gfa_stat_t st;
if (gfa_read_stats("graph.gfa", &st, GFA_STAT_DEGREES) != AK_OK) return 1;

printf("%lld segment(s), %lld link(s)\n", (long long)st.n_seg, (long long)st.n_link);
printf("segment length: mean %.2f, sd %.2f, %llu..%llu\n",
       st.seg_mean, st.seg_sd,
       (unsigned long long)st.seg_min, (unsigned long long)st.seg_max);

// -1 means no segment has a degree at all, which is not the same as 0
if (st.max_in >= 0)
    printf("in degree %d..%d, out degree %d..%d\n",
           st.min_in, st.max_in, st.min_out, st.max_out);

if (st.n_undefined)
    printf("warning: %lld id(s) are named but never defined\n", (long long)st.n_undefined);
```
