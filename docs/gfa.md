# `gfa` - (r)GFA graph model, reader and traversal

Source: [`src/lib/gfa.c`](../src/lib/gfa.c) &middot; Header: [`include/akhal/gfa.h`](../include/akhal/gfa.h)

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

- [Reading and releasing](#reading-and-releasing) - [`gfa_read`](#gfa_read), [`gfa_write`](#gfa_write), [`gfa_destroy`](#gfa_destroy)
- [Lookup and accessors](#lookup-and-accessors) - [`gfa_idx`](#gfa_idx), [`gfa_get`](#gfa_get), [counts and element accessors](#counts-and-element-accessors)
- [Traversal](#traversal) - [`gfa_arcs`](#gfa_arcs), [`gfa_has_arc`](#gfa_has_arc), [`gfa_path_segs`](#gfa_path_segs)
- [Fragmented paths](#fragmented-paths) - [`gfa_path_merge`](#gfa_path_merge), [`gfa_merge_segs`](#gfa_merge_segs), [`gfa_merge_destroy`](#gfa_merge_destroy)
- [Ranks](#ranks) - [`gfa_rank_paths`](#gfa_rank_paths), [`gfa_rank_mark`](#gfa_rank_mark)
- [Rewriting the path block](#rewriting-the-path-block) - [`gfa_clear_paths`](#gfa_clear_paths), [`gfa_add_path`](#gfa_add_path)
- [Ordering](#ordering) - [`gfa_toposort`](#gfa_toposort)

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

### `gfa_destroy`

```c
void gfa_destroy(gfa_t *g);
```

Releases a graph and everything it owns - segment sequences, links, the CSR
arrays, path names and the id index. Safe to call with `NULL`, which makes it
usable on every error path.

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
           (unsigned long long)s->id, s->len, s->rank, s->start, s->end);
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

## Fragmented paths

Builders such as `vg` split one reference across consecutive `P` lines
(`chr22[0]`, `chr22[21]`, ...). These three functions stitch them back
together. `gfa_merge_t` is CSR like everything else: chain `c` owns the path
indices `frag[off[c] .. off[c+1])`.

### `gfa_path_merge`

```c
gfa_merge_t *gfa_path_merge(const gfa_t *g, const char *key);
```

Selects fragments by name - vg's `name[start]` and `name:start-end` decoration
is stripped, and a PanSN name like `GRCh38#0#chr22[0]` is also found by its
bare contig name - then chains them through the links. `key` of `NULL` selects
every path and groups each base name separately. Requires
`GFA_LINKS | GFA_PATHS`.

```c
// Both flags are required: names select the fragments, links order them.
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

gfa_merge_t *m = gfa_path_merge(g, "chr22");   // NULL would group every path
if (!m) { gfa_destroy(g); return 1; }

for (int32_t c = 0; c < m->n; c++)
    printf("%s\tfrom %d P line(s)\n", m->name[c], m->off[c + 1] - m->off[c]);

gfa_merge_destroy(m);
gfa_destroy(g);
```

### `gfa_merge_segs`

```c
int64_t gfa_merge_segs(const gfa_t *g, const gfa_merge_t *m, int32_t c, uint32_t **segs);
```

Flattens chain `c` into one array of segment indices, dropping `GFA_NIL`
entries as it goes. **The caller frees `*segs`.** Returns the count, or a
negative `AK_E*` code.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
gfa_merge_t *m = gfa_path_merge(g, "chr22");
if (!m) { gfa_destroy(g); return 1; }

uint32_t *segs;
int64_t ns = gfa_merge_segs(g, m, 0, &segs);   // chain 0, already GFA_NIL-free
if (ns > 0) {
    for (int64_t i = 0; i < ns; i++) {
        const gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
        if (s->seq) fwrite(s->seq, 1, s->len, stdout);
    }
    free(segs);   // the array is yours; the graph still owns the sequences
}

gfa_merge_destroy(m);
gfa_destroy(g);
```

### `gfa_merge_destroy`

```c
void gfa_merge_destroy(gfa_merge_t *m);
```

Releases a chain set and the names it owns. Safe with `NULL`. It does not touch
the graph, so destroy order between the two does not matter.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

gfa_merge_t *m = gfa_path_merge(g, NULL);
if (m) {
    printf("%d chain(s)\n", m->n);
    gfa_merge_destroy(m);   // independent of the graph's lifetime
}

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
segment's `ref_name`, so this clears `ref_name` before freeing them - no
segment is left pointing at a freed name. Segments and ranks are untouched.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

gfa_clear_paths(g);   // also resets every seg->ref_name to NULL

printf("%d path(s) left\n", gfa_n_path(g));   // 0
gfa_destroy(g);
```

### `gfa_add_path`

```c
int gfa_add_path(gfa_t *g, const char *name, const uint32_t *segs, const char *ori, int64_t n);
```

Appends one path, laying it out as the reader would: each segment's reference
`start`/`end` is recomputed along it and its `ref_name` repointed at the new
name. `ori` may be `NULL` to treat every step as `'+'`, and `GFA_NIL` entries
are skipped. Pair it with `gfa_clear_paths()` to replace the path block
outright - which, with `gfa_rank_mark()`, is how an external reference becomes
the graph's backbone.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

// consolidate a vg-fragmented reference into a single P line
gfa_merge_t *m = gfa_path_merge(g, "chr22");
if (!m) {
    gfa_destroy(g);
    return 1;
}
uint32_t *segs;
int64_t ns = gfa_merge_segs(g, m, 0, &segs);   // flatten before clearing:
if (ns > 0) {                                  // it reads the very paths
    gfa_clear_paths(g);                        // clearing would free
    gfa_add_path(g, m->name[0], segs, NULL, ns);
    free(segs);
}

gfa_merge_destroy(m);
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
