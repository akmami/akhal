# `vg` - reader for vg's native `.vg` format

Source: [`src/lib/vg.c`](../src/lib/vg.c) &middot; Header: [`include/akhal/vg.h`](../include/akhal/vg.h)

A `.vg` file is a gzip/BGZF-compressed stream of length-delimited Protobuf
`Graph` messages. On the wire it is a repeating pattern of
`[varint count][count x (varint length + message bytes)]`, and each `Graph`
carries repeated `Node`, `Edge` and `Path` submessages. There is no index and no
global header: the file is a log of graph chunks, and reading it means
concatenating every chunk into one graph. That is what `vg view -g` does, and
what this module does.

Rather than depend on protobuf or libvgio, the module decodes the handful of
message types needed for GFA conversion directly from the wire format, with
field numbers taken from `vg.proto`: `Node { sequence = 1, id = 3 }`,
`Edge { from = 1, to = 2, from_start = 3, to_end = 4, overlap = 5 }`,
`Path { name = 1, mapping = 2, is_circular = 3 }` and, inside a `Mapping`,
`Position { node_id = 1, is_reverse = 4 }`. Every other field is skipped by wire
type. The only external dependency is zlib, for decompression.

Decoding is deliberately lenient. A tag it cannot parse, an unknown wire type or
a truncated message ends the scan rather than failing it, so a damaged file
yields the prefix that did decode. Messages larger than 1 GiB are refused
outright, which stops a non-`.vg` input from turning a garbage length into a
huge allocation.

```c
#include "akhal/vg.h"      // the graph structs, vg_read() and vg_graph_destroy()
#include "akhal/error.h"   // ak_log() for reporting what a lenient read produced
```

## Contents

- [The graph model](#the-graph-model) - `vg_node_t`, `vg_edge_t`, `vg_step_t`, `vg_path_t`, `vg_graph_t`
- [Reading and releasing](#reading-and-releasing) - [`vg_read`](#vg_read), [`vg_graph_destroy`](#vg_graph_destroy)
- [Converting to GFA](#converting-to-gfa) - walking nodes, edges and paths, and the edge orientation convention

## The graph model

Five plain structs, all of them public and all of them owned by the graph. There
are no accessor functions: read the fields directly.

`vg_node_t` - one sequence node.

| Field | Type | Notes |
| --- | --- | --- |
| `id` | `int64_t` | node id, positive and nonzero; `0` if the message had no id field |
| `seq` | `char *` | owned, NUL-terminated sequence, or `NULL` when the message carried none |
| `seq_len` | `uint32_t` | length of `seq`, `0` when `seq` is `NULL` |

`vg_edge_t` - one edge, by node id rather than by index.

| Field | Type | Notes |
| --- | --- | --- |
| `from`, `to` | `int64_t` | endpoint node **ids** |
| `from_start` | `int` | the edge leaves the 5' (start) side of `from` |
| `to_end` | `int` | the edge enters the 3' (end) side of `to` |
| `overlap` | `int32_t` | overlap length in bp, usually `0` |

`vg_step_t` - one visit inside a path.

| Field | Type | Notes |
| --- | --- | --- |
| `node_id` | `int64_t` | the visited node's id |
| `is_reverse` | `int` | visited in reverse-complement orientation |

`vg_path_t` - one named walk.

| Field | Type | Notes |
| --- | --- | --- |
| `name` | `char *` | owned; never `NULL` after a successful read - an unnamed path gets `""` |
| `step` | `vg_step_t *` | owned array of visits, in order |
| `n_step`, `m_step` | `int32_t` | used and allocated step counts |
| `is_circular` | `int` | decoded from the message; the GFA writer ignores it |

`vg_graph_t` - the three arrays, each with a used and an allocated count.

| Field | Type | Notes |
| --- | --- | --- |
| `node`, `n_node`, `m_node` | `vg_node_t *`, `int32_t` | nodes in file order |
| `edge`, `n_edge`, `m_edge` | `vg_edge_t *`, `int32_t` | edges in file order |
| `path`, `n_path`, `m_path` | `vg_path_t *`, `int32_t` | paths in file order |

The `m_*` fields are the reader's own growth bookkeeping; loop over the `n_*`
counts. Nodes are appended in the order the file lists them and are never sorted
or de-duplicated, and there is no id-to-index table, so looking a node up by id
means scanning `g->node`. Build your own map if you need repeated lookups.

## Reading and releasing

### `vg_read`

```c
vg_graph_t *vg_read(const char *fn);
```

Opens `fn` through zlib - which transparently accepts gzip, BGZF and raw
uncompressed input - and accumulates every `Graph` chunk into one freshly
allocated graph. Returns `NULL` only when the file cannot be opened or an
allocation fails; both cases are logged under the `vg` tag.

Decode problems are not failures. A truncated group header, a truncated message,
an implausible message size or an unparseable tag logs a warning (or, for a
malformed blob, nothing at all) and stops the scan, and the graph collected so
far is returned. A file that is not a `.vg` at all can therefore come back as a
valid but empty graph, so check the counts.

```c
vg_graph_t *g = vg_read("graph.vg");
if (!g) return 1;   // only a failed open or OOM; the reason is already logged

// Decoding is lenient: a truncated or non-vg input returns the prefix that did
// parse, so an empty graph is the signal that the input was not what you think.
if (g->n_node == 0) {
    ak_log(AK_LOG_ERROR, "vg", "no nodes decoded; is this really a .vg file?");
    vg_graph_destroy(g);
    return 1;
}

printf("%d nodes, %d edges, %d paths\n", g->n_node, g->n_edge, g->n_path);
printf("first node: id %lld, %u bp\n",
       (long long)g->node[0].id, g->node[0].seq_len);

vg_graph_destroy(g);
```

### `vg_graph_destroy`

```c
void vg_graph_destroy(vg_graph_t *g);
```

Releases the graph and everything it owns - every node sequence, the node and
edge arrays, and every path's name and step array. Safe to call with `NULL`,
which makes it usable on every error path. Nothing the graph hands out is
allocated separately, so there is nothing else to free and nothing borrowed from
it stays valid afterwards.

```c
vg_graph_t *g = vg_read("graph.vg");
if (!g) return 1;

FILE *out = fopen("out.gfa", "w");
if (!out) {
    vg_graph_destroy(g);   // release before bailing out
    return 1;
}

// Sequences point into the graph, so copy anything you want to outlive it.
fprintf(out, "H\tVN:Z:1.0\n");
if (g->n_node && g->node[0].seq)
    fprintf(out, "S\t%lld\t%s\n", (long long)g->node[0].id, g->node[0].seq);

fclose(out);
vg_graph_destroy(g);
```

## Converting to GFA

The three arrays map onto GFA1 one for one: nodes become `S` lines, edges become
`L` lines, paths become `P` lines. Node ids carry over unchanged, so no
renumbering is needed.

The part worth getting right is edge orientation. `vg` records an edge by the
*sides* it touches, while GFA records it by the *orientations* of the two
segments. The default edge leaves the 3' end of `from` and enters the 5' start
of `to`, which is GFA's `+ +`; each flag inverts its own end:

| `from_start` | `to_end` | GFA `L` line | Meaning |
| --- | --- | --- | --- |
| 0 | 0 | `L from + to + ...` | end of `from` to start of `to` (the ordinary case) |
| 1 | 0 | `L from - to + ...` | start of `from` to start of `to` |
| 0 | 1 | `L from + to - ...` | end of `from` to end of `to` |
| 1 | 1 | `L from - to - ...` | start of `from` to end of `to` |

So the from-orientation is `from_start ? '-' : '+'` and the to-orientation is
`to_end ? '-' : '+'`. The trap is the asymmetry in the names: it is `to_end`,
not `to_start`, that flips the second character, because the flags name the side
touched while the GFA characters name the strand traversed.

Path steps need no such translation - `is_reverse` maps straight onto `-`.

```c
vg_graph_t *g = vg_read("graph.vg");
if (!g) return 1;

printf("H\tVN:Z:1.0\n");

// S lines: seq is NULL when the Node message carried no sequence field.
for (int32_t i = 0; i < g->n_node; i++)
    printf("S\t%lld\t%s\n", (long long)g->node[i].id,
           g->node[i].seq ? g->node[i].seq : "*");

// L lines: each flag inverts its own end; overlap is usually 0, giving "0M".
for (int32_t i = 0; i < g->n_edge; i++) {
    const vg_edge_t *e = &g->edge[i];
    printf("L\t%lld\t%c\t%lld\t%c\t%dM\n",
           (long long)e->from, e->from_start ? '-' : '+',
           (long long)e->to,   e->to_end     ? '-' : '+',
           e->overlap);
}

// P lines: steps are already in visit order, and name is never NULL.
for (int32_t i = 0; i < g->n_path; i++) {
    const vg_path_t *p = &g->path[i];
    printf("P\t%s\t", p->name);
    for (int32_t s = 0; s < p->n_step; s++)
        printf("%s%lld%c", s ? "," : "", (long long)p->step[s].node_id,
               p->step[s].is_reverse ? '-' : '+');
    printf("\t*\n");   // no overlaps recorded per step
}

vg_graph_destroy(g);
```

This is exactly what [`src/cli/cmd_vg2gfa.c`](../src/cli/cmd_vg2gfa.c) does. Note
that the `L` lines reference node ids that an `S` line may never define: edges
and paths are not checked against the node set during decoding.

---

[Back to the library index](README.md)
