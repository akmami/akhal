# `dot` - Graphviz output, for looking at a graph

Source: [`src/lib/dot.c`](../src/lib/dot.c) &middot; Header: [`include/akhal/dot.h`](../include/akhal/dot.h)

Every other module here computes on a graph. This one draws it: a Graphviz
digraph with a box per segment and an arrow per link, laid out left to right.

It is meant for graphs small enough to look at - a bubble, a test case, a region
pulled out of something bigger. Nothing stops you handing it a whole pangenome,
but Graphviz will accept that and spend a very long time on it, so
[`dot_write`](#dot_write) says so first.

This module backs [`akhal gfa2dot`](../README.md#10-gfa2dot).

```c
#include "akhal/dot.h"     // dot_write() and the DOT_* constants
#include "akhal/gfa.h"     // gfa_read() and the GFA_* flags (dot.h pulls this in)
#include "akhal/error.h"   // AK_OK / AK_E* return codes, ak_log()
```

## Contents

- [What it draws](#what-it-draws) - [`dot_write`](#dot_write)

## What it draws

GFA is bidirected and DOT is not, so a link is drawn as one arrow with the
strands it joins written on it - `+/-` where it flips - and nothing written
where both ends are forward. A clean variation graph therefore comes out as
plain arrows, and a strand flip is visible precisely because it is the only
thing labelled. An overlap is labelled where there is one, so a graph of blunt
joins stays quiet too.

Nothing is lost in the drawing, but it is a directed reading of an undirected
thing: the arrow direction is the direction the `L` line was written in, not a
claim about the only way through.

A node shows its id, and its bases when the segment is at most `DOT_SEQ_MAX`
(20) of them - which in a variation graph is most of them, and is exactly what
you want to read off a test case. Longer segments show a length instead.
`DOT_IDS` drops the sequences whatever their length.

On a file that ranks itself - `g->has_sr`, the same test [`sort`](gfa.md#gfa_write)
uses - the rank-0 backbone is shaded apart from everything hanging off it, so an
rGFA's structure reads at a glance. A plain GFA is drawn unshaded rather than
being coloured by ranks the reader merely derived.

| Constant | Value | Means |
| --- | --- | --- |
| `DOT_SEQ_MAX` | 20 | bases at or under which a segment is labelled with its sequence |
| `DOT_BIG` | 10000 | nodes past which the size is reported before writing |
| `DOT_IDS` | flag | label with ids alone, however short the bases |

### `dot_write`

```c
int dot_write(const gfa_t *g, FILE *out, int flags);
```

Writes the whole digraph and returns `AK_OK`, or `AK_EIO` if the stream went
bad. `flags` is `DOT_IDS` or 0.

No read flags are required beyond the segments themselves, but without
`GFA_LINKS` there is nothing to draw between them, so a graph read without it
comes out as a field of unconnected boxes. `GFA_PATHS` is not used at all -
paths are not drawn.

```c
// links are what makes it a picture; paths are not drawn, so leave them out
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS);
if (!g) return 1;

FILE *out = fopen("graph.dot", "w");
if (!out) {
    gfa_destroy(g);
    return 1;
}

int rc = dot_write(g, out, 0);   // DOT_IDS for a graph that is getting crowded

fclose(out);
gfa_destroy(g);
if (rc != AK_OK) return 1;
```

Then, outside the program:

```sh
dot -Tpng graph.dot -o graph.png     # or -Tsvg, which scales better to read
```

---

[Back to the library index](README.md)
