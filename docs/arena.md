# `arena` - the append-only byte store behind the graph strings

Source: [`src/lib/arena.c`](../src/lib/arena.c) &middot; Header: [`include/akhal/arena.h`](../include/akhal/arena.h)

A graph holds one string per segment, and a large one holds a lot of them: a
chr22 graph carries 17.9 M segment sequences, a whole-genome one hundreds of
millions. Held as individual `malloc()`s they cost an allocator header apiece
and one `free()` each to release - and that teardown is not free. It walks
every chunk in address order that the allocator chose, which on a heap that has
just been streamed end to end is a cache miss per value and, under memory
pressure, a page fault per value.

The arena replaces that with a list of large blocks. Values are copied in and
handed back as pointers into a block; releasing the arena frees the blocks. On
the same chr22 graph that is **73 `free()`s instead of 17,872,082**, and 562 MB
of sequence rather than the 858 MB the same bytes occupy as separate chunks.

```c
#include "akhal/arena.h"
#include "akhal/error.h"   // AK_E* return codes, ak_log()
#include "akhal/gfa.h"     // only for the gfa_t example at the end
```

## Contents

- [Why blocks and not one buffer](#why-blocks-and-not-one-buffer)
- [Block sizes](#block-sizes) - ordinary growth, and values larger than a block
- [Putting and releasing](#putting-and-releasing) - [`ak_arena_put`](#ak_arena_put), [`ak_arena_destroy`](#ak_arena_destroy), [`ak_arena_bytes`](#ak_arena_bytes)

## Why blocks and not one buffer

The obvious shape is one flat buffer that grows. It does not work, because
growing it means `realloc`, `realloc` may move it, and every pointer already
handed out becomes garbage. You would have to store offsets instead of
pointers and resolve them at every use - which is a different, more invasive
design.

A block never moves once it has been cut. Only the *array of block pointers*
is `realloc`'d, and nothing outside the arena points into that. So callers keep
plain `const char *` and the arena keeps its promise:

> every pointer `ak_arena_put()` ever returned stays valid until
> `ak_arena_destroy()`.

Zero-initialize to create one - `ak_arena_t a = {0};`, or as a field of a
`calloc`'d parent, which is how `gfa_t` and `vg_graph_t` hold theirs. There is
no init function.

## Block sizes

Two rules decide how big each new block is.

**Ordinary growth.** The first block is `AK_ARENA_MIN` (64 KiB) and each later
one is twice the last, up to `AK_ARENA_MAX` (8 MiB). A ten-segment test graph
therefore never allocates megabytes, and a whole-chromosome one stops paying
per-block overhead after a handful of blocks.

**Values larger than a block.** The size keeps doubling until the value itself
fits, however large it is. A segment that a compaction folded into a megabase -
or a linear reference held as one 250 Mbp `S` line - still lands in one piece
rather than being refused. That growth only *sticks* while it stays under the
ceiling: past it the oversized block is a one-off, so a single huge value does
not leave every later block huge as well.

A value never straddles two blocks. When the current block has less room than
the value needs, whatever is left of it is abandoned and a new block is cut -
at most `AK_ARENA_MAX` wasted across the whole arena.

## Putting and releasing

### `ak_arena_put`

```c
const char *ak_arena_put(ak_arena_t *a, const char *s, size_t len);
```

Copies `len` bytes in and NUL-terminates the copy, returning a pointer to it or
`NULL` on allocation failure (logged). The terminator is added on top of `len`,
so the stored value is `len + 1` bytes and can be handed to anything expecting
a C string - while the caller's own length field stays the authority on how
many bytes it holds.

`s` may be `NULL` when `len` is 0, which stores an empty string.

```c
ak_arena_t a = {0};          // no init call; zeroing is enough

const char *x = ak_arena_put(&a, "ACGT", 4);
const char *y = ak_arena_put(&a, "TT", 2);
if (!x || !y) {              // the reason has already been logged
    ak_arena_destroy(&a);
    return 1;
}

// x is still valid: adding y cannot have moved it
printf("%s %s\n", x, y);     // ACGT TT

ak_arena_destroy(&a);
```

### `ak_arena_destroy`

```c
void ak_arena_destroy(ak_arena_t *a);
```

Frees every block and zeroes the arena, so a second call is a no-op and an
arena embedded in a `calloc`'d struct is safe to destroy whether or not it was
ever used. Every pointer the arena returned dangles afterwards.

```c
ak_arena_t a = {0};
for (int i = 0; i < 1000000; i++) ak_arena_put(&a, "ACGT", 4);

// One free() per block, not per value. That is the whole point: the same
// million values as separate malloc()s would be a million free()s here.
ak_arena_destroy(&a);
ak_arena_destroy(&a);        // safe; the first call zeroed it
```

### `ak_arena_bytes`

```c
size_t ak_arena_bytes(const ak_arena_t *a);
```

Total size of the blocks, allocated rather than used, for reporting. `NULL`
reads as 0.

```c
gfa_t *g = gfa_read("graph.gfa", 0);
if (!g) return 1;

printf("%d segment(s) in %d arena block(s), %.1f MB\n",
       gfa_n_seg(g), g->strs.n_blk, ak_arena_bytes(&g->strs) / 1e6);

gfa_destroy(g);
```

---

[Back to the library index](README.md)
