# libakhal - library reference

`akhal` is split the way `htslib` and `samtools` are: a reusable library
(`libakhal`, under [`src/lib/`](../src/lib)) that owns the file formats and data
structures, and a thin command-line front end (under [`src/cli/`](../src/cli))
that assembles those pieces into subcommands.

These pages document the library. One page per module, one section per public
function, each with a short snippet showing the data preparation the call needs.
For the command-line tool, see the [main README](../README.md).

## Modules

### Graphs

| Module | What it gives you |
| --- | --- |
| [`gfa`](gfa.md) | The (r)GFA graph model - segments, links, paths - plus the reader and writer every command uses, traversal over the CSR adjacency, fragmented-path merging, `SR` ranking, path-block rewriting and topological sort |
| [`vg`](vg.md) | Reader for vg's native `.vg` format (compressed Protobuf), decoded straight from the wire format with only zlib |

### Sequence and alignment formats

| Module | What it gives you |
| --- | --- |
| [`fasta`](fasta.md) | A whole FASTA file loaded into memory and indexed by sequence name |
| [`gaf`](gaf.md) | GAF alignment records, in both streaming and batch flavours, plus traversal of a GAF path string |
| [`sam`](sam.md) | The CIGAR operation alphabet, CIGAR expansion, and SAM header/record output |
| [`vcf`](vcf.md) | A minimal streaming VCF reader - the columns the toolkit needs, nothing more |

### Analysis

| Module | What it gives you |
| --- | --- |
| [`call`](call.md) | Labels a reference backbone through the graph, then reports everything the graph carries beyond it as variants, and writes them as VCF |
| [`diff`](diff.md) | Compares two graphs that disagree on segment ids: segments matched by sequence content, links through the labelling that produces, paths by the bases they spell |
| [`rgfa`](rgfa.md) | Works the rGFA stable-sequence tags out of a graph's own P lines: a backbone at rank 0, everything that detours off it one rank deeper and offset from where it left |
| [`compact`](compact.md) | Folds every run of non-branching segments into one node, keeping the bases, the links between runs and the paths exactly as they were |
| [`annot`](annot.md) | A node-origin store: which variant or which sample sequence explains each node, queryable and round-trippable through a compact binary file |

### Infrastructure

| Module | What it gives you |
| --- | --- |
| [`io`](io.md) | The line-oriented input handle every text reader pulls lines through |
| [`kstr`](kstr.md) | A growable, always-NUL-terminated string buffer |
| [`util`](util.md) | Complement / reverse-complement, suffix testing, and mean / variance / stddev |
| [`error`](error.md) | The `AK_E*` return codes and the `ak_log()` diagnostics every module reports through |

## Building against the library

`make` produces `lib/libakhal.a` alongside the `akhal` executable. Add
`include/` to the header search path and link the archive:

```sh
gcc -O2 -Wall -Iinclude my_tool.c lib/libakhal.a -lm -lz -o my_tool
```

`-lz` is needed because the `vg` reader decompresses `.vg` files; `-lm` because
`util` uses `sqrt`. Include the headers you actually use, or pull in the whole
public API at once:

```c
#include "akhal/akhal.h"   // umbrella: every public header
```

## Conventions

These hold across every module, so each page assumes them rather than repeating
them.

**Errors.** Code under `src/lib/` never calls `exit()`. An int-returning
function signals failure with a negative `AK_E*` code, a pointer-returning one
with `NULL`, and either may explain itself through `ak_log()`. Only the CLI
layer decides whether to abort. `AK_OK` is 0 and every error is negative, so
`if (rc < 0)` is the idiomatic test. See [`error`](error.md) for the full list.

**Ownership.** A `*_read()` / `*_init()` / `*_open()` function returns something
you release with the matching `*_destroy()` / `*_close()`. Every destroy is safe
to call with `NULL`, which is what makes the `goto done` error paths in the
commands work. Where a function hands you a buffer to free yourself - such as
[`gfa_merge_segs`](gfa.md#gfa_merge_segs) or
[`sam_rg_prefix`](sam.md#sam_rg_prefix) - its page says so explicitly.

**Array + dict.** The graph, the FASTA store and the annotation store all use
the same layout: records in one contiguous array, plus a hash table mapping the
external key to an array index. Anything cross-referencing a record does so by
index, not by pointer, so the arrays stay reallocatable.

**Reading only what you asked for.** [`gfa_read`](gfa.md#gfa_read) takes flags
(`GFA_LINKS`, `GFA_PATHS`, `GFA_VALIDATE`) and skips the work you did not
request. A function that needs one says so and returns an error when it is
missing, rather than quietly returning empty results.

**The file wins on ranks.** rGFA's `SR:i:` tag marks the reference backbone as
rank 0. A file that carries those tags is authoritative and the reader never
overwrites them; `g->has_sr` says whether it saw any, and only when it did not
are ranks derived from the paths. Overwriting them is something a caller asks
for explicitly, through [`gfa_rank_paths`](gfa.md#gfa_rank_paths) or
[`gfa_rank_mark`](gfa.md#gfa_rank_mark).

**Snippets.** Every example on these pages is a statement body: paste it into a
function. It assumes `<stdio.h>`, `<stdlib.h>` and `<string.h>` are included
alongside the `akhal/` headers listed at the top of that page. They are all
compile-checked against the real headers - `bash docs/.check_snippets.sh`
extracts every `c` block, wraps each in a function and builds the page as one
translation unit, so an example that drifts out of step with the API fails
loudly rather than quietly misleading someone.

## How the modules compose

The commands are mostly short assemblies of the pieces above. Three worked
shapes, with the module that owns each step:

**Extract a merged reference** - `akhal extract path`

```text
gfa_read(GFA_LINKS|GFA_PATHS)   gfa    load the graph
gfa_path_merge(g, "chr22")      gfa    chain the fragmented P lines
gfa_merge_segs(g, m, c, &segs)  gfa    flatten one chain to segment indices
  -> write FASTA                       spell the sequences out
```

**Call variants from a graph** - `akhal extract vcf`

```text
gfa_read(GFA_LINKS|GFA_PATHS)   gfa    load the graph
call_ref_path(g, "chr22")       call   label the backbone and its coordinates
  or call_ref_fasta(g, fa, ...) call   ... or trace an external reference
call_variants(g, ref)           call   collect the detours off it
call_write_vcf(c, ref, "out")   call   emit VCF 4.2
```

**Re-root a graph on a different reference** - `akhal rank --fasta`

```text
gfa_read(GFA_LINKS|GFA_PATHS)   gfa    load the graph
call_ref_fasta(g, fa, NULL)     call   trace the reference, labelling on/walk
gfa_rank_mark(g, ref->on)       gfa    rank 0 on that backbone, 1 elsewhere
gfa_clear_paths(g)              gfa    drop the old P lines
gfa_add_path(g, name, walk...)  gfa    install the traced walk in their place
gfa_write(g, out)               gfa    emit the re-ranked graph
```

**Compare two graphs that number their nodes differently** - `akhal compare`

```text
gfa_read(GFA_LINKS|GFA_PATHS)   gfa    load both graphs
diff_map(a, b)                  diff   match segments by sequence, label them
diff_graphs(a, b)               diff   ... and compare links and paths through it
diff_identical(d)               diff   the whole verdict, for an exit status
```

**Label a plain GFA as an rGFA** - `akhal gfa2rgfa`

```text
gfa_read(GFA_LINKS|GFA_PATHS)   gfa    load the graph
rgfa_build(g, "chr22", &st)     rgfa   chain the P lines, then walk them all,
                                       writing SN/SO/SR onto every segment
gfa_write_rgfa(g, out)          gfa    emit the graph, tags and all
```

**Trace where each node came from** - `akhal annotate`

```text
gfa_read(GFA_LINKS|GFA_PATHS)   gfa    load the graph
annot_init()                    annot  empty store
annot_backbone(a, g, ref)       annot  mark the reference path
annot_build_vcf(a, g, bb, vcf)  annot  match variants onto bubble sides
annot_build_fasta(a, g, fa)     annot  walk sample sequences through the graph
annot_write(a, "out.annot")     annot  persist for later lookup
```

---

[Back to the main README](../README.md)
