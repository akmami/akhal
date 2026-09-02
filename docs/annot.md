# `annot` - per-node origin annotations

Source: [`src/lib/annot.c`](../src/lib/annot.c) &middot; Header: [`include/akhal/annot.h`](../include/akhal/annot.h)

Answers one question for every graph node: where does it come from? An
annotated node maps to a single info string - `SNP chr1:12345 A>G rs99`,
`REF chr1 0-1200`, `SEQ sampleA 4501` - and the store keeps those strings
packed end to end in one NUL-separated buffer. A record holds only the node id,
the byte offset of its string, its length and a kind, so the whole store is one
allocation per component and a lookup is a hash probe plus a pointer into the
buffer. No per-node string allocation happens at all.

Builders fill the store from a graph plus a VCF (alt-path matching against the
backbone) or a FASTA (greedy walk of each sequence through the graph). Both are
optional; nodes nothing annotates simply have no record. The store round-trips
through a compact binary file, so the annotation is computed once and queried
later without the graph in hand.

```c
#include "akhal/annot.h"   // annot_t, annot_rec_t, ANNOT_* kinds, the annot_* API
#include "akhal/gfa.h"     // gfa_read() and the GFA_* flags (annot.h pulls this in)
#include "akhal/fasta.h"   // fasta_read(), for annot_build_fasta()
#include "akhal/error.h"   // AK_OK / AK_E* return codes, ak_log(), ak_strerror()
```

## Contents

- [Kinds](#kinds)
- [Layout](#layout)
- [Lifecycle](#lifecycle) - [`annot_init`](#annot_init), [`annot_destroy`](#annot_destroy)
- [Writing](#writing) - [`annot_set`](#annot_set), [`annot_add`](#annot_add)
- [Querying](#querying) - [`annot_get`](#annot_get), [`annot_kind_name`](#annot_kind_name), [counts and record accessors](#counts-and-record-accessors)
- [Builders](#builders) - [`annot_backbone`](#annot_backbone), [`annot_build_vcf`](#annot_build_vcf), [`annot_build_fasta`](#annot_build_fasta)
- [File format](#file-format) - [`annot_write`](#annot_write), [`annot_read`](#annot_read)

## Kinds

| Kind | Value | Meaning |
| --- | --- | --- |
| `ANNOT_UNKNOWN` | 0 | no record for this node; never stored, only returned |
| `ANNOT_BACKBONE` | 1 | the node lies on the backbone reference path |
| `ANNOT_INFO` | 2 | the node carries variant or origin information |

`ANNOT_UNKNOWN` is what `annot_get()` reports for a node the store has never
heard of. A record always carries `ANNOT_BACKBONE` or `ANNOT_INFO`, and
`ANNOT_INFO` is the default a freshly created record starts with.

## Layout

```c
typedef struct {
    uint64_t id;     // node (segment) id
    uint64_t off;    // byte offset of the info string in the shared buffer
    uint32_t len;    // info string length, excluding the NUL terminator
    uint8_t  kind;   // ANNOT_BACKBONE or ANNOT_INFO
} annot_rec_t;

typedef struct {
    annot_rec_t *rec;      // records, one per annotated node
    int64_t      n, m;     // used / allocated record count
    char        *buf;      // packed NUL-separated info strings
    uint64_t     buf_l;    // used bytes in buf (including NUL terminators)
    uint64_t     buf_m;    // allocated bytes in buf
    void        *idx;      // opaque id -> record-index hash table
} annot_t;
```

Records are keyed by the *segment id as written in the GFA*, not by array
index, so a store stays meaningful without the graph it was built from.
`buf` is append-only: `annot_set()` and `annot_add()` always write a fresh
copy and leave the previous bytes stranded in the buffer, which is the price of
never having to move a string once it is placed. Rewriting the same node many
times therefore grows the file even though the record count does not.

A record with `len == 0` is an empty annotation - it exists (so `annot_get()`
reports its kind) but has no string, and `off` is 0 rather than a real offset.

## Lifecycle

### `annot_init`

```c
annot_t *annot_init(void);
```

Allocates an empty store: no records, no buffer, and an empty id index.
Returns `NULL` on allocation failure, which is already logged.

```c
annot_t *a = annot_init();
if (!a) return 1;   // the reason is already logged

// Nodes are addressed by GFA segment id, so a store needs no graph to be filled.
if (annot_set(a, 7, ANNOT_INFO, "SEQ sampleA 0") != AK_OK) {
    annot_destroy(a);
    return 1;
}

printf("%lld record(s)\n", (long long)annot_n(a));

annot_destroy(a);
```

### `annot_destroy`

```c
void annot_destroy(annot_t *a);
```

Releases the record array, the shared buffer and the id index. Safe with
`NULL`, so it works on every error path. Every `const char *` handed out by
`annot_get()` or `annot_info_at()` points into `buf` and dies with the store.

```c
annot_t *a = annot_read("graph.annot");
if (!a) return 1;

const char *info;
if (annot_get(a, 7, &info) != ANNOT_UNKNOWN && info) {
    // Copy anything that has to outlive the store; the string is borrowed.
    char *keep = strdup(info);
    annot_destroy(a);
    free(keep);
    return 0;
}

annot_destroy(a);
```

## Writing

### `annot_set`

```c
int annot_set(annot_t *a, uint64_t id, int kind, const char *info);
```

Sets, or replaces, a node's annotation. The string is copied into the shared
buffer; `NULL` or `""` stores an empty annotation rather than removing the
record. `kind` is always written, so this is also how you change an existing
record's kind. Returns `AK_OK` or `AK_ENOMEM`. There is no delete operation.

```c
annot_t *a = annot_init();
if (!a) return 1;

int rc = annot_set(a, 12, ANNOT_BACKBONE, "REF chr1 0-1200");
// Setting the same node again replaces the string and the kind outright; the
// old bytes stay in the buffer, they are simply no longer referenced.
if (rc == AK_OK) rc = annot_set(a, 12, ANNOT_INFO, "SNP chr1:400 A>G");
if (rc != AK_OK) {
    ak_log(AK_LOG_ERROR, NULL, "%s", ak_strerror(rc));
    annot_destroy(a);
    return 1;
}

annot_destroy(a);
```

### `annot_add`

```c
int annot_add(annot_t *a, uint64_t id, int kind, const char *info);
```

Like `annot_set()`, except that a node which already has a non-empty string
keeps it and the new one is appended after `"; "`. Overlapping variants
therefore accumulate instead of overwriting each other. `kind` is only applied
when the record is created - an existing record keeps the kind it had, so
adding `ANNOT_INFO` to a backbone node does not demote it.

```c
annot_t *a = annot_init();
if (!a) return 1;

int rc = annot_add(a, 12, ANNOT_INFO, "SNP chr1:400 A>G");
if (rc == AK_OK) rc = annot_add(a, 12, ANNOT_INFO, "SNP chr1:402 C>T");
if (rc != AK_OK) { annot_destroy(a); return 1; }

const char *info;
annot_get(a, 12, &info);
printf("%s\n", info);   // "SNP chr1:400 A>G; SNP chr1:402 C>T"

annot_destroy(a);
```

## Querying

### `annot_get`

```c
int annot_get(const annot_t *a, uint64_t id, const char **info);
```

Looks up one node. Returns its kind, or `ANNOT_UNKNOWN` when no record exists.
`info` may be `NULL` if you only want the kind; otherwise it is set to the
string, or to `NULL` for an unknown node *and* for a record whose annotation is
empty - so test the return value, not the pointer, to decide whether a record
exists. The string points into the shared buffer and is invalidated by the next
`annot_set()` or `annot_add()`, which may reallocate it.

```c
annot_t *a = annot_read("graph.annot");
if (!a) return 1;

const char *info;
int kind = annot_get(a, 12, &info);
// A record can exist with no string, so the kind is the existence test.
if (kind == ANNOT_UNKNOWN)
    printf("node 12 is unannotated\n");
else
    printf("node 12 is %s: %s\n", annot_kind_name(kind), info ? info : "(empty)");

annot_destroy(a);
```

### `annot_kind_name`

```c
const char *annot_kind_name(int kind);
```

Maps a kind to a short label for output. `ANNOT_BACKBONE` is `"backbone"`,
`ANNOT_INFO` is `"annot"` (not "info"), and anything else - including
`ANNOT_UNKNOWN` and a corrupt value read from a file - is `"unknown"`. Never
returns `NULL`.

```c
annot_t *a = annot_init();
if (!a) return 1;
if (annot_set(a, 3, ANNOT_BACKBONE, "REF chr1 0-100") != AK_OK) {
    annot_destroy(a);
    return 1;
}

// Safe on any int, so a kind read straight out of a file needs no validation.
for (int k = 0; k < 4; k++)
    printf("%d -> %s\n", k, annot_kind_name(k));

annot_destroy(a);
```

### Counts and record accessors

```c
int64_t            annot_n(const annot_t *a);
const annot_rec_t *annot_at(const annot_t *a, int64_t i);
const char        *annot_info_at(const annot_t *a, int64_t i);
```

Inline, unchecked accessors for iterating the record array in insertion order.
They do no bounds checking, so keep `i` below `annot_n()`. `annot_info_at()`
returns `NULL` for an empty annotation, and both borrow - the store owns the
record and the string.

```c
annot_t *a = annot_read("graph.annot");
if (!a) return 1;

// Iteration order is the order records were created, not sorted by id.
for (int64_t i = 0; i < annot_n(a); i++) {
    const annot_rec_t *r = annot_at(a, i);
    const char *info = annot_info_at(a, i);
    printf("%llu\t%s\t%s\n", (unsigned long long)r->id,
           annot_kind_name(r->kind), info ? info : ".");
}

annot_destroy(a);
```

## Builders

The two builders expect to run in order: `annot_backbone()` first, so the
reference nodes are known, then `annot_build_vcf()` and/or
`annot_build_fasta()`, which both treat backbone membership as fixed.

### `annot_backbone`

```c
int32_t annot_backbone(annot_t *a, const gfa_t *g, const char *ref_path);
```

Marks every segment of one `P` line as `ANNOT_BACKBONE` with info
`REF <path> <start>-<end>`, where the coordinates come from the segment's
`start`/`end` fields as the reader laid them out. Returns the path index, which
`annot_build_vcf()` wants, or a negative `AK_E*` code (`AK_EINVAL` when the
graph has no paths or no path of that name), in which case nothing is marked.

Path names are matched exactly against `g->path[]`, with no fragment chaining -
a reference split by `vg` into `chr22[0]`, `chr22[21]`, ... must be named by
one of those literal fragment names, and only that fragment is marked.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
annot_t *a = annot_init();
if (!a) { gfa_destroy(g); return 1; }

// NULL takes the graph's first P line; the name must match a P line exactly.
int32_t ref = annot_backbone(a, g, "chr22");
if (ref < 0) {
    ak_log(AK_LOG_ERROR, NULL, "%s", ak_strerror(ref));
    annot_destroy(a);
    gfa_destroy(g);
    return 1;
}

annot_destroy(a);
gfa_destroy(g);
```

### `annot_build_vcf`

```c
int64_t annot_build_vcf(annot_t *a, const gfa_t *g, int32_t ref_path, const char *vcf_fn);
```

Annotates alternative-allele nodes from a VCF. For each record whose `CHROM`
equals the backbone path name, the shared REF/ALT prefix is stripped and the
branch point is found as the backbone segment whose reference `end` equals the
variant's 0-based position; each non-backbone successor is then walked through
its unique non-backbone chain, and a walk whose concatenated sequence equals
the allele exactly is annotated `<TYPE> <chrom>:<pos> <ref>><alt> [<id>]`.
Requires `GFA_LINKS | GFA_PATHS`.

Returns the number of alleles matched, or a negative `AK_E*` code. Matching is
strict and silent about misses: pure deletions have no node to annotate, a
bubble whose interior branches is ambiguous and is skipped, a walk longer than
4096 nodes is abandoned, and every such case is logged at DEBUG only. A return
of 0 is a normal result for a VCF that does not line up with the graph.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
annot_t *a = annot_init();
if (!a) { gfa_destroy(g); return 1; }

// The backbone has to be marked first: it supplies both the CHROM name to
// match on and the membership flags that stop the walk.
int32_t ref = annot_backbone(a, g, NULL);
int64_t hit = ref >= 0 ? annot_build_vcf(a, g, ref, "calls.vcf") : ref;
if (hit < 0)
    ak_log(AK_LOG_ERROR, NULL, "%s", ak_strerror((int)hit));
else
    printf("%lld allele(s) annotated\n", (long long)hit);

annot_destroy(a);
gfa_destroy(g);
```

### `annot_build_fasta`

```c
int64_t annot_build_fasta(annot_t *a, const gfa_t *g, const fasta_t *fa);
```

Annotates nodes by walking each FASTA record through the graph. A start node is
located (a source node whose sequence prefixes the record, falling back to any
matching node) and the walk greedily takes the out-link whose sequence
continues the record, honouring link overlaps. Every visited node that is not
already `ANNOT_BACKBONE` gets `SEQ <name> <offset>` appended through
`annot_add()`, so a node shared by several samples collects one entry per
sample.

Requires `GFA_LINKS`; `GFA_PATHS` only matters if you marked a backbone first.
Returns the number of records traced end to end - a stalled walk is logged at
WARN, keeps the annotations gathered so far, and is not counted, so a return
below `fasta_n(fa)` means partial coverage rather than failure.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS);
if (!g) return 1;
fasta_t *fa = fasta_read("samples.fa");
if (!fa) { gfa_destroy(g); return 1; }
annot_t *a = annot_init();
if (!a) { fasta_destroy(fa); gfa_destroy(g); return 1; }

int64_t done = annot_build_fasta(a, g, fa);
// Fewer than fasta_n(fa) means some walks stalled, not that the call failed.
if (done >= 0 && done < fasta_n(fa))
    ak_log(AK_LOG_WARN, NULL, "%lld of %lld sequence(s) fully traced",
           (long long)done, (long long)fasta_n(fa));

annot_destroy(a);
fasta_destroy(fa);
gfa_destroy(g);
```

## File format

A store serialises to a flat binary file in **native endianness**, with no
padding or checksum:

```
magic "AKANNOT1"                                     8 bytes
uint64 n_rec
uint64 buf_len
n_rec x { uint64 id, uint64 off, uint32 len, uint8 kind }   21 bytes each
buf_len bytes of NUL-separated info strings
```

The record fields are written one at a time rather than as a struct, so the
on-disk record is 21 bytes with no alignment padding. The file is not portable
between machines of different endianness.

### `annot_write`

```c
int annot_write(const annot_t *a, const char *fn);
```

Writes the store to `fn`, truncating an existing file. Returns `AK_OK`,
`AK_EOPEN` if the file cannot be created, or `AK_EIO` if any write or the final
`fclose()` fails. Orphaned bytes left in the buffer by repeated
`annot_set()`/`annot_add()` calls are written out as well, so the file can be
larger than the live annotations need.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
annot_t *a = annot_init();
if (!a) { gfa_destroy(g); return 1; }

int32_t ref = annot_backbone(a, g, NULL);
if (ref >= 0) annot_build_vcf(a, g, ref, "calls.vcf");

// The graph is not needed to query the result later, only to build it.
int rc = annot_write(a, "graph.annot");
if (rc != AK_OK) ak_log(AK_LOG_ERROR, NULL, "%s", ak_strerror(rc));

annot_destroy(a);
gfa_destroy(g);
```

### `annot_read`

```c
annot_t *annot_read(const char *fn);
```

Loads a file written by `annot_write()` and rebuilds the id index as it goes.
Returns `NULL` on a missing file, a wrong magic, or a truncated one; the reason
is logged. Sanity checks are cheap but real: every record's string must lie
wholly inside the buffer and the buffer must end on a NUL, and a duplicate id
keeps its first record.

```c
// No graph is involved: the store is keyed by segment id, not by index.
annot_t *a = annot_read("graph.annot");
if (!a) return 1;   // missing, wrong magic, or truncated - already logged

int64_t n_bb = 0;
for (int64_t i = 0; i < annot_n(a); i++)
    if (annot_at(a, i)->kind == ANNOT_BACKBONE) n_bb++;

printf("%lld of %lld record(s) are backbone\n",
       (long long)n_bb, (long long)annot_n(a));

annot_destroy(a);
```

---

[Back to the library index](README.md)
