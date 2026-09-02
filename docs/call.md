# `call` - variant discovery against a graph backbone

Source: [`src/lib/call.c`](../src/lib/call.c) &middot; Header: [`include/akhal/call.h`](../include/akhal/call.h)

Reads a graph as a reference plus its differences. One walk through the graph
is declared the backbone - a `P` line, or an external sequence traced through
the nodes - and every node on that walk is stamped with its offset on the
resulting reference sequence. Everything the backbone does not cover is, by
definition, a variation.

Each detour that leaves the backbone at one node and rejoins it at a later one
spells an alternate allele over the reference span between them; a link that
skips forward along the backbone itself spells a deletion of the bases it jumps
over. Detours sharing a REF span are grouped into one multi-allelic record, so
the result is an ordered array that maps row for row onto a VCF.

```c
#include "akhal/call.h"    // call_ref_t, call_var_t, call_t and the call_* API
#include "akhal/gfa.h"     // gfa_read() and the GFA_* flags (call.h pulls this in)
#include "akhal/fasta.h"   // fasta_read(), for a backbone that comes from a FASTA
#include "akhal/error.h"   // AK_OK / AK_E* return codes, ak_log(), ak_strerror()
```

## Contents

- [Coordinates](#coordinates)
- [Search caps](#search-caps)
- [The backbone](#the-backbone) - [`call_ref_path`](#call_ref_path), [`call_ref_fasta`](#call_ref_fasta), [`call_ref_destroy`](#call_ref_destroy)
- [Variants](#variants) - [`call_variants`](#call_variants), [`call_destroy`](#call_destroy), [`call_n`](#call_n), [`call_at`](#call_at)
- [Output](#output) - [`call_write_vcf`](#call_write_vcf)

## Coordinates

A backbone is a labelled reference: the concatenated sequence, plus two
parallel arrays indexed by *segment index* (not segment id).

```c
typedef struct {
    char    *name;    // owned: reference name, used as the VCF CHROM column
    char    *seq;     // owned: backbone sequence, 5' -> 3'
    int64_t  len;     // length of seq
    uint8_t *on;      // length n_seg: 1 when the segment lies on the backbone
    int64_t *pos;     // length n_seg: 0-based start offset, -1 when off it
    int32_t  n_seg;   // segment count of the graph this was built from
    uint32_t *walk;   // owned: the backbone's segments in the order it visits
    int64_t  n_walk;  // length of walk; a repeated segment appears each time
} call_ref_t;
```

A variant record is the grouped form of one branch point's alleles.

```c
typedef struct {
    int64_t  pos;     // 0-based offset of the first REF base
    int64_t  end;     // 0-based, exclusive: one past the last REF base
    char    *ref;     // owned: REF allele
    char    *alt;     // owned: ALT alleles, comma-separated
    char    *type;    // owned: SNP/MNP/INS/DEL/COMPLEX tags, one per ALT
} call_var_t;
```

`pos` and `end` are **0-based and half-open**, so the REF allele is exactly
`ref->seq[pos .. end)` and `strlen(v->ref) == v->end - v->pos`. VCF POS is
1-based, so [`call_write_vcf`](#call_write_vcf) writes `pos + 1`; the `END`
value it puts in INFO is `end` unchanged, which is already the 1-based
*inclusive* end of the same span.

Records where either side of the difference is empty (a pure insertion or
deletion) carry the preceding reference base on both alleles, as VCF requires.
`pos` is then one base earlier than the difference itself. An allele that would
need an anchor before offset 0 has nowhere to take it from, so `call_variants()`
logs it at debug level and drops it.

`n_seg` records the graph the backbone was built from. `call_variants()`
compares it against the graph it is given and refuses to run if they differ.

`on` and `pos` answer "is this segment on the backbone, and where"; `walk`
answers "what order does the backbone visit them in". The pair is what lets an
external reference be installed as the graph's own backbone - feed `on` to
[`gfa_rank_mark`](gfa.md#gfa_rank_mark) and `walk` to
[`gfa_add_path`](gfa.md#gfa_add_path), which is exactly what `akhal rank
--fasta` does.

## Search caps

A tangled region can spell an unbounded number of walks, so the search is
bounded. Exceeding a cap is not an error: the alleles found so far are kept,
the site is counted, and a summary is logged at WARN.

| Constant | Value | Limits |
| --- | --- | --- |
| `CALL_MAX_ALT` | 16 | alternate alleles collected at one branch point |
| `CALL_MAX_WALK` | 1024 | nodes visited along one alternate walk |
| `CALL_MAX_LEN` | 1000000 | bases spelled by one alternate allele |

`CALL_MAX_WALK` doubles as the cycle guard: a node already on the current walk
is not re-entered, so a cyclic graph terminates rather than spinning.

## The backbone

### `call_ref_path`

```c
call_ref_t *call_ref_path(const gfa_t *g, const char *path_name);
```

Builds the backbone from the graph's `P` lines. Because builders such as `vg`
split one reference over several consecutive `P` lines, the candidates are the
chains [`gfa_path_merge`](gfa.md#gfa_path_merge) stitches out of them, not the
raw paths - which is why this needs `GFA_LINKS | GFA_PATHS` even though it is a
path operation. Only the `GFA_PATHS` half is checked here; the missing
`GFA_LINKS` is reported by the merge and returns `NULL`.

`path_name` may name a whole reference (`chr22`) or one of its fragments. If
the name selects fragments but no chain ends up carrying that exact name, the
**longest** chain is used rather than failing. A segment the backbone visits
more than once keeps its first coordinate, though `walk` still lists it at
every visit - the chain flattened in path order is exactly what `walk` holds.

```c
// Both flags are needed: P lines name the fragments, L lines put them in order.
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

call_ref_t *r = call_ref_path(g, "chr22");   // NULL takes the graph's first path
if (!r) { gfa_destroy(g); return 1; }

printf("backbone %s, %lld bp\n", r->name, (long long)r->len);

// on[]/pos[] are indexed by segment index, so translate the file's id first
int32_t i = gfa_idx(g, 5);
if (i >= 0 && r->on[i])
    printf("segment 5 sits at offset %lld\n", (long long)r->pos[i]);

call_ref_destroy(r);
gfa_destroy(g);
```

### `call_ref_fasta`

```c
call_ref_t *call_ref_fasta(const gfa_t *g, const fasta_t *fa, const char *seq_name);
```

Builds the backbone by tracing an external sequence through the graph. A start
node is located (a source node whose sequence prefixes the record, falling back
to any matching node) and the walk greedily takes the out-link whose sequence
continues the record, honouring link overlaps. Requires `GFA_LINKS`; `P` lines
are not used at all.

`walk` records every node the trace visits, in order, so it can be written back
out as a `P` line (see [`gfa_add_path`](gfa.md#gfa_add_path)); a node the walk
returns to appears once per visit, while `on`/`pos` keep only its first.

A walk that stalls is logged at WARN and the backbone is **truncated there** -
`r->len` and `r->seq` are cut back to the covered prefix, so nothing past that
point yields variants. Compare `r->len` against the FASTA record's length if
that matters to you.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS);   // no P lines involved
if (!g) return 1;
fasta_t *fa = fasta_read("chr22.fa");
if (!fa) { gfa_destroy(g); return 1; }

call_ref_t *r = call_ref_fasta(g, fa, "chr22");   // NULL takes the first record
if (!r) { fasta_destroy(fa); gfa_destroy(g); return 1; }

// A stalled walk still returns a backbone; it is just shorter than the record.
const fasta_rec_t *fr = fasta_get(fa, "chr22");
if (fr && r->len < fr->len)
    ak_log(AK_LOG_WARN, NULL, "traced %lld of %lld bp",
           (long long)r->len, (long long)fr->len);

call_ref_destroy(r);
fasta_destroy(fa);
gfa_destroy(g);
```

### `call_ref_destroy`

```c
void call_ref_destroy(call_ref_t *r);
```

Releases a backbone and everything it owns - the name, the sequence, both
per-segment arrays and the walk. Safe with `NULL`, so it works on every error
path. It holds no reference to the graph, so destroy order between the two does
not matter.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;

call_ref_t *r = call_ref_path(g, NULL);
if (r) {
    // r->seq is NUL-terminated, so it is safe to hand to string functions
    printf("%zu bp of backbone sequence\n", strlen(r->seq));
    call_ref_destroy(r);   // independent of the graph's lifetime
}

gfa_destroy(g);
```

## Variants

### `call_variants`

```c
call_t *call_variants(const gfa_t *g, const call_ref_t *ref);
```

Walks every backbone node in reference order and collects what the graph
carries beside it. A link into an off-backbone node starts a bounded walk that
ends where it rejoins, and the bases it spells become an alternate allele; a
link that lands *further along* the backbone than the node it left is a plain
deletion of the span it skips. Alleles sharing a REF span are merged into one
multi-allelic record, duplicates are dropped, and the result is sorted by
`(pos, end)`.

Requires `GFA_LINKS`, and requires `ref` to come from the same graph -
`ref->n_seg != gfa_n_seg(g)` is refused outright. Returns `NULL` on failure,
with the reason logged.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
call_ref_t *r = call_ref_path(g, NULL);
if (!r) { gfa_destroy(g); return 1; }

// The backbone must have been built from this same graph: indices are shared.
call_t *c = call_variants(g, r);
if (!c) { call_ref_destroy(r); gfa_destroy(g); return 1; }

// Records are already ordered by (pos, end), the order a VCF wants.
printf("%lld record(s) over %s\n", (long long)call_n(c), r->name);

call_destroy(c);
call_ref_destroy(r);
gfa_destroy(g);
```

### `call_destroy`

```c
void call_destroy(call_t *c);
```

Releases a variant set together with the `ref`, `alt` and `type` strings of
every record. Safe with `NULL`. Any `call_var_t` pointer you kept is dangling
afterwards, so copy out what you need first.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
call_ref_t *r = call_ref_path(g, NULL);
call_t *c = r ? call_variants(g, r) : NULL;

if (c) {
    // Copy before destroying: the record's strings die with the set.
    char *first = call_n(c) > 0 ? strdup(call_at(c, 0)->alt) : NULL;
    call_destroy(c);
    free(first);
}

call_ref_destroy(r);   // NULL-safe, so no extra branch is needed
gfa_destroy(g);
```

### `call_n`

```c
int64_t call_n(const call_t *c);   // static inline in the header
```

Number of records in the set. Inline and unchecked, so `c` must not be `NULL`.
A set with zero records is a normal result, not an error - it means the graph
holds nothing the backbone does not already spell.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
call_ref_t *r = call_ref_path(g, NULL);
if (!r) { gfa_destroy(g); return 1; }
call_t *c = call_variants(g, r);
if (!c) { call_ref_destroy(r); gfa_destroy(g); return 1; }

// Zero is a legitimate answer; call_variants() only returns NULL on failure.
if (call_n(c) == 0)
    ak_log(AK_LOG_INFO, NULL, "%s carries no variation", r->name);

call_destroy(c);
call_ref_destroy(r);
gfa_destroy(g);
```

### `call_at`

```c
const call_var_t *call_at(const call_t *c, int64_t i);   // static inline in the header
```

The record at array index `i`, borrowed - the set still owns its strings. No
bounds checking, so keep `i` below `call_n()`.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
call_ref_t *r = call_ref_path(g, NULL);
if (!r) { gfa_destroy(g); return 1; }
call_t *c = call_variants(g, r);
if (!c) { call_ref_destroy(r); gfa_destroy(g); return 1; }

for (int64_t i = 0; i < call_n(c); i++) {
    const call_var_t *v = call_at(c, i);
    // pos is 0-based; alt and type are comma-separated and line up one to one
    printf("%s\t%lld\t%s\t%s\t%s\n", r->name, (long long)v->pos,
           v->ref, v->alt, v->type);
}

call_destroy(c);
call_ref_destroy(r);
gfa_destroy(g);
```

## Output

### `call_write_vcf`

```c
int call_write_vcf(const call_t *c, const call_ref_t *ref, const char *fn);
```

Writes the set as VCF 4.2: one row per record, POS converted to 1-based,
`ID`/`QUAL`/`FILTER` left as `.`, and INFO carrying `END` and `TYPE`. The
backbone supplies the `CHROM` column and the `##contig` header line, so `ref`
must be the one the variants were called against.

Returns `AK_OK`, `AK_EOPEN` if the file cannot be created, or `AK_EIO` if a
write or the final `fclose()` fails. The header is emitted even for an empty
set, so the output is always a valid VCF.

```c
gfa_t *g = gfa_read("graph.gfa", GFA_LINKS | GFA_PATHS);
if (!g) return 1;
call_ref_t *r = call_ref_path(g, NULL);
if (!r) { gfa_destroy(g); return 1; }
call_t *c = call_variants(g, r);
if (!c) { call_ref_destroy(r); gfa_destroy(g); return 1; }

// Pass the same backbone: it supplies CHROM and the ##contig length.
int rc = call_write_vcf(c, r, "out.vcf");
if (rc != AK_OK)
    ak_log(AK_LOG_ERROR, NULL, "%s", ak_strerror(rc));

call_destroy(c);
call_ref_destroy(r);
gfa_destroy(g);
```

---

[Back to the library index](README.md)
