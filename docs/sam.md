# `sam` - CIGAR alphabet and SAM output

Source: [`src/lib/sam.c`](../src/lib/sam.c) &middot; Header: [`include/akhal/sam.h`](../include/akhal/sam.h)

The SAM-side conventions used by the toolkit: the CIGAR operation alphabet,
expansion of a CIGAR string into a per-base op array, and writing a SAM header
and alignment records. Commands assemble the fields and hand them here.

The working representation is the *per-base op array*, not the run-length
string. `sam_cigar_expand()` turns `3=1I2X` into `===IXX`; a command walking a
graph path appends and rewrites single ops as it goes; `sam_write_record()`
collapses the array back into runs on the way out. Nothing here parses SAM - the
module only writes it.

```c
#include "akhal/sam.h"   // ops, flags, the record struct and the two writers
```

## Contents

- [Operations, predicates and flags](#operations-predicates-and-flags) - the CIGAR alphabet, `CIGAR_QUE`/`CIGAR_REF`, the `SAM_F*` bits, `SAM_MAX_CIGAR`
- [CIGAR expansion](#cigar-expansion) - [`sam_cigar_expand`](#sam_cigar_expand)
- [The alignment record](#the-alignment-record) - the `sam_rec_t` fields and the in-place rewrite
- [Writing](#writing) - [`sam_write_header`](#sam_write_header), [`sam_write_record`](#sam_write_record)
- [Read groups](#read-groups) - [`sam_rg_prefix`](#sam_rg_prefix)

## Operations, predicates and flags

The op alphabet, with the two properties that matter when you walk an
alignment. "Query" is the read; "reference" is what `pos` indexes into.

| Macro | Char | Consumes query | Consumes reference |
| --- | --- | --- | --- |
| `CIGAR_ALIGNMENT_MATCH` | `M` | yes | yes |
| `CIGAR_INSERTION` | `I` | yes | no |
| `CIGAR_DELETION` | `D` | no | yes |
| `CIGAR_SKIPPED` | `N` | no | yes |
| `CIGAR_SOFT_CLIP` | `S` | yes | no |
| `CIGAR_HARD_CLIP` | `H` | no | no |
| `CIGAR_PADDING` | `P` | no | no |
| `CIGAR_SEQUENCE_MATCH` | `=` | yes | yes |
| `CIGAR_SEQUENCE_MISMATCH` | `X` | yes | yes |

Two predicates read that table for you:

| Macro | True for |
| --- | --- |
| `CIGAR_QUE(x)` | `M`, `I`, `S`, `=`, `X` |
| `CIGAR_REF(x)` | `M`, `D`, `N`, `=`, `X` |

Both are function-like macros that expand their argument up to five times, so
do not pass an expression with side effects such as `CIGAR_QUE(ops[i++])`.
`CIGAR_QUE()` is false for `H` and `CIGAR_REF()` is false for `P`, so a loop
that advances a query cursor on one and a reference cursor on the other covers
the whole alphabet without a `default` case.

The FLAG bits the module defines:

| Macro | Value | Meaning |
| --- | --- | --- |
| `SAM_FMUNMAP` | `0x4` | unmapped |
| `SAM_FREVERSE` | `0x10` | reverse strand |
| `SAM_FSECONDARY` | `0x100` | secondary alignment |
| `SAM_FSUPPLEMENTARY` | `0x800` | supplementary alignment |

`SAM_FMUNMAP` is `0x4`, which is the SAM spec's *segment unmapped* bit, not the
*next segment unmapped* bit (`0x8`) that the `M` in the name suggests. Only
`SAM_FREVERSE` is set anywhere in the toolkit today; the rest are here for
callers.

`SAM_MAX_CIGAR` (`1048576`) is the upper bound on both a per-base op array and
the CIGAR string rendered from one. `sam_write_record()` renders into a
`SAM_MAX_CIGAR`-byte buffer *without* a bounds check, and a fully alternating op
array costs two characters per op, so an `n_ops` above half of `SAM_MAX_CIGAR`
is only safe when the ops come in runs. The commands size their scratch at
exactly `SAM_MAX_CIGAR` and refuse any alignment that fills it.

## CIGAR expansion

### `sam_cigar_expand`

```c
int sam_cigar_expand(const char *cigar, char *ops, int max_ops, int rev);
```

Expands a run-length CIGAR into one op character per base, reversing the result
when `rev` is non-zero (for an alignment reported on the opposite strand).
Returns the number of ops written, or `-1` when they do not fit. It does not
NUL-terminate `ops`, and it does not validate the alphabet: every non-digit is
copied through as an op character, so junk in is junk out.

The capacity test is `j + num >= max_ops`, so the buffer needs one byte of slack
- a CIGAR that would exactly fill `max_ops` is reported as an overflow.

```c
// The commands allocate this once and reuse it across every record.
char *ops = (char *)malloc(SAM_MAX_CIGAR);
if (!ops) return 1;

int n = sam_cigar_expand("3=1I2X", ops, SAM_MAX_CIGAR, 0);
if (n < 0) { free(ops); return 1; }   // did not fit
ops[n] = '\0';                        // expansion leaves the buffer unterminated

int ref = 0, que = 0;
for (int i = 0; i < n; i++) {
    if (CIGAR_REF(ops[i])) ref++;     // bases spanned on the reference
    if (CIGAR_QUE(ops[i])) que++;     // bases consumed from the read
}
printf("%s: %d ops, %d ref bp, %d query bp\n", ops, n, ref, que);

free(ops);
```

## The alignment record

`sam_rec_t` is a plain field bundle. `sam_write_record()` takes a non-const
pointer because it edits `ops` in place.

| Field | Type | Notes |
| --- | --- | --- |
| `qname` | `const char *` | QNAME, written verbatim |
| `flag` | `int` | FLAG, an OR of the `SAM_F*` bits |
| `rname` | `const char *` | RNAME, written verbatim |
| `pos` | `int` | POS, **1-based** |
| `mapq` | `int` | MAPQ |
| `ops` | `char *` | per-base op array, **mutated in place** (see below) |
| `n_ops` | `int` | length of `ops`; `0` renders the CIGAR as `*` |
| `simplify` | `int` | non-zero renders `=` and `X` as `M` |
| `seq` | `const char *` | SEQ; `NULL` renders as `*` |
| `nm` | `int` | `NM:i` |
| `as` | `double` | `AS:f`, `%.2f` |
| `dv` | `double` | `dv:f`, `%.6f` |
| `id` | `double` | `id:f`, `%.6f` |
| `rg` | `const char *` | `RG:Z`; `NULL` falls back to `akhal.0` |

The struct has no initialiser, so `memset()` it before filling it in. RNEXT is
always `*`, PNEXT and TLEN are always `0`, and QUAL is always `*`. The four
optional tags and `RG:Z` are always emitted, including when their values are
zero.

Two rewrites happen to `ops`, in this order, before the runs are collapsed:

1. **Clip folding.** From each end the writer skips any existing `S`, then turns
   the run of `I` and `X` that follows into `S`. This writes to the caller's
   buffer, and it happens whether or not `simplify` is set. Pass scratch you own
   rather than a shared or string-literal array, and re-expand the CIGAR if you
   still need the original ops afterwards.
2. **Simplification.** With `simplify` set, every remaining `=` and `X` becomes
   `M`. Ops already turned into `S` by step 1 stay `S`.

Only `I` and `X` fold into clips - a leading `D` or `N` is left alone, so a
record whose ops start on the reference still reports that op.

## Writing

### `sam_write_header`

```c
void sam_write_header(FILE *out, char **names, int n, const uint64_t *lens, const char *pg);
```

Writes `@HD VN:1.6 SO:unsorted GO:query`, one `@SQ` per reference sequence, a
`@PG`, and a single default `@RG` with id `<pg>.0` and `PL`/`PU`/`SM` all
`UNKNOWN`. `@SQ` lines come out in canonical chromosome order: a leading `chr`
is ignored, autosomes `1..22` rank first, then `X`, `Y`, `M`/`MT`, and every
other contig ranks last and is ordered by name. The caller's arrays are never
reordered - the sort runs over a private index.

There is no return value and no error path. If the temporary index cannot be
allocated the function returns straight after `@HD`, leaving a header with no
`@SQ`, `@PG` or `@RG`.

```c
// names is char ** rather than const char **, but both arrays are only read.
char *names[3] = { "chr10", "scaffold_7", "chr2" };
uint64_t lens[3] = { 133797422ULL, 4000ULL, 242193529ULL };

sam_write_header(stdout, names, 3, lens, "akhal");
// @SQ order is chr2, chr10, scaffold_7: canonical contigs by number, everything
// else after them by name. names[] itself is left as you built it.
printf("names[0] is still %s\n", names[0]);

// The @RG id is "<pg>.0", and a record that leaves sam_rec_t::rg NULL writes
// RG:Z:akhal.0 - so any pg other than "akhal" needs an explicit rg per record.
uint64_t total = 0;
for (int i = 0; i < 3; i++) total += lens[i];
printf("%llu bp of reference declared\n", (unsigned long long)total);
```

### `sam_write_record`

```c
void sam_write_record(FILE *out, sam_rec_t *r);
```

Renders one alignment line and writes it, newline included, after applying the
two `ops` rewrites described above. It returns nothing and reports nothing, so a
failed write is only visible through `ferror()` on the stream.

```c
char ops[8];
memcpy(ops, "IX==X=I", 7);      // per-base ops, not a run-length CIGAR

sam_rec_t r;
memset(&r, 0, sizeof(r));       // no initialiser: zero it before filling it in
r.qname = "read1";
r.flag  = SAM_FREVERSE;
r.rname = "chr1";
r.pos   = 1001;                 // 1-based
r.mapq  = 60;
r.ops   = ops;                  // writable scratch: the call rewrites it
r.n_ops = 7;
r.simplify = 1;                 // '=' and 'X' render as 'M'
r.seq   = "ACGTACG";
sam_write_record(stdout, &r);   // CIGAR 2S4M1S

printf("ops afterwards: %.7s\n", ops);   // SSMMMMS, no longer IX==X=I
```

## Read groups

### `sam_rg_prefix`

```c
char *sam_rg_prefix(const char *name);
```

Derives a read-group id from a read name: the leading component before the first
`.`, `/` or `_`. A leading `@` or `>` is skipped first, so a raw SAM or FASTA
line can be passed instead of a bare name. **The result is newly allocated and
the caller frees it.** `NULL` is returned only on allocation failure; a name with
no separator yields a copy of the whole string.

The separators are not equal. `.` and `/` are searched first and the earlier of
the two wins; `_` is consulted only when the argument contains neither anywhere.
On a whole SAM line that means a dot inside a later field - in an `AS:f` tag,
say - is preferred over an underscore in the read name.

```c
char *rg = sam_rg_prefix("@HG002.chr1_1234");   // a leading '@' or '>' is skipped
if (!rg) return 1;
printf("read group: %s\n", rg);                 // HG002
free(rg);                                       // the string is yours

// '_' applies only when there is no '.' and no '/' anywhere in the argument,
// so handing over a whole SAM line can cut inside its optional tags instead.
char *from_line = sam_rg_prefix("s_1\t0\tchr1\t1\t60\t5M\t*\t0\t0\tACGTA\t*\tAS:f:1.00");
if (!from_line) return 1;
printf("cut after %d bytes\n", (int)strlen(from_line));   // at the dot in 1.00
free(from_line);
```

---

[Back to the library index](README.md)
