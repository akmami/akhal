# `gaf` - GAF alignment records, readers and path traversal

Source: [`src/lib/gaf.c`](../src/lib/gaf.c) &middot; Header: [`include/akhal/gaf.h`](../include/akhal/gaf.h)

The reader for GAF (Graph Alignment Format), the format aligners such as
`minigraph`, `minimap2` and `GraphAligner` emit when the target is a graph
rather than a linear reference. A record carries the twelve mandatory columns
plus the handful of optional tags the toolkit uses, and the target is a *path
string* like `>12<7>3` rather than a single contig name.

Two consumption styles are offered over the same parser, streaming and batch,
and a small helper walks the path string one oriented node at a time.

```c
#include "akhal/gaf.h"
#include "akhal/error.h"   // AK_E* return codes, ak_strerror(), ak_log()
```

## Contents

- [The record](#the-record) - mandatory fields and optional tags
- [Streaming or batch](#streaming-or-batch) - which one to pick
- [Record lifetime](#record-lifetime) - [`gaf_rec_init`](#gaf_rec_init), [`gaf_rec_clear`](#gaf_rec_clear)
- [Streaming](#streaming) - [`gaf_open`](#gaf_open), [`gaf_read1`](#gaf_read1), [`gaf_close`](#gaf_close)
- [Batch](#batch) - [`gaf_slurp`](#gaf_slurp), [`gaf_free`](#gaf_free)
- [Path strings](#path-strings) - [`gaf_path_next`](#gaf_path_next)

## The record

The twelve mandatory columns, in file order:

| Field | Type | Notes |
| --- | --- | --- |
| `qname` | `char *` | owned; query/read name |
| `qlen` | `int64_t` | query length |
| `qstart` | `int64_t` | query start, 0-based |
| `qend` | `int64_t` | query end |
| `strand` | `char` | `'+'` or `'-'`; left as `'\0'` when the column was empty |
| `path` | `char *` | owned; path string such as `>1>2<3` |
| `plen` | `int64_t` | path length in bp |
| `pstart` | `int64_t` | path start |
| `pend` | `int64_t` | path end |
| `matches` | `int64_t` | matching bases |
| `block_len` | `int64_t` | alignment block length |
| `mapq` | `int` | mapping quality; 255 means unavailable |

Everything after column 12 is scanned as optional tags. Five are recognised and
the rest are ignored. Each recognised tag has a companion `has_*` flag which is
the only way to tell "the tag was absent" from "the tag was present and zero" -
the value fields are zeroed by `gaf_rec_clear()`, so a missing `NM` and an
`NM:i:0` both leave `nm == 0`.

| Tag | Field | Flag | Notes |
| --- | --- | --- | --- |
| `NM:i` | `nm` (`int`) | `has_nm` | mismatches + indels |
| `AS:f` | `as` (`double`) | `has_as` | alignment score; `AS:i` is also accepted |
| `dv:f` | `dv` (`double`) | `has_dv` | divergence |
| `id:f` | `id` (`double`) | `has_id` | identity |
| `cg:Z` | `cigar` (`char *`) | *none* | difference CIGAR; owned, `NULL` when absent |

`cigar` has no flag because `NULL` already says the tag was absent. A tag whose
type letter does not match the table - `dv:i`, say - is skipped as if it were
unrecognised, leaving the flag at 0.

## Streaming or batch

Both styles use the same line parser and produce identical records; they differ
only in who owns the memory and for how long.

- **Streaming** - `gaf_open()`, then `gaf_read1()` into one reusable record,
  then `gaf_close()`. One record's worth of memory is live at a time, so the
  file size does not matter, but each record is destroyed by the next read.
  Use it for a single pass that consumes each alignment as it arrives.
- **Batch** - `gaf_slurp()` returns the whole file as an array, released with
  `gaf_free()`. Every record is independently owned and stays valid, so the
  array can be sorted, indexed or walked repeatedly. Use it when a single
  forward pass is not enough. The cost is that the entire file is resident.

## Record lifetime

### `gaf_rec_init`

```c
void gaf_rec_init(gaf_rec_t *r);
```

Zeroes a record. It allocates nothing, and must be called before the record is
first handed to `gaf_read1()` - the first read clears its destination, so an
uninitialized record would have `free()` called on junk pointers.

```c
gaf_rec_t rec;
gaf_rec_init(&rec);          // zeroes the struct; nothing is allocated yet

gaf_reader_t *r = gaf_open("aln.gaf");
if (!r) return 1;

// One record is refilled for the whole file: nothing may keep a pointer into
// it across a read.
while (gaf_read1(r, &rec) == 1)
    printf("%s\t%s\n", rec.qname, rec.path);

gaf_close(r);
gaf_rec_clear(&rec);         // the last record read still owns its strings
```

### `gaf_rec_clear`

```c
void gaf_rec_clear(gaf_rec_t *r);
```

Frees `qname`, `path` and `cigar`, then zeroes the struct, which also resets
every `has_*` flag. Safe on a merely initialized record, and safe to call
twice. `gaf_close()` does not do it for you: the record is caller-owned, so
clear it once when the loop is over.

```c
gaf_rec_t rec;
gaf_rec_init(&rec);
gaf_reader_t *r = gaf_open("aln.gaf");
if (!r) return 1;

long n_scored = 0;
while (gaf_read1(r, &rec) == 1) {
    // has_as is the only way to distinguish "no AS tag" from "AS:f:0"
    if (rec.has_as && rec.as > 100.0) n_scored++;
}
printf("%ld alignment(s) scoring over 100\n", n_scored);

gaf_close(r);
gaf_rec_clear(&rec);         // frees whatever the last accepted line left behind
```

## Streaming

### `gaf_open`

```c
gaf_reader_t *gaf_open(const char *fn);
```

Opens a `.gaf` file for streaming. Returns `NULL` when the file cannot be opened
or the reader cannot be allocated; the reason is logged. `gaf_reader_t` is
opaque - it owns only the file handle and the reusable line buffer.

```c
gaf_reader_t *r = gaf_open("aln.gaf");
if (!r) return 1;            // the reason has already been logged

gaf_rec_t rec;
gaf_rec_init(&rec);

// GAF has no header, so the first successful read is the first alignment.
if (gaf_read1(r, &rec) == 1)
    printf("%s: %lld-%lld of %lld on %s\n", rec.qname,
           (long long)rec.qstart, (long long)rec.qend,
           (long long)rec.qlen, rec.path);

gaf_rec_clear(&rec);
gaf_close(r);
```

### `gaf_read1`

```c
int gaf_read1(gaf_reader_t *r, gaf_rec_t *rec);
```

Clears `rec` and refills it from the next well-formed line. Returns 1 on
success, 0 at end of file, `AK_EINVAL` if either argument is `NULL`, or
`AK_ENOMEM` on an allocation failure. Blank lines are skipped, and a line with
fewer than twelve fields is logged at `AK_LOG_WARN` and skipped as well, so a
return of 0 always means the file is exhausted.

Test the value that ended the loop to tell the cases apart: 0 is a clean end of
file and anything negative is a failure worth reporting through
`ak_strerror()`.

```c
gaf_reader_t *r = gaf_open("aln.gaf");
if (!r) return 1;

gaf_rec_t rec;
gaf_rec_init(&rec);

int rc;
while ((rc = gaf_read1(r, &rec)) == 1) {
    // mapq 255 is the format's "unavailable", not a high-quality mapping
    if (rec.mapq != 255 && rec.mapq >= 30)
        printf("%s\t%c\t%s\t%lld/%lld\n", rec.qname, rec.strand, rec.path,
               (long long)rec.matches, (long long)rec.block_len);
}
if (rc < 0) ak_log(AK_LOG_ERROR, "gaf", "%s", ak_strerror(rc));  // rc == 0 is EOF

gaf_rec_clear(&rec);
gaf_close(r);
```

### `gaf_close`

```c
void gaf_close(gaf_reader_t *r);
```

Closes the file, releases the line buffer and frees the reader. Safe with
`NULL`. It does not touch any record, so the destroy order between reader and
record does not matter.

```c
gaf_reader_t *r = gaf_open("aln.gaf");
if (!r) return 1;

gaf_rec_t rec;
gaf_rec_init(&rec);

int64_t bases = 0;
while (gaf_read1(r, &rec) == 1) bases += rec.matches;

gaf_close(r);                // file and line buffer go away here
gaf_rec_clear(&rec);         // the record is separate, and still yours
printf("%lld matching base(s)\n", (long long)bases);
```

## Batch

### `gaf_slurp`

```c
gaf_rec_t *gaf_slurp(const char *fn, int64_t *n);
```

Reads the whole file into one newly allocated array and writes the record count
through `n`. Each element owns its own strings, so the array survives
independently of any reader. Returns `NULL` on failure, with `*n` left at 0.

An empty file also returns `NULL` with `*n == 0`, so `NULL` alone does not mean
an error occurred - check `*n`, or rely on the diagnostic that a real failure
logs.

```c
int64_t n = 0;
gaf_rec_t *recs = gaf_slurp("aln.gaf", &n);
if (!recs) return 1;   // failure, or simply an empty file; n is 0 either way

// Unlike the streaming record, these stay valid together, so multiple passes
// and random access are fine.
int64_t mapped = 0;
for (int64_t i = 0; i < n; i++)
    if (recs[i].mapq != 255) mapped++;

printf("%lld/%lld record(s) carry a mapping quality\n",
       (long long)mapped, (long long)n);

gaf_free(recs, n);
```

### `gaf_free`

```c
void gaf_free(gaf_rec_t *recs, int64_t n);
```

Clears the first `n` records, then frees the array. Safe with `NULL`. Pass the
same `n` that `gaf_slurp()` reported: a smaller count leaks the tail's strings
and a larger one walks off the end.

```c
int64_t n = 0;
gaf_rec_t *recs = gaf_slurp("aln.gaf", &n);
if (!recs) return 1;

// A second pass over the same array; the batch reader exists for exactly this.
for (int64_t i = 0; i < n; i++)
    if (recs[i].cigar)   // cg:Z is absent unless the aligner emitted it
        printf("%s\t%s\n", recs[i].qname, recs[i].cigar);

gaf_free(recs, n);       // n must be the count gaf_slurp() reported
```

## Path strings

### `gaf_path_next`

```c
int gaf_path_next(const char *p, uint64_t *id, char *strand);
```

Decodes the oriented node at `p` in a GAF path string, writing the node id and
its orientation character (`>` forward, `<` reverse) and returning the number of
bytes it consumed. It returns 0 at the terminating NUL, which is what ends the
loop; a non-zero return is always at least 1, so the pointer always advances.

The parse is deliberately unchecked. Whatever byte sits at `p` becomes the
orientation, and a run of digits - possibly empty, in which case the id is 0 -
becomes the id. A stable path (`>`/`<` plus a decimal id) is decoded exactly;
anything else is decoded, not diagnosed.

```c
const char *p = ">12<7>3";   // as it appears in gaf_rec_t.path
uint64_t id;
char strand;
int n;

// The return value is the byte count, so the advance is p += n, not p++.
// The loop ends when gaf_path_next() reports 0 at the NUL terminator.
while ((n = gaf_path_next(p, &id, &strand)) > 0) {
    printf("node %llu %s\n", (unsigned long long)id,
           strand == '>' ? "forward" : "reverse");
    p += n;
}
```

---

[Back to the library index](README.md)
