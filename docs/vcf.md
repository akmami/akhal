# `vcf` - minimal streaming VCF reader

Source: [`src/lib/vcf.c`](../src/lib/vcf.c) &middot; Header: [`include/akhal/vcf.h`](../include/akhal/vcf.h)

A reader for the five VCF columns the toolkit actually uses. `vcf_open()` opens
a file and `vcf_read1()` pulls one data record at a time into a caller-owned
`vcf_rec_t`, so memory stays flat no matter how large the variant file is. The
shape is deliberately the same as the GAF streaming reader.

Header lines and blank lines are consumed by the reader and never surface.
Malformed data lines are logged through `ak_log()` and skipped, so a read only
fails on an allocation failure.

```c
#include "akhal/vcf.h"
#include "akhal/error.h"   // AK_E* return codes, ak_strerror(), ak_log()
```

## Contents

- [The record](#the-record) - fields, and which columns are parsed
- [Reuse contract](#reuse-contract) - who owns the strings, and for how long
- [Record lifetime](#record-lifetime) - [`vcf_rec_init`](#vcf_rec_init), [`vcf_rec_clear`](#vcf_rec_clear)
- [Streaming](#streaming) - [`vcf_open`](#vcf_open), [`vcf_read1`](#vcf_read1), [`vcf_close`](#vcf_close)

## The record

`vcf_rec_t` holds one data line. Every `char *` in it is owned by the record.

| Field | Type | Column | Notes |
| --- | --- | --- | --- |
| `chrom` | `char *` | CHROM | reference/contig name, copied verbatim |
| `pos` | `int64_t` | POS | 1-based position of the first REF base |
| `id` | `char *` | ID | `NULL` when the column was `.` |
| `ref` | `char *` | REF | reference allele, copied verbatim |
| `alt` | `char *` | ALT | comma-separated alternates; `NULL` when the column was `.` |

Columns 6 onwards - QUAL, FILTER, INFO, FORMAT and every sample column - are
not parsed and are not retained. A line is accepted on its first five
tab-separated fields alone; anything after them is ignored.

Only `id` and `alt` collapse a lone `.` to `NULL`. `chrom` and `ref` are
duplicated as written, so a `.` there arrives as the one-character string
`"."`. A line with fewer than five fields, or with a POS that does not parse to
a positive integer, is logged at `AK_LOG_WARN` and skipped.

## Reuse contract

One record serves the whole file. `vcf_read1()` clears the destination before
it fills it, which frees the strings from the previous call, so:

- Nothing may hold a pointer into the record across a `vcf_read1()` call. Copy
  `chrom`, `id`, `ref` or `alt` if you need them to outlive the iteration.
- The reader never frees the record itself. After the loop the last record read
  still owns its strings, so the caller must call `vcf_rec_clear()` once at the
  end. `vcf_close()` does not do it.
- A fresh record must be zeroed with `vcf_rec_init()` before its first use, or
  the first `vcf_read1()` will call `free()` on uninitialized pointers.

## Record lifetime

### `vcf_rec_init`

```c
void vcf_rec_init(vcf_rec_t *r);
```

Zeroes a record. It allocates nothing and must be called exactly once, before
the record is first passed to `vcf_read1()`. `vcf_rec_clear()` leaves the record
zeroed too, so a cleared record is immediately reusable without re-initializing.

```c
vcf_rec_t rec;
vcf_rec_init(&rec);          // zeroes the struct; nothing is allocated yet

vcf_reader_t *r = vcf_open("calls.vcf");
if (!r) return 1;

int rc;
while ((rc = vcf_read1(r, &rec)) == 1)
    printf("%s\t%lld\n", rec.chrom, (long long)rec.pos);

vcf_close(r);
vcf_rec_clear(&rec);         // the last record read still owns its strings
if (rc < 0) return 1;
```

### `vcf_rec_clear`

```c
void vcf_rec_clear(vcf_rec_t *r);
```

Frees `chrom`, `id`, `ref` and `alt`, then zeroes the struct. `free(NULL)` is
harmless, so it is safe on a record that was only initialized. Call it once
after the read loop; calling it twice in a row is also safe.

```c
vcf_rec_t rec;
vcf_rec_init(&rec);
vcf_reader_t *r = vcf_open("calls.vcf");
if (!r) return 1;

long n_snv = 0;
while (vcf_read1(r, &rec) == 1) {
    // alt is NULL when the ALT column was '.', so guard before touching it
    if (rec.alt && strlen(rec.ref) == 1 && strlen(rec.alt) == 1) n_snv++;
}
printf("%ld single-base substitution(s)\n", n_snv);

vcf_close(r);
vcf_rec_clear(&rec);         // frees whatever the last accepted line left behind
```

## Streaming

### `vcf_open`

```c
vcf_reader_t *vcf_open(const char *fn);
```

Opens a `.vcf` file for streaming. Returns `NULL` when the file cannot be opened
or the reader cannot be allocated; the reason is already logged, so callers can
just bail out. `vcf_reader_t` is opaque - the file handle, the line buffer and
the line counter used in diagnostics are all private.

```c
vcf_reader_t *r = vcf_open("calls.vcf");
if (!r) return 1;            // the reason has already been logged

vcf_rec_t rec;
vcf_rec_init(&rec);

// The '##'/'#CHROM' header is consumed inside the reader, so the first
// successful read already yields a data record.
if (vcf_read1(r, &rec) == 1)
    printf("first variant: %s:%lld %s -> %s\n", rec.chrom, (long long)rec.pos,
           rec.ref, rec.alt ? rec.alt : ".");

vcf_rec_clear(&rec);
vcf_close(r);
```

### `vcf_read1`

```c
int vcf_read1(vcf_reader_t *r, vcf_rec_t *rec);
```

Clears `rec` and refills it from the next usable data line. Returns 1 when a
record was parsed, 0 at end of file, or a negative `AK_E*` code - in practice
`AK_ENOMEM` - when an allocation fails. Header and malformed lines are skipped
internally, so 0 means the file is exhausted, never that a line was rejected.

The two terminating cases are told apart by the value that ended the loop: 0 is
a clean end of file, anything negative is a real failure and can be rendered
with `ak_strerror()`. On `AK_ENOMEM` the record may be left partly filled, so
clear it as usual.

```c
vcf_reader_t *r = vcf_open("calls.vcf");
if (!r) return 1;

vcf_rec_t rec;
vcf_rec_init(&rec);

int rc;
while ((rc = vcf_read1(r, &rec)) == 1) {
    // rec is refilled in place each iteration: copy anything you want to keep
    printf("%s\t%lld\t%s\t%s\t%s\n", rec.chrom, (long long)rec.pos,
           rec.id ? rec.id : ".", rec.ref, rec.alt ? rec.alt : ".");
}
// rc == 0 is a clean end of file; rc < 0 is the only real failure
if (rc < 0) ak_log(AK_LOG_ERROR, "vcf", "%s", ak_strerror(rc));

vcf_rec_clear(&rec);
vcf_close(r);
```

### `vcf_close`

```c
void vcf_close(vcf_reader_t *r);
```

Closes the file, releases the reader's line buffer and frees the reader. Safe
with `NULL`, so it works on every error path. It deliberately does not touch any
record: records are caller-owned and are cleared separately.

```c
vcf_reader_t *r = vcf_open("calls.vcf");
if (!r) return 1;

vcf_rec_t rec;
vcf_rec_init(&rec);

long n = 0;
while (vcf_read1(r, &rec) == 1) n++;

vcf_close(r);                // file and line buffer go away here
vcf_rec_clear(&rec);         // the record is separate, and still yours
printf("%ld record(s)\n", n);
```

---

[Back to the library index](README.md)
