# `fasta` - whole-file FASTA store indexed by name

Source: [`src/lib/fasta.c`](../src/lib/fasta.c) &middot; Header: [`include/akhal/fasta.h`](../include/akhal/fasta.h)

A FASTA file loaded into memory in one call. Records live in a contiguous
array and a hash table maps each sequence name to its array index, the same
"array + dict" layout the graph uses. This backs both the read store used by
`gaf2sam` and the reference used by `sampoke`.

There is no streaming interface and no random-access index on disk: the entire
file is resident for as long as the `fasta_t` lives, which is what makes the
by-name lookup O(1) and lets callers hand out borrowed `seq` pointers.

```c
#include "akhal/fasta.h"
#include "akhal/error.h"   // AK_E* return codes, ak_log()
```

## Contents

- [Layout](#layout) - the record array and the name index
- [Loading and releasing](#loading-and-releasing) - [`fasta_read`](#fasta_read), [`fasta_destroy`](#fasta_destroy)
- [Lookup and iteration](#lookup-and-iteration) - [`fasta_get`](#fasta_get), [`fasta_n`](#fasta_n)

## Layout

`fasta_rec_t` is one sequence:

| Field | Type | Notes |
| --- | --- | --- |
| `name` | `char *` | owned; the header text up to the first space or tab |
| `seq` | `char *` | owned; the sequence lines concatenated, NUL-terminated, case preserved |
| `len` | `int64_t` | length of `seq` in bytes |

`fasta_t` is the set:

| Field | Type | Notes |
| --- | --- | --- |
| `rec` | `fasta_rec_t *` | the records, in file order |
| `n` | `int64_t` | number of records in use |
| `m` | `int64_t` | allocated capacity; an implementation detail |
| `idx` | `void *` | opaque name to index hash table |

`rec` and everything it points at is public and may be read directly; `idx` is
opaque and is only reachable through `fasta_get()`. Both `name` and `seq` are
borrowed by the caller - they stay valid until `fasta_destroy()` and must not be
freed individually.

The name is the *first whitespace-delimited token* of the header, so a header
line of `>chr1 dna:chromosome GRCh38` is stored and indexed as `chr1`; the rest
of the description is discarded. Sequence lines are concatenated with their
newlines stripped and no case normalisation, so soft-masked lowercase bases
survive. Any sequence line appearing before the first `>` header is ignored.

Duplicate names are accepted: every record is kept in `rec`, a warning is
logged, and the index keeps pointing at the *first* one. `fasta_get()` will
therefore never return the later duplicates - reach them by scanning `rec`.

## Loading and releasing

### `fasta_read`

```c
fasta_t *fasta_read(const char *fn);
```

Reads the whole file into a freshly allocated set and builds the name index.
Returns `NULL` only when the file cannot be opened or an allocation fails; both
are logged. A file with no records is not an error - it yields a valid set with
`fasta_n() == 0`.

```c
fasta_t *fa = fasta_read("reads.fa");
if (!fa) return 1;   // unreadable file or out of memory; already logged

// The whole file is resident now: each record's seq is a private,
// NUL-terminated copy of that record's concatenated sequence lines.
printf("%lld sequence(s)\n", (long long)fasta_n(fa));

const fasta_rec_t *r = fasta_get(fa, "read1");
if (r) printf("read1 is %lld bp\n", (long long)r->len);

fasta_destroy(fa);
```

### `fasta_destroy`

```c
void fasta_destroy(fasta_t *fa);
```

Frees every record's name and sequence, the record array, the hash table and
the set itself. Safe to call with `NULL`, which makes it usable on every error
path. Every borrowed `name` and `seq` pointer dangles afterwards.

```c
fasta_t *fa = fasta_read("reads.fa");
if (!fa) return 1;

FILE *out = fopen("lengths.tsv", "w");
if (!out) {
    fasta_destroy(fa);   // release the store before bailing out
    return 1;
}

for (int64_t i = 0; i < fasta_n(fa); i++)
    fprintf(out, "%s\t%lld\n", fa->rec[i].name, (long long)fa->rec[i].len);

fclose(out);
fasta_destroy(fa);
```

## Lookup and iteration

### `fasta_get`

```c
const fasta_rec_t *fasta_get(const fasta_t *fa, const char *name);
```

Looks a record up by name through the hash table, in O(1). Returns `NULL` when
the name is absent. The returned record is borrowed and stays valid until
`fasta_destroy()`. `fa` must be non-`NULL`; the lookup dereferences it.

```c
fasta_t *fa = fasta_read("ref.fa");
if (!fa) return 1;

// The key is the first token of the header, so ">chr1 dna:chromosome" is
// found as "chr1" - passing the full header line would not match.
const fasta_rec_t *r = fasta_get(fa, "chr1");
if (!r) { fasta_destroy(fa); return 1; }

// seq is borrowed and NUL-terminated; len avoids a strlen() over the whole
// chromosome. Case is exactly as the file wrote it.
int64_t n = r->len < 60 ? r->len : 60;
fwrite(r->seq, 1, (size_t)n, stdout);
putchar('\n');

fasta_destroy(fa);
```

### `fasta_n`

```c
int64_t fasta_n(const fasta_t *fa);
```

The number of loaded records, as an inline accessor over `fa->n`. It is the
only valid bound for a loop over `fa->rec`; `fa->m` is the allocated capacity
and is usually larger.

```c
fasta_t *fa = fasta_read("reads.fa");
if (!fa) return 1;

// Iterating the array visits records in file order, including any duplicate
// names that fasta_get() cannot reach because the index kept the first one.
int64_t total = 0;
for (int64_t i = 0; i < fasta_n(fa); i++)
    total += fa->rec[i].len;

printf("%lld record(s), %lld bp total\n",
       (long long)fasta_n(fa), (long long)total);

fasta_destroy(fa);
```

---

[Back to the library index](README.md)
