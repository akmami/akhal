# `util` - sequence, string and statistics helpers

Source: [`src/lib/util.c`](../src/lib/util.c) &middot; Header: [`include/akhal/util.h`](../include/akhal/util.h)

Small, dependency-free helpers shared by the library and the CLI. Nothing here
allocates, nothing here logs, and nothing here can fail: every function either
returns a value or writes through a buffer you already own. That is deliberate,
so these can be called from the middle of a parser without an error path.

Two groups: a pair of DNA routines used wherever a reverse strand has to be
materialized, and the mean/variance/stddev trio the `stats` command reports
with.

```c
#include "akhal/util.h"   // ak_complement, ak_revcomp, ak_ends_with, ak_str2int, ak_mean/ak_variance/ak_stddev
```

## Contents

- [Sequences](#sequences) - [`ak_complement`](#ak_complement), [`ak_revcomp`](#ak_revcomp)
- [Strings](#strings) - [`ak_ends_with`](#ak_ends_with), [`ak_str2int`](#ak_str2int)
- [Summary statistics](#summary-statistics) - [`ak_mean`](#ak_mean), [`ak_variance`](#ak_variance), [`ak_stddev`](#ak_stddev)

## Sequences

### `ak_complement`

```c
char ak_complement(char base);
```

Complements one base. The lookup masks the character with `0xDF`, so case does
not matter on input - but the result is always uppercase, which means
soft-masking is lost. Anything that is not A, C, G or T after that mask,
including IUPAC ambiguity codes and gap characters, comes back as `'N'`.

```c
// Case-insensitive in, uppercase out: 'c' complements to 'G', not 'g'.
char bases[] = "AcGtN-";

for (size_t i = 0; i < strlen(bases); i++)
    printf("%c -> %c\n", bases[i], ak_complement(bases[i]));

// Everything unrecognised collapses to 'N', so an ambiguity code is not
// preserved and cannot be round-tripped back.
if (ak_complement('R') != 'N') return 1;
if (ak_complement('-') != 'N') return 1;
```

### `ak_revcomp`

```c
void ak_revcomp(char *seq, size_t len);
```

Reverse-complements the first `len` bytes of `seq` **in place**. There is no
return value and no copy: the buffer must be writable, so a string literal will
not do. `len` excludes any terminator and the function never touches `seq[len]`,
so a NUL-terminated buffer stays terminated. A `NULL` pointer or a `len` of 0 is
a no-op, and an odd length has its middle base complemented in place.

```c
// Writable storage: char seq[] copies the literal, const char *seq would not.
char seq[] = "AACCGGTT";
size_t len = strlen(seq);

ak_revcomp(seq, len);
printf("%s\n", seq);   // AACCGGTT - this one is its own reverse complement

// Any prefix or slice works the same way, because len is just a byte count.
ak_revcomp(seq, 4);
printf("%s\n", seq);

ak_revcomp(NULL, 10);  // no-op, safe to call on a segment with no sequence
ak_revcomp(seq, 0);
```

## Strings

### `ak_ends_with`

```c
int ak_ends_with(const char *str, const char *suffix);
```

Returns 1 when `str` ends with `suffix`, otherwise 0. The comparison is exact
and case-sensitive, so `".GFA"` does not match `".gfa"`. Both arguments are
passed straight to `strlen()`, so neither may be `NULL`. A suffix longer than
the string is simply not a match; the empty suffix always matches.

```c
const char *fn = "graph.rgfa";

// Used for format sniffing by extension - test each spelling you accept.
if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa"))
    printf("%s looks like a GFA file\n", fn);

// Case-sensitive, and an over-long suffix is a clean 0 rather than a read
// past the start of the string.
if (ak_ends_with(fn, ".RGFA")) return 1;
if (ak_ends_with("x", ".gfa")) return 1;
if (!ak_ends_with(fn, ""))     return 1;
```

### `ak_str2int`

```c
int ak_str2int(const char *str, int *out);
```

Returns 1 and writes the parsed value through `out` when `str` is a valid
base-10 int, otherwise returns 0 and leaves `out` untouched. The whole string
must be consumed - trailing junk, a bare sign, whitespace, and `""` are all
rejected the same as `NULL` - and a value outside `int` range is rejected as
overflow rather than silently truncated.

```c
int wrap_len;

// A clean base-10 int is accepted and written through out.
if (ak_str2int("120", &wrap_len))
    printf("wrap length %d\n", wrap_len);

// Trailing junk, empty input and out-of-range values are all a clean 0 -
// out is left untouched, so a caller-supplied default survives.
if (ak_str2int("120x", &wrap_len)) return 1;
if (ak_str2int("", &wrap_len))     return 1;
if (ak_str2int("99999999999999999999", &wrap_len)) return 1;
```

## Summary statistics

These three are meant to be used together, in order: `ak_mean()` over the array,
that mean into `ak_variance()`, that variance into `ak_stddev()`. The mean is
passed in rather than recomputed so a caller that already has it does not walk
the array twice.

All of them take `const size_t *`, not `double *` - they are written for the
counts and lengths the rest of the library deals in. `n == 0` is safe in both
array functions and yields `0.0`, so an empty graph needs no special case.

### `ak_mean`

```c
double ak_mean(const size_t *a, size_t n);
```

Arithmetic mean of `n` values, accumulated in `double`. Returns `0.0` when `n`
is 0, which is a stand-in for "undefined" rather than a real mean - check `n`
yourself if the distinction matters.

```c
size_t lens[] = { 100, 250, 175, 400, 75 };
size_t n = sizeof(lens) / sizeof(lens[0]);

// Values are size_t and are cast to double one at a time as they are summed.
double mean = ak_mean(lens, n);
printf("mean length %.2f over %lu value(s)\n", mean, (unsigned long)n);

// n == 0 short-circuits before touching the pointer, so this is safe even
// with a NULL array - but 0.0 is not distinguishable from a genuine mean of 0.
printf("empty: %.2f\n", ak_mean(NULL, 0));
```

### `ak_variance`

```c
double ak_variance(const size_t *a, size_t n, double mean);
```

**Population** variance: the sum of squared deviations divided by `n`, not by
`n - 1`. It trusts the `mean` you pass and never recomputes it, so passing the
mean of a different array silently produces a meaningless number. Returns `0.0`
when `n` is 0.

```c
size_t lens[] = { 100, 250, 175, 400, 75 };
size_t n = sizeof(lens) / sizeof(lens[0]);

// The mean must be the mean of this same array and this same n.
double mean = ak_mean(lens, n);
double var  = ak_variance(lens, n, mean);

printf("mean %.2f, variance %.2f\n", mean, var);

// Dividing by n rather than n - 1 means a single value has variance 0.0
// instead of a division by zero.
printf("one value: %.2f\n", ak_variance(lens, 1, ak_mean(lens, 1)));
```

### `ak_stddev`

```c
double ak_stddev(double variance);
```

`sqrt()` of the value it is given, and nothing else - it never sees the array.
It performs no domain check, so a negative argument would produce NaN;
`ak_variance()` never returns one.

```c
size_t depths[] = { 3, 5, 4, 9, 4, 6 };
size_t n = sizeof(depths) / sizeof(depths[0]);

// The three calls chain: array -> mean -> variance -> standard deviation.
double mean = ak_mean(depths, n);
double sd   = ak_stddev(ak_variance(depths, n, mean));

printf("depth %.2f +/- %.2f\n", mean, sd);

// Since it is only sqrt(), it is equally usable on a variance computed
// elsewhere - but feed it nothing negative.
printf("sd of 6.25 is %.2f\n", ak_stddev(6.25));
```

---

[Back to the library index](README.md)
