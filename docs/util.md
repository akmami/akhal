# `util` - sequence, string and statistics helpers

Source: [`src/lib/util.c`](../src/lib/util.c) &middot; Header: [`include/akhal/util.h`](../include/akhal/util.h)

Small, dependency-free helpers shared by the library and the CLI. Nothing here
allocates, nothing here logs, and nothing here can fail: every function either
returns a value or writes through a buffer you already own. That is deliberate,
so these can be called from the middle of a parser without an error path.

Four groups: a pair of DNA routines used wherever a reverse strand has to be
materialized, two string helpers, number formatting with thousands separators,
and the running-distribution accumulator both the graph and the alignment
`stats` report with.

```c
#include "akhal/util.h"   // ak_complement, ak_revcomp, ak_ends_with, ak_str2int, ak_format_*, ak_dist_t
```

## Contents

- [Sequences](#sequences) - [`ak_complement`](#ak_complement), [`ak_revcomp`](#ak_revcomp)
- [Strings](#strings) - [`ak_ends_with`](#ak_ends_with), [`ak_str2int`](#ak_str2int)
- [Number formatting](#number-formatting) - [`ak_format_u64`](#ak_format_u64), [`ak_format_i64`](#ak_format_i64), [`ak_format_f64`](#ak_format_f64), [`ak_format_bytes`](#ak_format_bytes)
- [Process measurements](#process-measurements) - [`ak_realtime`](#ak_realtime), [`ak_peak_rss`](#ak_peak_rss), [`ak_rss`](#ak_rss)
- [Summary statistics](#summary-statistics) - [`ak_dist_t`](#ak_dist_t), [`ak_dist_add`](#ak_dist_add), [`ak_dist_variance`](#ak_dist_variance), [`ak_dist_sd`](#ak_dist_sd)

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

## Number formatting

Thousands separators for output a person reads: `543020649` is a figure,
`543,020,649` is a number. All three write into a buffer the caller owns and
return it, so a call can sit inline in a `printf()` argument list; `AK_NUM_LEN`
(48) is enough for any value. A single buffer serves one call at a time - two
formatted numbers in the same `printf()` need two buffers, since the arguments
are all evaluated before anything is printed.

### `ak_format_u64`

```c
char *ak_format_u64(char *buf, uint64_t v);
```

`543020649` -> `543,020,649`; `18446744073709551615` -> `18,446,744,073,709,551,615`.

### `ak_format_i64`

```c
char *ak_format_i64(char *buf, int64_t v);
```

The same with a sign: `-1234567` -> `-1,234,567`. `INT64_MIN` is handled.

### `ak_format_f64`

```c
char *ak_format_f64(char *buf, double v, int prec);
```

Commas in the integer part, `prec` digits after the point (0..20), rounded the
way `printf("%.*f")` rounds: `12345.678` at 2 -> `12,345.68`, `999.999` at 2
-> `1,000.00`. NaN and infinities come out as `printf()` writes them.

```c
char b[AK_NUM_LEN], b2[AK_NUM_LEN];
printf("%s segments, %s bp\n", ak_format_i64(b, st.n_seg), ak_format_u64(b2, st.n_bp));
printf("mean length %s\n", ak_format_f64(b, st.seg_mean, 2));
```

### `ak_format_bytes`

```c
char *ak_format_bytes(char *buf, uint64_t bytes);
```

A byte count the way a person reads it, in powers of 1024 with one decimal:
`1536` -> `1.5 KB`, `2199023255552` -> `2.0 TB`. Below the first step up it
prints whole bytes, since a fraction of a byte is noise.

## Process measurements

Two readings for code that reports what it cost - what
[`gfa_read`](gfa.md#gfa_read) uses under `GFA_VERBOSE`. Neither can fail in a
way a caller must handle: where the value is unavailable they return zero.

### `ak_realtime`

```c
double ak_realtime(void);
```

Seconds from a monotonic clock, which no clock adjustment can move backwards.
The origin is arbitrary, so only the difference between two readings means
anything.

```c
double t0 = ak_realtime();
do_the_work();
printf("%.2f s\n", ak_realtime() - t0);
```

### `ak_peak_rss`

```c
size_t ak_peak_rss(void);
```

The high-water mark of this process's resident set, in bytes - not what is
held right now, so it never falls, and so it answers "how much did this need"
rather than "how much is left". `getrusage` reports it in kilobytes on Linux
and bytes on the BSDs; this returns bytes on both.

### `ak_rss`

```c
size_t ak_rss(void);
```

What the process holds right now, in bytes - `/proc/self/statm` on Linux,
`task_info` on macOS, `0` anywhere else. Unlike the peak this can fall, but
only as far as the allocator lets it: a freed block it keeps for reuse stays
resident, so this measures the process rather than any structure in it. To
attribute cost to a structure, count it - which is what
[`gfa_read`](gfa.md#gfa_read) does under `GFA_FOOTPRINT`.

## Summary statistics

One accumulator, `ak_dist_t`, fed one value at a time. It keeps the count, the
running mean, the sum of squared deviations (Welford's update, so a single
pass gives the variance without a second walk and without the cancellation a
naive sum of squares suffers on large values) and the extremes. Nothing is
collected into an array, which is what lets `gfa_read_stats()` and the GAF
side of `akhal stats` summarize files far larger than memory.

Values are `double`: the library feeds it counts and lengths, which are exact
up to 2^53, and the alignment stats feed it ratios.

### `ak_dist_t`

```c
typedef struct {
    int64_t n;           // values seen
    double  mean;        // running mean
    double  m2;          // sum of squared deviations from the running mean
    double  min, max;    // extremes; meaningless while n == 0
} ak_dist_t;
```

Zero-initialize it (`ak_dist_t d = {0};`) and it is ready. `mean` is always
current, so it can be read directly; `min` and `max` hold whatever the first
value was once `n > 0`, and are unspecified before that - check `n` first.

### `ak_dist_add`

```c
void ak_dist_add(ak_dist_t *d, double x);
```

Folds one value in. Order does not matter for the result beyond floating
point rounding.

```c
ak_dist_t lens = {0};
for (int32_t i = 0; i < g->n_seg; i++) ak_dist_add(&lens, (double)g->seg[i].len);

printf("%lld segments, mean %.2f, min %llu, max %llu\n",
       (long long)lens.n, lens.mean,
       (unsigned long long)lens.min, (unsigned long long)lens.max);
```

### `ak_dist_variance`

```c
double ak_dist_variance(const ak_dist_t *d);
```

**Population** variance: the sum of squared deviations divided by `n`, not by
`n - 1`. Returns `0.0` on an empty accumulator, and a single value has
variance `0.0` rather than a division by zero.

### `ak_dist_sd`

```c
double ak_dist_sd(const ak_dist_t *d);
```

`sqrt()` of `ak_dist_variance()`, so `0.0` on an empty accumulator and never
NaN.

```c
ak_dist_t depth = {0};
double vals[] = { 3, 5, 4, 9, 4, 6 };
for (size_t i = 0; i < 6; i++) ak_dist_add(&depth, vals[i]);

printf("depth %.2f +/- %.2f\n", depth.mean, ak_dist_sd(&depth));
```

---

[Back to the library index](README.md)
