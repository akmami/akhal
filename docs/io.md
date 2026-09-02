# `io` - line-oriented input for every text reader

Source: [`src/lib/io.c`](../src/lib/io.c) &middot; Header: [`include/akhal/io.h`](../include/akhal/io.h)

`ak_file` is an opaque handle around one open stream. It exists so that every
text format reader in the library - GFA, GAF, FASTA, SAM, VCF - pulls its lines
through one function, `ak_getline()`. A buffered or compressed backend can then
be slotted in behind the handle without touching a single parser.

The handle carries the stream; the *line buffer* belongs to the caller. You
declare one [`kstring_t`](kstr.md), hand its address to `ak_getline()` on every
iteration, and free it once when the loop ends. The buffer grows to fit the
longest line it has seen and is then reused, so a scan over a large file does a
handful of allocations in total rather than one per line.

```c
#include "akhal/io.h"     // ak_file, ak_open/ak_getline/ak_rewind/ak_close
#include "akhal/kstr.h"   // kstring_t, KS_INIT, ks_free() for the line buffer
#include "akhal/error.h"  // AK_OK / AK_EOF / AK_E* codes, ak_log(), ak_strerror()
```

## Contents

- [The line contract](#the-line-contract)
- [Opening and closing](#opening-and-closing) - [`ak_open`](#ak_open), [`ak_close`](#ak_close)
- [Reading](#reading) - [`ak_getline`](#ak_getline), [`ak_rewind`](#ak_rewind)

## The line contract

`ak_getline()` returns a `long`, and the whole contract hangs off its sign:

| Return | Meaning |
| --- | --- |
| `>= 0` | the length of the line, terminator excluded |
| `AK_EOF` | no more input |
| `AK_EINVAL` | the handle was `NULL` or already closed |
| `AK_ENOMEM` | the buffer could not be grown |

So the loop condition is `>= 0` and never a truth test: a blank line is a
perfectly good line of length 0. On success `ks->s` holds the line, `ks->l`
holds the same value the call returned, and a trailing `"\r\n"` or `"\n"` has
been removed. A final line with no newline at all is still returned normally;
the `AK_EOF` comes on the call after it.

Each call resets `ks->l` to 0 first, so lines never accumulate - if you need to
keep a line past the next iteration, copy it or take it with
[`ks_release()`](kstr.md#ks_release).

## Opening and closing

### `ak_open`

```c
ak_file *ak_open(const char *fn);
```

Opens `fn` for reading in text mode and wraps it in a handle. Returns `NULL` on
a `NULL` path, an unopenable file or an allocation failure; the reason is
already written to stderr through `ak_log()`, so a caller only has to bail out.
There is no stdin special case - pass a real path.

```c
const char *fn = "graph.gfa";

ak_file *f = ak_open(fn);
if (!f) return 1;   // the reason is already logged

kstring_t ks = KS_INIT;          // the caller owns the line buffer
long len = ak_getline(f, &ks);

if (len >= 0)
    printf("first line (%ld bytes): %s\n", len, ks.s);
else
    ak_log(AK_LOG_WARN, "io", "%s is empty", fn);

ks_free(&ks);
ak_close(f);
```

### `ak_close`

```c
void ak_close(ak_file *f);
```

Closes the stream and frees the handle. Safe to call with `NULL`, which makes it
usable on every error path. It does not touch your `kstring_t`; that is a
separate `ks_free()`.

```c
ak_file *f = ak_open("graph.gfa");
if (!f) return 1;

kstring_t ks = KS_INIT;
long len = ak_getline(f, &ks);

// AK_EOF just means the file was empty; anything else negative is a failure
// worth reporting before unwinding.
if (len < 0 && len != AK_EOF)
    ak_log(AK_LOG_ERROR, "io", "%s", ak_strerror((int)len));

ks_free(&ks);   // buffer first, handle second - they are independent
ak_close(f);
return len < 0 && len != AK_EOF ? 1 : 0;
```

## Reading

### `ak_getline`

```c
long ak_getline(ak_file *f, kstring_t *ks);
```

Reads the next line into `ks`, stripping the trailing `'\n'` and a preceding
`'\r'`, and returns its length. `ks` must be a valid `kstring_t` - `KS_INIT` is
enough - and is grown as needed; the same buffer is meant to be passed to every
call. Returns `AK_EOF` once the stream is exhausted.

```c
ak_file *f = ak_open("graph.gfa");
if (!f) return 1;

kstring_t ks = KS_INIT;   // one buffer, reused for every line
long len, n = 0;

// >= 0, not a truth test: a blank line legitimately returns 0. Every negative
// value leaves the loop - AK_EOF at the end of input, AK_ENOMEM if the buffer
// could not grow.
while ((len = ak_getline(f, &ks)) >= 0) {
    if (len == 0 || ks.s[0] == '#') continue;   // skip blanks and comments
    n++;
}

printf("%lld data line(s)\n", (long long)n);
ks_free(&ks);   // once, after the loop - not inside it
ak_close(f);
```

### `ak_rewind`

```c
int ak_rewind(ak_file *f);
```

Puts the stream back at byte 0 so the same handle can be read again, which is
how two-pass readers size their arrays before filling them. Returns `AK_OK`, or
`AK_EINVAL` for a `NULL` handle. Like `rewind()` it reports nothing else, so a
seek failure on an unseekable stream is not surfaced.

```c
ak_file *f = ak_open("graph.gfa");
if (!f) return 1;

kstring_t ks = KS_INIT;
size_t n_seg = 0;

while (ak_getline(f, &ks) >= 0)          // pass one: count the S lines
    if (ks.l > 0 && ks.s[0] == 'S') n_seg++;

if (ak_rewind(f) != AK_OK) {             // back to the top for pass two
    ks_free(&ks);
    ak_close(f);
    return 1;
}

while (ak_getline(f, &ks) >= 0) { /* pass two: parse into the sized arrays */ }

printf("%lu segment line(s)\n", (unsigned long)n_seg);
ks_free(&ks);
ak_close(f);
```

---

[Back to the library index](README.md)
