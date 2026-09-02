# `kstr` - growable string buffer

Source: [`src/lib/kstr.c`](../src/lib/kstr.c) &middot; Header: [`include/akhal/kstr.h`](../include/akhal/kstr.h)

A minimal growable byte buffer in the spirit of htslib's `kstring_t`. It is
three fields and no hidden state, which is the point: readers can hand it
around, grow it, and hand ownership of the finished string to a caller without
any allocator bookkeeping of their own.

Only `ks_resize()` has a real implementation, in `src/lib/kstr.c`. `ks_clear()`,
`ks_free()` and `ks_release()` are `static inline` in the header and compile
down to a few stores, so they cost nothing to call in a loop. `KS_INIT` is a
macro, not a function.

```c
#include "akhal/kstr.h"   // kstring_t, KS_INIT, ks_resize/ks_clear/ks_free/ks_release
#include "akhal/error.h"  // AK_OK / AK_ENOMEM, the return values of ks_resize()
```

## Contents

- [The struct](#the-struct) - [`KS_INIT`](#ks_init)
- [Growing and reusing](#growing-and-reusing) - [`ks_resize`](#ks_resize), [`ks_clear`](#ks_clear)
- [Releasing](#releasing) - [`ks_free`](#ks_free), [`ks_release`](#ks_release)

## The struct

```c
typedef struct {
    size_t l;   // length
    size_t m;   // capacity
    char  *s;   // string
} kstring_t;
```

- `l` is the used length in bytes, **excluding** the NUL terminator.
- `m` is the allocated size of `s` in bytes, terminator included. `m` must
  therefore be at least `l + 1` for any buffer you have written to.
- `s` is the buffer itself, and is `NULL` until the first successful
  `ks_resize()`.

The invariant the rest of the library relies on: **`s` is always
NUL-terminated once written to.** Nothing enforces it for you. `ks_resize()`
only changes the allocation, so code that writes bytes by hand is responsible
for setting `s[l] = '\0'` afterwards. Producers that fill the buffer for you -
[`ak_getline()`](io.md#ak_getline) is the main one - maintain it themselves.

### `KS_INIT`

```c
#define KS_INIT { 0, 0, NULL }
```

The zero initializer. It is a brace initializer, so it belongs in a
declaration; an already-declared struct is zeroed with `memset()` instead, or
reset without freeing by [`ks_clear()`](#ks_clear).

```c
kstring_t ks = KS_INIT;   // l = 0, m = 0, s = NULL

// KS_INIT expands to a brace list, so it only works in a declaration. A struct
// declared elsewhere is initialized by zeroing it, which is the same state.
kstring_t other;
memset(&other, 0, sizeof(other));

// Nothing is allocated yet: s stays NULL until the first ks_resize().
printf("l=%lu m=%lu s=%s\n", (unsigned long)ks.l, (unsigned long)ks.m,
       ks.s ? ks.s : "(null)");

ks_free(&ks);     // free(NULL) is a no-op, so this is safe on an unused buffer
ks_free(&other);
```

## Growing and reusing

### `ks_resize`

```c
int ks_resize(kstring_t *ks, size_t size);
```

Ensures the buffer holds at least `size` bytes, growing by doubling from a
floor of 16 and returning immediately when `m` is already large enough.
Returns `AK_OK`, or `AK_ENOMEM` when `realloc()` fails - in which case the old
buffer is left intact. `size` is the **total** capacity you want, not extra
bytes, and it has to include room for the terminator.

```c
kstring_t ks = KS_INIT;

// Total capacity, terminator included: 6 for a five-character string.
if (ks_resize(&ks, 6) != AK_OK) return 1;

memcpy(ks.s, "ACGTN", 5);
ks.l = 5;
ks.s[ks.l] = '\0';   // ks_resize() grows the allocation; it never terminates

// m comes back rounded up to a power of two (16 here), so repeated growth in
// a loop is amortized and existing contents are preserved across a resize.
printf("%s (l=%lu, m=%lu)\n", ks.s, (unsigned long)ks.l, (unsigned long)ks.m);

ks_free(&ks);
```

### `ks_clear`

```c
void ks_clear(kstring_t *ks);   // static inline in kstr.h
```

Resets `l` to 0 and writes a terminator at `s[0]`, keeping the allocation. This
is what makes one buffer usable for a whole file: clear, refill, repeat, and the
capacity converges on the longest item. Safe when `s` is still `NULL`.

```c
kstring_t ks = KS_INIT;
const char *names[] = { "chr1", "chr22_alt" };

for (int i = 0; i < 2; i++) {
    ks_clear(&ks);   // length back to 0, the allocation stays
    size_t n = strlen(names[i]);
    if (ks_resize(&ks, n + 1) != AK_OK) { ks_free(&ks); return 1; }
    memcpy(ks.s, names[i], n + 1);   // copy the terminator along with the text
    ks.l = n;
    printf("%s\n", ks.s);
}

// The second name grew the buffer; the capacity is still there to be reused.
printf("still %lu bytes allocated\n", (unsigned long)ks.m);
ks_free(&ks);
```

## Releasing

### `ks_free`

```c
void ks_free(kstring_t *ks);   // static inline in kstr.h
```

Frees `s` and zeroes all three fields, leaving the struct in its `KS_INIT`
state. It takes the struct by pointer and does not free the struct itself, so
it works on a stack `kstring_t`, and calling it twice is harmless.

```c
kstring_t ks = KS_INIT;
if (ks_resize(&ks, 32) != AK_OK) return 1;
memcpy(ks.s, "S\t1\tACGT", 9);
ks.l = 8;

ks_free(&ks);   // frees s, then sets s = NULL and l = m = 0

// Back to the KS_INIT state, so a second free is a no-op and the same struct
// can simply be used again.
ks_free(&ks);
if (ks_resize(&ks, 8) != AK_OK) return 1;
ks.s[0] = '\0';

ks_free(&ks);
```

### `ks_release`

```c
char *ks_release(kstring_t *ks);   // static inline in kstr.h
```

Detaches the allocation and hands it to you as a plain C string, leaving the
`kstring_t` empty. **Ownership transfers**: you must `free()` the returned
pointer, and `ks_free()` on that buffer will no longer do it for you. Returns
`NULL` if the buffer was never allocated. This is how a parser keeps one field
out of a line it is otherwise about to overwrite.

```c
kstring_t ks = KS_INIT;
if (ks_resize(&ks, 16) != AK_OK) return 1;
memcpy(ks.s, "chr22", 6);
ks.l = 5;

// The pointer moves out; ks is left as if freshly KS_INIT'd and is safe to
// keep filling. It will never free this allocation again.
char *name = ks_release(&ks);
if (!name) return 1;   // NULL only when s was still unallocated

printf("kept: %s (ks.l is now %lu)\n", name, (unsigned long)ks.l);
free(name);            // the caller frees, not ks_free()

ks_free(&ks);          // releases whatever ks holds now; here, nothing
```

---

[Back to the library index](README.md)
