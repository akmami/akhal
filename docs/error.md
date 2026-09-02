# `error` - return codes and logging

Source: [`src/lib/error.c`](../src/lib/error.c) &middot; Header: [`include/akhal/error.h`](../include/akhal/error.h)

The convention the whole codebase follows, stated once here:

> **Code under `src/lib/` never calls `exit()`.** A function signals failure by
> returning a negative `AK_E*` code (int-returning functions) or `NULL`
> (pointer-returning functions), and may emit a diagnostic through `ak_log()`.
> Only the CLI layer under `src/cli/` decides whether to abort.

That is why every library entry point is checkable and every failure is
recoverable by the caller. Success is `0` and every error is negative, so the
test is always `if (rc < 0)` - never `if (rc)`, which would treat `AK_OK` and a
positive count identically.

Logging is a second, independent mechanism: a process-wide severity threshold
and one variadic printer to stderr. It reports; it does not decide anything.

```c
#include "akhal/error.h"   // AK_OK / AK_E* codes, AK_LOG_* levels, ak_strerror(), ak_log()
```

## Contents

- [Return codes](#return-codes)
- [Logging levels](#logging-levels)
- [Describing a code](#describing-a-code) - [`ak_strerror`](#ak_strerror)
- [Logging](#logging) - [`ak_log_set_level`](#ak_log_set_level), [`ak_log_level`](#ak_log_level), [`ak_log`](#ak_log)

## Return codes

| Code | Value | Meaning |
| --- | --- | --- |
| `AK_OK` | `0` | success |
| `AK_EOPEN` | `-1` | could not open a file |
| `AK_ENOMEM` | `-2` | allocation failed |
| `AK_EFORMAT` | `-3` | malformed or unexpected input |
| `AK_EINVAL` | `-4` | invalid argument from the caller |
| `AK_EOF` | `-5` | end of input - not necessarily an error |
| `AK_EIO` | `-6` | read/write failure on an already open file |

`AK_EOF` is the one to watch: it is negative like the rest, but for a reader
such as [`ak_getline()`](io.md#ak_getline) it is the normal way a loop ends, not
a failure to report. The values are stable and new codes are appended at the
end, so it is safe to store one.

## Logging levels

| Level | Value | Used for |
| --- | --- | --- |
| `AK_LOG_ERROR` | `0` | a failure the caller is about to be told about |
| `AK_LOG_WARN` | `1` | something survivable, such as a skipped malformed record |
| `AK_LOG_INFO` | `2` | progress and counts |
| `AK_LOG_DEBUG` | `3` | per-record detail |

The threshold defaults to `AK_LOG_INFO`, and a message prints when its level is
at or below it - so raising the threshold shows *more*. Setting it to
`AK_LOG_WARN` suppresses INFO and DEBUG. The CLI configures it once at startup.

## Describing a code

### `ak_strerror`

```c
const char *ak_strerror(int code);
```

Maps a code to a fixed English description. The returned pointer is static
storage: never `NULL`, never to be freed, and valid for the life of the
process. Any value outside the enum maps to `"unknown error"` rather than
failing.

```c
int rc = AK_EFORMAT;   // as if returned by a reader

// The idiom is a sign test, because 0 is success and readers may return
// positive counts.
if (rc < 0) {
    ak_log(AK_LOG_ERROR, "gfa", "read failed: %s", ak_strerror(rc));
    printf("%s\n", ak_strerror(rc));
}

// Static storage, so it can be held or printed without copying, and an
// out-of-range value is still safe to pass.
const char *msg = ak_strerror(-99);
printf("%s\n", msg);   // "unknown error"
return rc < 0 ? 1 : 0;
```

## Logging

### `ak_log_set_level`

```c
void ak_log_set_level(int level);
```

Sets the process-wide threshold. It is a single global with no locking, so set
it once during startup rather than toggling it from worker code. There is no
"restore previous" helper - read the old value with
[`ak_log_level()`](#ak_log_level) if you need to put it back.

```c
int verbosity = 2;   // how many times -v was given on the command line

ak_log_set_level(verbosity >= 2 ? AK_LOG_DEBUG :
                 verbosity == 1 ? AK_LOG_INFO  : AK_LOG_WARN);

ak_log(AK_LOG_DEBUG, "cli", "verbosity %d", verbosity);   // printed

// Raising the threshold shows more, lowering it shows less; AK_LOG_ERROR is
// as quiet as it gets, since errors are not suppressible.
int saved = ak_log_level();
ak_log_set_level(AK_LOG_ERROR);
ak_log(AK_LOG_WARN, "cli", "this one is dropped");
ak_log_set_level(saved);   // no helper does this for you
```

### `ak_log_level`

```c
int ak_log_level(void);
```

Returns the current threshold. It exists so a caller can skip building a
message that would be thrown away: `ak_log()` does test the threshold before it
formats anything, but your own arguments are evaluated before the call is even
made.

```c
size_t ids[] = { 12, 7, 44, 9 };
size_t n = sizeof(ids) / sizeof(ids[0]);

// The guard is about the cost of assembling the argument, not about ak_log()
// itself - it already drops the message before formatting.
if (ak_log_level() >= AK_LOG_DEBUG) {
    char buf[128] = "";
    size_t off = 0;
    for (size_t i = 0; i < n; i++) {
        int k = snprintf(buf + off, sizeof(buf) - off, "%lu ", (unsigned long)ids[i]);
        if (k < 0 || (size_t)k >= sizeof(buf) - off) break;   // truncated
        off += (size_t)k;
    }
    ak_log(AK_LOG_DEBUG, "gfa", "visited: %s", buf);
}
```

### `ak_log`

```c
void ak_log(int level, const char *ctx, const char *fmt, ...);
```

Writes one line to stderr as `[LEVEL] (ctx) message`, adding the newline for
you, and returns nothing - there is no way for logging to fail a caller. `ctx`
is a short subsystem tag such as `"gfa"`, `"gaf"` or `"io"`, and may be `NULL`,
in which case the parenthesised part is omitted. Under GCC and clang the format
string is checked like `printf`'s, so argument types have to match.

```c
const char *fn = "graph.gfa";
int rc = AK_EOPEN;

// "[ERROR] (io) could not open graph.gfa: could not open file"
ak_log(AK_LOG_ERROR, "io", "could not open %s: %s", fn, ak_strerror(rc));

// ctx is optional; NULL drops the "(io)" and nothing else changes.
ak_log(AK_LOG_INFO, NULL, "no subsystem tag on this one");

// Levels are clamped, not rejected: anything below AK_LOG_ERROR is tagged
// ERROR, while anything above AK_LOG_DEBUG is above every AK_LOG_* threshold
// and so never prints at all.
ak_log(AK_LOG_WARN, "gfa", "%lu record(s) skipped", (unsigned long)3);
ak_log(AK_LOG_DEBUG + 1, "gfa", "silently discarded");
return rc < 0 ? 1 : 0;
```

---

[Back to the library index](README.md)
