#include "akhal/error.h"

#include <stdio.h>
#include <stdarg.h>

static int g_log_level = AK_LOG_INFO;

// map an AK_E* code to a static human-readable string
const char *ak_strerror(int code) {
    switch (code) {
        case AK_OK:      return "success";
        case AK_EOPEN:   return "could not open file";
        case AK_ENOMEM:  return "out of memory";
        case AK_EFORMAT: return "malformed input";
        case AK_EINVAL:  return "invalid argument";
        case AK_EOF:     return "end of input";
        case AK_EIO:     return "read/write failure";
        default:         return "unknown error";
    }
}

// set the global verbosity threshold
void ak_log_set_level(int level) {
    g_log_level = level;
}

/**
 * @return The current verbosity threshold
 */
int  ak_log_level(void) { 
    return g_log_level;
}

/**
 * Emit a diagnostic to stderr, gated by the verbosity threshold
 * @param level Severity (AK_LOG_*)
 * @param ctx Optional subsystem tag; may be NULL
 * @param fmt printf-style format followed by its arguments
 */
void ak_log(int level, const char *ctx, const char *fmt, ...) {
    static const char *tag[] = { "ERROR", "WARN", "INFO", "DEBUG" };

    if (level > g_log_level) return;
    if (level < AK_LOG_ERROR) level = AK_LOG_ERROR;
    if (level > AK_LOG_DEBUG) level = AK_LOG_DEBUG;

    fprintf(stderr, "[%s]", tag[level]);
    if (ctx) fprintf(stderr, " (%s)", ctx);
    fputc(' ', stderr);

    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);

    fputc('\n', stderr);
}
