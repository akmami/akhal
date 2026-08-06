#ifndef AKHAL_ERROR_H
#define AKHAL_ERROR_H

/**
 * Error handling and logging for libakhal.
 *
 * Library rule: code under src/lib/ NEVER calls exit(). Functions signal
 * failure by returning a negative AK_E* code (int-returning functions) or
 * NULL (pointer-returning functions), and may emit a diagnostic through
 * ak_log(). Only the CLI layer (src/cli/) decides whether to abort.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Return codes. 0 is success; all errors are negative so callers can test
// `if (rc < 0)`. Keep these stable; append new codes at the end.
enum {
    AK_OK      =  0,
    AK_EOPEN   = -1,   // could not open a file
    AK_ENOMEM  = -2,   // allocation failed
    AK_EFORMAT = -3,   // malformed / unexpected input
    AK_EINVAL  = -4,   // invalid argument from the caller
    AK_EOF     = -5    // end of input (not necessarily an error)
};

/**
 * Human-readable description of an AK_E* return code.
 * @param code An AK_OK / AK_E* value.
 * @return A static string; never NULL.
 */
const char *ak_strerror(int code);

// Logging levels, ordered from most to least severe. Setting the threshold
// to AK_LOG_WARN suppresses INFO and DEBUG. The CLI configures it once.
enum {
    AK_LOG_ERROR = 0,
    AK_LOG_WARN  = 1,
    AK_LOG_INFO  = 2,
    AK_LOG_DEBUG = 3
};

/**
 * Set the global verbosity threshold (default: AK_LOG_INFO).
 * @param level One of the AK_LOG_* levels.
 */
void ak_log_set_level(int level);

/**
 * Current verbosity threshold, handy for guarding expensive debug formatting.
 * @return The active AK_LOG_* level.
 */
int ak_log_level(void);

/**
 * Emit a diagnostic to stderr if level is at or below the threshold.
 * @param level Severity, one of the AK_LOG_* levels.
 * @param ctx Short subsystem tag such as "gfa" or "gaf"; may be NULL.
 * @param fmt printf-style format string, followed by its arguments.
 */
void ak_log(int level, const char *ctx, const char *fmt, ...)
#if defined(__GNUC__)
    __attribute__((format(printf, 3, 4)))
#endif
    ;

#ifdef __cplusplus
}
#endif

#endif  // AKHAL_ERROR_H
