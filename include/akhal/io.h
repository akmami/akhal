#ifndef AKHAL_IO_H
#define AKHAL_IO_H

#include "akhal/kstr.h"

/**
 * Line-oriented input handle.
 *
 * ak_file owns its own line buffer, so callers no longer manage a
 * `char *line` + capacity by hand. It is the single choke point through which
 * every text format reader (GFA, GAF, FASTA, SAM) pulls lines, which is where
 * buffered or compressed backends can later be slotted in without touching
 * the parsers.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ak_file ak_file;

/**
 * Open a file for reading.
 * @param fn Path to open.
 * @return A handle, or NULL on failure (the reason is logged).
 */
ak_file *ak_open(const char *fn);

/**
 * Read the next line, with the trailing '\n'/'\r' stripped.
 * @param f Open handle.
 * @param ks Reusable buffer to receive the NUL-terminated line; grows as needed.
 * @return The line length (>= 0), or AK_EOF at end of input.
 */
long ak_getline(ak_file *f, kstring_t *ks);

/**
 * Rewind to the beginning of the stream.
 * @param f Open handle.
 * @return AK_OK, or a negative AK_E* code.
 */
int ak_rewind(ak_file *f);

/**
 * Close a handle and free it. Safe to call with NULL.
 * @param f Handle to close.
 */
void ak_close(ak_file *f);

#ifdef __cplusplus
}
#endif

#endif  // AKHAL_IO_H
