#ifndef AKHAL_KSTR_H
#define AKHAL_KSTR_H

#include <stddef.h>
#include <string.h>
#include <stdlib.h>

/**
 * A minimal growable string buffer, in the spirit of htslib's kstring_t.
 *
 * `s` is always NUL-terminated once written to, `l` is the used length
 * (excluding the terminator), and `m` is the allocated capacity. Initialize
 * with `kstring_t ks = KS_INIT;` (or memset to zero).
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t l;   // length
    size_t m;   // capacity
    char  *s;   // string
} kstring_t;

#define KS_INIT { 0, 0, NULL }

/**
 * Ensure the buffer holds at least `size` bytes.
 * @param ks Buffer to grow.
 * @param size Minimum capacity required.
 * @return AK_OK on success, AK_ENOMEM on allocation failure.
 */
int ks_resize(kstring_t *ks, size_t size);

/**
 * Reset length to 0 without freeing the buffer, so it can be reused.
 * @param ks Buffer to clear.
 */
static inline void ks_clear(kstring_t *ks) {
    ks->l = 0;
    if (ks->s) ks->s[0] = '\0';
}

/**
 * Free the buffer and zero the struct.
 * @param ks Buffer to release.
 */
static inline void ks_free(kstring_t *ks) {
    free(ks->s);
    ks->s = NULL;
    ks->l = ks->m = 0;
}

/**
 * Detach the owned C string, leaving the kstring empty.
 * @param ks Buffer to take ownership from.
 * @return The owned string (caller must free), or NULL if never allocated.
 */
static inline char *ks_release(kstring_t *ks) {
    char *s = ks->s;
    ks->s = NULL;
    ks->l = ks->m = 0;
    return s;
}

#ifdef __cplusplus
}
#endif

#endif  // AKHAL_KSTR_H
