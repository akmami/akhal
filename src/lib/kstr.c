#include "akhal/kstr.h"
#include "akhal/error.h"

/**
 * Grow a kstring so it holds at least `size` bytes, doubling capacity
 * @param ks Buffer to grow
 * @param size Minimum capacity required
 * @return AK_OK on success, AK_ENOMEM on allocation failure
 */
int ks_resize(kstring_t *ks, size_t size) {
    if (ks->m >= size) return AK_OK;

    // round up to the next power of two, starting at a small floor
    size_t m = ks->m ? ks->m : 16;
    while (m < size) {
        size_t next = m << 1;
        if (next < m) { m = size; break; }  // overflow guard
        m = next;
    }

    char *s = (char *)realloc(ks->s, m);
    if (!s) return AK_ENOMEM;
    ks->s = s;
    ks->m = m;
    return AK_OK;
}
