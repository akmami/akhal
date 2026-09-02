#include "akhal/kstr.h"
#include "akhal/error.h"

// grows capacity to the next power of two; does not NUL-terminate
int ks_resize(kstring_t *ks, size_t size) {
    if (ks->m >= size) return AK_OK;

    // round up to the next power of two, starting at a small floor
    size_t m = ks->m ? ks->m : 16;
    while (m < size) {
        size_t next = m << 1;
        if (next < m) {   // overflow guard
            m = size;
            break;
        }
        m = next;
    }

    char *s = (char *)realloc(ks->s, m);
    if (!s) return AK_ENOMEM;
    ks->s = s;
    ks->m = m;
    return AK_OK;
}
