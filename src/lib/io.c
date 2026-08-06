#include "akhal/io.h"
#include "akhal/error.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct ak_file {
    FILE *fp;
    char *name;
};

/**
 * Open a file for reading and wrap it in an ak_file.
 * @param fn Path to open.
 * @return A handle, or NULL on failure (logged).
 */
ak_file *ak_open(const char *fn) {
    if (!fn) return NULL;

    FILE *fp = fopen(fn, "r");
    if (!fp) {
        ak_log(AK_LOG_ERROR, "io", "could not open %s", fn);
        return NULL;
    }

    ak_file *f = (ak_file *)calloc(1, sizeof(ak_file));
    if (!f) {
        fclose(fp);
        ak_log(AK_LOG_ERROR, "io", "out of memory opening %s", fn);
        return NULL;
    }
    f->fp = fp;
    f->name = NULL;  // reserved for future diagnostics
    return f;
}

/**
 * Read one line into a growable buffer, stripping the trailing newline.
 * @param f Open handle.
 * @param ks Reusable destination buffer (grows as needed).
 * @return Line length (>= 0), or AK_EOF / negative AK_E* code.
 */
long ak_getline(ak_file *f, kstring_t *ks) {
    if (!f || !f->fp) return AK_EINVAL;

    ks->l = 0;
    if (ks_resize(ks, 256) != AK_OK) return AK_ENOMEM;
    ks->s[0] = '\0';

    int any = 0;
    for (;;) {
        // read into the free tail of the buffer
        size_t space = ks->m - ks->l;
        if (space < 2) {
            if (ks_resize(ks, ks->m + 1) != AK_OK) return AK_ENOMEM;
            space = ks->m - ks->l;
        }

        if (!fgets(ks->s + ks->l, (int)space, f->fp)) break;
        any = 1;

        size_t chunk = strlen(ks->s + ks->l);
        ks->l += chunk;

        if (ks->l > 0 && ks->s[ks->l - 1] == '\n') break;   // full line read
        if (chunk + 1 < space) break;                        // EOF, no newline
    }

    if (!any && ks->l == 0) return AK_EOF;

    // strip a trailing "\r\n" or "\n"
    if (ks->l > 0 && ks->s[ks->l - 1] == '\n') ks->l--;
    if (ks->l > 0 && ks->s[ks->l - 1] == '\r') ks->l--;
    ks->s[ks->l] = '\0';

    return (long)ks->l;
}

/** Rewind the stream to the beginning. */
int ak_rewind(ak_file *f) {
    if (!f || !f->fp) return AK_EINVAL;
    rewind(f->fp);
    return AK_OK;
}

/** Close a handle and free it. Safe with NULL. */
void ak_close(ak_file *f) {
    if (!f) return;
    if (f->fp) fclose(f->fp);
    free(f->name);
    free(f);
}
