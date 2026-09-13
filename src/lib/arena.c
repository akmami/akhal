#include "akhal/arena.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

// Cut a new block with room for at least `need` bytes and make it the current
// one, abandoning whatever was left of the previous block.
//
// Two growth rules meet here. Ordinarily each block is twice the last, from
// AK_ARENA_MIN up to AK_ARENA_MAX, so a ten-segment graph never allocates
// megabytes and a whole-chromosome one stops paying per-block overhead. On top
// of that the size doubles until the value itself fits, however large it is -
// a segment that a compaction folded into a megabase still lands in one piece.
// That second growth only sticks while it stays under the ceiling: past it the
// oversized block is a one-off, so a single huge value does not leave every
// later block huge as well.
static int arena_block(ak_arena_t *a, size_t need) {
    size_t cap = a->blk_size ? a->blk_size : AK_ARENA_MIN;
    if (a->n_blk && cap < AK_ARENA_MAX) {
        cap <<= 1;
    }
    while (cap < need) {
        if (cap > (SIZE_MAX >> 1)) {   // one more doubling would wrap
            cap = need;
            break;
        }
        cap <<= 1;
    }

    if (a->n_blk == a->m_blk) {
        int32_t m = a->m_blk ? a->m_blk << 1 : 16;
        char **p = (char **)realloc(a->blk, (size_t)m * sizeof(*p));
        if (!p) return AK_ENOMEM;
        a->blk = p;
        a->m_blk = m;
    }

    char *b = (char *)malloc(cap);
    if (!b) return AK_ENOMEM;

    a->blk[a->n_blk++] = b;
    a->used   = 0;
    a->cap    = cap;
    a->total += cap;
    if (cap <= AK_ARENA_MAX) {
        a->blk_size = cap;
    } else if (a->blk_size == 0) {
        a->blk_size = AK_ARENA_MIN;    // the very first block was an outlier
    }
    return AK_OK;
}

// copy a value in, NUL-terminated; see akhal/arena.h
const char *ak_arena_put(ak_arena_t *a, const char *s, size_t len) {
    if (!a) return NULL;

    size_t need = len + 1;             // every value keeps a terminator
    if (need < len) {                  // len was SIZE_MAX
        ak_log(AK_LOG_ERROR, "arena", "value too long");
        return NULL;
    }

    // a->used never exceeds a->cap, so the subtraction cannot wrap
    if (a->n_blk == 0 || need > a->cap - a->used) {
        if (arena_block(a, need) != AK_OK) {
            ak_log(AK_LOG_ERROR, "arena", "out of memory");
            return NULL;
        }
    }

    char *dst = a->blk[a->n_blk - 1] + a->used;
    if (len) {
        memcpy(dst, s, len);
    }
    dst[len] = '\0';
    a->used += need;
    return dst;
}

// free every block; see akhal/arena.h
void ak_arena_destroy(ak_arena_t *a) {
    if (!a) return;
    for (int32_t i = 0; i < a->n_blk; i++) free(a->blk[i]);
    free(a->blk);
    memset(a, 0, sizeof(*a));          // so a second call is a no-op
}

size_t ak_arena_bytes(const ak_arena_t *a) {
    return a ? a->total : 0;
}
