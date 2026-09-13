#ifndef AKHAL_ARENA_H
#define AKHAL_ARENA_H

#include <stddef.h>
#include <stdint.h>

/**
 * Append-only byte arena.
 *
 * Values are copied into large blocks and handed back as pointers into them.
 * A block never moves once it has been cut, so every pointer the arena ever
 * returned stays valid however much more is added - which is exactly what a
 * single flat, realloc'd buffer cannot promise, and the reason this is a list
 * of blocks rather than one growing array. Releasing the arena frees the
 * blocks, not the values: one free() per few megabytes rather than one per
 * value.
 *
 * That is the whole point at graph scale. A chr22 graph carries 17.9 M segment
 * sequences; held as individual malloc()s they cost 17.9 M free()s to release
 * and an allocator header apiece, and on a whole-genome graph the teardown
 * alone walks hundreds of millions of chunks scattered across a cold heap.
 *
 * Every value is NUL-terminated, so an arena string can be handed to anything
 * expecting a C string, and the stored length is still the authority on how
 * many bytes it holds.
 *
 * Zero-initialize to create one - `ak_arena_t a = {0};`, or a calloc'd parent
 * struct - and release it with ak_arena_destroy(). There is no init function.
 */

#ifdef __cplusplus
extern "C" {
#endif

// The first block is cut at this size...
#define AK_ARENA_MIN ((size_t)64 << 10)   /* 64 KiB */
// ...and each later one at twice the last, up to this ceiling, so a small
// graph never allocates more than it needs and a large one stops paying
// per-block overhead.
#define AK_ARENA_MAX ((size_t)8 << 20)    /* 8 MiB */

typedef struct {
    char   **blk;        // the blocks; only the last one has room left
    int32_t  n_blk;      // blocks in use
    int32_t  m_blk;      // blocks the blk array has space for
    size_t   used;       // bytes taken in the last block
    size_t   cap;        // size of the last block
    size_t   blk_size;   // size the last ordinary block was cut at
    size_t   total;      // bytes across every block
} ak_arena_t;

/**
 * Copy a value into the arena and NUL-terminate the copy
 *
 * A value larger than a whole block is not refused: the block size doubles
 * until the value fits, so a segment that a compaction folded into a megabase
 * still lands in one piece
 * @param a Arena to append to
 * @param s Bytes to copy; may be NULL when len is 0
 * @param len Number of bytes, excluding any terminator
 * @return A pointer into the arena, valid until ak_arena_destroy(), or NULL
 *         on allocation failure (logged)
 */
const char *ak_arena_put(ak_arena_t *a, const char *s, size_t len);

/**
 * Free every block and zero the arena. Safe on a zeroed arena and safe to
 * call twice; every pointer ak_arena_put() returned dangles afterwards
 * @param a Arena to release
 */
void ak_arena_destroy(ak_arena_t *a);

/**
 * Bytes the arena holds across all of its blocks, allocated rather than used
 * @param a Arena to measure; NULL reads as 0
 * @return The total block size in bytes
 */
size_t ak_arena_bytes(const ak_arena_t *a);

#ifdef __cplusplus
}
#endif

#endif
