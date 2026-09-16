#ifndef AKHAL_GFA_H
#define AKHAL_GFA_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

#include "akhal/arena.h"

/**
 * In-memory model of an (r)GFA assembly/pangenome graph.
 *
 * Storage follows the "array + dict" design: every node (segment) and every
 * edge (link) lives in a contiguous array, and a hash table maps the external
 * segment id to its array index. Anything that cross-references a segment does
 * so by index, which keeps the structure cache-friendly and makes it cheap to
 * attach parallel per-node/per-edge arrays later.
 *
 * Both edges and paths are stored in CSR form (a flat array plus per-owner
 * offsets), so a node's out-edges and a path's ordered segments are each a
 * contiguous slice. A single segment may belong to many paths, which is why
 * paths are their own arrays rather than a pointer hung off a segment.
 *
 * Callers select how much work to do with the GFA_* flags.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Sentinel index for a path entry whose segment id was not found.
#define GFA_NIL UINT32_MAX

// Nodes (S lines)

typedef struct {
    uint64_t    id;          // segment id as it appears in the file
    const char *seq;         // borrowed from the graph's arena, NULL if empty; gfa_seg_set_seq() is the only way to set it
    uint32_t    len;         // sequence length (cached strlen)
    int32_t     rank;        // SR tag value, or -1 if the tag is absent
    int32_t     start;       // reference offset (SO tag or path layout), -1 unplaced
    int32_t     ref_path;    // index of the owning path in path[], -1 for none
} gfa_seg_t;

/**
 * End of a segment's span on its stable sequence: start + len
 *
 * Derived rather than stored - at 663 M segments a fourth 4-byte field is the
 * difference between a 40- and a 48-byte record. An unplaced segment carries
 * start == -1 and ends nowhere, which this reports as -1 rather than len-1
 * 
 * @param s Segment to measure
 * @return The end offset, or -1 when the segment has no offset
 */
static inline int32_t gfa_seg_end(const gfa_seg_t *s) {
    return s->start < 0 ? -1 : s->start + (int32_t)s->len;
}

// Edges (L lines)

typedef struct {
    uint32_t v;              // source segment index
    uint32_t w;              // destination segment index
    uint32_t overlap;        // overlap length in bp (the leading M count)
    char     from_orient;    // orientation of v on the L line ('+' or '-')
    char     to_orient;      // orientation of w on the L line ('+' or '-')
} gfa_link_t;

// Graph

typedef struct {
    // nodes
    gfa_seg_t  *seg;
    int32_t     n_seg, m_seg;

    // Per-segment degrees, built when GFA_DEGREES is set: in_degree[i] and
    // out_degree[i] count the links entering and leaving segment index i.
    // Both have length n_seg, and both are NULL if not built.
    int32_t    *in_degree, *out_degree;

    // edges
    gfa_link_t *link;
    int32_t     n_link, m_link;

    // Out-adjacency in CSR form (built when GFA_LINKS is set): the out-links
    // of segment index v are the link indices arc[arc_off[v] .. arc_off[v+1]).
    // Both are NULL if not built.
    uint32_t   *arc;         // length n_link
    int32_t    *arc_off;     // length n_seg + 1

    // paths (P lines)
    char      **path;        // owned path names, length n_path
    uint64_t   *path_len;    // total sequence length of each path (GFA_PATHS)
    int32_t     n_path, m_path;

    // Path membership in CSR form: the segments of path k are
    // path_seg[path_off[k] .. path_off[k+1]), with matching orientation chars
    // in path_ori. Entries may be GFA_NIL for unresolved ids.
    int32_t    *path_off;    // length n_path + 1
    uint32_t   *path_seg;    // length n_path_seg
    char       *path_ori;    // length n_path_seg ('+'/'-')
    int32_t     m_path_seg;
    uint64_t    n_path_seg;  // total steps across all paths

    // Backing store for every segment sequence. One allocation per few
    // megabytes instead of one per segment, so a graph with millions of
    // segments is released in a handful of free()s rather than millions.
    ak_arena_t  strs;

    void       *idx;         // opaque id -> index hash table
    int         flags;       // the GFA_* flags this graph was read with
    int         has_sr;      // 1 when the file itself carried SR:i: tags
} gfa_t;

// Read flags

// S lines - segments
#define GFA_SEGS       0x01  // parse S lines: seg[], lengths, tags and the id index (implied by every flag below except GFA_PATH_NAMES)
#define GFA_SEQ        0x02  // also copy the bases into the graph's arena

// L lines - links
#define GFA_LINKS      0x04  // record edges (link[])
#define GFA_ARCS       0x08  // also build the CSR out-adjacency (implies GFA_LINKS)
#define GFA_DEGREES    0x10  // also build in_degree/out_degree (implies GFA_LINKS)

// P lines - paths
#define GFA_PATH_NAMES 0x20  // record path names and step counts only; the steps are counted with one pass over the commas, never parsed
#define GFA_PATHS      0x40  // resolve every step: membership (CSR), orientation and reference layout (implies GFA_PATH_NAMES)

// checks
#define GFA_VALIDATE   0x80  // warn about L lines naming unknown segments and, when GFA_SEQ is also set, about overlap mismatches; path steps are checked under GFA_PATHS

// What a caller wanting the whole graph asks for
#define GFA_ALL        (GFA_SEGS | GFA_SEQ | GFA_LINKS | GFA_ARCS | GFA_DEGREES | GFA_PATH_NAMES | GFA_PATHS)

/**
 * Read an (r)GFA file into a freshly allocated graph
 *
 * Soft problems found under GFA_VALIDATE are reported via ak_log() but do not
 * fail the read. The result must be released with gfa_destroy()
 *
 * Ranks: a file's own SR:i: tags are authoritative and are never touched. Only
 * when the file carried none (g->has_sr stays 0) and GFA_PATHS was requested
 * are ranks derived, exactly as gfa_rank_paths() would - so a plain GFA with
 * P lines comes back with a rank-0 backbone, and one without them comes back
 * entirely rank 1. Read without GFA_PATHS there is nothing to derive from, and
 * ranks are left absent (-1)
 * 
 * @param fn Path to the .gfa / .rgfa file
 * @param flags Bitwise OR of the GFA_* read flags; 0 reads an empty graph
 * @return The graph, or NULL on a fatal error (unreadable file, OOM)
 */
gfa_t *gfa_read(const char *fn, int flags);

// Statistics without a graph

// What gfa_read_stats() should work out beyond the streaming totals. Each
// costs memory in proportion to the file (see below); neither is needed for
// the counts and the length / overlap distributions.
#define GFA_STAT_DEGREES 0x1   // in- and out-degree distributions, and n_undefined
#define GFA_STAT_RANKS   0x2   // n_rank0 from the P lines when the file has no SR tags

/**
 * What a summary pass over a file can say about it. Counts are of lines, so
 * n_seg is S lines and n_link L lines, whether or not every id they name is
 * defined elsewhere in the file.
 *
 * The distributions are populations, not samples: the standard deviations
 * divide by n. The degree distributions cover only segments with a non-zero
 * degree, so a graph whose segments all stand alone reports -1 for the four
 * extremes rather than 0.
 *
 * A figure that was not asked for is -1: n_rank0 when the file has no SR
 * tags and GFA_STAT_RANKS was off, and n_undefined and the degree extremes
 * when GFA_STAT_DEGREES was off.
 */
typedef struct {
    int64_t  n_seg;          // S lines
    int64_t  n_link;         // L lines
    int64_t  n_path;         // P lines
    uint64_t n_bp;           // total sequence length: the sum of every segment's length
    int64_t  n_rank0;        // segments at rank 0; -1 when not derived
    int      has_sr;         // whether those SR tags were the file's own

    double   seg_mean;       // segment length
    double   seg_sd;
    uint64_t seg_min, seg_max;

    double   ov_mean;        // link overlap
    double   ov_sd;

    double   in_mean, in_sd;     // degrees, over segments that have any
    double   out_mean, out_sd;
    int32_t  min_in, max_in;     // -1 for none
    int32_t  min_out, max_out;

    int64_t  n_undefined;    // distinct ids an L or P line names that no S line defines; -1 when not checked
} gfa_stat_t;

/**
 * Summarize a GFA in one streaming pass, without building a graph.
 *
 * The counts and the length and overlap distributions are running totals and
 * cost nothing per line. The rest are properties of the whole file - a
 * degree is how many L lines name a segment, "undefined" is named but never
 * defined, "rank 0" without SR tags is named by some P line - and are
 * answered from flat arrays of ids: appended as the lines stream past, sorted in
 * place afterwards, and read off as run lengths and merges. An array costs
 * 4 bytes per entry, or 8 once an id exceeds 32 bits, regardless of how the
 * ids are numbered.
 *
 * GFA_STAT_DEGREES keeps three arrays: the S ids, the L sources and the L
 * targets - about 4 * (n_seg + 2 * n_link) bytes, so 8 GB on a graph of 700
 * million segments and links. GFA_STAT_RANKS keeps the S ids and every P
 * step, which on a file carrying all haplotypes as P lines is far larger than
 * the graph itself; a file that ranks itself with SR tags needs nothing, and
 * answers n_rank0 with either setting. Without either flag the pass holds
 * nothing but the line buffer.
 *
 * @param fn Path to the .gfa / .rgfa file
 * @param st Filled in on success
 * @param flags Bitwise OR of GFA_STAT_DEGREES and GFA_STAT_RANKS, or 0
 * @return AK_OK, or a negative AK_E* code (unreadable file, OOM)
 */
int gfa_read_stats(const char *fn, gfa_stat_t *st, int flags);

/**
 * Write a graph back out as GFA: an H line, one S per segment (with SR:i:
 * where a rank is set), one L per link, and one P per path
 * 
 * @param g Graph to emit
 * @param out Destination stream
 * @return AK_OK, or AK_EIO if the stream went bad
 */
int gfa_write(const gfa_t *g, FILE *out);

/**
 * Give a segment its sequence, copying it into the graph's arena
 *
 * This is the only supported way to set gfa_seg_t.seq. The field points into
 * storage the graph owns and releases in one go, so assigning a malloc'd
 * buffer to it directly would be freed twice, and free()ing it would corrupt
 * the arena. An empty sequence leaves the segment with seq == NULL and len 0,
 * which is what the reader records for an S line carrying none
 * 
 * @param g Graph owning the segment
 * @param s Segment to set, obtained from gfa_seg_at()
 * @param seq Bytes to copy; may be NULL when len is 0
 * @param len Sequence length in bases
 * @return AK_OK, AK_EINVAL on a NULL graph or segment, or AK_ENOMEM
 */
int gfa_seg_set_seq(gfa_t *g, gfa_seg_t *s, const char *seq, size_t len);

/**
 * Write a graph as rGFA: the same lines gfa_write() emits, plus the stable
 * sequence each segment sits on - SN:Z: from `ref_name` and SO:i: from `start`
 * alongside the SR:i: rank.
 *
 * Each tag is emitted only where the segment carries it, so one left without a
 * name or an offset (a NULL `ref_name`, a negative `start`) simply comes out
 * with the tags it does have. See rgfa_build(), which works those out
 * 
 * @param g Graph to emit
 * @param out Destination stream
 * @return AK_OK, or AK_EIO if the stream went bad
 */
int gfa_write_rgfa(const gfa_t *g, FILE *out);

/**
 * Release a graph and everything it owns. Safe to call with NULL
 * 
 * @param g Graph to destroy
 */
void gfa_destroy(gfa_t *g);

// Accessors / traversal

/**
 * Look up a segment's array index by id. O(1)
 * 
 * @param g Graph to query
 * @param id Segment id
 * @return The array index, or -1 if absent
 */
int32_t gfa_idx(const gfa_t *g, uint64_t id);

/**
 * Look up a segment by id
 * 
 * @param g Graph to query
 * @param id Segment id
 * @return Pointer to the segment, or NULL if absent
 */
gfa_seg_t *gfa_get(const gfa_t *g, uint64_t id);

/**
 * The stable sequence a segment sits on - its SN name
 *
 * Held as an index into path[] rather than a borrowed pointer, so nothing
 * dangles when the path block is rewritten and the record stays 40 bytes
 * 
 * @param g Graph owning the segment
 * @param s Segment to ask about
 * @return The path name, or NULL when the segment belongs to none
 */
static inline const char *gfa_seg_ref(const gfa_t *g, const gfa_seg_t *s) {
    return s->ref_path < 0 ? NULL : g->path[s->ref_path];
}

/** 
 * @return Number of segments (nodes) in the graph
 */
static inline int32_t gfa_n_seg(const gfa_t *g) { 
    return g->n_seg; 
}

/** 
 * @return Number of links (edges) in the graph
 */
static inline int32_t gfa_n_link(const gfa_t *g) { 
    return g->n_link; 
}

/** 
 * @return Number of paths in the graph
 */
static inline int32_t gfa_n_path(const gfa_t *g) { 
    return g->n_path; 
}

/** 
 * @return Segment at array index i
 */
static inline gfa_seg_t *gfa_seg_at(const gfa_t *g, int32_t i) { 
    return &g->seg[i]; 
}

/** 
 * @return Link at array index i
 */
static inline gfa_link_t *gfa_link_at(const gfa_t *g, int32_t i) { 
    return &g->link[i]; 
}

/** 
 * @return Name of path k (borrowed)
 */
static inline const char *gfa_path_name(const gfa_t *g, int32_t k) { 
    return g->path[k]; 
}

/** 
 * @return Total sequence length of path k
 */
static inline uint64_t gfa_path_len(const gfa_t *g, int32_t k) { 
    return g->path_len[k]; 
}

/**
 * Out-edge traversal for a segment
 * 
 * @param g Graph to query (must have been read with GFA_LINKS)
 * @param v Segment index whose out-edges are wanted
 * @param arcs Set to an array of link indices leaving v; feed gfa_link_at()
 * @return Number of out-edges, or 0 if v has none or adjacency was not built
 */
int gfa_arcs(const gfa_t *g, int32_t v, const uint32_t **arcs);

/**
 * Test for a directed link v -> w. Requires GFA_LINKS. O(out-degree of v)
 * 
 * @param g Graph to query
 * @param v Source segment index
 * @param w Destination segment index
 * @return 1 if the link exists, else 0
 */
int gfa_has_arc(const gfa_t *g, int32_t v, int32_t w);

/**
 * Test for an oriented join v(ov) -> w(ow), the way a path walks it. An L
 * line states one strand of a join and implies the other, so this is
 * satisfied by "L v ov w ow" and equally by its complement "L w !ow v !ov".
 * Requires GFA_ARCS. O(out-degree of v + out-degree of w)
 * 
 * @param g Graph to query
 * @param v Source segment index
 * @param ov Orientation v is walked in, '+' or '-'
 * @param w Destination segment index
 * @param ow Orientation w is walked in, '+' or '-'
 * @return 1 if the join exists on either strand, else 0
 */
int gfa_has_link(const gfa_t *g, int32_t v, char ov, int32_t w, char ow);

/**
 * Ordered segments of a path. The matching orientation chars are in
 * g->path_ori at the same offset; an entry may be GFA_NIL for an unresolved id
 * 
 * @param g Graph to query (must have been read with GFA_PATHS)
 * @param k Path index
 * @param segs Set to an array of segment indices; feed gfa_seg_at()
 * @return Number of segments in the path, or 0 (with *segs NULL) if none or
 *         the steps were not built (GFA_PATH_NAMES alone)
 */
int gfa_path_segs(const gfa_t *g, int32_t k, const uint32_t **segs);

/**
 * Number of steps in path k, known from GFA_PATH_NAMES on (no step arrays needed)
 * 
 * @param g Graph to query
 * @param k Path index
 * @return The step count, or 0 if k is out of range or no path flag was set
 */
static inline int64_t gfa_path_n_steps(const gfa_t *g, int32_t k) {
    if (!g->path_off || k < 0 || k >= g->n_path) return 0;
    return (int64_t)g->path_off[k + 1] - g->path_off[k];
}

// Ranks

/**
 * Rank the segments against the graph's own paths.
 *
 * rGFA uses SR:i: to say how far a segment sits from the reference: 0 is the
 * reference backbone, anything higher came from a sample. Every segment that
 * any P line visits is stamped rank 0 and every other segment rank 1, so a
 * graph with no paths at all comes back entirely rank 1.
 *
 * Every P line counts as backbone, so a reference spread over several P lines
 * is ranked no differently from one that arrives whole.
 *
 * Existing SR values are overwritten. gfa_read() calls this itself when the
 * file carried no SR tags, so calling it explicitly is how you re-rank a graph
 * whose tags you want to replace. Requires GFA_PATHS
 * 
 * @param g Graph to rank, modified in place
 * @return Number of segments left at rank 0, or a negative AK_E* code
 */
int64_t gfa_rank_paths(gfa_t *g);

/**
 * Rank the segments against a caller-supplied backbone: rank 0 wherever `on`
 * is set, rank 1 everywhere else. This is the general form of
 * gfa_rank_paths(), for a backbone that did not come from the P lines - a
 * traced reference sequence, for instance (see call_ref_fasta)
 * 
 * @param g Graph to rank, modified in place
 * @param on Flags of length gfa_n_seg(g); non-zero marks a backbone segment
 * @return Number of segments set to rank 0, or a negative AK_E* code
 */
int64_t gfa_rank_mark(gfa_t *g, const uint8_t *on);

// Rewriting the path block

/**
 * Drop every path in the graph.
 *
 * Segments refer to their path by index, so this resets every ref_path to -1
 * rather than leaving one pointing at a path that no longer exists.
 * The segments themselves, and their ranks, are untouched
 * 
 * @param g Graph to modify
 */
void gfa_clear_paths(gfa_t *g);

/**
 * Append one path to the graph, laying it out as the reader would: each
 * segment's reference start/end is recomputed along the path and its ref_name
 * repointed at the new name.
 *
 * Together with gfa_clear_paths() and gfa_rank_mark() this is what makes an
 * externally supplied reference the graph's backbone - mark the ranks, then
 * install the walk that produced them in place of the old P lines.
 * Requires GFA_PATHS
 * 
 * @param g Graph to modify
 * @param name Name for the new path; copied
 * @param segs Ordered segment indices; GFA_NIL entries are skipped
 * @param ori Matching orientation chars, or NULL to treat every step as '+'
 * @param n Number of entries in segs
 * @return AK_OK, or a negative AK_E* code
 */
int gfa_add_path(gfa_t *g, const char *name, const uint32_t *segs, const char *ori, int64_t n);

// Ordering

/**
 * Topologically order the segments (Kahn's algorithm on the directed graph
 * given by the links). Ties in the ready set are broken by node sequence
 * content, alphabetically, so the ordering does not depend on the input's node
 * numbering (a NULL/empty sequence sorts first). Any nodes that remain inside
 * cycles are appended after the acyclic prefix, also by sequence, so `order`
 * is always a full permutation of 0..n_seg-1.
 *
 * Requires the graph was read with GFA_LINKS
 * 
 * @param g Graph to order
 * @param order Caller-allocated array of length n_seg; filled with segment
 *              indices in topological order
 * @return The number of nodes placed before any cycle (n_seg if the graph is
 *         acyclic), or a negative AK_E* code on error
 */
int gfa_toposort(const gfa_t *g, int32_t *order);

/**
 * Whether the graph carries per-segment degree arrays (read with GFA_DEGREES)
 * 
 * @param g Graph to query
 * @return Non-zero when in_degree and out_degree may be indexed
 */
static inline int gfa_has_degrees(const gfa_t *g) {
    return g->in_degree != NULL && g->out_degree != NULL;
}

/**
 * Whether the graph carries segments (read with GFA_SEGS or any flag implying
 * it). Read with GFA_PATH_NAMES alone it does not, and seg[] must not be used
 * 
 * @param g Graph to query
 * @return Non-zero when seg[] and the segment accessors may be used
 */
static inline int gfa_has_segs(const gfa_t *g) {
    return (g->flags & GFA_SEGS) != 0;
}

#ifdef __cplusplus
}
#endif

#endif
