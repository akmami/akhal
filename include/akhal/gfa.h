#ifndef AKHAL_GFA_H
#define AKHAL_GFA_H

#include <stdint.h>
#include <stddef.h>

/**
 * In-memory model of an (r)GFA assembly graph.
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
 * The single reader gfa_read() replaces the previously duplicated per-command
 * parsers. Callers select how much work to do with the GFA_* flags.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Sentinel index for a path entry whose segment id was not found.
#define GFA_NIL UINT32_MAX

// ---- Nodes (S lines) ----

typedef struct {
    uint64_t id;             // segment id as it appears in the file
    char    *seq;            // owned sequence, or NULL if empty
    uint32_t len;            // sequence length (cached strlen)
    int32_t  rank;           // SR tag value, or -1 if the tag is absent
    int32_t  start;          // reference offset (SO tag or path layout)
    int32_t  end;            // start + len
    int32_t  in_degree;      // populated when GFA_LINKS is set
    int32_t  out_degree;     // populated when GFA_LINKS is set
    const char *ref_name;    // borrowed: an owning path name (NULL if none)
} gfa_seg_t;

// ---- Edges (L lines) ----

typedef struct {
    uint32_t v;              // source segment index
    uint32_t w;              // destination segment index
    uint32_t overlap;        // overlap length in bp (the leading M count)
    char     from_orient;    // orientation of v on the L line ('+' or '-')
    char     to_orient;      // orientation of w on the L line ('+' or '-')
} gfa_link_t;

// ---- Graph ----

typedef struct {
    // nodes
    gfa_seg_t  *seg;
    int32_t     n_seg, m_seg;

    // edges
    gfa_link_t *link;
    int32_t     n_link, m_link;

    // Out-adjacency in CSR form (built when GFA_LINKS is set): the out-links
    // of segment index v are the link indices arc[arc_off[v] .. arc_off[v+1]).
    // Both are NULL if not built.
    uint32_t   *arc;         // length n_link
    int32_t    *arc_off;     // length n_seg + 1

    // paths (P lines), built when GFA_PATHS is set
    char      **path;        // owned path names, length n_path
    uint64_t   *path_len;    // total sequence length of each path
    int32_t     n_path, m_path;

    // Path membership in CSR form: the segments of path k are
    // path_seg[path_off[k] .. path_off[k+1]), with matching orientation chars
    // in path_ori. Entries may be GFA_NIL for unresolved ids.
    int32_t    *path_off;    // length n_path + 1
    uint32_t   *path_seg;    // length n_path_seg
    char       *path_ori;    // length n_path_seg ('+'/'-')
    int32_t     m_path_seg;
    uint64_t    n_path_seg;  // total segment occurrences across all paths

    void       *idx;         // opaque id -> index hash table
    int         flags;       // the GFA_* flags this graph was read with
} gfa_t;

// ---- Read flags ----

#define GFA_LINKS    0x1     // record edges, degrees, and out-adjacency
#define GFA_PATHS    0x2     // build path membership (CSR) and layout
#define GFA_VALIDATE 0x4     // check overlap consistency + integrity

/**
 * Read an (r)GFA file into a freshly allocated graph.
 *
 * Soft problems found under GFA_VALIDATE are reported via ak_log() but do not
 * fail the read. The result must be released with gfa_destroy().
 * @param fn Path to the .gfa / .rgfa file.
 * @param flags Bitwise OR of GFA_LINKS, GFA_PATHS, GFA_VALIDATE (may be 0).
 * @return The graph, or NULL on a fatal error (unreadable file, OOM).
 */
gfa_t *gfa_read(const char *fn, int flags);

/**
 * Release a graph and everything it owns. Safe to call with NULL.
 * @param g Graph to destroy.
 */
void gfa_destroy(gfa_t *g);

// ---- Accessors / traversal ----

/**
 * Look up a segment's array index by id. O(1).
 * @param g Graph to query.
 * @param id Segment id.
 * @return The array index, or -1 if absent.
 */
int32_t gfa_idx(const gfa_t *g, uint64_t id);

/**
 * Look up a segment by id.
 * @param g Graph to query.
 * @param id Segment id.
 * @return Pointer to the segment, or NULL if absent.
 */
gfa_seg_t *gfa_get(const gfa_t *g, uint64_t id);

/** @return Number of segments (nodes) in the graph. */
static inline int32_t     gfa_n_seg(const gfa_t *g)  { return g->n_seg; }
/** @return Number of links (edges) in the graph. */
static inline int32_t     gfa_n_link(const gfa_t *g) { return g->n_link; }
/** @return Number of paths in the graph. */
static inline int32_t     gfa_n_path(const gfa_t *g) { return g->n_path; }
/** @return Segment at array index i. */
static inline gfa_seg_t  *gfa_seg_at(const gfa_t *g, int32_t i)  { return &g->seg[i]; }
/** @return Link at array index i. */
static inline gfa_link_t *gfa_link_at(const gfa_t *g, int32_t i) { return &g->link[i]; }
/** @return Name of path k (borrowed). */
static inline const char *gfa_path_name(const gfa_t *g, int32_t k) { return g->path[k]; }
/** @return Total sequence length of path k. */
static inline uint64_t    gfa_path_len(const gfa_t *g, int32_t k)  { return g->path_len[k]; }

/**
 * Out-edge traversal for a segment.
 * @param g Graph to query (must have been read with GFA_LINKS).
 * @param v Segment index whose out-edges are wanted.
 * @param arcs Set to an array of link indices leaving v; feed gfa_link_at().
 * @return Number of out-edges, or 0 if v has none or adjacency was not built.
 */
int gfa_arcs(const gfa_t *g, int32_t v, const uint32_t **arcs);

/**
 * Test for a directed link v -> w. Requires GFA_LINKS. O(out-degree of v).
 * @param g Graph to query.
 * @param v Source segment index.
 * @param w Destination segment index.
 * @return 1 if the link exists, else 0.
 */
int gfa_has_arc(const gfa_t *g, int32_t v, int32_t w);

/**
 * Ordered segments of a path. The matching orientation chars are in
 * g->path_ori at the same offset; an entry may be GFA_NIL for an unresolved id.
 * @param g Graph to query (must have been read with GFA_PATHS).
 * @param k Path index.
 * @param segs Set to an array of segment indices; feed gfa_seg_at().
 * @return Number of segments in the path, or 0 if none / paths not built.
 */
int gfa_path_segs(const gfa_t *g, int32_t k, const uint32_t **segs);

// ---- Ordering ----

/**
 * Topologically order the segments (Kahn's algorithm on the directed graph
 * given by the links). Ties in the ready set are broken by node sequence
 * content, alphabetically, so the ordering does not depend on the input's node
 * numbering (a NULL/empty sequence sorts first). Any nodes that remain inside
 * cycles are appended after the acyclic prefix, also by sequence, so `order`
 * is always a full permutation of 0..n_seg-1.
 *
 * Requires the graph was read with GFA_LINKS.
 * @param g Graph to order.
 * @param order Caller-allocated array of length n_seg; filled with segment
 *              indices in topological order.
 * @return The number of nodes placed before any cycle (n_seg if the graph is
 *         acyclic), or a negative AK_E* code on error.
 */
int gfa_toposort(const gfa_t *g, int32_t *order);

#ifdef __cplusplus
}
#endif

#endif  // AKHAL_GFA_H
