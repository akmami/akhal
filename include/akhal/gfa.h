#ifndef AKHAL_GFA_H
#define AKHAL_GFA_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

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

// Nodes (S lines)

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
    int         has_sr;      // 1 when the file itself carried SR:i: tags
} gfa_t;

// Read flags

#define GFA_LINKS    0x1     // record edges, degrees, and out-adjacency
#define GFA_PATHS    0x2     // build path membership (CSR) and layout
#define GFA_VALIDATE 0x4     // check overlap consistency + integrity

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
 * @param fn Path to the .gfa / .rgfa file
 * @param flags Bitwise OR of GFA_LINKS, GFA_PATHS, GFA_VALIDATE (may be 0)
 * @return The graph, or NULL on a fatal error (unreadable file, OOM)
 */
gfa_t *gfa_read(const char *fn, int flags);

/**
 * Write a graph back out as GFA: an H line, one S per segment (with SR:i:
 * where a rank is set), one L per link, and one P per path
 * @param g Graph to emit
 * @param out Destination stream
 * @return AK_OK, or AK_EIO if the stream went bad
 */
int gfa_write(const gfa_t *g, FILE *out);

/**
 * Write a graph as rGFA: the same lines gfa_write() emits, plus the stable
 * sequence each segment sits on - SN:Z: from `ref_name` and SO:i: from `start`
 * alongside the SR:i: rank.
 *
 * Each tag is emitted only where the segment carries it, so one left without a
 * name or an offset (a NULL `ref_name`, a negative `start`) simply comes out
 * with the tags it does have. See rgfa_build(), which works those out
 * @param g Graph to emit
 * @param out Destination stream
 * @return AK_OK, or AK_EIO if the stream went bad
 */
int gfa_write_rgfa(const gfa_t *g, FILE *out);

/**
 * Release a graph and everything it owns. Safe to call with NULL
 * @param g Graph to destroy
 */
void gfa_destroy(gfa_t *g);

// Accessors / traversal

/**
 * Look up a segment's array index by id. O(1)
 * @param g Graph to query
 * @param id Segment id
 * @return The array index, or -1 if absent
 */
int32_t gfa_idx(const gfa_t *g, uint64_t id);

/**
 * Look up a segment by id
 * @param g Graph to query
 * @param id Segment id
 * @return Pointer to the segment, or NULL if absent
 */
gfa_seg_t *gfa_get(const gfa_t *g, uint64_t id);

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
 * @param g Graph to query (must have been read with GFA_LINKS)
 * @param v Segment index whose out-edges are wanted
 * @param arcs Set to an array of link indices leaving v; feed gfa_link_at()
 * @return Number of out-edges, or 0 if v has none or adjacency was not built
 */
int gfa_arcs(const gfa_t *g, int32_t v, const uint32_t **arcs);

/**
 * Test for a directed link v -> w. Requires GFA_LINKS. O(out-degree of v)
 * @param g Graph to query
 * @param v Source segment index
 * @param w Destination segment index
 * @return 1 if the link exists, else 0
 */
int gfa_has_arc(const gfa_t *g, int32_t v, int32_t w);

/**
 * Ordered segments of a path. The matching orientation chars are in
 * g->path_ori at the same offset; an entry may be GFA_NIL for an unresolved id
 * @param g Graph to query (must have been read with GFA_PATHS)
 * @param k Path index
 * @param segs Set to an array of segment indices; feed gfa_seg_at()
 * @return Number of segments in the path, or 0 if none / paths not built
 */
int gfa_path_segs(const gfa_t *g, int32_t k, const uint32_t **segs);

// Fragmented paths

/**
 * Chains of P-line fragments that spell one longer path.
 *
 * Graph builders such as vg emit a reference as several consecutive P lines
 * rather than one. A chain groups those fragments back together, in CSR form 
 * like everything else here: chain c owns the path indices frag[off[c] .. off[c+1]).
 */
typedef struct {
    char   **name;       // owned chain names, length n
    int32_t *frag;       // path indices, grouped by chain
    int32_t *off;        // length n + 1: chain c is frag[off[c] .. off[c+1])
    int32_t  n;          // number of chains
} gfa_merge_t;

/**
 * Group a graph's P-line fragments into chains.
 *
 * Fragments are selected by name: a path belongs to `key` when its name equals
 * `key` once a region suffix is stripped - "chr22:1000-2000" in the usual form,
 * or "chr22[1000]" as vg spells it - or when the last '#'-delimited field of
 * that base matches, so a PanSN name like "GRCh38#0#chr22:1000-2000" is found
 * by "chr22". Passing NULL selects every path and groups each base name separately.
 *
 * Selected fragments are ordered by the start offset in their name (those
 * without one keep file order and sort last), then chained: fragment A is
 * followed by B when a link joins A's last segment to B's first with matching
 * orientations, B has no other predecessor, and the link does not close a
 * cycle. Fragments that nothing joins simply end up in a chain of their own,
 * so every selected path appears in the result exactly once.
 *
 * A chain of several fragments is named after their shared base ("chr22"), or
 * "<base>_1", "<base>_2", ... when a base yields more than one such chain. A
 * chain of a single fragment keeps that path's original name.
 *
 * Chaining only follows the forward strand; a fragment stored reverse-
 * complemented relative to its neighbours is left unmerged. Requires the graph
 * to be read with GFA_LINKS | GFA_PATHS
 * @param g Graph to group
 * @param key Path name to select, or NULL for every path
 * @return The chains (release with gfa_merge_destroy), or NULL on failure
 */
gfa_merge_t *gfa_path_merge(const gfa_t *g, const char *key);

/**
 * Release a chain set and everything it owns. Safe to call with NULL
 * @param m Chain set to destroy
 */
void gfa_merge_destroy(gfa_merge_t *m);

/**
 * Segments of chain c, in order, as a flat list across its fragments.
 * Unresolved (GFA_NIL) entries are dropped
 * @param g Graph the chains came from
 * @param m Chain set from gfa_path_merge()
 * @param c Chain index
 * @param segs Set to a freshly allocated array of segment indices; the caller
 *             frees it. Set to NULL when the chain spells nothing
 * @return Number of segments written, or a negative AK_E* code
 */
int64_t gfa_merge_segs(const gfa_t *g, const gfa_merge_t *m, int32_t c, uint32_t **segs);

// Ranks

/**
 * Rank the segments against the graph's own paths.
 *
 * rGFA uses SR:i: to say how far a segment sits from the reference: 0 is the
 * reference backbone, anything higher came from a sample. Every segment that
 * any P line visits is stamped rank 0 and every other segment rank 1, so a
 * graph with no paths at all comes back entirely rank 1.
 *
 * Whether a reference arrives as one P line or as several vg-style fragments
 * makes no difference here, since every path counts as backbone either way;
 * gfa_path_merge() plus gfa_clear_paths()/gfa_add_path() is how you also
 * consolidate those fragments into one P line.
 *
 * Existing SR values are overwritten. gfa_read() calls this itself when the
 * file carried no SR tags, so calling it explicitly is how you re-rank a graph
 * whose tags you want to replace. Requires GFA_PATHS
 * @param g Graph to rank, modified in place
 * @return Number of segments left at rank 0, or a negative AK_E* code
 */
int64_t gfa_rank_paths(gfa_t *g);

/**
 * Rank the segments against a caller-supplied backbone: rank 0 wherever `on`
 * is set, rank 1 everywhere else. This is the general form of
 * gfa_rank_paths(), for a backbone that did not come from the P lines - a
 * traced reference sequence, for instance (see call_ref_fasta)
 * @param g Graph to rank, modified in place
 * @param on Flags of length gfa_n_seg(g); non-zero marks a backbone segment
 * @return Number of segments set to rank 0, or a negative AK_E* code
 */
int64_t gfa_rank_mark(gfa_t *g, const uint8_t *on);

// Rewriting the path block

/**
 * Drop every path in the graph.
 *
 * Path names are owned by the graph but borrowed by each segment's ref_name,
 * so this clears ref_name first - no segment is left pointing at a freed name.
 * The segments themselves, and their ranks, are untouched
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
 * install the walk that produced them in place of the old P lines. It is also
 * how a reference split across several fragments is consolidated back into one
 * P line. Requires GFA_PATHS
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
 * @param g Graph to order
 * @param order Caller-allocated array of length n_seg; filled with segment
 *              indices in topological order
 * @return The number of nodes placed before any cycle (n_seg if the graph is
 *         acyclic), or a negative AK_E* code on error
 */
int gfa_toposort(const gfa_t *g, int32_t *order);

#ifdef __cplusplus
}
#endif

#endif
