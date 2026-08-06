#ifndef AKHAL_VG_H
#define AKHAL_VG_H

#include <stdint.h>

/**
 * Reader for vg's native ".vg" format.
 *
 * A .vg file is a gzip/BGZF-compressed stream of length-delimited Protobuf
 * `Graph` messages (each holding repeated Node/Edge/Path), grouped as
 * [varint count][count x (varint length + message bytes)]. Rather than depend
 * on protobuf/libvgio, this module decodes the handful of message types needed
 * for GFA conversion directly from the wire format (field numbers taken from
 * vg.proto), and accumulates them into one in-memory graph — mirroring what
 * `vg view -g` does. Only zlib is required, for decompression.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int64_t  id;        // node id (positive, nonzero)
    char    *seq;       // owned sequence, or NULL if absent
    uint32_t seq_len;   // sequence length
} vg_node_t;

typedef struct {
    int64_t from, to;   // endpoint node ids
    int     from_start; // edge leaves the 5' (start) side of `from`
    int     to_end;     // edge enters the 3' (end) side of `to`
    int32_t overlap;    // overlap length in bp
} vg_edge_t;

typedef struct {
    int64_t node_id;    // visited node
    int     is_reverse; // visited in reverse-complement orientation
} vg_step_t;

typedef struct {
    char      *name;        // owned path name
    vg_step_t *step;        // ordered visits
    int32_t    n_step, m_step;
    int        is_circular;
} vg_path_t;

typedef struct {
    vg_node_t *node; int32_t n_node, m_node;
    vg_edge_t *edge; int32_t n_edge, m_edge;
    vg_path_t *path; int32_t n_path, m_path;
} vg_graph_t;

/**
 * Read a .vg file into a single accumulated graph.
 * @param fn Path to the .vg file (gzip/BGZF-compressed or raw).
 * @return The graph (release with vg_graph_destroy), or NULL on error.
 */
vg_graph_t *vg_read(const char *fn);

/**
 * Release a graph and everything it owns. Safe with NULL.
 * @param g Graph to destroy.
 */
void vg_graph_destroy(vg_graph_t *g);

#ifdef __cplusplus
}
#endif

#endif  // AKHAL_VG_H
