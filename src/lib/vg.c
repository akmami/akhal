#include "akhal/vg.h"
#include "akhal/error.h"

#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

// Refuse absurd message sizes (protobuf caps vg graphs at ~64 MB); this also
// stops a non-vg input from triggering a huge allocation off a garbage length.
#define VG_MAX_MSG (1u << 30)

// growable graph

static int reserve_node(vg_graph_t *g) {
    if (g->n_node < g->m_node) return AK_OK;
    int32_t m = g->m_node ? g->m_node << 1 : 4096;
    vg_node_t *p = (vg_node_t *)realloc(g->node, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    g->node = p;
    g->m_node = m;
    return AK_OK;
}

static int reserve_edge(vg_graph_t *g) {
    if (g->n_edge < g->m_edge) return AK_OK;
    int32_t m = g->m_edge ? g->m_edge << 1 : 4096;
    vg_edge_t *p = (vg_edge_t *)realloc(g->edge, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    g->edge = p;
    g->m_edge = m;
    return AK_OK;
}

static int reserve_path(vg_graph_t *g) {
    if (g->n_path < g->m_path) return AK_OK;
    int32_t m = g->m_path ? g->m_path << 1 : 64;
    vg_path_t *p = (vg_path_t *)realloc(g->path, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    g->path = p;
    g->m_path = m;
    return AK_OK;
}

static int reserve_step(vg_path_t *pt) {
    if (pt->n_step < pt->m_step) return AK_OK;
    int32_t m = pt->m_step ? pt->m_step << 1 : 16;
    vg_step_t *p = (vg_step_t *)realloc(pt->step, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    pt->step = p;
    pt->m_step = m;
    return AK_OK;
}

// protobuf decoding over an in-memory blob

typedef struct {
    const unsigned char *p, *end;
} pbuf;

// base-128 varint: 7 payload bits per byte, high bit means another byte follows
static int pb_varint(pbuf *b, uint64_t *out) {
    uint64_t v = 0;
    int shift = 0;
    while (b->p < b->end) {
        unsigned char c = *b->p++;
        v |= (uint64_t)(c & 0x7f) << shift;
        if (!(c & 0x80)) {
            *out = v;
            return 1;
        }
        shift += 7;
        if (shift >= 64) return 0;
    }
    return 0;
}

// wire type 2: a varint byte count followed by that many bytes
static int pb_bytes(pbuf *b, const unsigned char **start, uint64_t *len) {
    uint64_t l;
    if (!pb_varint(b, &l)) return 0;
    if ((uint64_t)(b->end - b->p) < l) return 0;
    *start = b->p;
    *len = l;
    b->p += l;
    return 1;
}

// wire types: 0 varint, 1 fixed64, 2 length-delimited, 5 fixed32
static int pb_skip(pbuf *b, int wire) {
    uint64_t tmp;
    const unsigned char *s;
    switch (wire) {
        case 0: return pb_varint(b, &tmp);
        case 1:
            if (b->end - b->p < 8) return 0;
            b->p += 8;
            return 1;
        case 2: return pb_bytes(b, &s, &tmp);
        case 5:
            if (b->end - b->p < 4) return 0;
            b->p += 4;
            return 1;
        default: return 0;   // start/end group (3,4) unused; anything else is bad
    }
}

// Position { int64 node_id = 1; bool is_reverse = 4; ... }
static void parse_position(const unsigned char *buf, uint64_t len, int64_t *node_id, int *is_rev) {
    pbuf b = { buf, buf + len };
    uint64_t tag, v;
    while (b.p < b.end) {
        if (!pb_varint(&b, &tag)) return;
        int field = (int)(tag >> 3), wire = (int)(tag & 7);
        if (field == 1 && wire == 0) {
            if (!pb_varint(&b, &v)) return;
            *node_id = (int64_t)v;
        } else if (field == 4 && wire == 0) {
            if (!pb_varint(&b, &v)) return;
            *is_rev = (v != 0);
        } else if (!pb_skip(&b, wire)) return;
    }
}

// Mapping { Position position = 1; ... } -> append one step to the path
static int parse_mapping(vg_path_t *pt, const unsigned char *buf, uint64_t len) {
    pbuf b = { buf, buf + len };
    int64_t node_id = 0;
    int is_rev = 0, have_pos = 0;
    uint64_t tag;
    while (b.p < b.end) {
        if (!pb_varint(&b, &tag)) break;
        int field = (int)(tag >> 3), wire = (int)(tag & 7);
        if (field == 1 && wire == 2) {
            const unsigned char *s;
            uint64_t l;
            if (!pb_bytes(&b, &s, &l)) break;
            parse_position(s, l, &node_id, &is_rev);
            have_pos = 1;
        } else if (!pb_skip(&b, wire)) break;
    }
    if (have_pos) {
        if (reserve_step(pt) != AK_OK) return AK_ENOMEM;
        pt->step[pt->n_step].node_id = node_id;
        pt->step[pt->n_step].is_reverse = is_rev;
        pt->n_step++;
    }
    return AK_OK;
}

// Path { string name = 1; repeated Mapping mapping = 2; bool is_circular = 3; }
static int parse_path(vg_graph_t *g, const unsigned char *buf, uint64_t len) {
    if (reserve_path(g) != AK_OK) return AK_ENOMEM;
    vg_path_t *pt = &g->path[g->n_path];
    memset(pt, 0, sizeof(*pt));

    pbuf b = { buf, buf + len };
    uint64_t tag, v;
    while (b.p < b.end) {
        if (!pb_varint(&b, &tag)) break;
        int field = (int)(tag >> 3), wire = (int)(tag & 7);
        if (field == 1 && wire == 2) {
            const unsigned char *s;
            uint64_t l;
            if (!pb_bytes(&b, &s, &l)) break;
            const char *q = ak_arena_put(&g->strs, (const char *)s, (size_t)l);
            if (!q) {
                free(pt->step);
                return AK_ENOMEM;
            }
            pt->name = q;
        } else if (field == 2 && wire == 2) {
            const unsigned char *s;
            uint64_t l;
            if (!pb_bytes(&b, &s, &l)) break;
            if (parse_mapping(pt, s, l) == AK_ENOMEM) {
                free(pt->step);
                return AK_ENOMEM;
            }
        } else if (field == 3 && wire == 0) {
            if (!pb_varint(&b, &v)) break;
            pt->is_circular = (v != 0);
        } else if (!pb_skip(&b, wire)) break;
    }
    if (!pt->name) {
        pt->name = "";   // an unnamed path, without an allocation of its own
    }
    g->n_path++;
    return AK_OK;
}

// Node { string sequence = 1; int64 id = 3; ... }
static int parse_node(vg_graph_t *g, const unsigned char *buf, uint64_t len) {
    if (reserve_node(g) != AK_OK) return AK_ENOMEM;
    vg_node_t *n = &g->node[g->n_node];
    memset(n, 0, sizeof(*n));

    pbuf b = { buf, buf + len };
    uint64_t tag, v;
    while (b.p < b.end) {
        if (!pb_varint(&b, &tag)) break;
        int field = (int)(tag >> 3), wire = (int)(tag & 7);
        if (field == 1 && wire == 2) {
            const unsigned char *s;
            uint64_t l;
            if (!pb_bytes(&b, &s, &l)) break;
            // a repeated sequence field simply overwrites; the earlier copy
            // stays in the arena, which costs a few bytes and cannot happen
            // in a file vg wrote
            const char *q = ak_arena_put(&g->strs, (const char *)s, (size_t)l);
            if (!q) return AK_ENOMEM;
            n->seq = q;
            n->seq_len = (uint32_t)l;
        } else if (field == 3 && wire == 0) {
            if (!pb_varint(&b, &v)) break;
            n->id = (int64_t)v;
        } else if (!pb_skip(&b, wire)) break;
    }
    g->n_node++;
    return AK_OK;
}

// Edge { int64 from=1; int64 to=2; bool from_start=3; bool to_end=4; int32 overlap=5; }
static int parse_edge(vg_graph_t *g, const unsigned char *buf, uint64_t len) {
    if (reserve_edge(g) != AK_OK) return AK_ENOMEM;
    vg_edge_t *e = &g->edge[g->n_edge];
    memset(e, 0, sizeof(*e));

    pbuf b = { buf, buf + len };
    uint64_t tag, v;
    while (b.p < b.end) {
        if (!pb_varint(&b, &tag)) break;
        int field = (int)(tag >> 3), wire = (int)(tag & 7);
        if (wire == 0) {
            if (!pb_varint(&b, &v)) break;
            switch (field) {
                case 1: e->from = (int64_t)v; break;
                case 2: e->to = (int64_t)v; break;
                case 3: e->from_start = (v != 0); break;
                case 4: e->to_end = (v != 0); break;
                case 5: e->overlap = (int32_t)v; break;
                default: break;
            }
        } else if (!pb_skip(&b, wire)) break;
    }
    g->n_edge++;
    return AK_OK;
}

// Graph { repeated Node node=1; repeated Edge edge=2; repeated Path path=3; }
static int parse_graph_blob(vg_graph_t *g, const unsigned char *buf, uint64_t len) {
    pbuf b = { buf, buf + len };
    uint64_t tag;
    while (b.p < b.end) {
        if (!pb_varint(&b, &tag)) return AK_OK;   // malformed / tag blob: stop leniently
        int field = (int)(tag >> 3), wire = (int)(tag & 7);
        if (wire == 2) {
            const unsigned char *s;
            uint64_t l;
            if (!pb_bytes(&b, &s, &l)) return AK_OK;
            int rc = AK_OK;
            if (field == 1) {
                rc = parse_node(g, s, l);
            } else if (field == 2) {
                rc = parse_edge(g, s, l);
            } else if (field == 3) {
                rc = parse_path(g, s, l);
            }
            if (rc == AK_ENOMEM) return AK_ENOMEM;
        } else if (!pb_skip(&b, wire)) {
            return AK_OK;   // invalid wire type (e.g. a type-tag message): skip blob
        }
    }
    return AK_OK;
}

// Hand back the slack the doubling left behind. Each array grows by doubling,
// so on arrival it holds up to twice what it needs - on a whole-genome graph
// that is tens of gigabytes of allocated, untouched pages. Shrinking is cheap:
// glibc remaps rather than copies, about 5 ms per GB. A failed shrink is not
// an error, since the oversized block is still perfectly valid.
static void shrink(void **p, size_t n, size_t esz) {
    if (!*p || n == 0) return;
    void *q = realloc(*p, n * esz);
    if (q) *p = q;
}

// buffered gzip input

typedef struct {
    gzFile        gz;
    unsigned char buf[1 << 16];
    int           len, pos;
} vgin;

static int vin_fill(vgin *r) {
    r->len = gzread(r->gz, r->buf, (unsigned)sizeof(r->buf));
    r->pos = 0;
    return r->len;
}

// all-or-nothing: 0 if the stream ends before n bytes are available
static int vin_read(vgin *r, unsigned char *dst, uint64_t n) {
    while (n) {
        if (r->pos >= r->len) {
            if (vin_fill(r) <= 0) return 0;
        }
        uint64_t avail = (uint64_t)(r->len - r->pos);
        uint64_t take = n < avail ? n : avail;
        memcpy(dst, r->buf + r->pos, take);
        r->pos += (int)take;
        dst += take;
        n -= take;
    }
    return 1;
}

// 0 only at a clean message boundary, -1 on a truncated or overlong varint
static int vin_varint(vgin *r, uint64_t *out) {
    uint64_t v = 0;
    int shift = 0, first = 1;
    for (;;) {
        if (r->pos >= r->len) {
            if (vin_fill(r) <= 0) return first ? 0 : -1;
        }
        unsigned char c = r->buf[r->pos++];
        first = 0;
        v |= (uint64_t)(c & 0x7f) << shift;
        if (!(c & 0x80)) {
            *out = v;
            return 1;
        }
        shift += 7;
        if (shift >= 64) return -1;
    }
}

// public API

// read a .vg file into an accumulated graph; see akhal/vg.h
vg_graph_t *vg_read(const char *fn) {
    gzFile gz = gzopen(fn, "rb");
    if (!gz) {
        ak_log(AK_LOG_ERROR, "vg", "could not open %s", fn);
        return NULL;
    }
    // zlib reads the file 8 KiB at a time by default, which turns a large .vg
    // into millions of syscalls - punishing on a network filesystem
    gzbuffer(gz, 1 << 22);   /* 4 MiB */

    vg_graph_t *g = (vg_graph_t *)calloc(1, sizeof(vg_graph_t));
    if (!g) {
        gzclose(gz);
        ak_log(AK_LOG_ERROR, "vg", "out of memory");
        return NULL;
    }

    vgin r;
    r.gz = gz;
    r.len = 0;
    r.pos = 0;

    unsigned char *msg = NULL;
    uint64_t msgcap = 0;
    int oom = 0;

    for (;;) {
        uint64_t count;
        int v = vin_varint(&r, &count);
        if (v == 0) break;   // clean EOF at a group boundary
        if (v < 0) {
            ak_log(AK_LOG_WARN, "vg", "truncated group header");
            break;
        }

        int stop = 0;
        for (uint64_t i = 0; i < count && !stop; i++) {
            uint64_t len;
            if (vin_varint(&r, &len) <= 0) {
                ak_log(AK_LOG_WARN, "vg", "truncated message length");
                stop = 1;
                break;
            }
            if (len > VG_MAX_MSG) {
                ak_log(AK_LOG_WARN, "vg", "implausible message size; stopping");
                stop = 1;
                break;
            }

            if (len > msgcap) {
                uint64_t nc = msgcap ? msgcap : 4096;
                while (nc < len) nc <<= 1;
                unsigned char *tmp = (unsigned char *)realloc(msg, nc);
                if (!tmp) {
                    oom = 1;
                    stop = 1;
                    break;
                }
                msg = tmp;
                msgcap = nc;
            }
            if (!vin_read(&r, msg, len)) {
                ak_log(AK_LOG_WARN, "vg", "truncated message body");
                stop = 1;
                break;
            }
            if (parse_graph_blob(g, msg, len) == AK_ENOMEM) {
                oom = 1;
                stop = 1;
                break;
            }
        }
        if (stop) break;
    }

    free(msg);
    gzclose(gz);

    if (oom) {
        vg_graph_destroy(g);
        ak_log(AK_LOG_ERROR, "vg", "out of memory");
        return NULL;
    }

    // the graph is complete; give back what the doubling over-reserved
    shrink((void **)&g->node, (size_t)g->n_node, sizeof(*g->node));
    g->m_node = g->n_node;
    shrink((void **)&g->edge, (size_t)g->n_edge, sizeof(*g->edge));
    g->m_edge = g->n_edge;
    shrink((void **)&g->path, (size_t)g->n_path, sizeof(*g->path));
    g->m_path = g->n_path;

    return g;
}

// free a vg graph; see akhal/vg.h
void vg_graph_destroy(vg_graph_t *g) {
    if (!g) return;
    ak_arena_destroy(&g->strs);   // every sequence and path name, in one go
    free(g->node);
    free(g->edge);
    for (int32_t i = 0; i < g->n_path; i++) free(g->path[i].step);
    free(g->path);
    free(g);
}
