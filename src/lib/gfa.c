#include "akhal/gfa.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/error.h"

#include "khashl.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// id -> array index
KHASHL_MAP_INIT(KH_LOCAL, idxmap_t, idxmap, uint64_t, uint32_t,
                kh_hash_uint64, kh_eq_generic)

// growth helpers

/** Ensure room for one more segment. @return AK_OK or AK_ENOMEM. */
static int reserve_seg(gfa_t *g) {
    if (g->n_seg < g->m_seg) return AK_OK;
    int32_t m = g->m_seg ? g->m_seg << 1 : 1024;
    gfa_seg_t *p = (gfa_seg_t *)realloc(g->seg, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    g->seg = p;
    g->m_seg = m;
    return AK_OK;
}

/** Ensure room for one more link. @return AK_OK or AK_ENOMEM. */
static int reserve_link(gfa_t *g) {
    if (g->n_link < g->m_link) return AK_OK;
    int32_t m = g->m_link ? g->m_link << 1 : 1024;
    gfa_link_t *p = (gfa_link_t *)realloc(g->link, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    g->link = p;
    g->m_link = m;
    return AK_OK;
}

/** Ensure room for one more path (name, length, offset). @return AK_OK or AK_ENOMEM. */
static int reserve_path(gfa_t *g) {
    if (g->n_path < g->m_path) return AK_OK;
    int32_t m = g->m_path ? g->m_path << 1 : 16;
    char **np = (char **)realloc(g->path, (size_t)m * sizeof(*np));
    if (!np) return AK_ENOMEM;
    g->path = np;
    uint64_t *nl = (uint64_t *)realloc(g->path_len, (size_t)m * sizeof(*nl));
    if (!nl) return AK_ENOMEM;
    g->path_len = nl;
    int32_t *no = (int32_t *)realloc(g->path_off, ((size_t)m + 1) * sizeof(*no));
    if (!no) return AK_ENOMEM;
    g->path_off = no;
    g->m_path = m;
    return AK_OK;
}

/** Grow the flat path-membership arrays to hold one more entry. @return AK_OK or AK_ENOMEM. */
static int reserve_pathseg(gfa_t *g) {
    if ((int64_t)g->n_path_seg < g->m_path_seg) return AK_OK;
    int32_t m = g->m_path_seg ? g->m_path_seg << 1 : 4096;
    uint32_t *ns = (uint32_t *)realloc(g->path_seg, (size_t)m * sizeof(*ns));
    if (!ns) return AK_ENOMEM;
    g->path_seg = ns;
    char *no = (char *)realloc(g->path_ori, (size_t)m * sizeof(*no));
    if (!no) return AK_ENOMEM;
    g->path_ori = no;
    g->m_path_seg = m;
    return AK_OK;
}

// line handlers

/**
 * Split "TAG:TYPE:VALUE" in place; VALUE keeps any embedded ':'.
 * @param token Token to split (modified in place).
 * @param tag Set to the TAG part.
 * @param type Set to the TYPE part.
 * @param value Set to the VALUE part.
 * @return 1 on success, 0 if the token is not a TAG:TYPE:VALUE triple.
 */
static int split_tag(char *token, char **tag, char **type, char **value) {
    char *c1 = strchr(token, ':');
    if (!c1) return 0;
    char *c2 = strchr(c1 + 1, ':');
    if (!c2) return 0;
    *c1 = '\0';
    *c2 = '\0';
    *tag   = token;
    *type  = c1 + 1;
    *value = c2 + 1;
    return 1;
}

/**
 * Parse one S (segment) line into a new node and index it.
 * @param g Graph being built.
 * @param line The S line (modified in place by tokenizing).
 * @param h The id -> index map.
 * @return AK_OK, AK_EFORMAT, or AK_ENOMEM.
 */
static int handle_S(gfa_t *g, char *line, idxmap_t *h) {
    if (reserve_seg(g) != AK_OK) return AK_ENOMEM;
    gfa_seg_t *s = &g->seg[g->n_seg];
    memset(s, 0, sizeof(*s));
    s->rank = -1;   // -1 until an SR tag says otherwise

    char *save;
    char *tok = strtok_r(line, "\t", &save);   // 'S'
    tok = strtok_r(NULL, "\t", &save);         // id
    if (!tok) { ak_log(AK_LOG_WARN, "gfa", "S line without id"); return AK_EFORMAT; }
    s->id = strtoull(tok, NULL, 10);

    tok = strtok_r(NULL, "\t", &save);         // sequence
    if (!tok || tok[0] == '\0') {
        ak_log(AK_LOG_WARN, "gfa", "segment %llu has empty sequence", (unsigned long long)s->id);
        s->seq = NULL;
        s->len = 0;
    } else {
        size_t n = strlen(tok);
        s->seq = (char *)malloc(n + 1);
        if (!s->seq) return AK_ENOMEM;
        memcpy(s->seq, tok, n + 1);
        s->len = (uint32_t)n;
    }
    s->start = 0;
    s->end = (int32_t)s->len;

    // optional tags: SN:Z:name  SO:i:offset  SR:i:rank
    while ((tok = strtok_r(NULL, "\t", &save)) != NULL) {
        char *tag, *type, *val;
        if (!split_tag(tok, &tag, &type, &val)) continue;
        if (!strcmp(tag, "SO") && !strcmp(type, "i")) {
            s->start = atoi(val);
            s->end   = s->start + (int32_t)s->len;
        } else if (!strcmp(tag, "SR") && !strcmp(type, "i")) {
            s->rank = atoi(val);
        }
        // SN is handled via path names; segment->ref_name is set there.
    }

    int absent;
    khint_t k = idxmap_put(h, s->id, &absent);
    if (absent) kh_val(h, k) = (uint32_t)g->n_seg;
    else ak_log(AK_LOG_WARN, "gfa", "duplicate segment id %llu", (unsigned long long)s->id);

    g->n_seg++;
    return AK_OK;
}

/**
 * Parse one L (link) line: optional validation, then record the edge.
 * @param g Graph being built.
 * @param line The L line.
 * @param h The id -> index map.
 * @param flags Active GFA_* flags (GFA_LINKS / GFA_VALIDATE).
 * @return AK_OK, AK_EFORMAT, or AK_ENOMEM.
 */
static int handle_L(gfa_t *g, char *line, idxmap_t *h, int flags) {
    uint64_t id1, id2;
    char st1, st2;
    size_t overlap = 0;
    if (sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM",
               &id1, &st1, &id2, &st2, &overlap) < 4) {
        ak_log(AK_LOG_WARN, "gfa", "malformed L line");
        return AK_EFORMAT;
    }

    khint_t k1 = idxmap_get(h, id1);
    khint_t k2 = idxmap_get(h, id2);
    int have1 = (k1 < kh_end(h)), have2 = (k2 < kh_end(h));

    if (flags & GFA_VALIDATE) {
        if (!have1) ak_log(AK_LOG_WARN, "gfa", "L references unknown segment %llu", (unsigned long long)id1);
        if (!have2) ak_log(AK_LOG_WARN, "gfa", "L references unknown segment %llu", (unsigned long long)id2);
        if (have1 && have2 && overlap > 0) {
            gfa_seg_t *a = &g->seg[kh_val(h, k1)];
            gfa_seg_t *b = &g->seg[kh_val(h, k2)];
            if (a->seq && b->seq &&
                overlap < a->len && overlap < b->len &&
                strncmp(a->seq, b->seq + (b->len - overlap), overlap) != 0) {
                ak_log(AK_LOG_WARN, "gfa", "overlap mismatch %lu -> %lu (len %lu)", (unsigned long)id1, (unsigned long)id2, (unsigned long)overlap);
            }
        }
    }

    if (!(flags & GFA_LINKS)) return AK_OK;
    if (!have1 || !have2) return AK_OK;

    if (reserve_link(g) != AK_OK) return AK_ENOMEM;
    uint32_t v = kh_val(h, k1), w = kh_val(h, k2);
    g->link[g->n_link].v = v;
    g->link[g->n_link].w = w;
    g->link[g->n_link].overlap = (uint32_t)overlap;
    g->link[g->n_link].from_orient = (st1 == '-') ? '-' : '+';
    g->link[g->n_link].to_orient   = (st2 == '-') ? '-' : '+';
    g->n_link++;

    g->seg[v].out_degree++;
    g->seg[w].in_degree++;
    return AK_OK;
}

/**
 * Parse one P (path) line: record the ordered, oriented segment membership and
 * lay out reference coordinates. Without GFA_PATHS only n_path_seg (the
 * occurrence count) is maintained.
 * @param g Graph being built.
 * @param line The P line.
 * @param h The id -> index map.
 * @param flags Active GFA_* flags.
 * @return AK_OK, AK_EFORMAT, or AK_ENOMEM.
 */
static int handle_P(gfa_t *g, char *line, idxmap_t *h, int flags) {
    char *save;
    char *tok = strtok_r(line, "\t", &save);        // 'P'
    tok = strtok_r(NULL, "\t", &save);              // path name
    if (!tok) {
        ak_log(AK_LOG_WARN, "gfa", "P line without name");
        return AK_EFORMAT;
    }

    char *name = NULL;
    int32_t pi = -1;
    if (flags & GFA_PATHS) {
        if (reserve_path(g) != AK_OK) return AK_ENOMEM;
        pi = g->n_path;
        name = strdup(tok);
        if (!name) return AK_ENOMEM;
        g->path[pi] = name;
        g->path_len[pi] = 0;
        g->path_off[pi] = (int32_t)g->n_path_seg;   // start of this path's slice
        g->n_path++;
    }

    tok = strtok_r(NULL, "\t", &save);              // comma list of segments
    if (!tok) {
        if (flags & GFA_PATHS) {
            g->path_off[pi + 1] = (int32_t)g->n_path_seg;   // empty slice
            ak_log(AK_LOG_WARN, "gfa", "path %s has no segments", name ? name : "?");
        }
        return AK_OK;
    }

    char *sp;
    char *seg_tok = strtok_r(tok, ",", &sp);
    int32_t ref_pos = 0;

    while (seg_tok) {
        // orientation suffix ('+'/'-'), stripped before the id is parsed
        size_t len = strlen(seg_tok);
        char ori = '+';
        if (len && (seg_tok[len - 1] == '+' || seg_tok[len - 1] == '-')) {
            ori = seg_tok[len - 1];
            seg_tok[len - 1] = '\0';
        }

        uint64_t sid = strtoull(seg_tok, NULL, 10);
        khint_t k = idxmap_get(h, sid);
        uint32_t si = (k < kh_end(h)) ? kh_val(h, k) : GFA_NIL;

        if (si == GFA_NIL && (flags & (GFA_PATHS | GFA_VALIDATE)))
            ak_log(AK_LOG_WARN, "gfa", "segment %llu in path %s not found", (unsigned long long)sid, name ? name : "?");

        if (flags & GFA_PATHS) {
            if (reserve_pathseg(g) != AK_OK) return AK_ENOMEM;
            g->path_seg[g->n_path_seg] = si;
            g->path_ori[g->n_path_seg] = ori;
            if (si != GFA_NIL) {
                gfa_seg_t *cur = &g->seg[si];
                cur->ref_name = name;
                cur->start = ref_pos;
                ref_pos += (int32_t)cur->len;
                cur->end = ref_pos;
                g->path_len[pi] += cur->len;
            }
        }
        g->n_path_seg++;
        seg_tok = strtok_r(NULL, ",", &sp);
    }

    if (flags & GFA_PATHS)
        g->path_off[pi + 1] = (int32_t)g->n_path_seg;   // end of this path's slice
    return AK_OK;
}

// CSR out-adjacency

/**
 * Build the out-adjacency (arc_off / arc) from the link array.
 * @param g Graph whose links have been read.
 * @return AK_OK or AK_ENOMEM.
 */
static int build_arcs(gfa_t *g) {
    if (g->n_seg <= 0 || g->n_link <= 0) return AK_OK;

    g->arc_off = (int32_t *)calloc((size_t)g->n_seg + 1, sizeof(int32_t));
    g->arc     = (uint32_t *)malloc((size_t)g->n_link * sizeof(uint32_t));
    if (!g->arc_off || !g->arc) return AK_ENOMEM;

    for (int32_t k = 0; k < g->n_link; k++)
        g->arc_off[g->link[k].v + 1]++;
    for (int32_t i = 0; i < g->n_seg; i++)
        g->arc_off[i + 1] += g->arc_off[i];

    int32_t *cursor = (int32_t *)malloc((size_t)g->n_seg * sizeof(int32_t));
    if (!cursor) return AK_ENOMEM;
    for (int32_t i = 0; i < g->n_seg; i++) cursor[i] = g->arc_off[i];
    for (int32_t k = 0; k < g->n_link; k++) {
        uint32_t v = g->link[k].v;
        g->arc[cursor[v]++] = (uint32_t)k;
    }
    free(cursor);
    return AK_OK;
}

// public API

/** Read an (r)GFA into a graph; see akhal/gfa.h. */
gfa_t *gfa_read(const char *fn, int flags) {
    ak_file *f = ak_open(fn);
    if (!f) return NULL;

    gfa_t *g = (gfa_t *)calloc(1, sizeof(gfa_t));
    if (!g) {
        ak_close(f); ak_log(AK_LOG_ERROR, "gfa", "out of memory");
        return NULL;
    }

    idxmap_t *h = idxmap_init();
    if (!h) {
        free(g);
        ak_close(f);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory");
        return NULL;
    }
    g->idx = h;
    g->flags = flags;

    kstring_t ks = KS_INIT;
    int rc = AK_OK;
    long len;

    while ((len = ak_getline(f, &ks)) >= 0) {
        if (len == 0) continue;
        switch (ks.s[0]) {
            case 'S': rc = handle_S(g, ks.s, h); break;
            case 'L': rc = handle_L(g, ks.s, h, flags); break;
            case 'P': rc = handle_P(g, ks.s, h, flags); break;
            default:  rc = AK_OK; break;   // ignore H, W, comments, etc.
        }
        if (rc == AK_ENOMEM) break;   // only allocation failures are fatal
        rc = AK_OK;
    }

    ks_free(&ks);
    ak_close(f);

    if (rc == AK_ENOMEM) {
        gfa_destroy(g);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory");
        return NULL;
    }

    if (flags & GFA_LINKS) {
        if (build_arcs(g) != AK_OK) {
            gfa_destroy(g);
            ak_log(AK_LOG_ERROR, "gfa", "out of memory building adjacency");
            return NULL;
        }
    }

    return g;
}

/** 
 * @brief Free a graph and everything it owns; see akhal/gfa.h.
 */
void gfa_destroy(gfa_t *g) {
    if (!g) return;
    for (int32_t i = 0; i < g->n_seg; i++) free(g->seg[i].seq);
    free(g->seg);
    free(g->link);
    free(g->arc);
    free(g->arc_off);
    if (g->path) {
        for (int32_t i = 0; i < g->n_path; i++) free(g->path[i]);
        free(g->path);
    }
    free(g->path_len);
    free(g->path_off);
    free(g->path_seg);
    free(g->path_ori);
    if (g->idx) idxmap_destroy((idxmap_t *)g->idx);
    free(g);
}

/** 
 * @brief Segment index for an id, or -1 if absent. 
 */
int32_t gfa_idx(const gfa_t *g, uint64_t id) {
    idxmap_t *h = (idxmap_t *)g->idx;
    khint_t k = idxmap_get(h, id);
    return (k < kh_end(h)) ? (int32_t)kh_val(h, k) : -1;
}

/** 
 * @brief Segment for an id, or NULL if absent. 
 */
gfa_seg_t *gfa_get(const gfa_t *g, uint64_t id) {
    int32_t i = gfa_idx(g, id);
    return (i < 0) ? NULL : &g->seg[i];
}

/** 
 * @brief Out-edges of segment v; see akhal/gfa.h. 
 */
int gfa_arcs(const gfa_t *g, int32_t v, const uint32_t **arcs) {
    if (!g->arc_off || v < 0 || v >= g->n_seg) { *arcs = NULL; return 0; }
    int32_t beg = g->arc_off[v], end = g->arc_off[v + 1];
    *arcs = g->arc + beg;
    return (int)(end - beg);
}

/** 
 * @return 1 if a link v -> w exists, else 0. 
 */
int gfa_has_arc(const gfa_t *g, int32_t v, int32_t w) {
    const uint32_t *a;
    int n = gfa_arcs(g, v, &a);
    for (int k = 0; k < n; k++)
        if ((int32_t)g->link[a[k]].w == w) return 1;
    return 0;
}

/** 
 * @brief Ordered segments of path k; see akhal/gfa.h. 
 */
int gfa_path_segs(const gfa_t *g, int32_t k, const uint32_t **segs) {
    if (!g->path_off || k < 0 || k >= g->n_path) { *segs = NULL; return 0; }
    int32_t beg = g->path_off[k], end = g->path_off[k + 1];
    *segs = g->path_seg + beg;
    return (int)(end - beg);
}

// topological sort

/**
 * @return 1 if segment a's sequence sorts before b's, alphabetically.
 * Ordering by sequence content (rather than id) makes the result independent
 * of the input's node numbering: two graphs identical in topology and content
 * sort the same way. A NULL/empty sequence sorts first.
 */
static int seq_lt(const gfa_t *g, int32_t a, int32_t b) {
    const char *sa = g->seg[a].seq ? g->seg[a].seq : "";
    const char *sb = g->seg[b].seq ? g->seg[b].seq : "";
    return strcmp(sa, sb) < 0;
}

/** 
 * @brief Sift the last heap element up to restore the min-heap order.
 */
static void heap_push(const gfa_t *g, int32_t *heap, int *hn, int32_t v) {
    int i = (*hn)++;
    heap[i] = v;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (!seq_lt(g, heap[i], heap[p])) break;
        int32_t t = heap[i]; heap[i] = heap[p]; heap[p] = t;
        i = p;
    }
}

/** 
 * @brief Pop and return the alphabetically-smallest node from the heap. 
 */
static int32_t heap_pop(const gfa_t *g, int32_t *heap, int *hn) {
    int32_t top = heap[0];
    int n = --(*hn);
    heap[0] = heap[n];
    int i = 0;
    for (;;) {
        int l = 2 * i + 1, r = 2 * i + 2, m = i;
        if (l < n && seq_lt(g, heap[l], heap[m])) m = l;
        if (r < n && seq_lt(g, heap[r], heap[m])) m = r;
        if (m == i) break;
        int32_t t = heap[i]; heap[i] = heap[m]; heap[m] = t;
        i = m;
    }
    return top;
}

/** 
 * @brief Topological order with alphabetical id tie-break; see akhal/gfa.h. 
 */
int gfa_toposort(const gfa_t *g, int32_t *order) {
    int32_t n = g->n_seg;
    if (n == 0) return 0;
    if (!(g->flags & GFA_LINKS)) {
        ak_log(AK_LOG_ERROR, "gfa", "toposort requires the graph to be read with GFA_LINKS");
        return AK_EINVAL;
    }

    int32_t *indeg = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t *heap  = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!indeg || !heap) { free(indeg); free(heap); return AK_ENOMEM; }

    for (int32_t i = 0; i < n; i++) indeg[i] = g->seg[i].in_degree;

    int hn = 0;
    for (int32_t i = 0; i < n; i++)
        if (indeg[i] == 0) heap_push(g, heap, &hn, i);

    int32_t placed = 0;
    while (hn > 0) {
        int32_t u = heap_pop(g, heap, &hn);
        order[placed++] = u;
        const uint32_t *arcs;
        int na = gfa_arcs(g, u, &arcs);
        for (int k = 0; k < na; k++) {
            uint32_t w = g->link[arcs[k]].w;
            if (--indeg[w] == 0) heap_push(g, heap, &hn, (int32_t)w);
        }
    }

    int32_t acyclic = placed;
    if (placed < n) {
        // Remaining nodes are inside cycles; append them alphabetically.
        for (int32_t i = 0; i < n; i++)
            if (indeg[i] > 0) heap_push(g, heap, &hn, i);
        while (hn > 0) order[placed++] = heap_pop(g, heap, &hn);
    }

    free(indeg);
    free(heap);
    return acyclic;
}
