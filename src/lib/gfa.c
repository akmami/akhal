#include "akhal/gfa.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/error.h"

#include "khashl.h"

#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

// id -> array index
KHASHL_MAP_INIT(KH_LOCAL, idxmap_t, idxmap, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

// growth helpers

static int reserve_seg(gfa_t *g) {
    if (g->n_seg < g->m_seg) return AK_OK;
    int32_t m_seg = g->m_seg ? g->m_seg << 1 : 1024;
    gfa_seg_t *seg = (gfa_seg_t *)realloc(g->seg, (size_t)m_seg * sizeof(gfa_seg_t));
    if (!seg) return AK_ENOMEM;
    g->seg = seg;
    g->m_seg = m_seg;
    return AK_OK;
}

static inline int reserve_link(gfa_t *g) {
    if (g->n_link < g->m_link) return AK_OK;
    int32_t m_link = g->m_link ? g->m_link << 1 : 1024;
    gfa_link_t *link = (gfa_link_t *)realloc(g->link, (size_t)m_link * sizeof(gfa_link_t));
    if (!link) return AK_ENOMEM;
    g->link = link;
    g->m_link = m_link;
    return AK_OK;
}

static int reserve_path(gfa_t *g) {
    if (g->n_path < g->m_path) return AK_OK;
    int32_t m_path = g->m_path ? g->m_path << 1 : 16;
    char **path = (char **)realloc(g->path, (size_t)m_path * sizeof(char *));
    if (!path) return AK_ENOMEM;
    g->path = path;
    uint64_t *path_len = (uint64_t *)realloc(g->path_len, (size_t)m_path * sizeof(uint64_t));
    if (!path_len) return AK_ENOMEM;
    g->path_len = path_len;
    int32_t *path_off = (int32_t *)realloc(g->path_off, ((size_t)m_path + 1) * sizeof(int32_t));
    if (!path_off) return AK_ENOMEM;
    g->path_off = path_off;
    g->m_path = m_path;
    return AK_OK;
}

static int reserve_pathseg(gfa_t *g) {
    if ((int64_t)g->n_path_seg < g->m_path_seg) return AK_OK;
    int32_t m_path_seg = g->m_path_seg ? g->m_path_seg << 1 : 4096;
    uint32_t *path_seg = (uint32_t *)realloc(g->path_seg, (size_t)m_path_seg * sizeof(uint32_t));
    if (!path_seg) return AK_ENOMEM;
    g->path_seg = path_seg;
    char *path_ori = (char *)realloc(g->path_ori, (size_t)m_path_seg * sizeof(char));
    if (!path_ori) return AK_ENOMEM;
    g->path_ori = path_ori;
    g->m_path_seg = m_path_seg;
    return AK_OK;
}

// Hand back the slack the doubling left behind.
//
// Each array grows by doubling, so on arrival it holds up to twice what it
// needs - on a whole-genome graph that is tens of gigabytes of allocated,
// untouched pages. One realloc down to the exact count reclaims it, and
// shrinking is cheap: glibc remaps rather than copies, about 5 ms per GB.
// A failed shrink is not an error, since the oversized block is still valid.
static void shrink(void **ptr, size_t n, size_t esz) {
    if (!*ptr || n == 0) return;
    void *new_ptr = realloc(*ptr, n * esz);
    if (new_ptr) *ptr = new_ptr;
}

// line handlers

// splits "TAG:TYPE:VALUE" in place; VALUE keeps any embedded ':'
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

// parses one S line into a new node and indexes it; tokenizes line in place
static int handle_S(gfa_t *g, char *line, idxmap_t *h) {
    if (reserve_seg(g) != AK_OK) return AK_ENOMEM;
    gfa_seg_t *s = &g->seg[g->n_seg];
    memset(s, 0, sizeof(*s));
    s->rank = -1;       // -1 until an SR tag says otherwise
    s->ref_path = -1;   // memset made it 0, which would mean path 0

    char *save;
    char *tok = strtok_r(line, "\t", &save);   // 'S'
    tok = strtok_r(NULL, "\t", &save);         // id
    if (!tok) {
        ak_log(AK_LOG_WARN, "gfa", "S line without id");
        return AK_EFORMAT;
    }
    s->id = strtoull(tok, NULL, 10);

    tok = strtok_r(NULL, "\t", &save);         // sequence
    if (!tok || tok[0] == '\0') {
        ak_log(AK_LOG_WARN, "gfa", "segment %llu has empty sequence", (unsigned long long)s->id);
        s->seq = NULL;
        s->len = 0;
    } else if (g->flags & GFA_SEQ) {
        if (gfa_seg_set_seq(g, s, tok, strlen(tok)) != AK_OK) return AK_ENOMEM;
    } else {
        // We store the length is still cheap and still needed, so it's recorded anyway.
        s->seq = NULL;
        s->len = (uint32_t)strlen(tok);
    }
    s->start = 0;

    // optional tags: SN:Z:name  SO:i:offset  SR:i:rank
    while ((tok = strtok_r(NULL, "\t", &save)) != NULL) {
        char *tag, *type, *val;
        if (!split_tag(tok, &tag, &type, &val)) continue;
        if (!strcmp(tag, "SO") && !strcmp(type, "i")) {
            s->start = atoi(val);
        } else if (!strcmp(tag, "SR") && !strcmp(type, "i")) {
            s->rank = atoi(val);
            g->has_sr = 1;   // the file ranks itself; nothing may overwrite it
        }
        // SN is handled via path names; segment->ref_path is set there
    }

    int absent;
    khint_t k = idxmap_put(h, s->id, &absent);
    if (absent) {
        kh_val(h, k) = (uint32_t)g->n_seg;
    } else {
        ak_log(AK_LOG_WARN, "gfa", "duplicate segment id %llu", (unsigned long long)s->id);
    }

    g->n_seg++;
    return AK_OK;
}

// parses one L line: optional validation, then records the edge
static int handle_L(gfa_t *g, char *line, idxmap_t *h, int flags) {
    uint64_t id1, id2;
    char st1, st2;
    size_t overlap = 0;
    if (sscanf(line, "L\t%" SCNu64 "\t%c\t%" SCNu64 "\t%c\t%zuM", &id1, &st1, &id2, &st2, &overlap) < 4) {
        ak_log(AK_LOG_WARN, "gfa", "malformed L line");
        return AK_EFORMAT;
    }

    khint_t k1 = idxmap_get(h, id1);
    khint_t k2 = idxmap_get(h, id2);
    int have1 = (k1 < kh_end(h)), have2 = (k2 < kh_end(h));

    if (flags & GFA_VALIDATE) {
        if (!have1) {
            ak_log(AK_LOG_WARN, "gfa", "L references unknown segment %llu", (unsigned long long)id1);
        }
        if (!have2) {
            ak_log(AK_LOG_WARN, "gfa", "L references unknown segment %llu", (unsigned long long)id2);
        }
        if (have1 && have2 && overlap > 0) {
            gfa_seg_t *a = &g->seg[kh_val(h, k1)];
            gfa_seg_t *b = &g->seg[kh_val(h, k2)];
            if (a->seq && b->seq && overlap < a->len && overlap < b->len && strncmp(a->seq, b->seq + (b->len - overlap), overlap) != 0) {
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

    return AK_OK;
}

// parses one P line. Only called under GFA_PATH_NAMES
static int handle_P(gfa_t *g, char *line, idxmap_t *h, int flags) {
    char *save;
    strtok_r(line, "\t", &save);                    // 'P'
    char *tok = strtok_r(NULL, "\t", &save);        // path name
    if (!tok) {
        ak_log(AK_LOG_WARN, "gfa", "P line without name");
        return AK_EFORMAT;
    }

    if (reserve_path(g) != AK_OK) return AK_ENOMEM;
    int32_t pi = g->n_path;
    char *name = strdup(tok);
    if (!name) return AK_ENOMEM;
    g->path[pi] = name;
    g->path_len[pi] = 0;
    g->path_off[pi] = (int32_t)g->n_path_seg;       // start of this path's slice
    g->n_path++;

    tok = strtok_r(NULL, "\t", &save);              // comma list of steps
    if (!tok || tok[0] == '\0') {
        g->path_off[pi + 1] = (int32_t)g->n_path_seg;   // empty slice
        ak_log(AK_LOG_WARN, "gfa", "path %s has no segments", name);
        return AK_OK;
    }

    if (!(flags & GFA_PATHS)) {
        // steps = commas + 1; nothing else about them is needed
        uint64_t n = 1;
        for (const char *p = strchr(tok, ','); p; p = strchr(p + 1, ',')) n++;
        g->n_path_seg += n;
        g->path_off[pi + 1] = (int32_t)g->n_path_seg;
        return AK_OK;
    }

    int32_t ref_pos = 0;
    const char *p = tok;
    for (;;) {
        // "<id>[+|-]" up to the next ',' or the end of the list
        char *end;
        uint64_t sid = strtoull(p, &end, 10);
        char ori = '+';
        if (end != p && (*end == '+' || *end == '-')) ori = *end++;
        while (*end == ' ' || *end == '\r') end++;   // stray trailing whitespace

        if (end == p || (*end != ',' && *end != '\0')) {
            ak_log(AK_LOG_WARN, "gfa", "malformed step '%.*s' in path %s", (int)strcspn(p, ","), p, name);
            const char *next = strchr(p, ',');
            if (!next) break;
            p = next + 1;
            continue;
        }

        khint_t k = idxmap_get(h, sid);
        uint32_t si = (k < kh_end(h)) ? kh_val(h, k) : GFA_NIL;
        if (si == GFA_NIL) {
            ak_log(AK_LOG_WARN, "gfa", "segment %llu in path %s not found", (unsigned long long)sid, name);
        }

        if (reserve_pathseg(g) != AK_OK) return AK_ENOMEM;
        g->path_seg[g->n_path_seg] = si;
        g->path_ori[g->n_path_seg] = ori;
        g->n_path_seg++;
        if (si != GFA_NIL) {
            gfa_seg_t *cur = &g->seg[si];
            cur->ref_path = pi;
            cur->start = ref_pos;
            ref_pos += (int32_t)cur->len;
            g->path_len[pi] += cur->len;
        }

        if (*end != ',') break;
        p = end + 1;
    }

    g->path_off[pi + 1] = (int32_t)g->n_path_seg;   // end of this path's slice
    return AK_OK;
}

// CSR out-adjacency

static int build_arcs(gfa_t *g) {
    if (g->n_seg <= 0 || g->n_link <= 0) return AK_OK;

    g->arc_off = (int32_t *)calloc((size_t)g->n_seg + 1, sizeof(int32_t));
    g->arc     = (uint32_t *)malloc((size_t)g->n_link * sizeof(uint32_t));
    if (!g->arc_off || !g->arc) return AK_ENOMEM;

    for (int32_t k = 0; k < g->n_link; k++) g->arc_off[g->link[k].v + 1]++;
    for (int32_t i = 0; i < g->n_seg; i++) g->arc_off[i + 1] += g->arc_off[i];

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

static int build_degrees(gfa_t *g) {
    if (g->n_seg <= 0) return AK_OK;

    g->in_degree  = (int32_t *)calloc((size_t)g->n_seg, sizeof(int32_t));
    g->out_degree = (int32_t *)calloc((size_t)g->n_seg, sizeof(int32_t));
    if (!g->in_degree || !g->out_degree) return AK_ENOMEM;

    for (int32_t k = 0; k < g->n_link; k++) {
        g->out_degree[g->link[k].v]++;
        g->in_degree[g->link[k].w]++;
    }
    return AK_OK;
}

// public API

// read an (r)GFA into a graph; see akhal/gfa.h
gfa_t *gfa_read(const char *fn, int flags) {
    ak_file *f = ak_open(fn);
    if (!f) return NULL;

    gfa_t *g = (gfa_t *)calloc(1, sizeof(gfa_t));
    if (!g) {
        ak_close(f);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory");
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
    // the overlap check compares the bases either side of a join, so asking to validate implies asking for the sequences
    if (flags & GFA_VALIDATE) flags |= GFA_SEQ;
    // the adjacency is an index over the edges, so it cannot be built without them
    if (flags & GFA_ARCS) flags |= GFA_LINKS;
    // degrees are counted over the edges, so the same holds for them
    if (flags & GFA_DEGREES) flags |= GFA_LINKS;
    // the step arrays hang off the per-path offsets, so resolving the steps implies recording the paths
    if (flags & GFA_PATHS) flags |= GFA_PATH_NAMES;
    // everything but a bare path listing resolves ids against seg[]
    if (flags & (GFA_LINKS | GFA_PATHS | GFA_SEQ | GFA_VALIDATE | GFA_ARCS | GFA_DEGREES)) flags |= GFA_SEGS;
    g->flags = flags;

    kstring_t ks = KS_INIT;
    int rc = AK_OK;
    long len;

    while ((len = ak_getline(f, &ks)) >= 0) {
        if (len == 0) continue;
        switch (ks.s[0]) {
            // no line type is tokenized unless a flag asked for what it
            // carries; GFA_SEGS is implied by everything but GFA_PATH_NAMES
            case 'S': rc = (flags & GFA_SEGS) ? handle_S(g, ks.s, h) : AK_OK; break;
            case 'L': rc = (flags & (GFA_LINKS | GFA_VALIDATE)) ? handle_L(g, ks.s, h, flags) : AK_OK; break;
            case 'P': rc = (flags & GFA_PATH_NAMES) ? handle_P(g, ks.s, h, flags) : AK_OK; break;
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

    if (flags & GFA_ARCS) {
        if (build_arcs(g) != AK_OK) {
            gfa_destroy(g);
            ak_log(AK_LOG_ERROR, "gfa", "out of memory building adjacency");
            return NULL;
        }
    }

    if (flags & GFA_DEGREES) {
        if (build_degrees(g) != AK_OK) {
            gfa_destroy(g);
            ak_log(AK_LOG_ERROR, "gfa", "out of memory building degrees");
            return NULL;
        }
    }

    // a file that ranks itself is authoritative; only fill the gap when it did
    // not, and only when we actually read the paths to derive a backbone from
    if (!g->has_sr && (flags & GFA_PATHS)) {
        gfa_rank_paths(g);
    }

    // nothing appends to the graph after this point except gfa_add_path(),
    // which grows the path block again on its own
    shrink((void **)&g->seg,  (size_t)g->n_seg,  sizeof(*g->seg));
    g->m_seg = g->n_seg;
    shrink((void **)&g->link, (size_t)g->n_link, sizeof(*g->link));
    g->m_link = g->n_link;
    shrink((void **)&g->path,      (size_t)g->n_path, sizeof(*g->path));
    shrink((void **)&g->path_len,  (size_t)g->n_path, sizeof(*g->path_len));
    shrink((void **)&g->path_off,  (size_t)g->n_path + 1, sizeof(*g->path_off));
    g->m_path = g->n_path;
    shrink((void **)&g->path_seg, (size_t)g->n_path_seg, sizeof(*g->path_seg));
    shrink((void **)&g->path_ori, (size_t)g->n_path_seg, sizeof(*g->path_ori));
    g->m_path_seg = (int32_t)g->n_path_seg;

    return g;
}

// emit a graph as GFA; see akhal/gfa.h
// shared by gfa_write() and gfa_write_rgfa(); `tags` adds SN and SO
static int write_graph(const gfa_t *g, FILE *out, int tags) {
    // Without GFA_SEQ every S line would come out as "*", which is a valid
    // GFA but a silently different graph. Fail instead.
    if (g->n_seg > 0 && !(g->flags & GFA_SEQ)) {
        ak_log(AK_LOG_ERROR, "gfa", "graph was read without GFA_SEQ; its segments carry no sequence to write");
        return AK_EINVAL;
    }
    // Likewise a graph that only kept path names would write every P line
    // with an empty step list.
    if (g->n_path > 0 && !(g->flags & GFA_PATHS)) {
        ak_log(AK_LOG_ERROR, "gfa", "graph was read without GFA_PATHS; its paths carry no steps to write");
        return AK_EINVAL;
    }

    fprintf(out, "H\tVN:Z:1.0\n");

    for (int32_t i = 0; i < g->n_seg; i++) {
        const gfa_seg_t *s = &g->seg[i];
        fprintf(out, "S\t%llu\t%s", (unsigned long long)s->id, s->seq ? s->seq : "*");
        const char *sn = gfa_seg_ref(g, s);
        if (tags && sn) {
            fprintf(out, "\tSN:Z:%s", sn);
        }
        if (tags && s->start >= 0) {
            fprintf(out, "\tSO:i:%d", s->start);
        }
        if (s->rank >= 0) {
            fprintf(out, "\tSR:i:%d", s->rank);
        }
        fputc('\n', out);
    }

    for (int32_t i = 0; i < g->n_link; i++) {
        const gfa_link_t *e = &g->link[i];
        fprintf(out, "L\t%llu\t%c\t%llu\t%c\t%uM\n", (unsigned long long)g->seg[e->v].id, e->from_orient, (unsigned long long)g->seg[e->w].id, e->to_orient, e->overlap);
    }

    for (int32_t k = 0; k < g->n_path; k++) {
        fprintf(out, "P\t%s\t", g->path[k]);
        const uint32_t *segs;
        int ns = gfa_path_segs(g, k, &segs);
        const char *ori = g->path_ori + g->path_off[k];
        int written = 0;
        for (int i = 0; i < ns; i++) {
            if (segs[i] == GFA_NIL) continue;
            fprintf(out, "%s%llu%c", written ? "," : "", (unsigned long long)g->seg[segs[i]].id, ori[i]);
            written = 1;
        }
        fprintf(out, "\t*\n");
    }

    if (ferror(out)) {
        ak_log(AK_LOG_ERROR, "gfa", "write failed");
        return AK_EIO;
    }
    return AK_OK;
}

// give a segment its sequence; see akhal/gfa.h
int gfa_seg_set_seq(gfa_t *g, gfa_seg_t *s, const char *seq, size_t len) {
    if (!g || !s) return AK_EINVAL;
    if (!seq || len == 0) {
        s->seq = NULL;
        s->len = 0;
        return AK_OK;
    }
    const char *p = ak_arena_put(&g->strs, seq, len);
    if (!p) return AK_ENOMEM;
    s->seq = p;
    s->len = (uint32_t)len;
    return AK_OK;
}

// emit GFA; see akhal/gfa.h
int gfa_write(const gfa_t *g, FILE *out) {
    return write_graph(g, out, 0);
}

// emit rGFA, stable-sequence tags and all; see akhal/gfa.h
int gfa_write_rgfa(const gfa_t *g, FILE *out) {
    return write_graph(g, out, 1);
}

// free a graph and everything it owns; see akhal/gfa.h
void gfa_destroy(gfa_t *g) {
    if (!g) return;
    ak_arena_destroy(&g->strs);   // every segment sequence, released in one go
    free(g->seg);
    free(g->in_degree);
    free(g->out_degree);
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
    if (g->idx) {
        idxmap_destroy((idxmap_t *)g->idx);
    }
    free(g);
}

// segment index for an id, or -1 if absent
int32_t gfa_idx(const gfa_t *g, uint64_t id) {
    idxmap_t *h = (idxmap_t *)g->idx;
    khint_t k = idxmap_get(h, id);
    return (k < kh_end(h)) ? (int32_t)kh_val(h, k) : -1;
}

// segment for an id, or NULL if absent
gfa_seg_t *gfa_get(const gfa_t *g, uint64_t id) {
    int32_t i = gfa_idx(g, id);
    return (i < 0) ? NULL : &g->seg[i];
}

// out-edges of segment v; see akhal/gfa.h
int gfa_arcs(const gfa_t *g, int32_t v, const uint32_t **arcs) {
    if (!g->arc_off || v < 0 || v >= g->n_seg) {
        *arcs = NULL;
        return 0;
    }
    int32_t beg = g->arc_off[v], end = g->arc_off[v + 1];
    *arcs = g->arc + beg;
    return (int)(end - beg);
}

// whether a link v -> w exists; see akhal/gfa.h
int gfa_has_arc(const gfa_t *g, int32_t v, int32_t w) {
    const uint32_t *a;
    int n = gfa_arcs(g, v, &a);
    for (int k = 0; k < n; k++)
        if ((int32_t)g->link[a[k]].w == w) return 1;
    return 0;
}

// ordered segments of path k; see akhal/gfa.h
int gfa_path_segs(const gfa_t *g, int32_t k, const uint32_t **segs) {
    // path_seg is NULL when only GFA_PATH_NAMES was set
    if (!g->path_off || !g->path_seg || k < 0 || k >= g->n_path) {
        *segs = NULL;
        return 0;
    }
    int32_t beg = g->path_off[k], end = g->path_off[k + 1];
    *segs = g->path_seg + beg;
    return (int)(end - beg);
}

// summary statistics, without a graph

// Welford's online mean and variance: one pass, no array
// ak_variance() divides by n, so this keeps the population form too
typedef struct {
    int64_t  n;
    double   mean, m2;
    uint64_t mn, mx;
} sdist_t;

static void sdist_add(sdist_t *d, uint64_t x) {
    if (!d->n || x < d->mn) d->mn = x;
    if (!d->n || x > d->mx) d->mx = x;
    d->n++;
    double delta = (double)x - d->mean;
    d->mean += delta / (double)d->n;
    d->m2   += delta * ((double)x - d->mean);
}

static double sdist_sd(const sdist_t *d) {
    return d->n ? sqrt(d->m2 / (double)d->n) : 0.0;
}

// What the summary keeps per segment id: its two degrees and three bits
// saying where the id was seen. Indexed by id - lo, so a chunk numbered from
// forty million costs no more than one numbered from one
#define SEEN_SEG  0x1
#define SEEN_PATH 0x2
#define SEEN_LINK 0x4

typedef struct {
    uint32_t *in, *out;
    uint8_t  *seen;
    uint64_t  lo, n;     // ids lo .. lo + n - 1 are addressable
} sidx_t;

static void sidx_free(sidx_t *x) {
    free(x->in);
    free(x->out);
    free(x->seen);
    x->in = NULL;
    x->out = NULL;
    x->seen = NULL;
    x->lo = x->n = 0;
}

// Make room for `id`, growing in whichever direction it lies. Returns 0 when
// the range asked for is too sparse to be worth an array - `budget` is the
// caller's idea of how many entries the file has earned so far
static int sidx_fit(sidx_t *x, uint64_t id, int64_t budget) {
    if (x->n && id >= x->lo && id < x->lo + x->n) return 1;

    uint64_t lo = x->n ? (id < x->lo ? id : x->lo) : id;
    uint64_t hi = x->n ? (id >= x->lo + x->n ? id : x->lo + x->n - 1) : id;
    uint64_t need = hi - lo + 1;

    uint64_t cap = (uint64_t)budget * 4 + (1u << 20);
    if (need > cap) return 0;

    // grow geometrically so a file numbered 1..N reallocates log N times
    uint64_t n = x->n ? x->n : (1u << 16);
    while (n < need) {
        uint64_t next = n << 1;
        if (next < n) return 0;
        n = next;
    }
    if (n > cap) n = need;
    // extend downwards only as far as asked; the common case grows upwards
    if (lo < x->lo || !x->n) {
        // keep the new low end, and let the array run up from there
        if (x->n && lo + n < x->lo + x->n) n = x->lo + x->n - lo;
    } else {
        lo = x->lo;
    }
    if (n > cap) return 0;

    uint32_t *ni = (uint32_t *)calloc((size_t)n, sizeof(*ni));
    uint32_t *no = (uint32_t *)calloc((size_t)n, sizeof(*no));
    uint8_t  *ns = (uint8_t *)calloc((size_t)n, 1);
    if (!ni || !no || !ns) {
        free(ni);
        free(no);
        free(ns);
        return -1;
    }
    if (x->n) {
        size_t at = (size_t)(x->lo - lo);
        memcpy(ni + at, x->in,   (size_t)x->n * sizeof(*ni));
        memcpy(no + at, x->out,  (size_t)x->n * sizeof(*no));
        memcpy(ns + at, x->seen, (size_t)x->n);
        free(x->in);
        free(x->out);
        free(x->seen);
    }
    x->in = ni;
    x->out = no;
    x->seen = ns;
    x->lo = lo;
    x->n = n;
    return 1;
}

// the last field of an S line that looks like SR:i:<n>, or -1 when there is none
static int32_t sr_tag(const char *line) {
    for (const char *p = line; (p = strstr(p, "SR:i:")) != NULL; p += 5) {
        if (p != line && p[-1] != '\t') continue;   // a tag starts a field
        return (int32_t)strtol(p + 5, NULL, 10);
    }
    return -1;
}

// fall back to the graph when the ids defeat the index; the answer is the same
static int stats_via_graph(const char *fn, gfa_stat_t *st) {
    gfa_t *g = gfa_read(fn, GFA_LINKS | GFA_PATHS | GFA_DEGREES);
    if (!g) return AK_EOPEN;

    sdist_t sl = {0}, ov = {0};
    for (int32_t i = 0; i < g->n_seg; i++) sdist_add(&sl, g->seg[i].len);
    for (int32_t i = 0; i < g->n_link; i++) sdist_add(&ov, g->link[i].overlap);

    st->n_seg  = g->n_seg;
    st->n_link = g->n_link;
    st->n_path = g->n_path;
    st->has_sr = g->has_sr;
    st->seg_mean = sl.mean;
    st->seg_sd   = sdist_sd(&sl);
    st->seg_min  = g->n_seg ? sl.mn : 0;
    st->seg_max  = g->n_seg ? sl.mx : 0;
    st->ov_mean  = ov.mean;
    st->ov_sd    = sdist_sd(&ov);

    st->min_in = st->max_in = st->min_out = st->max_out = -1;
    st->n_rank0 = 0;
    for (int32_t i = 0; i < g->n_seg; i++) {
        const gfa_seg_t *s = &g->seg[i];
        if (s->rank == 0) st->n_rank0++;
        if (g->in_degree[i]) {
            if (st->min_in < 0 || g->in_degree[i] < st->min_in) st->min_in = g->in_degree[i];
            if (g->in_degree[i] > st->max_in) st->max_in = g->in_degree[i];
        }
        if (g->out_degree[i]) {
            if (st->min_out < 0 || g->out_degree[i] < st->min_out) st->min_out = g->out_degree[i];
            if (g->out_degree[i] > st->max_out) st->max_out = g->out_degree[i];
        }
    }
    st->n_undefined = 0;   // gfa_read() drops those lines rather than counting them
    gfa_destroy(g);
    return AK_OK;
}

// summarize a file in one pass; see akhal/gfa.h
int gfa_read_stats(const char *fn, gfa_stat_t *st) {
    ak_file *f = ak_open(fn);
    if (!f) return AK_EOPEN;

    sdist_t sl = {0}, ov = {0};
    sidx_t x = {0};
    int64_t n_seg = 0, n_link = 0, n_path = 0, sr0 = 0;
    int has_sr = 0, sparse = 0, oom = 0;

    kstring_t ks = KS_INIT;
    long len;
    while ((len = ak_getline(f, &ks)) >= 0 && !sparse && !oom) {
        if (len == 0) continue;
        char *p = ks.s;

        if (p[0] == 'S' && p[1] == '\t') {
            p += 2;
            uint64_t id = strtoull(p, &p, 10);
            if (*p == '\t') p++;
            const char *seq = p;
            while (*p && *p != '\t') p++;
            n_seg++;
            sdist_add(&sl, (uint64_t)(p - seq));

            int32_t r = sr_tag(ks.s);
            if (r >= 0) {
                has_sr = 1;
                if (r == 0) sr0++;
            }
            int ok = sidx_fit(&x, id, n_seg);
            if (ok < 0) { oom = 1; break; }
            if (!ok)    { sparse = 1; break; }
            x.seen[id - x.lo] |= SEEN_SEG;

        } else if (p[0] == 'L' && p[1] == '\t') {
            p += 2;
            char *q;
            uint64_t a = strtoull(p, &q, 10);
            if (q == p || *q != '\t') continue;          // malformed; gfa_read warns and skips
            p = q + 1;
            if (!*p || p[1] != '\t') continue;
            p += 2;
            uint64_t b = strtoull(p, &q, 10);
            // an L line needs both ids and both orientations; the overlap may
            // be absent or a "*", which counts as 0 - the same four fields
            // gfa_read()'s sscanf insists on
            if (q == p || *q != '\t' || !q[1]) continue;
            p = q + 1;                                   // the to-orientation
            uint64_t o = 0;
            if (p[1] == '\t') {
                p += 2;
                o = strtoull(p, &q, 10);
                if (q == p) o = 0;
            }
            n_link++;
            sdist_add(&ov, o);

            int ok = sidx_fit(&x, a, n_seg);
            if (ok > 0) ok = sidx_fit(&x, b, n_seg);
            if (ok < 0) { oom = 1; break; }
            if (!ok)    { sparse = 1; break; }
            if (x.out[a - x.lo] != UINT32_MAX) x.out[a - x.lo]++;
            if (x.in[b - x.lo]  != UINT32_MAX) x.in[b - x.lo]++;
            x.seen[a - x.lo] |= SEEN_LINK;
            x.seen[b - x.lo] |= SEEN_LINK;

        } else if (p[0] == 'P' && p[1] == '\t') {
            p += 2;
            while (*p && *p != '\t') p++;                // past the name
            if (!*p) continue;
            p++;
            n_path++;
            while (*p && *p != '\t') {
                char *q;
                uint64_t id = strtoull(p, &q, 10);
                if (q == p) break;
                int ok = sidx_fit(&x, id, n_seg);
                if (ok < 0) { oom = 1; break; }
                if (!ok)    { sparse = 1; break; }
                x.seen[id - x.lo] |= SEEN_PATH;
                p = q;
                while (*p && *p != ',' && *p != '\t') p++;
                if (*p == ',') p++;
            }
        }
    }
    ks_free(&ks);
    ak_close(f);

    if (oom) {
        sidx_free(&x);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory summarizing %s", fn);
        return AK_ENOMEM;
    }
    if (sparse) {
        sidx_free(&x);
        ak_log(AK_LOG_INFO, "gfa", "segment ids are too scattered to index; summarizing through the graph instead");
        return stats_via_graph(fn, st);
    }

    memset(st, 0, sizeof(*st));
    st->n_seg  = n_seg;
    st->n_link = n_link;
    st->n_path = n_path;
    st->has_sr = has_sr;
    st->seg_mean = sl.mean;
    st->seg_sd   = sdist_sd(&sl);
    st->seg_min  = n_seg ? sl.mn : 0;
    st->seg_max  = n_seg ? sl.mx : 0;
    st->ov_mean  = ov.mean;
    st->ov_sd    = sdist_sd(&ov);
    st->min_in = st->max_in = st->min_out = st->max_out = -1;

    for (uint64_t i = 0; i < x.n; i++) {
        uint8_t s = x.seen[i];
        if (!s) continue;
        if (!(s & SEEN_SEG)) {
            if (s & (SEEN_LINK | SEEN_PATH)) st->n_undefined++;
            continue;   // only a real segment carries a degree or a rank
        }
        uint32_t di = x.in[i], dou = x.out[i];
        if (di) {
            if (st->min_in < 0 || (int32_t)di < st->min_in) st->min_in = (int32_t)di;
            if ((int32_t)di > st->max_in) st->max_in = (int32_t)di;
        }
        if (dou) {
            if (st->min_out < 0 || (int32_t)dou < st->min_out) st->min_out = (int32_t)dou;
            if ((int32_t)dou > st->max_out) st->max_out = (int32_t)dou;
        }
        if (!has_sr && (s & SEEN_PATH)) st->n_rank0++;
    }
    if (has_sr) st->n_rank0 = sr0;

    sidx_free(&x);
    return AK_OK;
}

// ranks

// rank against a caller-supplied backbone; see akhal/gfa.h
int64_t gfa_rank_mark(gfa_t *g, const uint8_t *on) {
    if (!on) return AK_EINVAL;

    int64_t n0 = 0;
    for (int32_t i = 0; i < g->n_seg; i++) {
        if (on[i]) {
            g->seg[i].rank = 0;
            n0++;
        } else {
            g->seg[i].rank = 1;
        }
    }
    return n0;
}

// rank against the graph's own paths; see akhal/gfa.h
int64_t gfa_rank_paths(gfa_t *g) {
    if (!(g->flags & GFA_PATHS)) {
        ak_log(AK_LOG_ERROR, "gfa", "ranking requires the graph to be read with GFA_PATHS");
        return AK_EINVAL;
    }
    if (g->n_seg <= 0) return 0;

    uint8_t *on = (uint8_t *)calloc((size_t)g->n_seg, 1);
    if (!on) {
        ak_log(AK_LOG_ERROR, "gfa", "out of memory");
        return AK_ENOMEM;
    }

    for (int32_t k = 0; k < g->n_path; k++) {
        const uint32_t *segs;
        int ns = gfa_path_segs(g, k, &segs);
        for (int i = 0; i < ns; i++) {
            if (segs[i] != GFA_NIL) {
                on[segs[i]] = 1;
            }
        }
    }

    int64_t n0 = gfa_rank_mark(g, on);
    free(on);

    if (g->n_path == 0) {
        ak_log(AK_LOG_DEBUG, "gfa", "no paths to rank against; all %d segment(s) left at rank 1", g->n_seg);
    } else {
        ak_log(AK_LOG_DEBUG, "gfa", "%lld segment(s) at rank 0 over %d path(s), %lld at rank 1", (long long)n0, g->n_path, (long long)((int64_t)g->n_seg - n0));
    }
    return n0;
}

// rewriting the path block

// drop every path; see akhal/gfa.h
void gfa_clear_paths(gfa_t *g) {
    // every segment refers to a path by index, so reset those before the
    // names they stand for go away
    for (int32_t i = 0; i < g->n_seg; i++) {
        g->seg[i].ref_path = -1;
    }
    if (g->path) {
        for (int32_t i = 0; i < g->n_path; i++) free(g->path[i]);
    }
    g->n_path = 0;
    g->n_path_seg = 0;
    if (g->path_off) {
        g->path_off[0] = 0;
    }
}

// append one laid-out path; see akhal/gfa.h
int gfa_add_path(gfa_t *g, const char *name, const uint32_t *segs, const char *ori, int64_t n) {
    if (!(g->flags & GFA_PATHS)) {
        ak_log(AK_LOG_ERROR, "gfa", "adding a path requires the graph to be read with GFA_PATHS");
        return AK_EINVAL;
    }
    if (!name || (n > 0 && !segs)) return AK_EINVAL;

    if (reserve_path(g) != AK_OK) return AK_ENOMEM;
    int32_t pi = g->n_path;
    char *owned = strdup(name);
    if (!owned) return AK_ENOMEM;

    g->path[pi] = owned;
    g->path_len[pi] = 0;
    g->path_off[pi] = (int32_t)g->n_path_seg;
    g->n_path++;

    int32_t ref_pos = 0;
    for (int64_t i = 0; i < n; i++) {
        if (segs[i] == GFA_NIL) continue;
        if (reserve_pathseg(g) != AK_OK) return AK_ENOMEM;
        g->path_seg[g->n_path_seg] = segs[i];
        g->path_ori[g->n_path_seg] = ori ? ori[i] : '+';

        gfa_seg_t *cur = &g->seg[segs[i]];
        cur->ref_path = pi;
        cur->start = ref_pos;
        ref_pos += (int32_t)cur->len;
        g->path_len[pi] += cur->len;
        g->n_path_seg++;
    }

    g->path_off[pi + 1] = (int32_t)g->n_path_seg;
    return AK_OK;
}

// topological sort

// ordering by sequence content rather than id keeps the result independent of
// the input's node numbering; a NULL/empty sequence sorts first
static int seq_lt(const gfa_t *g, int32_t a, int32_t b) {
    const char *sa = g->seg[a].seq ? g->seg[a].seq : "";
    const char *sb = g->seg[b].seq ? g->seg[b].seq : "";
    return strcmp(sa, sb) < 0;
}

// sift the last heap element up to restore the min-heap order
static void heap_push(const gfa_t *g, int32_t *heap, int *hn, int32_t v) {
    int i = (*hn)++;
    heap[i] = v;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (!seq_lt(g, heap[i], heap[p])) break;
        int32_t t = heap[i];
        heap[i] = heap[p];
        heap[p] = t;
        i = p;
    }
}

// pop and return the alphabetically-smallest node from the heap
static int32_t heap_pop(const gfa_t *g, int32_t *heap, int *hn) {
    int32_t top = heap[0];
    int n = --(*hn);
    heap[0] = heap[n];
    int i = 0;
    for (;;) {
        int l = 2 * i + 1, r = 2 * i + 2, m = i;
        if (l < n && seq_lt(g, heap[l], heap[m])) {
            m = l;
        }
        if (r < n && seq_lt(g, heap[r], heap[m])) {
            m = r;
        }
        if (m == i) break;
        int32_t t = heap[i];
        heap[i] = heap[m];
        heap[m] = t;
        i = m;
    }
    return top;
}

// topological order with alphabetical id tie-break; see akhal/gfa.h
int gfa_toposort(const gfa_t *g, int32_t *order) {
    if (!g->arc_off) {
        ak_log(AK_LOG_ERROR, "gfa", "toposort needs the CSR adjacency; read with GFA_ARCS");
        return AK_EINVAL;
    }
    int32_t n = g->n_seg;
    if (n == 0) return 0;
    if (!gfa_has_degrees(g)) {
        ak_log(AK_LOG_ERROR, "gfa", "toposort needs the in-degrees; read with GFA_DEGREES");
        return AK_EINVAL;
    }

    int32_t *indeg = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t *heap  = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!indeg || !heap) {
        free(indeg);
        free(heap);
        return AK_ENOMEM;
    }

    for (int32_t i = 0; i < n; i++) indeg[i] = g->in_degree[i];

    int hn = 0;
    for (int32_t i = 0; i < n; i++) {
        if (indeg[i] == 0) {
            heap_push(g, heap, &hn, i);
        }
    }

    int32_t placed = 0;
    while (hn > 0) {
        int32_t u = heap_pop(g, heap, &hn);
        order[placed++] = u;
        const uint32_t *arcs;
        int na = gfa_arcs(g, u, &arcs);
        for (int k = 0; k < na; k++) {
            uint32_t w = g->link[arcs[k]].w;
            if (--indeg[w] == 0) {
                heap_push(g, heap, &hn, (int32_t)w);
            }
        }
    }

    int32_t acyclic = placed;
    if (placed < n) {
        // remaining nodes are inside cycles; append them alphabetically
        for (int32_t i = 0; i < n; i++) {
            if (indeg[i] > 0) {
                heap_push(g, heap, &hn, i);
            }
        }
        while (hn > 0) order[placed++] = heap_pop(g, heap, &hn);
    }

    free(indeg);
    free(heap);
    return acyclic;
}
