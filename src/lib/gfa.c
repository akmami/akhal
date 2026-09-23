#include "akhal/gfa.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/error.h"

#include "khashl.h"

#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <stdio.h>

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

    if (g->flags & GFA_IDX) {
        int absent;
        khint_t k = idxmap_put(h, s->id, &absent);
        if (absent) {
            kh_val(h, k) = (uint32_t)g->n_seg;
        } else {
            ak_log(AK_LOG_WARN, "gfa", "duplicate segment id %llu", (unsigned long long)s->id);
        }
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
    if (g->n_seg <= 0) return AK_OK;

    g->arc_off = (int32_t *)calloc((size_t)g->n_seg + 1, sizeof(int32_t));
    if (!g->arc_off) return AK_ENOMEM;
    if (g->n_link <= 0) return AK_OK;

    g->arc     = (uint32_t *)malloc((size_t)g->n_link * sizeof(uint32_t));
    if (!g->arc) return AK_ENOMEM;

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

    // GFA_VALIDATE checks what is loaded: unknown ids always, overlaps only when the sequences are there to compare
    // (GFA_SEQ), so a caller can check a whole-genome file's references without holding its bases
    // the adjacency is an index over the edges, so it cannot be built without them
    if (flags & GFA_ARCS) flags |= GFA_LINKS;
    // degrees are counted over the edges, so the same holds for them
    if (flags & GFA_DEGREES) flags |= GFA_LINKS;
    // the step arrays hang off the per-path offsets, so resolving the steps implies recording the paths
    if (flags & GFA_PATHS) flags |= GFA_PATH_NAMES;
    // an L line resolves both its ends through the index, and so does a path step once it is
    // resolved rather than merely counted - those are what need it built while the file is read
    if (flags & (GFA_LINKS | GFA_VALIDATE | GFA_PATHS)) flags |= GFA_IDX;
    // everything but a bare path listing resolves ids against seg[]
    if (flags & (GFA_IDX | GFA_LINKS | GFA_PATHS | GFA_SEQ | GFA_VALIDATE | GFA_ARCS | GFA_DEGREES)) flags |= GFA_SEGS;
    g->flags = flags;

    // On a whole-genome graph the index rivals the segments for size, so it is built only when it was asked for - left out
    idxmap_t *h = NULL;
    if (flags & GFA_IDX) {
        h = idxmap_init();
        if (!h) {
            free(g);
            ak_close(f);
            ak_log(AK_LOG_ERROR, "gfa", "out of memory");
            return NULL;
        }
    }
    g->idx = h;

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
    // Without GFA_SEQ every S line would come out as "*", which is a valid GFA but a silently different graph. Fail instead
    if (g->n_seg > 0 && !(g->flags & GFA_SEQ)) {
        ak_log(AK_LOG_ERROR, "gfa", "graph was read without GFA_SEQ; its segments carry no sequence to write");
        return AK_EINVAL;
    }
    // Likewise a graph that only kept path names would write every P line with an empty step list
    if (g->n_path > 0 && !(g->flags & GFA_PATHS)) {
        ak_log(AK_LOG_ERROR, "gfa", "graph was read without GFA_PATHS; its paths carry no steps to write");
        return AK_EINVAL;
    }

    fprintf(out, "H\tVN:Z:1.0\n");

    for (int32_t i = 0; i < g->n_seg; i++) {
        const gfa_seg_t *s = &g->seg[i];
        fprintf(out, "S\t%llu\t%s", (unsigned long long)s->id, s->seq ? s->seq : "*");
        const char *sn = gfa_seg_ref(g, s);
        if (tags && sn)             fprintf(out, "\tSN:Z:%s", sn);
        if (tags && s->start >= 0)  fprintf(out, "\tSO:i:%d", s->start);
        if (s->rank >= 0)           fprintf(out, "\tSR:i:%d", s->rank);
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

// release the parts named by `what`; see akhal/gfa.h
void gfa_drop(gfa_t *g, int what) {
    if (!g) return;

    if (what & GFA_SEQ) {
        // the bases live in one arena, so they go in one call - but every
        // segment points into it, and a dangling seq is worse than no seq
        ak_arena_destroy(&g->strs);
        for (int32_t i = 0; i < g->n_seg; i++) g->seg[i].seq = NULL;
    }
    if (what & GFA_ARCS) {
        free(g->arc);
        free(g->arc_off);
        g->arc = NULL;
        g->arc_off = NULL;
    }
    if (what & GFA_DEGREES) {
        free(g->in_degree);
        free(g->out_degree);
        g->in_degree = NULL;
        g->out_degree = NULL;
    }
    if (what & GFA_LINKS) {
        free(g->link);
        g->link = NULL;
        g->n_link = g->m_link = 0;
        what |= GFA_ARCS;               // an index over edges that are gone
        free(g->arc);
        free(g->arc_off);
        g->arc = NULL;
        g->arc_off = NULL;
    }
    if (what & GFA_PATHS) {
        free(g->path_seg);
        free(g->path_ori);
        g->path_seg = NULL;
        g->path_ori = NULL;
        g->m_path_seg = 0;
    }
    if (what & GFA_PATH_NAMES) {
        if (g->path) {
            for (int32_t i = 0; i < g->n_path; i++) free(g->path[i]);
            free(g->path);
            g->path = NULL;
        }
        free(g->path_len);
        free(g->path_off);
        g->path_len = NULL;
        g->path_off = NULL;
        g->n_path = g->m_path = 0;
        g->n_path_seg = 0;
    }
    if (what & GFA_IDX) {
        idxmap_destroy((idxmap_t *)g->idx);
        g->idx = NULL;
    }

    // what the graph no longer carries, it was no longer read with
    g->flags &= ~what;
}

// segment index for an id, or -1 if absent
int32_t gfa_idx(const gfa_t *g, uint64_t id) {
    idxmap_t *h = (idxmap_t *)g->idx;
    if (!h) return -1;                  // dropped with gfa_drop(g, GFA_IDX)
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

// whether the oriented join v(ov) -> w(ow) is in the graph; see akhal/gfa.h
int gfa_has_link(const gfa_t *g, int32_t v, char ov, int32_t w, char ow) {
    const uint32_t *a;
    int n = gfa_arcs(g, v, &a);
    for (int k = 0; k < n; k++) {
        const gfa_link_t *e = &g->link[a[k]];
        if ((int32_t)e->w == w && e->from_orient == ov && e->to_orient == ow) return 1;
    }
    char cv = ov == '+' ? '-' : '+', cw = ow == '+' ? '-' : '+';
    n = gfa_arcs(g, w, &a);
    for (int k = 0; k < n; k++) {
        const gfa_link_t *e = &g->link[a[k]];
        if ((int32_t)e->w == v && e->from_orient == cw && e->to_orient == cv) return 1;
    }
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

// ranks

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

/**
 * @brief ordering by sequence content rather than id keeps the result independent of
 * the input's node numbering; a NULL/empty sequence sorts first
 */
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

// topological order with sequence content being tie-break; see akhal/gfa.h
int gfa_toposort(const gfa_t *g, int32_t *order) {
    if (!g->arc_off) {
        ak_log(AK_LOG_ERROR, "gfa", "toposort needs the CSR adjacency; read with GFA_ARCS");
        return AK_EINVAL;
    }
    int32_t n = g->n_seg;
    if (n == 0) return 0;

    // The sort consumes the in-degrees as it goes, so it works on a copy
    // either way; a graph read without GFA_DEGREES is counted here instead,
    // which spares the caller the graph's two arrays for the sake of this one
    int32_t *indeg = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    int32_t *heap  = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!indeg || !heap) {
        free(indeg);
        free(heap);
        return AK_ENOMEM;
    }

    if (gfa_has_degrees(g)) {
        for (int32_t i = 0; i < n; i++) indeg[i] = g->in_degree[i];
    } else {
        for (int32_t k = 0; k < g->n_link; k++) indeg[g->link[k].w]++;
    }

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
        // remaining nodes are inside cycles; append them based in sequence content (to make deterministic)
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
