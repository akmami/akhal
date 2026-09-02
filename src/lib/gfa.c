#include "akhal/gfa.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/error.h"

#include "khashl.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// id -> array index
KHASHL_MAP_INIT(KH_LOCAL, idxmap_t, idxmap, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

// growth helpers

static int reserve_seg(gfa_t *g) {
    if (g->n_seg < g->m_seg) return AK_OK;
    int32_t m = g->m_seg ? g->m_seg << 1 : 1024;
    gfa_seg_t *p = (gfa_seg_t *)realloc(g->seg, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    g->seg = p;
    g->m_seg = m;
    return AK_OK;
}

static int reserve_link(gfa_t *g) {
    if (g->n_link < g->m_link) return AK_OK;
    int32_t m = g->m_link ? g->m_link << 1 : 1024;
    gfa_link_t *p = (gfa_link_t *)realloc(g->link, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    g->link = p;
    g->m_link = m;
    return AK_OK;
}

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
    s->rank = -1;   // -1 until an SR tag says otherwise

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
            g->has_sr = 1;   // the file ranks itself; nothing may overwrite it
        }
        // SN is handled via path names; segment->ref_name is set there.
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
    if (sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM",
               &id1, &st1, &id2, &st2, &overlap) < 4) {
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

// parses one P line and lays out reference coordinates; without GFA_PATHS only
// the occurrence count is maintained
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

        if (si == GFA_NIL && (flags & (GFA_PATHS | GFA_VALIDATE))) {
            ak_log(AK_LOG_WARN, "gfa", "segment %llu in path %s not found", (unsigned long long)sid, name ? name : "?");
        }

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

    if (flags & GFA_PATHS) {
        g->path_off[pi + 1] = (int32_t)g->n_path_seg;   // end of this path's slice
    }
    return AK_OK;
}

// CSR out-adjacency

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

    // a file that ranks itself is authoritative; only fill the gap when it did
    // not, and only when we actually read the paths to derive a backbone from
    if (!g->has_sr && (flags & GFA_PATHS)) {
        gfa_rank_paths(g);
    }

    return g;
}

// emit a graph as GFA; see akhal/gfa.h
int gfa_write(const gfa_t *g, FILE *out) {
    fprintf(out, "H\tVN:Z:1.0\n");

    for (int32_t i = 0; i < g->n_seg; i++) {
        const gfa_seg_t *s = &g->seg[i];
        fprintf(out, "S\t%llu\t%s", (unsigned long long)s->id, s->seq ? s->seq : "*");
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

// free a graph and everything it owns; see akhal/gfa.h
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
    if (!g->path_off || k < 0 || k >= g->n_path) {
        *segs = NULL;
        return 0;
    }
    int32_t beg = g->path_off[k], end = g->path_off[k + 1];
    *segs = g->path_seg + beg;
    return (int)(end - beg);
}

// fragmented paths

// name length with vg's subpath decoration stripped: "chr22[1000]" -> "chr22"
static size_t frag_base_len(const char *name) {
    const char *b = strchr(name, '[');
    if (b) return (size_t)(b - name);
    // a ':' only decorates the name when digits follow it; PanSN uses '#'
    for (const char *p = strchr(name, ':'); p; p = strchr(p + 1, ':'))
        if (p[1] >= '0' && p[1] <= '9') return (size_t)(p - name);
    return strlen(name);
}

// start offset encoded in a path name, or -1 when it carries none
static int64_t frag_base_off(const char *name) {
    const char *b = strchr(name, '[');
    if (!b) {
        for (const char *p = strchr(name, ':'); p; p = strchr(p + 1, ':'))
            if (p[1] >= '0' && p[1] <= '9') {
                b = p;
                break;
            }
    }
    return b ? (int64_t)strtoll(b + 1, NULL, 10) : -1;
}

// matches a base name against a group key; a bare PanSN contig name matches too
static int frag_in_group(const char *name, size_t base_l, const char *key) {
    size_t kl = strlen(key);
    if (base_l == kl && !strncmp(name, key, kl)) return 1;
    // PanSN "sample#hap#contig": the bare contig name selects it too
    for (size_t i = base_l; i > 0; i--) {
        if (name[i - 1] != '#') continue;
        return base_l - i == kl && !strncmp(name + i, key, kl);
    }
    return 0;
}

// one selected P-line fragment, sorted into its group before chaining
typedef struct {
    int32_t     pi;         // path index
    int32_t     ord;        // selection order, a stable tie-break
    int64_t     off;        // start offset from the name, or -1 when absent
    const char *base;       // borrowed: first byte of the base name
    size_t      base_l;     // base name length
    uint32_t    head, tail; // first / last resolved segment index, or GFA_NIL
    char        head_ori, tail_ori;
    int32_t     succ;       // fragment that follows this one, or -1
    int32_t     claimed;    // 1 once another fragment claims this as successor
} frag_t;

// group fragments by base name, then order each group along the reference
static int frag_cmp(const void *A, const void *B) {
    const frag_t *a = (const frag_t *)A, *b = (const frag_t *)B;
    size_t n = a->base_l < b->base_l ? a->base_l : b->base_l;
    int c = strncmp(a->base, b->base, n);
    if (c) return c;
    if (a->base_l != b->base_l) return a->base_l < b->base_l ? -1 : 1;
    if (a->off != b->off) {
        if (a->off < 0) return 1;    // an undecorated name sorts last
        if (b->off < 0) return -1;
        return a->off < b->off ? -1 : 1;
    }
    return a->ord < b->ord ? -1 : (a->ord > b->ord);
}

// fills in a fragment's end segments; 0 when the path resolved no segment
static int frag_ends(const gfa_t *g, frag_t *f) {
    const uint32_t *segs;
    int ns = gfa_path_segs(g, f->pi, &segs);
    const char *ori = g->path_ori + g->path_off[f->pi];
    int seen = 0;
    for (int i = 0; i < ns; i++) {
        if (segs[i] == GFA_NIL) continue;
        if (!seen) {
            f->head = segs[i];
            f->head_ori = ori[i];
            seen = 1;
        }
        f->tail = segs[i];
        f->tail_ori = ori[i];
    }
    return seen;
}

// whether a link joins v to w with these orientations on both ends
static int frag_joined(const gfa_t *g, uint32_t v, char vo, uint32_t w, char wo) {
    const uint32_t *arcs;
    int na = gfa_arcs(g, (int32_t)v, &arcs);
    for (int i = 0; i < na; i++) {
        const gfa_link_t *e = &g->link[arcs[i]];
        if (e->w == w && e->from_orient == vo && e->to_orient == wo) return 1;
    }
    return 0;
}

// a merged chain takes the shared base name, numbered when a base yields more
// than one; a lone fragment keeps the name it already had
static int frag_name_chains(const gfa_t *g, gfa_merge_t *m) {
    for (int32_t c = 0; c < m->n; ) {
        const char *base = g->path[m->frag[m->off[c]]];
        size_t bl = frag_base_len(base);

        // extent of this base's chains, and whether merging happened at all
        int32_t e = c, n_multi = 0, base_taken = 0;
        while (e < m->n) {
            const char *nm = g->path[m->frag[m->off[e]]];
            if (frag_base_len(nm) != bl || strncmp(nm, base, bl) != 0) break;
            if (m->off[e + 1] - m->off[e] > 1) {
                n_multi++;
            } else if (strlen(nm) == bl) {
                base_taken = 1;   // an undecorated path
            }
            e++;
        }

        int32_t seen = 0;
        for (int32_t k = c; k < e; k++) {
            const char *nm = g->path[m->frag[m->off[k]]];
            if (m->off[k + 1] - m->off[k] == 1) {
                m->name[k] = strdup(nm);
            } else if ((m->name[k] = (char *)malloc(bl + 16)) != NULL) {
                if (n_multi == 1 && !base_taken) {
                    snprintf(m->name[k], bl + 16, "%.*s", (int)bl, base);
                } else {
                    snprintf(m->name[k], bl + 16, "%.*s_%d", (int)bl, base, ++seen);
                }
            }
            if (!m->name[k]) return AK_ENOMEM;
        }
        c = e;
    }
    return AK_OK;
}

// group P-line fragments into chains; see akhal/gfa.h
gfa_merge_t *gfa_path_merge(const gfa_t *g, const char *key) {
    if (!(g->flags & GFA_PATHS) || !g->path_off) {
        ak_log(AK_LOG_ERROR, "gfa", "path merging requires the graph to be read with GFA_PATHS");
        return NULL;
    }
    if (!(g->flags & GFA_LINKS)) {
        ak_log(AK_LOG_ERROR, "gfa", "path merging requires the graph to be read with GFA_LINKS");
        return NULL;
    }
    if (g->n_path == 0) {
        ak_log(AK_LOG_ERROR, "gfa", "graph has no P lines");
        return NULL;
    }

    frag_t *f = (frag_t *)malloc((size_t)g->n_path * sizeof(frag_t));
    if (!f) {
        ak_log(AK_LOG_ERROR, "gfa", "out of memory");
        return NULL;
    }

    // select the fragments of interest and note the segments they end on
    int32_t nf = 0;
    for (int32_t k = 0; k < g->n_path; k++) {
        const char *nm = g->path[k];
        size_t bl = frag_base_len(nm);
        if (key && !frag_in_group(nm, bl, key)) continue;

        frag_t *c = &f[nf];
        c->pi      = k;
        c->ord     = nf;
        c->off     = frag_base_off(nm);
        c->base    = nm;
        c->base_l  = bl;
        c->succ    = -1;
        c->claimed = 0;
        if (!frag_ends(g, c)) {   // no resolved segment; nothing to chain on
            c->head = c->tail = GFA_NIL;
            c->head_ori = c->tail_ori = '+';
        }
        nf++;
    }
    if (nf == 0) {
        free(f);
        ak_log(AK_LOG_ERROR, "gfa", "no path named '%s' in the graph", key);
        return NULL;
    }
    qsort(f, (size_t)nf, sizeof(frag_t), frag_cmp);

    // chain a fragment to the first later one its last segment links into
    for (int32_t i = 0; i < nf; i++) {
        if (f[i].tail == GFA_NIL) continue;
        for (int32_t j = 0; j < nf; j++) {
            if (j == i || f[j].claimed || f[j].head == GFA_NIL) continue;
            if (f[j].base_l != f[i].base_l || strncmp(f[j].base, f[i].base, f[i].base_l) != 0) continue;
            if (!frag_joined(g, f[i].tail, f[i].tail_ori, f[j].head, f[j].head_ori)) continue;

            // refuse a join that would close a cycle
            int32_t t = j, guard = 0;
            while (t >= 0 && t != i && guard++ <= nf) t = f[t].succ;
            if (t == i) continue;

            f[i].succ = j;
            f[j].claimed = 1;
            break;
        }
    }

    int32_t nc = 0;
    for (int32_t i = 0; i < nf; i++) if (!f[i].claimed) nc++;

    gfa_merge_t *m = (gfa_merge_t *)calloc(1, sizeof(*m));
    if (m) {
        m->name = (char **)calloc((size_t)nc, sizeof(char *));
        m->frag = (int32_t *)malloc((size_t)nf * sizeof(int32_t));
        m->off  = (int32_t *)malloc(((size_t)nc + 1) * sizeof(int32_t));
    }
    if (!m || !m->name || !m->frag || !m->off) {
        free(f);
        gfa_merge_destroy(m);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory");
        return NULL;
    }

    // one chain per fragment nothing else claimed, walked through succ
    int32_t used = 0;
    m->off[0] = 0;
    for (int32_t i = 0; i < nf; i++) {
        if (f[i].claimed) continue;
        for (int32_t t = i, guard = 0; t >= 0 && guard <= nf; t = f[t].succ, guard++)
            m->frag[used++] = f[t].pi;
        m->off[++m->n] = used;
    }
    free(f);

    if (frag_name_chains(g, m) != AK_OK) {
        gfa_merge_destroy(m);
        ak_log(AK_LOG_ERROR, "gfa", "out of memory");
        return NULL;
    }

    ak_log(AK_LOG_INFO, "gfa", "%d path fragment(s) grouped into %d chain(s)", nf, m->n);
    return m;
}

// free a chain set; see akhal/gfa.h
void gfa_merge_destroy(gfa_merge_t *m) {
    if (!m) return;
    if (m->name) {
        for (int32_t i = 0; i < m->n; i++) free(m->name[i]);
        free(m->name);
    }
    free(m->frag);
    free(m->off);
    free(m);
}

// flat segment list of one chain; see akhal/gfa.h
int64_t gfa_merge_segs(const gfa_t *g, const gfa_merge_t *m, int32_t c, uint32_t **segs) {
    *segs = NULL;
    if (c < 0 || c >= m->n) return AK_EINVAL;

    int64_t n = 0;
    for (int32_t i = m->off[c]; i < m->off[c + 1]; i++) {
        const uint32_t *ss;
        int ns = gfa_path_segs(g, m->frag[i], &ss);
        for (int t = 0; t < ns; t++) if (ss[t] != GFA_NIL) n++;
    }
    if (n == 0) return 0;

    uint32_t *out = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    if (!out) return AK_ENOMEM;

    int64_t k = 0;
    for (int32_t i = m->off[c]; i < m->off[c + 1]; i++) {
        const uint32_t *ss;
        int ns = gfa_path_segs(g, m->frag[i], &ss);
        for (int t = 0; t < ns; t++) if (ss[t] != GFA_NIL) out[k++] = ss[t];
    }
    *segs = out;
    return n;
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
    // ref_name borrows from path[], so drop those references before freeing
    for (int32_t i = 0; i < g->n_seg; i++) {
        g->seg[i].ref_name = NULL;
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
        cur->ref_name = owned;
        cur->start = ref_pos;
        ref_pos += (int32_t)cur->len;
        cur->end = ref_pos;
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
    int32_t n = g->n_seg;
    if (n == 0) return 0;
    if (!(g->flags & GFA_LINKS)) {
        ak_log(AK_LOG_ERROR, "gfa", "toposort requires the graph to be read with GFA_LINKS");
        return AK_EINVAL;
    }

    int32_t *indeg = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    int32_t *heap  = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    if (!indeg || !heap) {
        free(indeg);
        free(heap);
        return AK_ENOMEM;
    }

    for (int32_t i = 0; i < n; i++) indeg[i] = g->seg[i].in_degree;

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
