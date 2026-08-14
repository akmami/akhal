#include "akhal/annot.h"
#include "akhal/vcf.h"
#include "akhal/error.h"

#include "khashl.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// id -> record index
KHASHL_MAP_INIT(KH_LOCAL, annmap_t, annmap, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

static const char ANNOT_MAGIC[8] = "AKANNOT1";

// longest alternative-allele walk (in nodes) matched against a variant
#define ANNOT_MAX_WALK 4096

// growth helpers

/**
 * Ensure room for one more record. 
 * @return AK_OK or AK_ENOMEM
 */ 
static int reserve_rec(annot_t *a) {
    if (a->n < a->m) return AK_OK;
    int64_t m = a->m ? a->m << 1 : 1024;
    annot_rec_t *p = (annot_rec_t *)realloc(a->rec, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    a->rec = p;
    a->m = m;
    return AK_OK;
}

/**
 * Append a string to the shared info buffer, NUL-terminated
 * @param a Store whose buffer grows
 * @param s String to copy (may be NULL / empty)
 * @param off Set to the byte offset of the copy
 * @param len Set to the string length (excluding the NUL)
 * @return AK_OK or AK_ENOMEM
 */
static int buf_put(annot_t *a, const char *s, uint64_t *off, uint32_t *len) {
    size_t l = s ? strlen(s) : 0;
    if (l == 0) { *off = 0; *len = 0; return AK_OK; }
    if (a->buf_l + l + 1 > a->buf_m) {
        uint64_t m = a->buf_m ? a->buf_m : 4096;
        while (m < a->buf_l + l + 1) m <<= 1;
        char *p = (char *)realloc(a->buf, m);
        if (!p) return AK_ENOMEM;
        a->buf = p;
        a->buf_m = m;
    }
    memcpy(a->buf + a->buf_l, s, l + 1);
    *off = a->buf_l;
    *len = (uint32_t)l;
    a->buf_l += l + 1;
    return AK_OK;
}

/**
 * Find a node's record, or create an empty one
 * @param a Store to modify
 * @param id Node id
 * @param created Set to 1 if the record was just created, else 0
 * @return Record index, or a negative AK_E* code
 */
static int64_t rec_upsert(annot_t *a, uint64_t id, int *created) {
    annmap_t *h = (annmap_t *)a->idx;
    int absent;
    khint_t k = annmap_put(h, id, &absent);
    if (!absent) { *created = 0; return (int64_t)kh_val(h, k); }
    if (reserve_rec(a) != AK_OK) { annmap_del(h, k); return AK_ENOMEM; }
    kh_val(h, k) = (uint32_t)a->n;
    annot_rec_t *r = &a->rec[a->n];
    r->id = id;
    r->off = 0;
    r->len = 0;
    r->kind = ANNOT_INFO;
    *created = 1;
    return a->n++;
}

// public API: store

// allocate an empty store; see akhal/annot.h
annot_t *annot_init(void) {
    annot_t *a = (annot_t *)calloc(1, sizeof(*a));
    if (!a) { ak_log(AK_LOG_ERROR, "annot", "out of memory"); return NULL; }
    a->idx = annmap_init();
    if (!a->idx) {
        free(a);
        ak_log(AK_LOG_ERROR, "annot", "out of memory");
        return NULL;
    }
    return a;
}

// free a store and everything it owns; see akhal/annot.h
void annot_destroy(annot_t *a) {
    if (!a) return;
    free(a->rec);
    free(a->buf);
    if (a->idx) annmap_destroy((annmap_t *)a->idx);
    free(a);
}

// set or replace a node's annotation; see akhal/annot.h
int annot_set(annot_t *a, uint64_t id, int kind, const char *info) {
    int created;
    int64_t i = rec_upsert(a, id, &created);
    if (i < 0) return (int)i;
    annot_rec_t *r = &a->rec[i];
    if (buf_put(a, info, &r->off, &r->len) != AK_OK) return AK_ENOMEM;
    r->kind = (uint8_t)kind;
    return AK_OK;
}

// append to a node's annotation; see akhal/annot.h
int annot_add(annot_t *a, uint64_t id, int kind, const char *info) {
    int created;
    int64_t i = rec_upsert(a, id, &created);
    if (i < 0) return (int)i;
    annot_rec_t *r = &a->rec[i];
    if (created || r->len == 0) {
        if (buf_put(a, info, &r->off, &r->len) != AK_OK) return AK_ENOMEM;
        if (created) r->kind = (uint8_t)kind;
        return AK_OK;
    }
    if (!info || !*info) return AK_OK;

    // merge "old; new" into a fresh buffer entry (old bytes stay orphaned)
    size_t ol = r->len, nl = strlen(info);
    char *merged = (char *)malloc(ol + 2 + nl + 1);
    if (!merged) return AK_ENOMEM;
    memcpy(merged, a->buf + r->off, ol);
    memcpy(merged + ol, "; ", 2);
    memcpy(merged + ol + 2, info, nl + 1);
    int rc = buf_put(a, merged, &r->off, &r->len);
    free(merged);
    return rc;
}

// look up a node's annotation; see akhal/annot.h
int annot_get(const annot_t *a, uint64_t id, const char **info) {
    if (info) *info = NULL;
    annmap_t *h = (annmap_t *)a->idx;
    khint_t k = annmap_get(h, id);
    if (k == kh_end(h)) return ANNOT_UNKNOWN;
    const annot_rec_t *r = &a->rec[kh_val(h, k)];
    if (info && r->len) *info = a->buf + r->off;
    return r->kind;
}

// name of an annotation kind; see akhal/annot.h
const char *annot_kind_name(int kind) {
    switch (kind) {
        case ANNOT_BACKBONE: return "backbone";
        case ANNOT_INFO:     return "annot";
        default:             return "unknown";
    }
}

// builders

// mark the backbone reference path; see akhal/annot.h
int32_t annot_backbone(annot_t *a, const gfa_t *g, const char *ref_path) {
    if (!g->path || g->n_path == 0) {
        ak_log(AK_LOG_WARN, "annot", "graph has no P lines; no backbone to mark");
        return AK_EINVAL;
    }

    int32_t k = 0;
    if (ref_path) {
        for (k = 0; k < g->n_path; k++)
            if (!strcmp(g->path[k], ref_path)) break;
        if (k == g->n_path) {
            ak_log(AK_LOG_ERROR, "annot", "no path named '%s' in the graph", ref_path);
            return AK_EINVAL;
        }
    } else if (g->n_path > 1) {
        ak_log(AK_LOG_INFO, "annot", "%d paths present; using '%s' as backbone",
               g->n_path, g->path[0]);
    }

    const char *name = g->path[k];
    size_t cap = strlen(name) + 48;
    char *info = (char *)malloc(cap);
    if (!info) return AK_ENOMEM;

    const uint32_t *segs;
    int np = gfa_path_segs(g, k, &segs);
    for (int i = 0; i < np; i++) {
        if (segs[i] == GFA_NIL) continue;
        const gfa_seg_t *s = &g->seg[segs[i]];
        snprintf(info, cap, "REF %s %d-%d", name, s->start, s->end);
        if (annot_set(a, s->id, ANNOT_BACKBONE, info) != AK_OK) {
            free(info);
            return AK_ENOMEM;
        }
    }
    free(info);
    return k;
}

/**
 * Classify a variant from its prefix-stripped core alleles
 * @param rl Length of the remaining REF core
 * @param al Length of the remaining ALT core
 * @return A static type tag such as "SNP"
 */
static const char *variant_type(size_t rl, size_t al) {
    if (rl == 0) return "INS";
    if (al == 0) return "DEL";
    if (rl == 1 && al == 1) return "SNP";
    if (rl == al) return "MNP";
    return "COMPLEX";
}

/**
 * Try to match an alternate allele along a non-backbone walk starting at w0.
 * The walk follows the unique non-backbone successor of each node and
 * succeeds only if the concatenated node sequences equal the allele exactly
 * @param g Graph to walk
 * @param bb Backbone membership flags, one per segment index
 * @param w0 First candidate node (a non-backbone successor of the branch)
 * @param ac Alternate allele core (prefix-stripped, not NUL-terminated)
 * @param al Length of ac
 * @param nodes Out: segment indices of the walk, at most ANNOT_MAX_WALK
 * @return Number of nodes in the matched walk, or 0 if it does not match
 */
static int walk_alt(const gfa_t *g, const uint8_t *bb, uint32_t w0, const char *ac, size_t al, uint32_t *nodes) {
    size_t pos = 0;
    int nn = 0;
    uint32_t cur = w0;

    for (;;) {
        const gfa_seg_t *s = &g->seg[cur];
        if (!s->seq || s->len == 0) return 0;
        if ((size_t)s->len > al - pos) return 0;
        if (memcmp(s->seq, ac + pos, s->len) != 0) return 0;
        if (nn >= ANNOT_MAX_WALK) return 0;
        nodes[nn++] = cur;
        pos += s->len;
        if (pos == al) return nn;

        // continue through the unique non-backbone successor
        const uint32_t *arcs;
        int na = gfa_arcs(g, (int32_t)cur, &arcs);
        int64_t next = -1;
        for (int t = 0; t < na; t++) {
            uint32_t w = g->link[arcs[t]].w;
            if (bb[w]) continue;
            if (next >= 0) return 0;   // ambiguous branch inside the bubble
            next = (int64_t)w;
        }
        if (next < 0) return 0;
        cur = (uint32_t)next;
    }
}

// annotate alt-path nodes from a VCF; see akhal/annot.h
int64_t annot_build_vcf(annot_t *a, const gfa_t *g, int32_t ref_path, const char *vcf_fn) {
    if (ref_path < 0 || ref_path >= g->n_path) {
        ak_log(AK_LOG_ERROR, "annot", "invalid backbone path index %d", ref_path);
        return AK_EINVAL;
    }
    if (!g->arc_off) {
        ak_log(AK_LOG_ERROR, "annot", "VCF annotation requires the graph to be read with GFA_LINKS");
        return AK_EINVAL;
    }
    const char *refname = g->path[ref_path];

    // backbone membership + branch index: reference end coordinate -> segment
    uint8_t *bb = (uint8_t *)calloc((size_t)g->n_seg, 1);
    annmap_t *cmap = annmap_init();
    uint32_t *walk = (uint32_t *)malloc(ANNOT_MAX_WALK * sizeof(uint32_t));
    if (!bb || !cmap || !walk) {
        free(bb); free(walk);
        if (cmap) annmap_destroy(cmap);
        return AK_ENOMEM;
    }

    const uint32_t *segs;
    int np = gfa_path_segs(g, ref_path, &segs);
    for (int i = 0; i < np; i++) {
        if (segs[i] == GFA_NIL) continue;
        bb[segs[i]] = 1;
        int absent;
        khint_t k = annmap_put(cmap, (uint64_t)g->seg[segs[i]].end, &absent);
        if (absent) kh_val(cmap, k) = segs[i];
    }

    vcf_reader_t *r = vcf_open(vcf_fn);
    if (!r) {
        free(bb); free(walk); annmap_destroy(cmap);
        return AK_EOPEN;
    }

    vcf_rec_t rec;
    vcf_rec_init(&rec);
    int64_t n_var = 0, n_allele = 0, n_hit = 0;
    int rc = AK_OK, warned_chrom = 0;
    char *info = NULL;
    size_t info_m = 0;

    int vrc;
    while ((vrc = vcf_read1(r, &rec)) == 1) {
        n_var++;
        if (strcmp(rec.chrom, refname) != 0) {
            if (!warned_chrom) {
                ak_log(AK_LOG_WARN, "annot", "VCF contig '%s' does not match backbone '%s'; such records are skipped", rec.chrom, refname);
                warned_chrom = 1;
            }
            continue;
        }
        if (!rec.alt) continue;   // ALT '.'

        size_t rlen = strlen(rec.ref);
        const char *alt = rec.alt;
        while (alt && *alt) {
            const char *comma = strchr(alt, ',');
            size_t alen = comma ? (size_t)(comma - alt) : strlen(alt);
            n_allele++;

            // strip the shared REF/ALT prefix to find the branch coordinate
            size_t sp = 0;
            while (sp < rlen && sp < alen && rec.ref[sp] == alt[sp]) sp++;
            size_t rl = rlen - sp, al = alen - sp;
            uint64_t p0 = (uint64_t)rec.pos - 1 + sp;   // 0-based branch position

            if (al == 0) {   // pure deletion: an edge, not a node
                ak_log(AK_LOG_DEBUG, "annot", "%s:%lld %s>%.*s is a deletion; no node to annotate", rec.chrom, (long long)rec.pos, rec.ref, (int)alen, alt);
                alt = comma ? comma + 1 : NULL;
                continue;
            }

            khint_t k = annmap_get(cmap, p0);
            if (k == kh_end(cmap)) {
                ak_log(AK_LOG_DEBUG, "annot", "%s:%lld: no backbone branch node ends at offset %llu", rec.chrom, (long long)rec.pos, (unsigned long long)p0);
                alt = comma ? comma + 1 : NULL;
                continue;
            }
            uint32_t v = kh_val(cmap, k);

            // try each non-backbone successor of the branch node
            const uint32_t *arcs;
            int na = gfa_arcs(g, (int32_t)v, &arcs), nn = 0;
            for (int t = 0; t < na && nn == 0; t++) {
                uint32_t w = g->link[arcs[t]].w;
                if (bb[w]) continue;
                nn = walk_alt(g, bb, w, alt + sp, al, walk);
            }
            if (nn == 0) {
                ak_log(AK_LOG_DEBUG, "annot", "%s:%lld %s>%.*s: no alternative path matches", rec.chrom, (long long)rec.pos, rec.ref, (int)alen, alt);
                alt = comma ? comma + 1 : NULL;
                continue;
            }

            // "<TYPE> <chrom>:<pos> <ref>><alt> [<id>]"
            size_t need = strlen(rec.chrom) + rlen + alen +
                          (rec.id ? strlen(rec.id) : 0) + 48;
            if (need > info_m) {
                char *p = (char *)realloc(info, need);
                if (!p) { rc = AK_ENOMEM; break; }
                info = p;
                info_m = need;
            }
            snprintf(info, info_m, "%s %s:%lld %s>%.*s%s%s", variant_type(rl, al), rec.chrom, (long long)rec.pos, rec.ref, (int)alen, alt, rec.id ? " " : "", rec.id ? rec.id : "");

            for (int t = 0; t < nn; t++) {
                if (annot_add(a, g->seg[walk[t]].id, ANNOT_INFO, info) != AK_OK) {
                    rc = AK_ENOMEM;
                    break;
                }
            }
            if (rc != AK_OK) break;
            n_hit++;
            alt = comma ? comma + 1 : NULL;
        }
        if (rc != AK_OK) break;
    }
    if (vrc < 0) rc = vrc;

    vcf_rec_clear(&rec);
    vcf_close(r);
    free(info);
    free(bb);
    free(walk);
    annmap_destroy(cmap);

    if (rc != AK_OK) {
        ak_log(AK_LOG_ERROR, "annot", "VCF annotation failed: %s", ak_strerror(rc));
        return rc;
    }
    ak_log(AK_LOG_INFO, "annot", "%lld variants read, %lld of %lld alleles matched to alternative paths",
           (long long)n_var, (long long)n_hit, (long long)n_allele);
    return n_hit;
}

// annotate nodes by tracing FASTA sequences; see akhal/annot.h
int64_t annot_build_fasta(annot_t *a, const gfa_t *g, const fasta_t *fa) {
    if (!g->arc_off) {
        ak_log(AK_LOG_ERROR, "annot", "FASTA annotation requires the graph to be read with GFA_LINKS");
        return AK_EINVAL;
    }

    int64_t n_done = 0;
    for (int64_t fi = 0; fi < fasta_n(fa); fi++) {
        const fasta_rec_t *fr = &fa->rec[fi];
        const char *S = fr->seq;
        int64_t L = fr->len;
        if (L <= 0) continue;

        // start node: prefer a source whose sequence prefixes the record
        int32_t start = -1, fallback = -1;
        for (int32_t i = 0; i < g->n_seg; i++) {
            const gfa_seg_t *s = &g->seg[i];
            if (!s->seq || s->len == 0) continue;
            size_t cl = (int64_t)s->len < L ? s->len : (size_t)L;
            if (memcmp(s->seq, S, cl) != 0) continue;
            if (s->in_degree == 0) { start = i; break; }
            if (fallback < 0) fallback = i;
        }
        if (start < 0) start = fallback;
        if (start < 0) {
            ak_log(AK_LOG_WARN, "annot", "no starting node matches sequence '%s'", fr->name);
            continue;
        }

        size_t cap = strlen(fr->name) + 40;
        char *info = (char *)malloc(cap);
        if (!info) return AK_ENOMEM;

        int64_t pos = 0;          // next unconsumed base of S
        int32_t cur = start;
        uint32_t cur_ov = 0;      // overlap consumed on entry to cur

        for (;;) {
            const gfa_seg_t *s = &g->seg[cur];
            int64_t novel = (int64_t)s->len - cur_ov;

            if (annot_get(a, s->id, NULL) != ANNOT_BACKBONE) {
                snprintf(info, cap, "SEQ %s %lld", fr->name, (long long)(pos - cur_ov));
                if (annot_add(a, s->id, ANNOT_INFO, info) != AK_OK) {
                    free(info);
                    return AK_ENOMEM;
                }
            }
            pos += novel < L - pos ? novel : L - pos;
            if (pos >= L) { n_done++; break; }

            // successor whose sequence continues S at pos (past any overlap)
            const uint32_t *arcs;
            int na = gfa_arcs(g, cur, &arcs);
            int32_t next = -1;
            uint32_t next_ov = 0;
            for (int t = 0; t < na; t++) {
                const gfa_link_t *e = &g->link[arcs[t]];
                const gfa_seg_t *w = &g->seg[e->w];
                if (!w->seq || e->overlap >= w->len) continue;
                int64_t cl = (int64_t)(w->len - e->overlap);
                if (cl > L - pos) cl = L - pos;
                if (memcmp(w->seq + e->overlap, S + pos, (size_t)cl) != 0) continue;
                next = (int32_t)e->w;
                next_ov = e->overlap;
                break;
            }
            if (next < 0) {
                ak_log(AK_LOG_WARN, "annot", "walk of '%s' stalled at %lld/%lld bp (node %llu)",
                       fr->name, (long long)pos, (long long)L,
                       (unsigned long long)s->id);
                break;
            }
            cur = next;
            cur_ov = next_ov;
        }
        free(info);
    }

    ak_log(AK_LOG_INFO, "annot", "%lld of %lld sequences fully traced",
           (long long)n_done, (long long)fasta_n(fa));
    return n_done;
}

// file I/O

// write a store to a binary .annot file; see akhal/annot.h
int annot_write(const annot_t *a, const char *fn) {
    FILE *fp = fopen(fn, "wb");
    if (!fp) {
        ak_log(AK_LOG_ERROR, "annot", "cannot open '%s' for writing", fn);
        return AK_EOPEN;
    }

    int ok = fwrite(ANNOT_MAGIC, 1, 8, fp) == 8;
    uint64_t n = (uint64_t)a->n;
    ok = ok && fwrite(&n, sizeof(n), 1, fp) == 1;
    ok = ok && fwrite(&a->buf_l, sizeof(a->buf_l), 1, fp) == 1;
    for (int64_t i = 0; ok && i < a->n; i++) {
        const annot_rec_t *r = &a->rec[i];
        ok = fwrite(&r->id,   sizeof(r->id),   1, fp) == 1 &&
             fwrite(&r->off,  sizeof(r->off),  1, fp) == 1 &&
             fwrite(&r->len,  sizeof(r->len),  1, fp) == 1 &&
             fwrite(&r->kind, sizeof(r->kind), 1, fp) == 1;
    }
    if (ok && a->buf_l)
        ok = fwrite(a->buf, 1, a->buf_l, fp) == a->buf_l;

    if (fclose(fp) != 0) ok = 0;
    if (!ok) {
        ak_log(AK_LOG_ERROR, "annot", "write failed on '%s'", fn);
        return AK_EIO;
    }
    return AK_OK;
}

// load a store from a .annot file; see akhal/annot.h
annot_t *annot_read(const char *fn) {
    FILE *fp = fopen(fn, "rb");
    if (!fp) {
        ak_log(AK_LOG_ERROR, "annot", "cannot open '%s'", fn);
        return NULL;
    }

    char magic[8];
    uint64_t n = 0, buf_l = 0;
    if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, ANNOT_MAGIC, 8) != 0 ||
        fread(&n, sizeof(n), 1, fp) != 1 ||
        fread(&buf_l, sizeof(buf_l), 1, fp) != 1) {
        fclose(fp);
        ak_log(AK_LOG_ERROR, "annot", "'%s' is not an akhal annotation file", fn);
        return NULL;
    }

    annot_t *a = annot_init();
    if (!a) { fclose(fp); return NULL; }

    int ok = 1;
    if (n) {
        a->rec = (annot_rec_t *)malloc((size_t)n * sizeof(annot_rec_t));
        ok = a->rec != NULL;
        a->m = ok ? (int64_t)n : 0;
    }
    annmap_t *h = (annmap_t *)a->idx;
    for (uint64_t i = 0; ok && i < n; i++) {
        annot_rec_t *r = &a->rec[i];
        ok = fread(&r->id,   sizeof(r->id),   1, fp) == 1 &&
             fread(&r->off,  sizeof(r->off),  1, fp) == 1 &&
             fread(&r->len,  sizeof(r->len),  1, fp) == 1 &&
             fread(&r->kind, sizeof(r->kind), 1, fp) == 1;
        // the info string plus its NUL must fit inside the buffer
        if (ok && r->len && r->off + r->len >= buf_l) ok = 0;
        if (ok) {
            int absent;
            khint_t k = annmap_put(h, r->id, &absent);
            if (absent) kh_val(h, k) = (uint32_t)i;
        }
    }
    a->n = ok ? (int64_t)n : 0;

    if (ok && buf_l) {
        a->buf = (char *)malloc(buf_l);
        ok = a->buf && fread(a->buf, 1, buf_l, fp) == buf_l;
        if (ok) {
            a->buf_l = a->buf_m = buf_l;
            ok = a->buf[buf_l - 1] == '\0';   // buffer must end on a separator
        }
    }
    fclose(fp);

    if (!ok) {
        ak_log(AK_LOG_ERROR, "annot", "'%s' is truncated or corrupt", fn);
        annot_destroy(a);
        return NULL;
    }
    return a;
}
