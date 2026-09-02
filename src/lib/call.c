#include "akhal/call.h"
#include "akhal/kstr.h"
#include "akhal/error.h"
#include "akhal/version.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// growth helpers

static int reserve_var(call_t *c) {
    if (c->n < c->m) return AK_OK;
    int64_t m = c->m ? c->m << 1 : 1024;
    call_var_t *p = (call_var_t *)realloc(c->var, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    c->var = p;
    c->m = m;
    return AK_OK;
}

// keeps the buffer NUL-terminated
static int ks_add(kstring_t *ks, const char *s, size_t len) {
    if (ks_resize(ks, ks->l + len + 1) != AK_OK) return AK_ENOMEM;
    memcpy(ks->s + ks->l, s, len);
    ks->l += len;
    ks->s[ks->l] = '\0';
    return AK_OK;
}

// backbone

static call_ref_t *ref_alloc(const gfa_t *g, const char *name) {
    call_ref_t *r = (call_ref_t *)calloc(1, sizeof(*r));
    if (!r) {
        ak_log(AK_LOG_ERROR, "call", "out of memory");
        return NULL;
    }

    r->n_seg = gfa_n_seg(g);
    size_t n = (size_t)(r->n_seg > 0 ? r->n_seg : 1);
    r->name = strdup(name);
    r->on   = (uint8_t *)calloc(n, 1);
    r->pos  = (int64_t *)malloc(n * sizeof(int64_t));
    if (!r->name || !r->on || !r->pos) {
        call_ref_destroy(r);
        ak_log(AK_LOG_ERROR, "call", "out of memory");
        return NULL;
    }
    for (int32_t i = 0; i < r->n_seg; i++) r->pos[i] = -1;
    return r;
}

// the requested name wins; without one, the file's first path, else the longest
static int32_t ref_pick(const gfa_t *g, const gfa_merge_t *m, const char *path_name) {
    if (m->n == 1) return 0;

    if (!path_name) {
        // no preference: stay with the path the file listed first
        for (int32_t k = 0; k < m->n; k++)
            for (int32_t i = m->off[k]; i < m->off[k + 1]; i++)
                if (m->frag[i] == 0) return k;
        return 0;
    }

    // the name asked for wins outright; otherwise take the longest chain
    int32_t best = 0;
    int64_t best_len = -1;
    for (int32_t k = 0; k < m->n; k++) {
        if (!strcmp(m->name[k], path_name)) return k;
        int64_t len = 0;
        for (int32_t i = m->off[k]; i < m->off[k + 1]; i++)
            len += (int64_t)gfa_path_len(g, m->frag[i]);
        if (len > best_len) {
            best_len = len;
            best = k;
        }
    }
    return best;
}

// label the backbone from the P lines; see akhal/call.h
call_ref_t *call_ref_path(const gfa_t *g, const char *path_name) {
    if (!(g->flags & GFA_PATHS) || g->n_path == 0) {
        ak_log(AK_LOG_ERROR, "call", "graph has no P lines to use as a backbone");
        return NULL;
    }

    // one reference is often split over several P lines, so the backbone is a
    // chain of fragments; an unfragmented graph simply yields chains of one
    gfa_merge_t *m = gfa_path_merge(g, path_name);
    if (!m) return NULL;

    int32_t k = ref_pick(g, m, path_name);
    if (m->n > 1) {
        ak_log(AK_LOG_INFO, "call", "%d candidate backbone(s); using '%s'", m->n, m->name[k]);
    }

    uint32_t *segs = NULL;
    int64_t ns = gfa_merge_segs(g, m, k, &segs);
    if (ns < 0) {
        ak_log(AK_LOG_ERROR, "call", "cannot assemble backbone '%s': %s", m->name[k], ak_strerror((int)ns));
        gfa_merge_destroy(m);
        return NULL;
    }

    int64_t total = 0;
    for (int64_t i = 0; i < ns; i++) total += gfa_seg_at(g, (int32_t)segs[i])->len;

    call_ref_t *r = ref_alloc(g, m->name[k]);
    if (r) {
        r->seq = (char *)malloc((size_t)total + 1);
    }
    if (!r || !r->seq) {
        call_ref_destroy(r);
        free(segs);
        gfa_merge_destroy(m);
        ak_log(AK_LOG_ERROR, "call", "out of memory");
        return NULL;
    }

    for (int64_t i = 0; i < ns; i++) {
        const gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
        // a segment the backbone visits twice keeps its first coordinate
        if (!r->on[segs[i]]) {
            r->on[segs[i]] = 1;
            r->pos[segs[i]] = r->len;
        }
        if (s->seq && s->len) {
            memcpy(r->seq + r->len, s->seq, s->len);
        }
        r->len += s->len;
    }
    r->seq[r->len] = '\0';

    // the flattened chain is already the ordered walk; keep it
    r->walk = segs;
    r->n_walk = ns;

    ak_log(AK_LOG_INFO, "call", "backbone '%s': %lld bp over %lld node(s), stitched from %d P line(s)", r->name, (long long)r->len, (long long)ns, m->off[k + 1] - m->off[k]);
    gfa_merge_destroy(m);
    return r;
}

// append one visited node to the backbone walk, growing it as needed
static int walk_record(call_ref_t *r, int64_t *m, uint32_t v) {
    if (r->n_walk >= *m) {
        int64_t nm = *m ? *m << 1 : 1024;
        uint32_t *p = (uint32_t *)realloc(r->walk, (size_t)nm * sizeof(*p));
        if (!p) return AK_ENOMEM;
        r->walk = p;
        *m = nm;
    }
    r->walk[r->n_walk++] = v;
    return AK_OK;
}

// label the backbone by tracing a FASTA record; see akhal/call.h
call_ref_t *call_ref_fasta(const gfa_t *g, const fasta_t *fa, const char *seq_name) {
    if (!(g->flags & GFA_LINKS)) {
        ak_log(AK_LOG_ERROR, "call", "tracing a backbone requires the graph to be read with GFA_LINKS");
        return NULL;
    }

    const fasta_rec_t *fr = NULL;
    if (seq_name) {
        fr = fasta_get(fa, seq_name);
        if (!fr) {
            ak_log(AK_LOG_ERROR, "call", "no sequence named '%s' in the FASTA", seq_name);
            return NULL;
        }
    } else {
        if (fasta_n(fa) == 0) {
            ak_log(AK_LOG_ERROR, "call", "the FASTA holds no sequence");
            return NULL;
        }
        fr = &fa->rec[0];
        if (fasta_n(fa) > 1) {
            ak_log(AK_LOG_INFO, "call", "%lld sequences present; using '%s' as the backbone", (long long)fasta_n(fa), fr->name);
        }
    }
    if (fr->len <= 0) {
        ak_log(AK_LOG_ERROR, "call", "sequence '%s' is empty", fr->name);
        return NULL;
    }

    call_ref_t *r = ref_alloc(g, fr->name);
    if (!r) return NULL;
    r->seq = (char *)malloc((size_t)fr->len + 1);
    if (!r->seq) {
        call_ref_destroy(r);
        ak_log(AK_LOG_ERROR, "call", "out of memory");
        return NULL;
    }
    memcpy(r->seq, fr->seq, (size_t)fr->len);
    r->seq[fr->len] = '\0';
    r->len = fr->len;

    const char *S = r->seq;
    int64_t L = r->len;

    // start node: prefer a source whose sequence prefixes the record
    int32_t start = -1, fallback = -1;
    for (int32_t i = 0; i < gfa_n_seg(g); i++) {
        const gfa_seg_t *s = gfa_seg_at(g, i);
        if (!s->seq || s->len == 0) continue;
        size_t cl = (int64_t)s->len < L ? s->len : (size_t)L;
        if (memcmp(s->seq, S, cl) != 0) continue;
        if (s->in_degree == 0) {
            start = i;
            break;
        }
        if (fallback < 0) {
            fallback = i;
        }
    }
    if (start < 0) {
        start = fallback;
    }
    if (start < 0) {
        call_ref_destroy(r);
        ak_log(AK_LOG_ERROR, "call", "no starting node matches sequence '%s'", fr->name);
        return NULL;
    }

    int64_t p = 0;         // next unconsumed base of the record
    int32_t cur = start;
    uint32_t cur_ov = 0;   // overlap consumed on entry to cur
    int warned_ov = 0, n_node = 0;
    int64_t walk_m = 0;

    for (;;) {
        const gfa_seg_t *s = gfa_seg_at(g, cur);
        if (cur_ov && !warned_ov) {
            ak_log(AK_LOG_WARN, "call", "backbone enters node %llu over a %u bp overlap; variant coordinates assume blunt joins", (unsigned long long)s->id, cur_ov);
            warned_ov = 1;
        }
        if (!r->on[cur]) {
            r->on[cur] = 1;
            r->pos[cur] = p - cur_ov;
            n_node++;
        }
        // every visit is recorded, so a node the walk returns to appears twice
        if (walk_record(r, &walk_m, (uint32_t)cur) != AK_OK) {
            call_ref_destroy(r);
            ak_log(AK_LOG_ERROR, "call", "out of memory");
            return NULL;
        }

        int64_t novel = (int64_t)s->len - cur_ov;
        p += novel < L - p ? novel : L - p;
        if (p >= L) break;

        // successor whose sequence continues the record at p, past any overlap
        const uint32_t *arcs;
        int na = gfa_arcs(g, cur, &arcs);
        int32_t next = -1;
        uint32_t next_ov = 0;
        for (int t = 0; t < na; t++) {
            const gfa_link_t *e = gfa_link_at(g, (int32_t)arcs[t]);
            const gfa_seg_t *wn = gfa_seg_at(g, (int32_t)e->w);
            if (!wn->seq || e->overlap >= wn->len) continue;
            int64_t cl = (int64_t)(wn->len - e->overlap);
            if (cl > L - p) {
                cl = L - p;
            }
            if (memcmp(wn->seq + e->overlap, S + p, (size_t)cl) != 0) continue;
            next = (int32_t)e->w;
            next_ov = e->overlap;
            break;
        }
        if (next < 0) {
            // keep the covered prefix; nothing past it can be placed
            ak_log(AK_LOG_WARN, "call", "backbone walk of '%s' stalled at %lld/%lld bp (node %llu); variants past that point are not reported", fr->name, (long long)p, (long long)L, (unsigned long long)s->id);
            r->len = p;
            r->seq[p] = '\0';
            break;
        }
        cur = next;
        cur_ov = next_ov;
    }

    ak_log(AK_LOG_INFO, "call", "backbone '%s': %lld bp over %d node(s)", r->name, (long long)r->len, n_node);
    return r;
}

// free a backbone; see akhal/call.h
void call_ref_destroy(call_ref_t *r) {
    if (!r) return;
    free(r->name);
    free(r->seq);
    free(r->on);
    free(r->pos);
    free(r->walk);
    free(r);
}

// variants

// one detour off the backbone, before alleles sharing a REF span are grouped
typedef struct {
    int64_t  end;   // 0-based backbone offset where the detour rejoins
    char    *seq;   // owned: bases it spells, "" for a deletion
} cand_t;

// state carried through the walk of one branch point
typedef struct {
    const gfa_t      *g;
    const call_ref_t *ref;
    int64_t           from;     // 0-based offset the detours branch from
    uint8_t          *stack;    // n_seg flags: 1 while a node is on the walk
    kstring_t         buf;      // bases spelled so far
    cand_t           *cand;     // collected detours, at most CALL_MAX_ALT
    int               n_cand;
    int               capped;   // 1 once a cap cut the search short
} walk_t;

// order the detours of one branch point so a shared REF span sits together
static int cand_cmp(const void *A, const void *B) {
    const cand_t *a = (const cand_t *)A, *b = (const cand_t *)B;
    if (a->end != b->end) return a->end < b->end ? -1 : 1;
    return strcmp(a->seq, b->seq);
}

// drop the detours of a branch point, ready for the next one
static void cand_clear(walk_t *w) {
    for (int i = 0; i < w->n_cand; i++) free(w->cand[i].seq);
    w->n_cand = 0;
}

// records one detour, ignoring joins back up the reference and alleles already
// collected at this branch point
static int cand_add(walk_t *w, int64_t end) {
    if (end < w->from) return AK_OK;
    const char *seq = w->buf.l ? w->buf.s : "";

    for (int i = 0; i < w->n_cand; i++)
        if (w->cand[i].end == end && !strcmp(w->cand[i].seq, seq)) return AK_OK;
    if (w->n_cand >= CALL_MAX_ALT) {
        w->capped = 1;
        return AK_OK;
    }

    char *cp = strdup(seq);
    if (!cp) return AK_ENOMEM;
    w->cand[w->n_cand].end = end;
    w->cand[w->n_cand].seq = cp;
    w->n_cand++;
    return AK_OK;
}

// follows one detour off the backbone; a nested bubble contributes one allele
// per way through it
static int walk_detour(walk_t *w, uint32_t v, int depth) {
    if (depth >= CALL_MAX_WALK || w->stack[v]) {
        w->capped = 1;
        return AK_OK;
    }

    const gfa_seg_t *s = gfa_seg_at(w->g, (int32_t)v);
    size_t mark = w->buf.l;
    if (mark + s->len > CALL_MAX_LEN) {
        w->capped = 1;
        return AK_OK;
    }
    if (s->len && ks_add(&w->buf, s->seq, s->len) != AK_OK) return AK_ENOMEM;
    w->stack[v] = 1;

    const uint32_t *arcs;
    int na = gfa_arcs(w->g, (int32_t)v, &arcs);
    int rc = AK_OK;
    for (int i = 0; i < na && rc == AK_OK; i++) {
        uint32_t u = gfa_link_at(w->g, (int32_t)arcs[i])->w;
        rc = w->ref->on[u] ? cand_add(w, w->ref->pos[u]) : walk_detour(w, u, depth + 1);
    }

    // unwind, so the next branch spells from the same prefix
    w->stack[v] = 0;
    w->buf.l = mark;
    if (w->buf.s) {
        w->buf.s[mark] = '\0';
    }
    return rc;
}

// classifies an allele from the lengths of the differing cores
static const char *var_type(int64_t rl, size_t al) {
    if (rl == 0) return "INS";
    if (al == 0) return "DEL";
    if (rl == 1 && al == 1) return "SNP";
    if ((size_t)rl == al) return "MNP";
    return "COMPLEX";
}

// one record per REF span: alleles sharing a span become a single multi-allelic
// row, and a span with an empty side is anchored on the preceding base
static int emit_group(call_t *c, const call_ref_t *ref, walk_t *w) {
    int64_t p = w->from;
    qsort(w->cand, (size_t)w->n_cand, sizeof(cand_t), cand_cmp);

    for (int i = 0; i < w->n_cand; ) {
        int j = i;
        while (j < w->n_cand && w->cand[j].end == w->cand[i].end) j++;

        int64_t q = w->cand[i].end, rl = q - p;
        int anchor = (rl == 0);
        for (int k = i; k < j && !anchor; k++) {
            if (w->cand[k].seq[0] == '\0') {
                anchor = 1;
            }
        }
        if (anchor && p == 0) {   // no preceding base to anchor on
            ak_log(AK_LOG_DEBUG, "call", "%s:1: allele needs an anchor before the reference start; skipped", ref->name);
            i = j;
            continue;
        }
        int64_t beg = anchor ? p - 1 : p;

        char *rf = (char *)malloc((size_t)(q - beg) + 1);
        int rc = rf ? AK_OK : AK_ENOMEM;
        if (rc == AK_OK) {
            memcpy(rf, ref->seq + beg, (size_t)(q - beg));
            rf[q - beg] = '\0';
        }

        kstring_t alt = KS_INIT, type = KS_INIT;
        for (int k = i; k < j && rc == AK_OK; k++) {
            const char *as = w->cand[k].seq;
            size_t al = strlen(as);
            if (k > i) {
                rc = ks_add(&alt, ",", 1);
                if (rc == AK_OK) {
                    rc = ks_add(&type, ",", 1);
                }
            }
            if (rc == AK_OK && anchor) {
                rc = ks_add(&alt, ref->seq + p - 1, 1);
            }
            if (rc == AK_OK) {
                rc = ks_add(&alt, as, al);
            }
            if (rc == AK_OK) {
                const char *tag = var_type(rl, al);
                rc = ks_add(&type, tag, strlen(tag));
            }
        }
        if (rc == AK_OK) {
            rc = reserve_var(c);
        }
        if (rc != AK_OK) {
            free(rf);
            ks_free(&alt);
            ks_free(&type);
            return AK_ENOMEM;
        }

        call_var_t *v = &c->var[c->n++];
        v->pos  = beg;
        v->end  = q;
        v->ref  = rf;
        v->alt  = ks_release(&alt);
        v->type = ks_release(&type);
        i = j;
    }
    return AK_OK;
}

// order the output the way a VCF expects it
static int var_cmp(const void *A, const void *B) {
    const call_var_t *a = (const call_var_t *)A, *b = (const call_var_t *)B;
    if (a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
    if (a->end != b->end) return a->end < b->end ? -1 : 1;
    return 0;
}

// collect the graph's variations over a backbone; see akhal/call.h
call_t *call_variants(const gfa_t *g, const call_ref_t *ref) {
    if (!(g->flags & GFA_LINKS)) {
        ak_log(AK_LOG_ERROR, "call", "variant discovery requires the graph to be read with GFA_LINKS");
        return NULL;
    }
    if (ref->n_seg != gfa_n_seg(g)) {
        ak_log(AK_LOG_ERROR, "call", "backbone was built from a different graph");
        return NULL;
    }

    call_t *c = (call_t *)calloc(1, sizeof(*c));
    walk_t w;
    memset(&w, 0, sizeof(w));
    w.g = g;
    w.ref = ref;
    w.stack = (uint8_t *)calloc((size_t)(ref->n_seg > 0 ? ref->n_seg : 1), 1);
    w.cand  = (cand_t *)calloc(CALL_MAX_ALT, sizeof(cand_t));
    if (!c || !w.stack || !w.cand) {
        call_destroy(c);
        free(w.stack);
        free(w.cand);
        ak_log(AK_LOG_ERROR, "call", "out of memory");
        return NULL;
    }

    int rc = AK_OK;
    int64_t n_capped = 0;
    for (int32_t v = 0; v < gfa_n_seg(g) && rc == AK_OK; v++) {
        if (!ref->on[v]) continue;

        // the bubble opens right after this node, where the next reference base sits
        int64_t p = ref->pos[v] + gfa_seg_at(g, v)->len;
        if (p <= 0 || p >= ref->len) continue;

        w.from = p;
        w.capped = 0;
        ks_clear(&w.buf);

        const uint32_t *arcs;
        int na = gfa_arcs(g, v, &arcs);
        for (int i = 0; i < na && rc == AK_OK; i++) {
            uint32_t u = gfa_link_at(g, (int32_t)arcs[i])->w;
            if (!ref->on[u]) {
                rc = walk_detour(&w, u, 1);
            } else if (ref->pos[u] > p) {
                rc = cand_add(&w, ref->pos[u]);   // skips ahead: a plain deletion
            }
        }

        if (rc == AK_OK && w.n_cand > 0) {
            if (w.capped) {
                ak_log(AK_LOG_DEBUG, "call", "%s:%lld: search cut short at the allele/walk caps", ref->name, (long long)p);
                n_capped++;
            }
            rc = emit_group(c, ref, &w);
        }
        cand_clear(&w);
    }

    cand_clear(&w);
    ks_free(&w.buf);
    free(w.stack);
    free(w.cand);

    if (rc != AK_OK) {
        call_destroy(c);
        ak_log(AK_LOG_ERROR, "call", "variant discovery failed: %s", ak_strerror(rc));
        return NULL;
    }
    if (c->n > 1) {
        qsort(c->var, (size_t)c->n, sizeof(call_var_t), var_cmp);
    }

    if (n_capped) {
        ak_log(AK_LOG_WARN, "call", "%lld site(s) were too tangled to enumerate fully", (long long)n_capped);
    }
    ak_log(AK_LOG_INFO, "call", "%lld variant record(s) over '%s'", (long long)c->n, ref->name);
    return c;
}

// free a variant set; see akhal/call.h
void call_destroy(call_t *c) {
    if (!c) return;
    for (int64_t i = 0; i < c->n; i++) {
        free(c->var[i].ref);
        free(c->var[i].alt);
        free(c->var[i].type);
    }
    free(c->var);
    free(c);
}

// file I/O

// write a variant set as VCF 4.2; see akhal/call.h
int call_write_vcf(const call_t *c, const call_ref_t *ref, const char *fn) {
    FILE *fp = fopen(fn, "w");
    if (!fp) {
        ak_log(AK_LOG_ERROR, "call", "cannot open '%s' for writing", fn);
        return AK_EOPEN;
    }

    fprintf(fp, "##fileformat=VCFv4.2\n");
    fprintf(fp, "##source=akhal %s\n", AKHAL_VERSION_STR);
    fprintf(fp, "##contig=<ID=%s,length=%lld>\n", ref->name, (long long)ref->len);
    fprintf(fp, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the REF allele\">\n");
    fprintf(fp, "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of each ALT allele\">\n");
    fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    for (int64_t i = 0; i < c->n; i++) {
        const call_var_t *v = &c->var[i];
        fprintf(fp, "%s\t%lld\t.\t%s\t%s\t.\t.\tEND=%lld;TYPE=%s\n", ref->name, (long long)(v->pos + 1), v->ref, v->alt, (long long)v->end, v->type);
    }

    int ok = !ferror(fp);
    if (fclose(fp) != 0) {
        ok = 0;
    }
    if (!ok) {
        ak_log(AK_LOG_ERROR, "call", "write failed on '%s'", fn);
        return AK_EIO;
    }
    return AK_OK;
}
