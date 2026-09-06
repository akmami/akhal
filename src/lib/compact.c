#include "akhal/compact.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

// which end of a segment a link lands on
#define END_L 0
#define END_R 1

// everything the run-finding pass needs, one entry per segment
typedef struct {
    int32_t  deg[2];    // links touching the left end and the right end
    uint32_t succ;      // the segment a forward join leaves this one for
    uint32_t pred;      // the segment such a join arrives from
    uint8_t  term[2];   // a path starts or ends at this end
} node_t;

// Count what sits on each end, and note the forward joins.
//
// A link records the end it uses at each of its two segments: leaving `v` it
// uses v's right end when the orientation is '+' and its left end when '-',
// and arriving at `w` it uses w's left end on '+' and its right end on '-'.
// The two spellings of one forward join, `L u + v +` and `L v - u -`, both come
// out as u's right joined to v's left
static void scan_links(const gfa_t *g, node_t *n) {
    for (int32_t i = 0; i < gfa_n_link(g); i++) {
        const gfa_link_t *e = gfa_link_at(g, i);
        n[e->v].deg[e->from_orient == '+' ? END_R : END_L]++;
        n[e->w].deg[e->to_orient   == '+' ? END_L : END_R]++;

        // bases shared across a join would have to be dropped to concatenate
        if (e->overlap != 0) continue;

        if (e->from_orient == '+' && e->to_orient == '+') {
            n[e->v].succ = e->w;
            n[e->w].pred = e->v;
        } else if (e->from_orient == '-' && e->to_orient == '-') {
            n[e->w].succ = e->v;
            n[e->v].pred = e->w;
        }
        // a join that flips strand is left for the degree counts alone
    }
}

// Mark the ends where a path stops. The end a path dangles from is the one it
// does not carry on through: its first step leaves the left end of a '+' step
// unentered, its last step leaves the right end of a '+' step unused, and a '-'
// step is the same the other way round
static void scan_paths(const gfa_t *g, node_t *n) {
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        const uint32_t *segs;
        int ns = gfa_path_segs(g, k, &segs);
        const char *ori = g->path_ori + g->path_off[k];

        int32_t first = -1, last = -1;
        for (int i = 0; i < ns; i++) {
            if (segs[i] == GFA_NIL) continue;
            if (first < 0) first = i;
            last = i;
        }
        if (first < 0) continue;

        n[segs[first]].term[ori[first] == '+' ? END_L : END_R] = 1;
        n[segs[last]].term[ori[last]   == '+' ? END_R : END_L] = 1;
    }
}

// 1 when the two carry the same stable sequence and sit back to back on it.
// Only asked of a file that ranks itself; a plain GFA has nothing to preserve
static int tags_agree(const gfa_seg_t *u, const gfa_seg_t *v) {
    if (u->rank != v->rank) return 0;
    if ((u->ref_name == NULL) != (v->ref_name == NULL)) return 0;
    if (u->ref_name && strcmp(u->ref_name, v->ref_name) != 0) return 0;
    if (u->start < 0 || v->start < 0) return u->start < 0 && v->start < 0;
    return v->start == u->start + (int32_t)u->len;
}

// keep only the joins that may actually be taken; see the header for the rules
static void filter_joins(const gfa_t *g, node_t *n) {
    for (int32_t u = 0; u < gfa_n_seg(g); u++) {
        uint32_t v = n[u].succ;
        int ok = v != GFA_NIL &&
                 v != (uint32_t)u &&
                 n[u].deg[END_R] == 1 &&
                 n[v].deg[END_L] == 1 &&
                 n[v].pred == (uint32_t)u &&
                 !n[u].term[END_R] &&
                 !n[v].term[END_L] &&
                 (!g->has_sr || tags_agree(gfa_seg_at(g, u), gfa_seg_at(g, (int32_t)v)));
        if (!ok) {
            if (v != GFA_NIL) n[v].pred = GFA_NIL;
            n[u].succ = GFA_NIL;
        }
    }

    // A segment with two forward joins off one end kept only the last of them
    // in `succ`, while every target kept its `pred`. Those left pointing back
    // at a segment that no longer points forward at them begin runs of their
    // own, and this is what says so
    for (int32_t v = 0; v < gfa_n_seg(g); v++) {
        uint32_t u = n[v].pred;
        if (u != GFA_NIL && n[u].succ != (uint32_t)v) n[v].pred = GFA_NIL;
    }
}

// walk one run from its first segment, numbering the segments as it goes
static void take_run(const node_t *n, compact_t *c, uint32_t start, int32_t r) {
    c->first[r] = start;
    c->count[r] = 0;

    uint32_t s = start, last = GFA_NIL;
    while (s != GFA_NIL && c->run[s] == GFA_NIL) {
        c->run[s] = (uint32_t)r;
        c->pos[s] = c->count[r]++;
        last = s;
        s = n[s].succ;
    }
    // whatever the walk stopped at is not part of this run, whether the chain
    // simply ended or closed back on itself. Cutting the link here is what
    // keeps a cycle from being walked round for ever
    if (last != GFA_NIL) c->next[last] = GFA_NIL;
}

// work out the runs; see akhal/compact.h
compact_t *compact_runs(const gfa_t *g) {
    if (!(g->flags & GFA_LINKS)) {
        ak_log(AK_LOG_ERROR, "compact", "compaction requires the graph to be read with GFA_LINKS");
        return NULL;
    }

    int32_t ns = gfa_n_seg(g);
    size_t one = (size_t)(ns > 0 ? ns : 1);

    compact_t *c = (compact_t *)calloc(1, sizeof(compact_t));
    node_t *n = (node_t *)calloc(one, sizeof(node_t));
    if (c) {
        c->run   = (uint32_t *)malloc(one * sizeof(uint32_t));
        c->pos   = (int32_t *)malloc(one * sizeof(int32_t));
        c->next  = (uint32_t *)malloc(one * sizeof(uint32_t));
        c->first = (uint32_t *)malloc(one * sizeof(uint32_t));
        c->count = (int32_t *)malloc(one * sizeof(int32_t));
    }
    if (!c || !n || !c->run || !c->pos || !c->next || !c->first || !c->count) {
        ak_log(AK_LOG_ERROR, "compact", "out of memory");
        free(n);
        compact_destroy(c);
        return NULL;
    }
    c->n_seg = ns;

    for (int32_t i = 0; i < ns; i++) {
        n[i].succ = n[i].pred = GFA_NIL;
        c->run[i] = GFA_NIL;
    }

    scan_links(g, n);
    if (g->flags & GFA_PATHS) {
        scan_paths(g, n);
    }
    filter_joins(g, n);
    for (int32_t i = 0; i < ns; i++) c->next[i] = n[i].succ;

    // a run begins wherever nothing joins into the segment from behind
    int32_t r = 0;
    for (int32_t i = 0; i < ns; i++) {
        if (n[i].pred == GFA_NIL && c->run[i] == GFA_NIL) {
            take_run(n, c, (uint32_t)i, r++);
        }
    }
    // what is left is a cycle of joins with no beginning; cut it at the lowest
    // segment so it comes out as one node carrying a self-loop
    for (int32_t i = 0; i < ns; i++) {
        if (c->run[i] == GFA_NIL) {
            take_run(n, c, (uint32_t)i, r++);
        }
    }

    c->n_run = r;
    c->n_merged = ns - r;
    free(n);
    return c;
}

// free a run set; see akhal/compact.h
void compact_destroy(compact_t *c) {
    if (!c) return;
    free(c->run);
    free(c->pos);
    free(c->next);
    free(c->first);
    free(c->count);
    free(c);
}

// the id a run answers to: the one its first segment had
static uint64_t run_id(const gfa_t *g, const compact_t *c, uint32_t r) {
    return gfa_seg_at(g, (int32_t)c->first[r])->id;
}

// one S line per run: every segment's bases end to end, and the first
// segment's tags, which the join rules made sure still describe the whole run
static void write_segs(const gfa_t *g, const compact_t *c, FILE *out, int tags) {
    for (int32_t r = 0; r < c->n_run; r++) {
        const gfa_seg_t *head = gfa_seg_at(g, (int32_t)c->first[r]);
        fprintf(out, "S\t%llu\t", (unsigned long long)head->id);

        int any = 0;
        for (uint32_t s = c->first[r]; s != GFA_NIL; s = c->next[s]) {
            const gfa_seg_t *seg = gfa_seg_at(g, (int32_t)s);
            if (seg->seq && seg->len) {
                fwrite(seg->seq, 1, seg->len, out);
                any = 1;
            }
        }
        if (!any) {
            fputc('*', out);
        }

        if (tags && head->ref_name) {
            fprintf(out, "\tSN:Z:%s", head->ref_name);
        }
        if (tags && head->start >= 0) {
            fprintf(out, "\tSO:i:%d", head->start);
        }
        if (tags && head->rank >= 0) {
            fprintf(out, "\tSR:i:%d", head->rank);
        }
        fputc('\n', out);
    }
}

// every link except the ones that became the inside of a run. An end that
// survived belongs to the run's first or last segment, so the orientation the
// link carried still points the same way once the endpoints are remapped
static void write_links(const gfa_t *g, const compact_t *c, FILE *out) {
    for (int32_t i = 0; i < gfa_n_link(g); i++) {
        const gfa_link_t *e = gfa_link_at(g, i);
        uint32_t rv = c->run[e->v], rw = c->run[e->w];

        // the join inside a run, in either of the two ways it can be spelled;
        // nothing else linking the same run is one of them
        int inside = rv == rw &&
                     ((e->from_orient == '+' && e->to_orient == '+' && c->pos[e->w] == c->pos[e->v] + 1) ||
                      (e->from_orient == '-' && e->to_orient == '-' && c->pos[e->v] == c->pos[e->w] + 1));
        if (inside) continue;

        fprintf(out, "L\t%llu\t%c\t%llu\t%c\t%uM\n",
                (unsigned long long)run_id(g, c, rv), e->from_orient,
                (unsigned long long)run_id(g, c, rw), e->to_orient,
                e->overlap);
    }
}

// paths, with the steps that walked a run one by one collapsed into the single
// node that replaced it
static void write_paths(const gfa_t *g, const compact_t *c, FILE *out) {
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        fprintf(out, "P\t%s\t", gfa_path_name(g, k));

        const uint32_t *segs;
        int ns = gfa_path_segs(g, k, &segs);
        const char *ori = g->path_ori + g->path_off[k];

        int written = 0;
        uint32_t prev_run = GFA_NIL;
        int32_t prev_pos = 0;
        char prev_ori = '+';

        for (int i = 0; i < ns; i++) {
            if (segs[i] == GFA_NIL) continue;
            uint32_t r = c->run[segs[i]];

            if (written && r == prev_run) {
                // still walking the run the last step was on: one more step
                // along it folds into the node that replaced it
                int32_t want = ori[i] == '-' ? prev_pos - 1 : prev_pos + 1;
                if (ori[i] == prev_ori && c->pos[segs[i]] == want) {
                    prev_pos = c->pos[segs[i]];
                    continue;
                }

                // not carrying on. A cyclic run wrapping round to where it
                // began is a second traversal and gets a step of its own;
                // anything else means the P line did not follow the links,
                // and the record it spells is not this module's to fix
                int32_t end = c->count[r] - 1;
                int wrap = ori[i] == prev_ori &&
                           ((ori[i] == '+' && prev_pos == end && c->pos[segs[i]] == 0) ||
                            (ori[i] == '-' && prev_pos == 0 && c->pos[segs[i]] == end));
                if (!wrap) {
                    ak_log(AK_LOG_WARN, "compact", "path %s reaches segment %llu without a link to it; it is written as a fresh step", gfa_path_name(g, k), (unsigned long long)gfa_seg_at(g, (int32_t)segs[i])->id);
                }
            }

            fprintf(out, "%s%llu%c", written ? "," : "", (unsigned long long)run_id(g, c, r), ori[i]);
            written = 1;
            prev_run = r;
            prev_pos = c->pos[segs[i]];
            prev_ori = ori[i];
        }
        fprintf(out, "\t*\n");
    }
}

// write the compacted graph; see akhal/compact.h
int compact_write(const gfa_t *g, const compact_t *c, FILE *out) {
    fprintf(out, "H\tVN:Z:1.0\n");

    write_segs(g, c, out, g->has_sr);
    write_links(g, c, out);
    write_paths(g, c, out);

    if (ferror(out)) {
        ak_log(AK_LOG_ERROR, "compact", "write failed");
        return AK_EIO;
    }
    return AK_OK;
}
