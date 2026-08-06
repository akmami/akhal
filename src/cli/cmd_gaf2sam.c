/**
 * `akhal gaf2sam <r/GFA> <GAF> <reads.fa> <out.sam> [--simple]`
 *
 * Convert graph alignments (GAF) into linear SAM against the reference paths
 * of an (r)GFA. The graph, the read store and the GAF records all come from
 * libakhal; this file holds only the conversion algorithm, which projects each
 * graph alignment onto its reference path and builds the corresponding CIGAR.
 *
 * A segment is treated as "reference" when it belongs to a path
 * (ref_name != NULL); bases on non-reference (alt) segments become insertions,
 * and gaps between consecutive reference segments become deletions.
 */

#include "akhal/gfa.h"
#include "akhal/gaf.h"
#include "akhal/fasta.h"
#include "akhal/sam.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NODE 32768
#define OPS_CAP  SAM_MAX_CIGAR

// Per-run counters, kept together for the closing summary.
typedef struct {
    int reads, segments, strands, ranks, quals, others; // discards by reason
    int fwd, rc;
    int rank0, rank1, rankboth, rankinv;
} stats_t;

/**
 * Convert one GAF record and, on success, write a SAM line.
 * @param g Reference graph (read with GFA_PATHS).
 * @param rec GAF record to convert.
 * @param read The read's sequence, looked up from the FASTA store.
 * @param simplify If non-zero, emit '='/'X' as 'M' in the CIGAR.
 * @param sam Destination SAM stream.
 * @param ops Reusable per-base op scratch (>= OPS_CAP).
 * @param cigar_ops Reusable expanded-CIGAR scratch (>= OPS_CAP).
 * @param nodes Reusable node-index scratch (>= MAX_NODE).
 * @param st Counters, updated with the outcome.
 */
static void convert_one(const gfa_t *g, const gaf_rec_t *rec, const char *read,
                        int simplify, FILE *sam,
                        char *ops, char *cigar_ops, uint32_t *nodes,
                        stats_t *st) {
    // local, mutable copies of the coordinates we may rewrite for '-' reads
    int qlen = (int)rec->qlen, qstart = (int)rec->qstart, qend = (int)rec->qend;
    int plen = (int)rec->plen, pstart = (int)rec->pstart, pend = (int)rec->pend;

    // walk the path, collecting node indices
    const char *p = rec->path;
    uint64_t id = 0;
    char pstrand = '>';
    int node_size = 0, fwd = 0, rev = 0, is_ref = 0, is_alt = 0;
    const char *ref_name = NULL;
    int step;

    while ((step = gaf_path_next(p, &id, &pstrand))) {
        p += step;
        int32_t si = gfa_idx(g, id);
        if (si < 0) { st->segments++; return; }
        if (node_size >= MAX_NODE) { st->others++; return; }
        nodes[node_size++] = (uint32_t)si;
        if (pstrand == '>') fwd++; else rev++;
        const gfa_seg_t *s = gfa_seg_at(g, si);
        if (s->ref_name) { is_ref = 1; ref_name = s->ref_name; }
        else is_alt = 1;
    }

    if (fwd && rev) { st->strands++; return; }

    if (is_ref && !is_alt)      st->rank0++;
    else if (!is_ref && is_alt) st->rank1++;
    else if (is_ref && is_alt)  st->rankboth++;
    else { st->rankinv++; st->ranks++; return; }

    char out_strand;
    if (fwd) { st->fwd++; out_strand = '+'; }
    else     { st->rc++;  out_strand = '-'; }

    if (!is_ref) return;                 // discard alt-only alignments

    // reverse strand: flip node order, mirror coordinates, rc the read
    char *read_owned = NULL;
    if (rev) {
        for (int i = 0; i < node_size / 2; i++) {
            uint32_t t = nodes[i];
            nodes[i] = nodes[node_size - i - 1];
            nodes[node_size - i - 1] = t;
        }
        int pre = pstart, post = plen - pend, mid = plen - pre - post;
        pstart = post; pend = post + mid;
        int rpre = qstart, rpost = qlen - qend, rmid = qlen - rpre - rpost;
        qstart = rpost; qend = rpost + rmid;

        read_owned = strdup(read);
        if (!read_owned) { st->others++; return; }
        ak_revcomp(read_owned, strlen(read_owned));
        read = read_owned;
    }

    // build the per-base CIGAR op array
    int op_index = 0, overflow = 0;
#define PUSH(c) do { if (op_index < OPS_CAP) ops[op_index++] = (c); else overflow = 1; } while (0)

    const gfa_seg_t *seg = gfa_seg_at(g, nodes[0]);
    const gfa_seg_t *prev_ref = seg->ref_name ? seg : NULL;
    int p_length = (int)seg->len - pstart;
    int node_index = 0;

    for (int i = 0; i < qstart && !overflow; i++) PUSH(CIGAR_SOFT_CLIP);   // leading clip

    int n_ops = sam_cigar_expand(rec->cigar, cigar_ops, OPS_CAP, rev > 0);
    if (n_ops < 0) { st->others++; free(read_owned); return; }

    int reference_start = seg->start + pstart + 1;
    int ref_start_set = 0;
    int ci = 0;

    while (ci < n_ops && !overflow) {
        char op = cigar_ops[ci++];

        if (!ref_start_set && seg->ref_name && CIGAR_REF(op)) {
            if (op == CIGAR_SEQUENCE_MATCH) ref_start_set = 1;
            else reference_start++;
        }

        if (seg->ref_name)      PUSH(op);
        else if (CIGAR_QUE(op)) PUSH(CIGAR_INSERTION);

        if (CIGAR_REF(op)) {
            if (--p_length) continue;              // still inside this segment
            if (++node_index == node_size) break;  // consumed all nodes
            seg = gfa_seg_at(g, nodes[node_index]);
            p_length = (int)seg->len;
            if (seg->ref_name) {
                if (!ref_start_set) reference_start = seg->start + 1;
                if (prev_ref) {                    // gap between ref segments -> D
                    for (int x = prev_ref->end; x < seg->start && !overflow; x++)
                        PUSH(CIGAR_DELETION);
                }
                prev_ref = seg;
            }
        }
    }

    for (int i = qend; i < qlen && !overflow; i++) PUSH(CIGAR_SOFT_CLIP);   // trailing clip
#undef PUSH

    if (overflow) { st->others++; free(read_owned); return; }

    sam_rec_t sr;
    memset(&sr, 0, sizeof(sr));
    sr.qname    = rec->qname;
    sr.flag     = (out_strand == '-') ? SAM_FREVERSE : 0;
    sr.rname    = ref_name;
    sr.pos      = reference_start;
    sr.mapq     = (int)rec->mapq;
    sr.ops      = ops;
    sr.n_ops    = op_index;
    sr.simplify = simplify;
    sr.seq      = read;
    sr.nm       = rec->nm;
    sr.as       = rec->as;
    sr.dv       = rec->dv;
    sr.id       = rec->id;
    sr.rg       = "akhal.0";
    sam_write_record(sam, &sr);

    free(read_owned);
}

/** Print the gaf2sam usage line. @return 1 (a convenient failure code). */
static int usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal gaf2sam <r/GFA> <GAF> <reads.fa> <out.sam> [--simple]");
    return 1;
}

/** `gaf2sam` entry point; see cli.h. */
int cmd_gaf2sam(int argc, char **argv) {
    if (argc < 6) return usage();

    const char *gfa_fn = argv[2], *gaf_fn = argv[3];
    const char *fa_fn  = argv[4], *sam_fn = argv[5];
    if (!ak_ends_with(gfa_fn, ".gfa") && !ak_ends_with(gfa_fn, ".rgfa")) return usage();
    if (!ak_ends_with(gaf_fn, ".gaf")) return usage();
    if (!ak_ends_with(fa_fn, ".fa") && !ak_ends_with(fa_fn, ".fasta")) return usage();
    if (!ak_ends_with(sam_fn, ".sam")) return usage();

    int simplify = 0;
    if (argc >= 7) {
        if (strcmp(argv[6], "--simple") == 0) simplify = 1;
        else return usage();
    }

    gfa_t *g = gfa_read(gfa_fn, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;
    for (int32_t k = 0; k < gfa_n_link(g); k++) {
        if (gfa_link_at(g, k)->overlap != 0) {
            ak_log(AK_LOG_ERROR, "gaf2sam", "overlaps are non-zero; cannot convert");
            gfa_destroy(g);
            return 1;
        }
    }

    fasta_t *reads = fasta_read(fa_fn);
    if (!reads) { gfa_destroy(g); return 1; }
    ak_log(AK_LOG_INFO, NULL, "loaded %ld reads", (long)fasta_n(reads));

    gaf_reader_t *r = gaf_open(gaf_fn);
    if (!r) { fasta_destroy(reads); gfa_destroy(g); return 1; }

    FILE *sam = fopen(sam_fn, "w");
    if (!sam) {
        ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", sam_fn);
        gaf_close(r); fasta_destroy(reads); gfa_destroy(g);
        return 1;
    }

    sam_write_header(sam, g->path, gfa_n_path(g), g->path_len, "akhal");

    char *ops       = (char *)malloc(OPS_CAP);
    char *cigar_ops = (char *)malloc(OPS_CAP);
    uint32_t *nodes = (uint32_t *)malloc(MAX_NODE * sizeof(uint32_t));
    if (!ops || !cigar_ops || !nodes) {
        ak_log(AK_LOG_ERROR, "gaf2sam", "out of memory");
        free(ops); free(cigar_ops); free(nodes);
        fclose(sam); gaf_close(r); fasta_destroy(reads); gfa_destroy(g);
        return 1;
    }

    stats_t st;
    memset(&st, 0, sizeof(st));

    gaf_rec_t rec;
    gaf_rec_init(&rec);
    while (gaf_read1(r, &rec) == 1) {
        if (rec.mapq == 0) { st.quals++; continue; }
        const fasta_rec_t *rr = fasta_get(reads, rec.qname);
        if (!rr) { st.reads++; continue; }
        if (!rec.cigar) { st.others++; continue; }
        convert_one(g, &rec, rr->seq, simplify, sam, ops, cigar_ops, nodes, &st);
    }
    gaf_rec_clear(&rec);

    printf("[INFO] reads (fasta): %ld\n", (long)fasta_n(reads));
    printf("[INFO] aligned forward / reverse: %d / %d\n", st.fwd, st.rc);
    printf("[INFO] ranks only-0 / only-1 / both / invalid: %d / %d / %d / %d\n", st.rank0, st.rank1, st.rankboth, st.rankinv);
    printf("[INFO] discarded read/segment/strand/rank/qual/other: %d/%d/%d/%d/%d/%d\n", st.reads, st.segments, st.strands, st.ranks, st.quals, st.others);

    free(ops); free(cigar_ops); free(nodes);
    fclose(sam);
    gaf_close(r);
    fasta_destroy(reads);
    gfa_destroy(g);
    return 0;
}
