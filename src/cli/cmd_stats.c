#include "akhal/gfa.h"
#include "akhal/gaf.h"
#include "akhal/sam.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "khashl.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// read name -> index into the read array
KHASHL_MAP_INIT(KH_LOCAL, rdmap_t, rdmap, const char *, uint32_t, kh_hash_str, kh_eq_str)

// print the stats usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal stats <r/GFA|GAF> [--cigar]");
}

// graph statistics

// `stats` over an r/GFA graph
static int stats_gfa(const char *fn) {
    // the summary needs counters, two distributions and two degrees per
    // segment - none of which is a graph, so it does not build one
    gfa_stat_t st;
    if (gfa_read_stats(fn, &st) != AK_OK) return 1;

    if (st.n_undefined > 0) {
        ak_log(AK_LOG_WARN, "stats", "%lld id(s) named by an L or P line are defined by no S line",
               (long long)st.n_undefined);
    }

    printf("Segment count: %lld\n", (long long)st.n_seg);
    printf("Rank 0 segment count: %lld\n", (long long)st.n_rank0);
    printf("Rank 0< segment count: %lld\n", (long long)(st.n_seg - st.n_rank0));
    printf("Segment avg length: %f\n", st.seg_mean);
    printf("Segment std length: %f\n", st.seg_sd);
    printf("Segment min. length %lu\n", (unsigned long)st.seg_min);
    printf("Segment max. length %lu\n", (unsigned long)st.seg_max);
    printf("Link count: %lld\n", (long long)st.n_link);
    printf("Link overlapping avg length: %f\n", st.ov_mean);
    printf("Link overlapping std length: %f\n", st.ov_sd);
    printf("Minimum in degree: %d\n", st.min_in);
    printf("Maximum in degree: %d\n", st.max_in);
    printf("Minimum out degree: %d\n", st.min_out);
    printf("Maximum out degree: %d\n", st.max_out);
    return 0;
}

// alignment statistics

// A GAF is read in one streaming pass and can be arbitrarily large, so the
// distributions are accumulated as they arrive rather than collected into an
// array first the way the graph ones are. Welford's update keeps the running
// mean and the sum of squared deviations without a second pass.
typedef struct {
    int64_t n;
    double  mean, m2;
    double  min, max;
} dist_t;

static void dist_add(dist_t *d, double x) {
    if (d->n == 0) {
        d->min = d->max = x;
    } else if (x < d->min) {
        d->min = x;
    } else if (x > d->max) {
        d->max = x;
    }

    d->n++;
    double delta = x - d->mean;
    d->mean += delta / (double)d->n;
    d->m2 += delta * (x - d->mean);
}

// population standard deviation, the same definition ak_variance() uses
static double dist_stddev(const dist_t *d) {
    return d->n ? ak_stddev(d->m2 / (double)d->n) : 0.0;
}

// print a distribution the way the graph stats print theirs: the mean and the
// standard deviation as reals, and the extremes in the unit they were measured
// in, so counts and lengths do not come out with a fractional part
static void dist_print(const char *label, const dist_t *d, int integral) {
    printf("%s avg: %f\n", label, d->mean);
    printf("%s std: %f\n", label, dist_stddev(d));
    if (integral) {
        printf("%s min.: %lld\n", label, (long long)d->min);
        printf("%s max.: %lld\n", label, (long long)d->max);
    } else {
        printf("%s min.: %f\n", label, d->min);
        printf("%s max.: %f\n", label, d->max);
    }
}

// one read, with the best alignment it received
typedef struct {
    char   *name;
    int64_t n_aln;
    double  best_ratio;   // highest (qend - qstart) / qlen
    double  best_ident;   // highest matches / block_len
} read_t;

// the array + dict layout the graph and the FASTA store use: records in one
// array, a hash table from read name to array index, and the array owns the
// names the keys point at
typedef struct {
    read_t  *rec;
    int64_t  n, m;
    rdmap_t *h;
} reads_t;

static int reads_grow(reads_t *t) {
    if (t->n < t->m) return AK_OK;
    int64_t m = t->m ? t->m << 1 : 4096;
    read_t *p = (read_t *)realloc(t->rec, (size_t)m * sizeof(*p));
    if (!p) return AK_ENOMEM;
    t->rec = p;
    t->m = m;
    return AK_OK;
}

// record one alignment against its read, keeping only the best figures
static int reads_add(reads_t *t, const char *name, double ratio, double ident) {
    khint_t k = rdmap_get(t->h, name);
    if (k < kh_end(t->h)) {
        read_t *r = &t->rec[kh_val(t->h, k)];
        r->n_aln++;
        if (ratio > r->best_ratio) {
            r->best_ratio = ratio;
        }
        if (ident > r->best_ident) {
            r->best_ident = ident;
        }
        return AK_OK;
    }

    if (reads_grow(t) != AK_OK) return AK_ENOMEM;

    // the record's own copy becomes the key, since the caller's name belongs to
    // a reusable gaf_rec_t that the next read frees
    char *own = strdup(name);
    if (!own) return AK_ENOMEM;

    read_t *r = &t->rec[t->n];
    r->name       = own;
    r->n_aln      = 1;
    r->best_ratio = ratio;
    r->best_ident = ident;

    int absent;
    k = rdmap_put(t->h, own, &absent);
    kh_val(t->h, k) = (uint32_t)t->n;
    t->n++;
    return AK_OK;
}

static void reads_free(reads_t *t) {
    for (int64_t i = 0; i < t->n; i++) free(t->rec[i].name);
    free(t->rec);
    rdmap_destroy(t->h);
}

// A run of this many bases or more is reported separately, which is the
// threshold gaftools uses to separate structural events from small ones.
#define STATS_LARGE_RUN 50

// D/I/X/= run counts, each with the subset spanning STATS_LARGE_RUN or more
typedef struct {
    int64_t n_del, n_del_large;
    int64_t n_ins, n_ins_large;
    int64_t n_mis, n_mis_large;
    int64_t n_mat, n_mat_large;
    int64_t n_perfect;
    int64_t n_missing;   // primary alignments carrying no cg:Z tag
} cigar_stat_t;

// tally the runs in a difference CIGAR such as "120=1X33=2I"
static void cigar_count(const char *cg, cigar_stat_t *cs) {
    int64_t n_run = 0;
    char last = '\0';

    for (const char *p = cg; *p; ) {
        int64_t len = 0;
        while (*p >= '0' && *p <= '9') {
            len = len * 10 + (int64_t)(*p - '0');
            p++;
        }
        if (!*p) break;   // a trailing count with no operation

        char op = *p++;
        int large = (len >= STATS_LARGE_RUN);
        switch (op) {
            case CIGAR_DELETION:
                cs->n_del++;
                cs->n_del_large += large;
                break;
            case CIGAR_INSERTION:
                cs->n_ins++;
                cs->n_ins_large += large;
                break;
            case CIGAR_SEQUENCE_MISMATCH:
                cs->n_mis++;
                cs->n_mis_large += large;
                break;
            case CIGAR_SEQUENCE_MATCH:
                cs->n_mat++;
                cs->n_mat_large += large;
                break;
            default:
                break;
        }
        n_run++;
        last = op;
    }

    // a lone '=' run is the whole alignment matching end to end
    if (n_run == 1 && last == CIGAR_SEQUENCE_MATCH) {
        cs->n_perfect++;
    }
}

// oriented nodes in a path string; a stable path name carries no '>' or '<',
// and is one target rather than none
static int64_t path_nodes(const char *p) {
    int64_t n = 0;
    for (; *p; p++) {
        if (*p == '>' || *p == '<') {
            n++;
        }
    }
    return n ? n : 1;
}

// `stats` over a GAF alignment file
static int stats_gaf(const char *fn, int want_cigar) {
    gaf_reader_t *r = gaf_open(fn);
    if (!r) return 1;

    reads_t reads = { NULL, 0, 0, rdmap_init() };
    if (!reads.h) {
        gaf_close(r);
        ak_log(AK_LOG_ERROR, NULL, "out of memory");
        return 1;
    }

    int64_t n_aln = 0, n_primary = 0, n_secondary = 0, aligned_bases = 0;
    dist_t d_mapq = {0}, d_ident = {0}, d_ratio = {0};
    dist_t d_qlen = {0}, d_block = {0}, d_nodes = {0};
    cigar_stat_t cs = {0};

    gaf_rec_t rec;
    gaf_rec_init(&rec);

    int rc;
    while ((rc = gaf_read1(r, &rec)) == 1) {
        n_aln++;

        // only a tp:A tag that is not P demotes an alignment, so a file without
        // the tag is all primary; a non-positive mapping quality demotes too
        int primary = !(rec.has_tp && rec.tp != 'P' && rec.tp != 'p');
        if (!primary || rec.mapq <= 0) {
            n_secondary++;
            continue;
        }
        n_primary++;

        // an empty query or block would divide by zero, and says nothing anyway
        double ratio = rec.qlen > 0 ? (double)(rec.qend - rec.qstart) / (double)rec.qlen : 0.0;
        double ident = rec.block_len > 0 ? (double)rec.matches / (double)rec.block_len : 0.0;

        aligned_bases += rec.matches;
        dist_add(&d_mapq, (double)rec.mapq);
        dist_add(&d_ident, ident);
        dist_add(&d_ratio, ratio);
        dist_add(&d_qlen, (double)rec.qlen);
        dist_add(&d_block, (double)rec.block_len);
        dist_add(&d_nodes, (double)path_nodes(rec.path ? rec.path : ""));

        if (want_cigar) {
            if (rec.cigar) {
                cigar_count(rec.cigar, &cs);
            } else {
                cs.n_missing++;
            }
        }

        if (reads_add(&reads, rec.qname ? rec.qname : "", ratio, ident) != AK_OK) {
            rc = AK_ENOMEM;
            break;
        }
    }

    gaf_rec_clear(&rec);
    gaf_close(r);

    if (rc < 0) {
        reads_free(&reads);
        ak_log(AK_LOG_ERROR, NULL, "cannot read %s: %s", fn, ak_strerror(rc));
        return 1;
    }

    // the average of each read's best alignment, not of every alignment: a read
    // aligned in ten places contributes its best figure once
    double avg_best_ident = 0.0, avg_best_ratio = 0.0;
    for (int64_t i = 0; i < reads.n; i++) {
        avg_best_ident += reads.rec[i].best_ident;
        avg_best_ratio += reads.rec[i].best_ratio;
    }
    if (reads.n > 0) {
        avg_best_ident /= (double)reads.n;
        avg_best_ratio /= (double)reads.n;
    }

    printf("Total alignments: %lld\n", (long long)n_aln);
    printf("  Primary: %lld\n", (long long)n_primary);
    printf("  Secondary: %lld\n", (long long)n_secondary);
    printf("Reads with at least one alignment: %lld\n", (long long)reads.n);
    printf("Total aligned bases: %lld\n", (long long)aligned_bases);
    printf("Average highest sequence identity: %f\n", avg_best_ident);
    printf("Average highest map ratio: %f\n", avg_best_ratio);
    dist_print("Mapping quality", &d_mapq, 1);
    dist_print("Sequence identity", &d_ident, 0);
    dist_print("Map ratio", &d_ratio, 0);
    dist_print("Query length", &d_qlen, 1);
    dist_print("Alignment block length", &d_block, 1);
    dist_print("Path node count", &d_nodes, 1);

    if (want_cigar) {
        printf("Cigar string statistics:\n");
        printf("  Total deletion regions: %lld (%lld >=%dbps)\n", (long long)cs.n_del, (long long)cs.n_del_large, STATS_LARGE_RUN);
        printf("  Total insertion regions: %lld (%lld >=%dbps)\n", (long long)cs.n_ins, (long long)cs.n_ins_large, STATS_LARGE_RUN);
        printf("  Total substitution regions: %lld (%lld >=%dbps)\n", (long long)cs.n_mis, (long long)cs.n_mis_large, STATS_LARGE_RUN);
        printf("  Total match regions: %lld (%lld >=%dbps)\n", (long long)cs.n_mat, (long long)cs.n_mat_large, STATS_LARGE_RUN);
        printf("  Total perfect alignments (exact match): %lld\n", (long long)cs.n_perfect);
    }

    printf("\n* Figures cover primary alignments with a mapping quality above 0\n");

    if (want_cigar && cs.n_missing > 0) {
        ak_log(AK_LOG_WARN, NULL, "%lld primary alignment(s) carry no cg:Z tag and were left out of the CIGAR statistics", (long long)cs.n_missing);
    }

    reads_free(&reads);
    return 0;
}

// `stats` entry point; see cli.h
int cmd_stats(int argc, char **argv) {
    const char *fn = NULL;
    int want_cigar = 0;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--cigar")) {
            want_cigar = 1;
        } else if (argv[i][0] == '-') {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[i]);
            usage();
            return 1;
        } else if (!fn) {
            fn = argv[i];
        } else {
            usage();
            return 1;
        }
    }
    if (!fn) {
        usage();
        return 1;
    }

    if (ak_ends_with(fn, ".gaf")) {
        return stats_gaf(fn, want_cigar);
    }
    if (ak_ends_with(fn, ".gfa") || ak_ends_with(fn, ".rgfa")) {
        if (want_cigar) {
            ak_log(AK_LOG_WARN, NULL, "--cigar only applies to GAF input; ignoring it");
        }
        return stats_gfa(fn);
    }

    ak_log(AK_LOG_ERROR, NULL, "expected a .gfa, .rgfa or .gaf file: %s", fn);
    return 1;
}
