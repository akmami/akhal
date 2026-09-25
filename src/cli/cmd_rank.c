#include "akhal/gfa.h"
#include "akhal/call.h"
#include "akhal/fasta.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// print the rank usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal rank <in.gfa> [out.gfa] [--fasta <FILE>] [--ref <NAME>] [--verbose] [--footprint]");
}

// make an external reference the backbone: rank against its walk, then install
// that walk in place of whatever the graph called a path before
static int rank_from_fasta(gfa_t *g, const char *fa_fn, const char *ref_name, int64_t *n0) {
    fasta_t *fa = fasta_read(fa_fn);
    if (!fa) return AK_EOPEN;

    call_ref_t *ref = call_ref_fasta(g, fa, ref_name);
    fasta_destroy(fa);
    if (!ref) return AK_EINVAL;

    int rc = AK_OK;
    int64_t marked = gfa_rank_mark(g, ref->on);
    if (marked < 0) {
        rc = (int)marked;
    } else {
        *n0 = marked;
        gfa_clear_paths(g);
        rc = gfa_add_path(g, ref->name, ref->walk, NULL, ref->n_walk);
    }

    call_ref_destroy(ref);
    return rc;
}

// `rank` entry point; see cli.h
int cmd_rank(int argc, char **argv) {
    const char *in = NULL, *out_fn = NULL, *fa_fn = NULL, *ref_name = NULL;
    int report = cli_take_report(&argc, argv, 2);

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--fasta") && i + 1 < argc) {
            fa_fn = argv[++i];
        } else if (!strcmp(argv[i], "--ref") && i + 1 < argc) {
            ref_name = argv[++i];
        } else if (argv[i][0] == '-') {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[i]);
            usage();
            return 1;
        } else if (!in) {
            in = argv[i];
        } else if (!out_fn) {
            out_fn = argv[i];
        } else {
            usage();
            return 1;
        }
    }
    if (!in) {
        usage();
        return 1;
    }
    if (!ak_ends_with(in, ".gfa") && !ak_ends_with(in, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", in);
        return 1;
    }
    if (out_fn && !ak_ends_with(out_fn, ".gfa") && !ak_ends_with(out_fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .gfa/.rgfa file: %s", out_fn);
        return 1;
    }

    gfa_t *g = gfa_read(in, GFA_ALL | report);
    if (!g) return 1;

    int rc = AK_OK;
    int64_t n0 = 0;
    int32_t n_paths = 0;

    if (fa_fn) {
        rc = rank_from_fasta(g, fa_fn, ref_name, &n0);
        n_paths = gfa_n_path(g);
    } else if (gfa_n_path(g) == 0) {
        // with no P lines there is nothing to call a backbone, and ranking
        // against nothing would flatten an rGFA's own SR tags to all-rank-1
        ak_log(AK_LOG_ERROR, NULL, "%s has no P lines to rank against; supply a reference with --fasta", in);
        gfa_destroy(g);
        return 1;
    } else {
        if (ref_name) {
            ak_log(AK_LOG_WARN, NULL, "--ref only applies with --fasta; ignoring it");
        }
        // every P line is backbone, so the paths are ranked exactly as they
        // were read - nothing is rewritten
        n_paths = gfa_n_path(g);
        int64_t marked = gfa_rank_paths(g);
        if (marked < 0) {
            rc = (int)marked;
        } else {
            n0 = marked;
        }
    }
    if (rc != AK_OK) {
        ak_log(AK_LOG_ERROR, NULL, "cannot re-rank %s: %s", in, ak_strerror(rc));
        gfa_destroy(g);
        return 1;
    }

    FILE *out = stdout;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) {
            ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn);
            gfa_destroy(g);
            return 1;
        }
    }

    rc = gfa_write(g, out);
    if (out_fn) {
        fclose(out);
    }

    if (rc == AK_OK) {
        if (fa_fn) {
            ak_log(AK_LOG_INFO, NULL, "backbone taken from %s: %lld node(s) at rank 0, %lld at rank 1", fa_fn, (long long)n0, (long long)((int64_t)gfa_n_seg(g) - n0));
        } else {
            ak_log(AK_LOG_INFO, NULL, "backbone taken from %d P line(s): %lld node(s) at rank 0, %lld at rank 1", n_paths, (long long)n0, (long long)((int64_t)gfa_n_seg(g) - n0));
        }
    }

    gfa_destroy(g);
    return rc == AK_OK ? 0 : 1;
}
