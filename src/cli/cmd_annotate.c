#include "akhal/gfa.h"
#include "akhal/fasta.h"
#include "akhal/annot.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <string.h>

// print the annotate usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal annotate <r/GFA> <out.annot> [--vcf <FILE>] [--fasta <FILE>] [--ref <PATH-NAME>]");
}

// `annotate` entry point; see cli.h
int cmd_annotate(int argc, char **argv) {
    const char *gfa_fn = NULL, *out_fn = NULL;
    const char *vcf_fn = NULL, *fa_fn = NULL, *ref_name = NULL;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--vcf") && i + 1 < argc) {
            vcf_fn = argv[++i];
        } else if (!strcmp(argv[i], "--fasta") && i + 1 < argc) {
            fa_fn = argv[++i];
        } else if (!strcmp(argv[i], "--ref") && i + 1 < argc) {
            ref_name = argv[++i];
        } else if (argv[i][0] == '-') {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[i]);
            usage();
            return 1;
        } else if (!gfa_fn) {
            gfa_fn = argv[i];
        } else if (!out_fn) {
            out_fn = argv[i];
        } else {
            usage();
            return 1;
        }
    }
    if (!gfa_fn || !out_fn) { usage(); return 1; }
    if (!ak_ends_with(gfa_fn, ".gfa") && !ak_ends_with(gfa_fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", gfa_fn);
        return 1;
    }

    gfa_t *g = gfa_read(gfa_fn, GFA_LINKS | GFA_PATHS);
    if (!g) return 1;

    annot_t *an = annot_init();
    if (!an) { gfa_destroy(g); return 1; }
    int ret = 1;

    // backbone first; without one, VCF matching is impossible and nodes simply stay unknown
    int32_t bk = annot_backbone(an, g, ref_name);
    if (bk == AK_ENOMEM) goto done;
    if (bk < 0 && ref_name) goto done;   // an explicitly named path must exist

    if (vcf_fn) {
        if (bk < 0) {
            ak_log(AK_LOG_WARN, NULL, "--vcf ignored: no backbone path in the graph");
        } else if (annot_build_vcf(an, g, bk, vcf_fn) < 0) {
            goto done;
        }
    }

    if (fa_fn) {
        fasta_t *fa = fasta_read(fa_fn);
        if (!fa) goto done;
        int64_t rc = annot_build_fasta(an, g, fa);
        fasta_destroy(fa);
        if (rc < 0) goto done;
    }

    if (annot_write(an, out_fn) != AK_OK) goto done;

    int64_t n_bb = 0, n_info = 0;
    for (int64_t i = 0; i < annot_n(an); i++) {
        if (annot_at(an, i)->kind == ANNOT_BACKBONE) n_bb++;
        else n_info++;
    }
    printf("Nodes: %ld\n", (long)gfa_n_seg(g));
    printf("Backbone: %lld\n", (long long)n_bb);
    printf("Annotated: %lld\n", (long long)n_info);
    printf("Unknown: %lld\n", (long long)(gfa_n_seg(g) - n_bb - n_info));
    ret = 0;

done:
    annot_destroy(an);
    gfa_destroy(g);
    return ret;
}
