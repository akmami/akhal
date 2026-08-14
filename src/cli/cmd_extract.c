#include "akhal/gfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <string.h>

#define FASTA_WRAP 80

/**
 * Stream a sequence out as FASTA, wrapping at FASTA_WRAP columns
 * @param out Destination stream
 * @param seq Bases to write
 * @param len Number of bases
 * @param col Current column position, carried across calls (updated in place)
 */
static void emit_wrapped(FILE *out, const char *seq, size_t len, int *col) {
    size_t printed = 0;
    while (printed < len) {
        size_t space = (size_t)(FASTA_WRAP - *col);
        size_t take = (len - printed < space) ? (len - printed) : space;
        fwrite(seq + printed, 1, take, out);
        printed += take;
        *col += (int)take;
        if (*col >= FASTA_WRAP) { fputc('\n', out); *col = 0; }
    }
}

// print the extract usage line
static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal extract fa <r/GFA> <out.fa|.fasta>");
}

// `extract` entry point; see cli.h
int cmd_extract(int argc, char **argv) {
    if (argc < 5 || strcmp(argv[2], "fa") != 0) { usage(); return 1; }

    const char *in = argv[3], *out_fn = argv[4];
    if (!ak_ends_with(in, ".gfa") && !ak_ends_with(in, ".rgfa")) { usage(); return 1; }
    if (!ak_ends_with(out_fn, ".fa") && !ak_ends_with(out_fn, ".fasta")) { usage(); return 1; }

    FILE *out = fopen(out_fn, "w");
    if (!out) { ak_log(AK_LOG_ERROR, NULL, "cannot open output %s", out_fn); return 1; }

    gfa_t *g = gfa_read(in, GFA_PATHS);
    if (!g) { fclose(out); return 1; }

    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        fprintf(out, ">%s\n", gfa_path_name(g, k));

        const uint32_t *segs;
        int n = gfa_path_segs(g, k, &segs);
        int col = 0;
        for (int i = 0; i < n; i++) {
            if (segs[i] == GFA_NIL) continue;
            gfa_seg_t *s = gfa_seg_at(g, (int32_t)segs[i]);
            if (s->seq) emit_wrapped(out, s->seq, s->len, &col);
        }
        if (col) fputc('\n', out);
    }

    fclose(out);
    gfa_destroy(g);
    ak_log(AK_LOG_INFO, NULL, "extracted reference to %s", out_fn);
    return 0;
}
