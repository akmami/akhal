/**
 * `akhal sampoke <ref.fa> <in.sam> [out.sam]`
 *
 * Validate a SAM file (typically produced by `gaf2sam`) against a reference:
 * every '=' CIGAR position must match the reference base. Optionally writes a
 * filtered SAM containing only the lines that validate, annotated with read
 * groups derived from read-name prefixes.
 *
 * The reference is loaded through the FASTA module (no .fai required), and
 * CIGAR handling comes from the SAM module.
 */

#include "akhal/fasta.h"
#include "akhal/sam.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Check that each '=' in the CIGAR matches the reference base.
 * @param read Query sequence (SAM SEQ field).
 * @param read_len Length of read.
 * @param ref Reference sequence for this alignment's rname.
 * @param ref_len Length of ref.
 * @param pos 1-based reference start position.
 * @param cigar CIGAR string to walk.
 * @param ops Caller-provided scratch buffer of at least SAM_MAX_CIGAR bytes.
 * @return 1 if the alignment is consistent, 0 otherwise.
 */
static int validate_alignment(const char *read, int64_t read_len,
                              const char *ref, int64_t ref_len,
                              int pos, const char *cigar, char *ops) {
    int n = sam_cigar_expand(cigar, ops, SAM_MAX_CIGAR, 0);
    if (n < 0) { ak_log(AK_LOG_WARN, "sampoke", "CIGAR parse error"); return 0; }

    int64_t ref_i = (int64_t)pos - 1, read_i = 0;
    for (int i = 0; i < n; i++) {
        char op = ops[i];
        if (op == CIGAR_SEQUENCE_MATCH || op == CIGAR_SEQUENCE_MISMATCH ||
            op == CIGAR_ALIGNMENT_MATCH) {
            if (ref_i < 0 || ref_i >= ref_len || read_i >= read_len) return 0;
            if (op == CIGAR_SEQUENCE_MATCH && ref[ref_i] != 'N' && ref[ref_i] != read[read_i])
                return 0;
            read_i++; ref_i++;
        } else if (op == CIGAR_INSERTION) {
            read_i++;
        } else if (op == CIGAR_DELETION) {
            ref_i++;
        } else if (op == CIGAR_SOFT_CLIP || op == CIGAR_HARD_CLIP) {
            read_i++;
        }
        // N, P consume neither in this check
    }
    return 1;
}

/**
 * Collect the distinct read-group prefixes seen among alignment lines.
 * @param sam_fn Path to the SAM file to scan.
 * @param out Set to a newly allocated array of owned prefix strings.
 * @return Number of prefixes collected.
 */
static int collect_rg(const char *sam_fn, char ***out) {
    *out = NULL;
    ak_file *f = ak_open(sam_fn);
    if (!f) return 0;

    int size = 0, cap = 0;
    char **hdr = NULL;
    kstring_t line = KS_INIT;
    long len;

    while ((len = ak_getline(f, &line)) >= 0) {
        if (len == 0) continue;
        if (line.s[0] == '@') {
            if (line.s[1] == 'R' && line.s[2] == 'G') break;  // already grouped
            continue;
        }
        char *prefix = sam_rg_prefix(line.s);
        if (!prefix) continue;

        int exists = 0;
        for (int i = 0; i < size; i++)
            if (strcmp(prefix, hdr[i]) == 0) { exists = 1; break; }
        if (exists) { free(prefix); continue; }

        if (size == cap) {
            cap = cap ? cap << 1 : 128;
            char **t = (char **)realloc(hdr, (size_t)cap * sizeof(char *));
            if (!t) { free(prefix); break; }
            hdr = t;
        }
        hdr[size++] = prefix;
    }

    ks_free(&line);
    ak_close(f);
    *out = hdr;
    return size;
}

/** Write one @RG header line per collected read-group prefix. */
static void emit_rg_lines(FILE *out, char **rg, int rg_n) {
    for (int i = 0; i < rg_n; i++)
        fprintf(out, "@RG\tID:akhal.%d\tPL:PACBIO\tPU:%s\tSM:sample\n", i, rg[i]);
}

/**
 * Split a SAM line copy into tab-delimited fields.
 * @param s Line to tokenize in place.
 * @param field Destination array of field pointers.
 * @param max Capacity of field.
 * @return Number of fields found.
 */
static int split_fields(char *s, char **field, int max) {
    int n = 0;
    char *save, *tok = strtok_r(s, "\t", &save);
    while (tok && n < max) {
        field[n++] = tok;
        tok = strtok_r(NULL, "\t", &save);
    }
    return n;
}

/**
 * Validate every alignment in a SAM file against the reference, reporting a
 * correct/incorrect tally and optionally writing a filtered, RG-annotated SAM.
 * @param sam_fn Input SAM path.
 * @param out_fn Output SAM path, or NULL to only report counts.
 * @param ref Reference sequences.
 * @param rg Collected read-group prefixes.
 * @param rg_n Number of prefixes.
 * @return 0 on success, non-zero if a file could not be opened.
 */
static int check_sam(const char *sam_fn, const char *out_fn,
                     const fasta_t *ref, char **rg, int rg_n) {
    ak_file *f = ak_open(sam_fn);
    if (!f) return 1;

    FILE *out = NULL;
    if (out_fn) {
        out = fopen(out_fn, "w");
        if (!out) { ak_log(AK_LOG_ERROR, "sampoke", "cannot open output %s", out_fn); ak_close(f); return 1; }
    }

    char *ops  = (char *)malloc(SAM_MAX_CIGAR);
    if (!ops) { ak_log(AK_LOG_ERROR, "sampoke", "out of memory"); if (out) fclose(out); ak_close(f); return 1; }

    kstring_t line = KS_INIT, work = KS_INIT;
    char *field[16];
    int correct = 0, invalid = 0, print_rg = 1;
    long len;

    while ((len = ak_getline(f, &line)) >= 0) {
        if (len == 0) continue;

        if (line.s[0] == '@') {
            if (out) {
                if (line.s[1] == 'P' && line.s[2] == 'G' && rg_n && print_rg) {
                    emit_rg_lines(out, rg, rg_n);
                    print_rg = 0;
                }
                fprintf(out, "%s\n", line.s);
            }
            continue;
        }
        if (out && rg_n && print_rg) { emit_rg_lines(out, rg, rg_n); print_rg = 0; }

        // tokenize a copy so the original line stays intact for output
        if (ks_resize(&work, line.l + 1) != AK_OK) break;
        memcpy(work.s, line.s, line.l + 1);
        work.l = line.l;

        int nf = split_fields(work.s, field, 16);
        if (nf < 11) {
            ak_log(AK_LOG_WARN, "sampoke", "malformed SAM line");
            invalid++;
            continue;
        }

        const char *rname = field[2];
        int   pos   = atoi(field[3]);
        const char *cigar = field[5];
        const char *seq   = field[9];

        if (strcmp(rname, "*") == 0) {           // unmapped
            if (out) fprintf(out, "%s\n", line.s);
            continue;
        }

        const fasta_rec_t *chrom = fasta_get(ref, rname);
        if (!chrom) {
            ak_log(AK_LOG_WARN, "sampoke", "reference not found: %s", rname);
            invalid++;
            continue;
        }

        if (!validate_alignment(seq, (int64_t)strlen(seq), chrom->seq, chrom->len,
                                pos, cigar, ops)) {
            invalid++;
            continue;
        }

        correct++;
        if (!out) continue;

        if (rg_n) {
            char *prefix = sam_rg_prefix(field[0]);
            if (prefix) {
                for (int i = 0; i < rg_n; i++)
                    if (strcmp(prefix, rg[i]) == 0)
                        fprintf(out, "%s\tRG:Z:akhal.%d\n", line.s, i);
                free(prefix);
            } else {
                fprintf(out, "%s\n", line.s);
            }
        } else {
            fprintf(out, "%s\n", line.s);
        }
    }

    printf("[INFO] correct: %d, incorrect: %d\n", correct, invalid);

    free(ops);
    ks_free(&line);
    ks_free(&work);
    ak_close(f);
    if (out) fclose(out);
    return 0;
}

/** `sampoke` entry point; see cli.h. */
int cmd_sampoke(int argc, char **argv) {
    if (argc < 4) {
        ak_log(AK_LOG_ERROR, NULL, "usage: akhal sampoke <ref.fa> <in.sam> [out.sam]");
        return 1;
    }
    const char *fa_fn  = argv[2];
    const char *sam_fn = argv[3];
    const char *out_fn = (argc >= 5) ? argv[4] : NULL;

    if (!ak_ends_with(fa_fn, ".fa") && !ak_ends_with(fa_fn, ".fasta")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a FASTA reference: %s", fa_fn);
        return 1;
    }
    if (!ak_ends_with(sam_fn, ".sam")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .sam file: %s", sam_fn);
        return 1;
    }
    if (out_fn && !ak_ends_with(out_fn, ".sam")) {
        ak_log(AK_LOG_ERROR, NULL, "output must be a .sam file: %s", out_fn);
        return 1;
    }

    fasta_t *ref = fasta_read(fa_fn);
    if (!ref) return 1;

    char **rg = NULL;
    int rg_n = collect_rg(sam_fn, &rg);

    int rc = check_sam(sam_fn, out_fn, ref, rg, rg_n);

    for (int i = 0; i < rg_n; i++) free(rg[i]);
    free(rg);
    fasta_destroy(ref);
    return rc;
}
