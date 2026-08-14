#include "akhal/sam.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// expand a CIGAR string into a per-base op array; see akhal/sam.h
int sam_cigar_expand(const char *cigar, char *ops, int max_ops, int rev) {
    int i = 0, j = 0, num = 0;
    while (cigar[i]) {
        if ('0' <= cigar[i] && cigar[i] <= '9') {
            num = num * 10 + (cigar[i] - '0');
        } else {
            if (j + num >= max_ops) return -1;   // out of range
            for (int k = 0; k < num; k++) ops[j++] = cigar[i];
            num = 0;
        }
        i++;
    }
    if (j > max_ops) return -1;

    if (rev) {
        int k = 0;
        while (k < j / 2) {
            char t = ops[k];
            ops[k] = ops[j - k - 1];
            ops[j - k - 1] = t;
            k++;
        }
    }
    return j;   // number of ops
}

/**
 * Sort key giving canonical chromosome ordering for header emission
 * @param chr Reference sequence name (a leading "chr" is ignored)
 * @return 1..22 for autosomes, 23/24/25 for X/Y/M, 1000 otherwise
 */
static int chrom_rank(const char *chr) {
    if (strncmp(chr, "chr", 3) == 0) chr += 3;
    if (isdigit((unsigned char)chr[0])) {
        int n = atoi(chr);
        if (n >= 1 && n <= 22) return n;
    }
    if (strcmp(chr, "X") == 0) return 23;
    if (strcmp(chr, "Y") == 0) return 24;
    if (strcmp(chr, "M") == 0 || strcmp(chr, "MT") == 0) return 25;
    return 1000;   // non-canonical contigs go last
}

// write @HD, chromosome-ordered @SQ lines, @PG and @RG; see akhal/sam.h
void sam_write_header(FILE *out, char **names, int n,
                      const uint64_t *lens, const char *pg) {
    fprintf(out, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");

    int *printed = (int *)calloc(n > 0 ? (size_t)n : 1, sizeof(int));
    if (!printed) return;

    // selection sort by (chrom_rank, name) so we don't reorder the caller arrays
    for (int done = 0; done < n; done++) {
        int min_i = -1, min_rank = 0;
        for (int i = 0; i < n; i++) {
            if (printed[i]) continue;
            int r = chrom_rank(names[i]);
            if (min_i == -1 || r < min_rank || (r == min_rank && strcmp(names[i], names[min_i]) < 0)) {
                min_i = i;
                min_rank = r;
            }
        }
        if (min_i == -1) break;
        fprintf(out, "@SQ\tSN:%s\tLN:%lu\n", names[min_i], (unsigned long)lens[min_i]);
        printed[min_i] = 1;
    }
    free(printed);

    fprintf(out, "@PG\tID:%s\tPN:%s\tVN:1.0\n", pg, pg);
    fprintf(out, "@RG\tID:%s.0\tPL:%s\tPU:%s\tSM:%s\n", pg, "UNKNOWN", "UNKNOWN", "UNKNOWN");
}

// render and write one SAM alignment line; see akhal/sam.h
void sam_write_record(FILE *out, sam_rec_t *r) {
    char cigar_string[SAM_MAX_CIGAR];
    int  cigar_pos = 0;
    char *ops = r->ops;
    int   c_size = r->n_ops;

    if (c_size == 0) {
        strcpy(cigar_string, "*");
    } else {
        // fold leading/trailing insertion & mismatch runs into soft clips
        int i = 0;
        while (i < c_size && ops[i] == CIGAR_SOFT_CLIP) i++;
        while (i < c_size && (ops[i] == CIGAR_INSERTION || ops[i] == CIGAR_SEQUENCE_MISMATCH))
            ops[i++] = CIGAR_SOFT_CLIP;

        int j = c_size - 1;
        // fix vs original: test ops[j], not ops[i]. The original read ops[i] (where i may 
        // equal c_size), an out-of-bounds read that also disabled trailing-clip folding
        while (j >= 0 && ops[j] == CIGAR_SOFT_CLIP) j--;
        while (j >= 0 && (ops[j] == CIGAR_INSERTION || ops[j] == CIGAR_SEQUENCE_MISMATCH))
            ops[j--] = CIGAR_SOFT_CLIP;

        if (r->simplify) {
            for (int k = 0; k < c_size; k++)
                if (ops[k] == CIGAR_SEQUENCE_MATCH || ops[k] == CIGAR_SEQUENCE_MISMATCH)
                    ops[k] = CIGAR_ALIGNMENT_MATCH;
        }

        int count = 1;
        for (int k = 1; k < c_size; k++) {
            if (ops[k] == ops[k - 1]) {
                count++;
            } else {
                cigar_pos += sprintf(cigar_string + cigar_pos, "%d%c", count, ops[k - 1]);
                count = 1;
            }
        }
        cigar_pos += sprintf(cigar_string + cigar_pos, "%d%c", count, ops[c_size - 1]);
    }
    (void)cigar_pos;

    const char *rnext = "*";
    const char *seq   = r->seq ? r->seq : "*";
    const char *qual  = "*";

    fprintf(out,
            "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t"
            "NM:i:%d\tAS:f:%.2f\tdv:f:%.6f\tid:f:%.6f\tRG:Z:%s\n",
            r->qname, r->flag, r->rname, r->pos, r->mapq, cigar_string,
            rnext, 0, 0, seq, qual,
            r->nm, r->as, r->dv, r->id, r->rg ? r->rg : "akhal.0");
}

// derive a read-group prefix from a read name/line; see akhal/sam.h
char *sam_rg_prefix(const char *name) {
    if (name[0] == '@' || name[0] == '>') name++;

    const char *dot   = strchr(name, '.');
    const char *slash = strchr(name, '/');
    const char *uline = strchr(name, '_');
    const char *end   = NULL;
    if (dot && slash) end = dot < slash ? dot : slash;
    else if (dot)     end = dot;
    else if (slash)   end = slash;
    else if (uline)   end = uline;
    else              end = name + strlen(name);

    size_t len = (size_t)(end - name);
    char *prefix = (char *)malloc(len + 1);
    if (!prefix) return NULL;
    memcpy(prefix, name, len);
    prefix[len] = '\0';
    return prefix;
}
