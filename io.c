#include "io.h"

FILE *io_open(const char* file_path, char **line, int cap) {

    FILE* file = fopen(file_path, "r");
    if (!file) { fprintf(stderr, "[ERROR] Failed to open file %s\n", file_path); exit(EXIT_FAILURE); }
    
    *line = (char *)malloc(cap);

    if (!(*line)) { fprintf(stderr, "[ERROR] Failed to allocate memory."); exit(EXIT_FAILURE); }

    return file;
}

void io_close(FILE *file, char **str) {
    if (file)
        fclose(file);
    if (*str) {
        free(*str);
    }
    file = NULL;
}

int io_read(FILE *file, char **str, size_t *cap) {

    char *line = *str;
    size_t line_len = *cap;

    if (fgets(line, line_len, file)) {
        size_t len = strlen(line);

        // read the line and fit it to `line`
        while (len == line_len - 1 && line[len - 1] != '\n') {
            line_len *= 2;
            *cap = line_len;
            char *temp_line = (char *)realloc(line, line_len);
            if (!temp_line) { fprintf(stderr, "[ERROR] Memory reallocation failed.\n"); return 0; }
            line = temp_line;
            *str = temp_line;

            if (fgets(line + len, line_len - len, file) == NULL) { break; }
            len = strlen(line);
        }
        if (len && line[len - 1] == '\n') { line[len - 1] = '\0'; len--; }

        return len;
    }

    return 0;
}

void write_sam_record(FILE *out_sam, alignment *aln, char *ops, int c_size, const char *rname, int pos, int simplify) {
    int flag = 0;
    if (aln->strand == '-') flag |= 0x10;  // reverse complemented

    // Compose the CIGAR string from `ops`
    char cigar_string[MAX_CIGAR];
    int cigar_pos = 0;

    if (c_size == 0) {
        strcpy(cigar_string, "*");
    } else {
        int i = 0;
        while (i < c_size && ops[i] == CIGAR_SOFT_CLIP) {
            i++;
        }
        while (i < c_size && (ops[i] == CIGAR_INSERTION || ops[i] == CIGAR_SEQUENCE_MISMATCH)) {
            ops[i++] = CIGAR_SOFT_CLIP;
        }

        int j = c_size - 1;
        while (j >= 0 && ops[i] == CIGAR_SOFT_CLIP) {
            j--;
        }
        while (j >= 0 && (ops[j] == CIGAR_INSERTION || ops[j] == CIGAR_SEQUENCE_MISMATCH)) {
            ops[j--] = CIGAR_SOFT_CLIP;
        }

        if (simplify) {
            int k = 0; 
            while (k < c_size) {
                if (ops[k] == CIGAR_SEQUENCE_MATCH || ops[k] == CIGAR_SEQUENCE_MISMATCH)
                    ops[k] = CIGAR_ALIGNMENT_MATCH;
                k++;
            }
        }

        int count = 1;
        for (int i = 1; i < c_size; i++) {
            if (ops[i] == ops[i - 1]) {
                count++;
            } else {
                cigar_pos += sprintf(cigar_string + cigar_pos, "%d%c", count, ops[i-1]);
                count = 1;
            }
        }
        // Final op
        cigar_pos += sprintf(cigar_string + cigar_pos, "%d%c", count, ops[c_size-1]);
    }

    // Default values
    const char *rnext = "*";
    int pnext = 0;
    int tlen = 0;
    const char *seq = aln->read;
    const char *qual = "*";

    // POS is 1-based in SAM, so add 1
    int mapq = aln->qual;  // mapping quality (can be 255)

    fprintf(out_sam,
            "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t"
            "NM:i:%d\tAS:f:%.2f\tdv:f:%.6f\tid:f:%.6f\tRG:Z:%s.%d\n",
            aln->readName,
            flag,
            rname,
            pos,
            mapq,
            cigar_string,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            aln->xdi,
            aln->score,
            aln->divergence,
            aln->identity,
            "akhal",
            0
    );
}