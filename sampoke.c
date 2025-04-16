#include "sampoke.h"

void read_fasta(const char *fasta_path, struct ref_seq *seqs) {

    // Initialize chromosome array
    int chrom_count_cap = MAX_CHROM_SIZE;
    seqs->size = 0;
    seqs->chrs = (struct chr*)malloc(chrom_count_cap * sizeof(struct chr));
    int chrom_size_cap = MAX_LINE;

    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(fasta_path, &line, line_cap);
    int line_len = 0;
    int exist = 0;

    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == '>') {
            if (seqs->size && seqs->chrs[seqs->size].seq_size != 0) {
                seqs->size++;
                if (seqs->size == chrom_count_cap) {
                    chrom_count_cap *= 2;
                    struct chr *temp_chr = (struct chr*)realloc(seqs->chrs, chrom_count_cap);
                    if (temp_chr == NULL) {
                        fprintf(stderr, "[ERROR] Reallocation of chromosome array failed.\n");
                        exit(EXIT_FAILURE);
                    }
                    seqs->chrs = temp_chr; 
                }
            }
            seqs->chrs[seqs->size].seq_name = strdup(line+1);
            chrom_size_cap = MAX_LINE;
            seqs->chrs[seqs->size].seq_size = 0;
            seqs->chrs[seqs->size].seq = malloc(chrom_size_cap);
            exist = 1;
        } else {
            if (seqs->chrs[seqs->size].seq_size + line_len >= chrom_size_cap) {
                chrom_size_cap *= 2;
                char *temp_seq = realloc(seqs->chrs[seqs->size].seq, chrom_size_cap);
                if (temp_seq == NULL) {
                    fprintf(stderr, "[ERROR] Reallocation of chromosome sequence failed.\n");
                    exit(EXIT_FAILURE);
                }
                seqs->chrs[seqs->size].seq = temp_seq;
            }
            memcpy(seqs->chrs[seqs->size].seq + seqs->chrs[seqs->size].seq_size, line, line_len);
            seqs->chrs[seqs->size].seq_size += line_len;
        }
    }
    // process reference file
    if (exist) {
        seqs->size++;
    }

    io_close(file, &line);
}

int validate_alignment(const char *read, const char *ref, int pos, char *cigar) {

    char ops[MAX_CIGAR];
    int n_ops = parse_cigar(cigar, ops, MAX_CIGAR, 0);
    if (n_ops < 0) {
        fprintf(stderr, "[ERROR] CIGAR parse error.\n");
        return 0;
    }

    int ref_len_count = 0;

    for (int i = 0; i < n_ops; i++) {
        if (CIGAR_REF(ops[i])) ref_len_count++;
    }
    
    int ref_i = pos - 1;
    int read_i = 0;
    for (int i = 0; i < n_ops; i++) {
        if (ops[i] == CIGAR_SEQUENCE_MATCH || ops[i] == CIGAR_SEQUENCE_MISMATCH || ops[i] == CIGAR_ALIGNMENT_MATCH) {
            if (ops[i] == CIGAR_SEQUENCE_MATCH && ref[ref_i] != 'N' && ref[ref_i] != read[read_i]) {
                return 0;
            }
            read_i++; ref_i++;
        } else if (ops[i] == CIGAR_INSERTION) {
            read_i++;
        } else if (ops[i] == CIGAR_DELETION) {
            ref_i++;
        } else if (ops[i] == CIGAR_SOFT_CLIP || ops[i] == CIGAR_HARD_CLIP) {
            read_i++;
        }
        // ignoring other types (N, H, P)
    }
    
    return 1;
}

void check_sam(const char *sam_file, const char *out_file, const struct ref_seq *seqs) {
    FILE *file = fopen(sam_file, "r");
    if (!file) {
        fprintf(stderr, "[ERROR] SAM open failed- %s\n", sam_file);
        exit(1);
    }

    FILE *outfile = NULL;
    if (out_file != NULL) {
        outfile = fopen(out_file, "w");
        if (!outfile) {
            fprintf(stderr, "[ERROR] Failed to open output SAM file\n");
            fclose(file);
            exit(1);
        }
    }

    char line[MAX_LINE];
    int line_no = 0;
    int invalid_count = 0, correct_count = 0;
    while (fgets(line, sizeof(line), file)) {
        line_no++;
        if (line[0] == '@') { if (outfile != NULL) fputs(line, outfile); continue; }

        char qname[100], rname[100], cigar[131072], seq[MAX_LINE];
        int flag, pos, mapq, pnext, tlen;
        char rnext[100], qual[MAX_LINE];
        int fields = sscanf(line, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
                            qname, &flag, rname, &pos, &mapq, cigar,
                            rnext, &pnext, &tlen, seq, qual);
        
        if (fields < 11) {
            fprintf(stderr, "Line %d: Malformed SAM line\n", line_no);
            continue;
        }

        int chrom_index = -1;
        for (int i=0; i<seqs->size; i++) {
            if (strcmp(rname, seqs->chrs[i].seq_name) == 0) { chrom_index = i; break; }
        }
        if (chrom_index == -1) {
            fprintf(stderr, "Line %d: Reference not found: %s\n", line_no, rname);
            continue;
        }

        if (!validate_alignment(seq, seqs->chrs[chrom_index].seq, pos, cigar)) invalid_count++;
        else { correct_count++; if (outfile != NULL) fputs(line, outfile); }
    }

    printf("[INFO] correct: %d, incorrect: %d\n", correct_count, invalid_count);
    fclose(file);

    if (outfile != NULL) fclose(outfile);
}

int sampoke(int argc, char *argv[]) {
    
    if (argc < 4) {
        fprintf(stderr, "[ERROR] Usage: ./akhal sampoke [FASTA FILE] [SAM FILE] [OUT SAM]\n");
        return -1;
    }

    if (!ends_with(argv[2], ".fa") && !ends_with(argv[2], ".fasta")) {
        printf("[ERROR] Usage: ./akhal sampoke [FASTA FILE] [SAM FILE] [OUT SAM]\n");
        return -1;
    }

    if (!ends_with(argv[3], ".sam")) {
        printf("[ERROR] Usage: ./akhal sampoke [FASTA FILE] [SAM FILE] [OUT SAM]\n");
        return -1;
    }

    if (argc == 5 && !ends_with(argv[4], ".sam")) {
        printf("[ERROR] Usage: ./akhal sampoke [FASTA FILE] [SAM FILE] [OUT SAM]\n");
        return -1;
    }

    char *out_file = NULL;
    if (argc == 5) out_file = argv[4];

    struct ref_seq seqs;
    read_fasta(argv[2], &seqs);
    check_sam(argv[3], out_file, &seqs);

    free_ref_seq(&seqs);
    return 0;
}
