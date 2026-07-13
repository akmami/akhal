#include "sampoke.h"

void read_fasta(const char *fasta_path, struct ref_seq *seqs) {

    char *fai_path = malloc(strlen(fasta_path)+5);
    if (fai_path == NULL) {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    sprintf(fai_path, "%s.fai", fasta_path);

    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *fasta_fai_file = io_open(fai_path, &line, line_cap);
    int line_len = 0;
    int chrom_index = 0;

    while ((line_len = io_read(fasta_fai_file, &line, &line_cap))) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "[ERROR] Index file is empty.\n");
        exit(EXIT_FAILURE);
    }

    seqs->size = chrom_index;
    seqs->chrs = (struct chr *)malloc(chrom_index*sizeof(struct chr));
    if (seqs->chrs == NULL) {
        fprintf(stderr, "[ERROR] Couldn't allocate memory to ref sequences\n");
        exit(EXIT_FAILURE);
    }

    rewind(fasta_fai_file);
    chrom_index = 0;

    while ((line_len = io_read(fasta_fai_file, &line, &line_cap))) {
        char *name, *length;

        // assign name
        char *saveptr;
        name = strtok_r(line, "\t", &saveptr);
        uint64_t name_len = strlen(name);
        seqs->chrs[chrom_index].seq_name = (char *)malloc(name_len+1);
        memcpy(seqs->chrs[chrom_index].seq_name, name, name_len);
        seqs->chrs[chrom_index].seq_name[name_len] = '\0';

        // assign size and allocate in memory
        length = strtok_r(NULL, "\t", &saveptr);
        seqs->chrs[chrom_index].seq_size = strtol(length, NULL, 10);

        seqs->chrs[chrom_index].seq = (char *)malloc(seqs->chrs[chrom_index].seq_size+1);
        if (seqs->chrs[chrom_index].seq == NULL) {
            fprintf(stderr, "[ERROR] Couldn't allocate memory to chromosome string.\n");
            exit(EXIT_FAILURE);
        }
        seqs->chrs[chrom_index].seq[seqs->chrs[chrom_index].seq_size] = '\0';
        chrom_index++;
    }

    io_close(fasta_fai_file, &line);

    // Initialize chromosome array
    FILE *fasta_file = io_open(fasta_path, &line, line_cap);
    chrom_index = 0;
    uint64_t seq_size = 0;

    while ((line_len = io_read(fasta_file, &line, &line_cap))) {
        if (line[0] == '>') {
            if (seq_size != 0) {
                seq_size = 0;
                chrom_index++;
            }
        } else {
            memcpy(seqs->chrs[chrom_index].seq + seq_size, line, line_len);
            seq_size += line_len;
        }
    }

    io_close(fasta_file, &line);
}

int get_rg_headers(const char *file_path, char ***rg_headers) {
    
    int rg_headers_size = 0, rg_headers_capacity = 128;
    char **temp_rg_headers = (char **)malloc(sizeof(char *) * rg_headers_capacity);
    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;

    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] != '@') {
            char *prefix = parse_rg_prefix(line);
            if (prefix != NULL) {
                int exist = 0;
                for (int i=0; i<rg_headers_size; i++) {
                    if (strcmp(prefix, temp_rg_headers[i]) == 0) {
                        exist = 1;
                        break;
                    }
                }
                if (!exist) {
                    if (rg_headers_size == rg_headers_capacity) {
                        rg_headers_capacity *= 2;
                        char **temp = realloc(temp_rg_headers, rg_headers_capacity * sizeof(char *));
                        temp_rg_headers = temp;
                        temp_rg_headers[rg_headers_size++] = prefix;
                    }   
                    else temp_rg_headers[rg_headers_size++] = prefix;
                } else {
                    free(prefix);
                }
            }
        } else if (line[1] == 'R' && line[2] == 'G') {
            break; // no need to find RG tags, already exist
        }
    }

    *rg_headers = temp_rg_headers;
    io_close(file, &line);
    return rg_headers_size;
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

void check_sam(const char *sam_file, const char *out_file, const struct ref_seq *seqs, char **rg_headers, int rg_headers_size) {
    
    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(sam_file, &line, line_cap);
    int line_len = 0;

    FILE *outfile = NULL;
    if (out_file != NULL) {
        outfile = fopen(out_file, "w");
        if (!outfile) {
            fprintf(stderr, "[ERROR] Failed to open output SAM file\n");
            fclose(file);
            exit(1);
        }
    }

    int line_no = 0;
    int invalid_count = 0, correct_count = 0;
    int print_rg = 1;
    while ((line_len = io_read(file, &line, &line_cap))) {
        line_no++;
        if (line[0] == '@') { 
            if (outfile != NULL) {
                if (line[1] == 'P' && line[2] == 'G' && rg_headers_size && print_rg) {
                    for (int i = 0; i < rg_headers_size; i++) {
                        fprintf(outfile, "@RG\tID:%s.%d\tPL:%s\tPU:%s\tSM:%s\n", "akhal", i, "PACBIO", rg_headers[i], "sample");
                    }
                    print_rg = 0;
                }
                fprintf(outfile, "%s\n", line);
            }
            continue; 
        } else if (outfile != NULL && rg_headers_size && print_rg) {
            for (int i = 0; i < rg_headers_size; i++) {
                fprintf(outfile, "@RG\tID:%s.%d\tPL:%s\tPU:%s\tSM:%s\n", "akhal", i, "PACBIO", rg_headers[i], "sample");
            }
            print_rg = 0;
        }

        char qname[128], rname[128], cigar[MAX_CIGAR], seq[MAX_LINE];
        int flag, pos, mapq, pnext, tlen;
        char rnext[128], qual[MAX_LINE];
        int fields = sscanf(line, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
                            qname, &flag, rname, &pos, &mapq, cigar,
                            rnext, &pnext, &tlen, seq, qual);
        
        if (fields < 11) {
            fprintf(stderr, "Line %d, id: %s: Malformed SAM line\n", line_no, qname);
            invalid_count++;
            continue;
        }

        if (strcmp(rname, "*") == 0) {
            if (outfile != NULL) fprintf(outfile, "%s\n", line);
            continue; // unmapped read
        }

        int chrom_index = -1;
        for (int i=0; i<seqs->size; i++) {
            if (strcmp(rname, seqs->chrs[i].seq_name) == 0) { chrom_index = i; break; }
        }
        if (chrom_index == -1) {
            fprintf(stderr, "Line %d, id: %s: Reference not found: %s\n", line_no, qname, rname);
            invalid_count++;
            continue;
        }

        if (!validate_alignment(seq, seqs->chrs[chrom_index].seq, pos, cigar)) {
            invalid_count++;
        } else { 
            correct_count++; 
            if (rg_headers_size) {
                char *prefix = parse_rg_prefix(line);
                if (prefix != NULL) {
                    for (int i=0; i<rg_headers_size; i++) {
                        if (strcmp(prefix, rg_headers[i]) == 0) {
                            fprintf(outfile, "%s\tRG:Z:akhal.%d\n", line, i);
                        }
                    }
                } else {
                    fprintf(outfile, "%s\n", line);
                }
                continue;
            }
            if (outfile != NULL) 
                fprintf(outfile, "%s\n", line);
        }
    }

    printf("[INFO] correct: %d, incorrect: %d\n", correct_count, invalid_count);
    io_close(file, &line);

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

    char **rg_headers;
    int rg_headers_size = get_rg_headers(argv[3], &rg_headers);

    struct ref_seq seqs;
    read_fasta(argv[2], &seqs);
    check_sam(argv[3], out_file, &seqs, rg_headers, rg_headers_size);

    free_ref_seq(&seqs);
    for (int i = 0; i < rg_headers_size; i++) free(rg_headers[i]);
    free(rg_headers);

    return 0;
}
