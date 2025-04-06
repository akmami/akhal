#include "gaf2sam.h"

// KHASHL_MAP_INIT(SCOPE, HType, prefix, khkey_t, kh_val_t, __hash_fn, __hash_eq)
KHASHL_MAP_INIT(KH_LOCAL, map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

inline static int next_cigar(const char *cigar, char *op, int *size) {
    if (*cigar != '\0') {
        int index = 0, length = 0;
        while (*cigar >= '0' && *cigar <= '9') {
            length = length * 10 + (*cigar - '0');
            cigar++; index++;
        }
        *op = *cigar; *size = length;
        return index+1;
    }
    return 0;
}

inline static int next_node(const char *path, uint64_t *id, char *strand) {
    if (*path != '\0') {
        *strand = *path;
        path++;
        int index = 1; uint64_t i = 0;
        while (*path >= '0' && *path <= '9') {
            i = i * 10 + (*path - '0');
            path++; index++;
        }
        *id = i;
        return index;
    }
    return 0;
}

void print_alignment(const alignment *aln) {
    printf("Read Name: %s\n", aln->readName);
    printf("Read Length: %d\n", aln->readLen);
    printf("Read Start: %d\n", aln->readStart);
    printf("Read End: %d\n", aln->readEnd);
    printf("Strand: %c\n", aln->strand);
    printf("Path: %s\n", aln->path);
    printf("Path Length: %d\n", aln->pathLen);
    printf("Path Start: %d\n", aln->pathStart);
    printf("Path End: %d\n", aln->pathEnd);
    printf("Matches: %d\n", aln->matches);
    printf("Block Length: %d\n", aln->blockLen);
    printf("Quality: %d\n", aln->qual);
    printf("XDI: %d\n", aln->xdi);
    printf("Score: %.2f\n", aln->score);
    printf("Divergence: %.6f\n", aln->divergence);
    printf("Identity: %.6f\n", aln->identity);
    printf("CIGAR: %s\n", aln->cigar ? aln->cigar : "N/A");

    int index = 0;
    const char *cigar = aln->cigar;
    char op;
    int length;
    while ((index = next_cigar(cigar, &op, &length))) {
        printf("Operation: %c, Length: %d\n", op, length);
        cigar += index;
    }
}

void write_sam_record(FILE *out_sam, alignment *aln, char *ops, int c_size, const char *rname, int pos) {
    int flag = 0;
    if (aln->strand == '<') flag |= 0x10;  // reverse complemented

    // Compose the CIGAR string from `ops`
    char cigar_string[65536];
    int cigar_pos = 0;

    if (c_size == 0) {
        strcpy(cigar_string, "*");
    } else {
        int i = 0;
        while (i < c_size && (ops[i] == CIGAR_INSERTION || ops[i] == CIGAR_SEQUENCE_MISMATCH)) {
            ops[i++] = CIGAR_SOFT_CLIP;
        }

        int j = c_size - 1;
        while (j >= 0 && (ops[j] == CIGAR_INSERTION || ops[j] == CIGAR_SEQUENCE_MISMATCH)) {
            ops[j--] = CIGAR_SOFT_CLIP;
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
    const char *seq = "*";
    const char *qual = "*";

    // POS is 1-based in SAM, so add 1
    int mapq = aln->qual;  // mapping quality (can be 255)

    fprintf(out_sam,
            "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t"
            "NM:i:%d\tAS:f:%.2f\tdv:f:%.6f\tid:f:%.6f\n",
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
            aln->identity
    );
}

int gaf2sam_parse_gaf(const char* file_path, segment *segments, map32_t *h, FILE *out_sam) {

    size_t line_cap = 1048576;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;

    int fwd_aln_count = 0, rc_aln_count = 0, inv_aln_count = 0;
    int rank_0_aln_count = 0, rank_1_aln_count = 0, rank_both_aln_count = 0, rank_inv_aln_count = 0;
    
    while ((line_len = io_read(file, &line, &line_cap))) {

        alignment aln;
        memset(&aln, 0, sizeof(alignment));

        char *saveptr;
        char *token;
        int field = 0;

        token = strtok_r(line, "\t", &saveptr);
        while (token != NULL) {
            switch (field) {
                case 0: aln.readName = strdup(token); break;
                case 1: aln.readLen = atoi(token); break;
                case 2: aln.readStart = atoi(token); break;
                case 3: aln.readEnd = atoi(token); break;
                case 4: aln.strand = token[0]; break;
                case 5: aln.path = strdup(token); break;
                case 6: aln.pathLen = atoi(token); break;
                case 7: aln.pathStart = atoi(token); break;
                case 8: aln.pathEnd = atoi(token); break;
                case 9: aln.matches = atoi(token); break;
                case 10: aln.blockLen = atoi(token); break;
                case 11: aln.qual = atoi(token); break;
                default:
                    if (strncmp(token, "NM:i:", 5) == 0) {
                        aln.xdi = atoi(token + 5);
                    } else if (strncmp(token, "AS:f:", 5) == 0) {
                        aln.score = atof(token + 5);
                    } else if (strncmp(token, "dv:f:", 5) == 0) {
                        aln.divergence = atof(token + 5);
                    } else if (strncmp(token, "id:f:", 5) == 0) {
                        aln.identity = atof(token + 5);
                    } else if (strncmp(token, "cg:Z:", 5) == 0) {
                        aln.cigar = strdup(token + 5);
                    }
                    break;
            }
            field++;
            token = strtok_r(NULL, "\t", &saveptr);
        }

        const char *cigar = aln.cigar; char op; int c_length;
        int c_index;
        const char *path = aln.path; uint64_t id = 0; char strand = 'N';
        int p_index = next_node(path, &id, &strand);
        path += p_index;
        khint_t k = map32_get(h, id);
        if (k == kh_end(h)) fprintf(stderr, "[ERROR] Segment not found %lu\n", id);
        segment *seg = segments + kh_val(h, k);
        int p_length = strlen(seg->seq) - aln.pathStart;

        if (strand == '>') fwd_aln_count++;
        else if (strand == '<') rc_aln_count++;
        else inv_aln_count++;

        char ops[131072]; // 2^17
        memset(ops, 0, 131072); // note need to do
        
        int c_size = 0;
        int is_alt = 0, is_ref = 0;

        // soft clip for the prefix part of the read
        for(int i=0; i<aln.readStart; i++) {
            ops[c_size++] = CIGAR_SOFT_CLIP;
        }

        int reference_start = seg->start + aln.pathStart + 1;
        int isRefName = 1;
        char *ref_name = NULL;

        int found_ref_anchor = 0;

        if (seg->rank == 0) {
            if (seg->ref_name == NULL) fprintf(stderr, "[ERROR] Reference name is not set for segment (rank:0)\n");
            else ref_name = seg->ref_name;
            isRefName = 0;
            found_ref_anchor = 1;
        }

        while ((c_index = next_cigar(cigar, &op, &c_length)) && p_index) {
            cigar += c_index;
            while (c_length--) {
                if (seg->rank == 0) {
                    ops[c_size++] = op;
                    is_ref = 1;
                } else {
                    if (CIGAR_QUE(op)) ops[c_size++] = CIGAR_INSERTION;
                    is_alt = 1;
                }

                if (--p_length < 0) {
                    p_index = next_node(path, &id, &strand);
                    if (!p_index) break;
                    path += p_index;
                    khint_t k = map32_get(h, id);
                    if (k == kh_end(h)) fprintf(stderr, "[ERROR] Segment not found %lu\n", id);
                    seg = segments + kh_val(h, k);
                    p_length = strlen(seg->seq);

                    if (isRefName && seg->rank == 0) {
                        if (seg->ref_name == NULL) fprintf(stderr, "[ERROR] Reference name is not set for segment (rank:0)\n");
                        else ref_name = seg->ref_name;
                        isRefName = 0;
                    }
                }

                if (!found_ref_anchor && seg->rank == 0) {
                    if (CIGAR_REF(op)) {
                        if (op == CIGAR_SEQUENCE_MATCH) {
                            found_ref_anchor = 1;
                        } else {
                            reference_start++;
                        }
                    }
                }
            }
        }

        for(int i=aln.readEnd; i<aln.readLen; i++) {
            ops[c_size++] = CIGAR_SOFT_CLIP;
        }
        
        if (is_alt && !is_ref) rank_0_aln_count++;
        else if (!is_alt && is_ref) rank_1_aln_count++;
        else if (is_alt && is_ref) rank_both_aln_count++;
        else rank_inv_aln_count++;

        // if (is_alt && is_ref) {
        //     for (int i=0; i<c_size; i++) {
        //         printf("%c", ops[i]);
        //     }
        //     printf("\n");
        //     printf("c_size: %d\n", c_size);
        //     printf("%s\n", aln.path);
        //     printf("%s\n", aln.cigar);

        //     path = aln.path;
        //     while ((p_index = next_node(path, &id, &strand))) {
        //         path += p_index;
        //         khint_t k = map32_get(h, id);
        //         if (k == kh_end(h)) fprintf(stderr, "[ERROR] Segment not found %lu\n", id);
        //         seg = segments + kh_val(h, k);
        //         printf("%lu-%d,", seg->id, seg->rank);
        //     }
        //     printf("\n\n");
        // }

        // for now, only print if there is any match mapping on reference backbone
        if (is_ref)
            write_sam_record(out_sam, &aln, ops, c_size, ref_name, reference_start);
        
        if (aln.readName) free(aln.readName);
        if (aln.path) free(aln.path);
        if (aln.cigar) free(aln.cigar);

    }    
    fclose(file);
    free(line);

    printf("[INFO] Stats:\n");
    printf("[INFO] # of reads aligned in forward strands: %d\n", fwd_aln_count);
    printf("[INFO] # of reads aligned in reverse strands: %d\n", rc_aln_count);
    printf("[INFO] # of reads aligned provided in invalid format: %d\n", inv_aln_count);

    printf("[INFO] # of alignments that only mapped to segments with rank 0: %d\n", rank_0_aln_count);
    printf("[INFO] # of alignments that only mapped to segments with rank 1: %d\n", rank_1_aln_count);
    printf("[INFO] # of alignments that mapped to segments with ranks 0 and 1: %d\n", rank_both_aln_count);
    printf("[INFO] # of alignments that mapped to segments with ranks 1<: %d\n", rank_inv_aln_count);

    return 1;
}

int gaf2sam_parse_gfa(const char* file_path, segment **segments, int *size, map32_t *h, FILE *sam_out, char ***paths_ptr, int *path_size_ptr, uint64_t **path_lens_ptr) {

    size_t line_cap = 1048576;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;
    
    int segment_size = 0, segment_capacity = 1000000;
    segment *temp_segments = (segment *)malloc(segment_capacity * sizeof(segment));
    int path_size = 0, path_capacity = 128;
    char **paths = (char**)malloc(path_capacity * sizeof(char *));
    uint64_t *path_lens = (uint64_t *)calloc(path_capacity, sizeof(uint64_t));
    if (paths == NULL || path_lens == NULL) {
        fprintf(stderr, "[ERROR] Malloc failed.\n");
        free(line); free_segments(&temp_segments, segment_size); map32_destroy(h);
        exit(EXIT_FAILURE);
    }

    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == 'S') {
            if (segment_size == segment_capacity) {
                segment_capacity = segment_capacity * 2;
                segment *temp = realloc(temp_segments, segment_capacity * sizeof(segment));
                if (temp == NULL) {
                    fprintf(stderr, "[ERROR] Realloc failed.\n");
                    free(line); free_segments(&temp_segments, segment_size); map32_destroy(h);
                    exit(EXIT_FAILURE);
                }
                temp_segments = temp;
            }
            
            char *saveptr;
            char *token = strtok_r(line, "\t", &saveptr); // Skip 'S'
            // ID
            token = strtok_r(NULL, "\t", &saveptr); // ID
            temp_segments[segment_size].id = strtoull(token, NULL, 10);
            
            // Sequence
            token = strtok_r(NULL, "\t", &saveptr); // Sequence string
            temp_segments[segment_size].seq = malloc(strlen(token) + 1);
            if (!temp_segments[segment_size].seq) {
                fprintf(stderr, "[ERROR] Memory allocation failed\n");
                exit(EXIT_FAILURE);
            }
            strcpy(temp_segments[segment_size].seq, token);
            temp_segments[segment_size].rank = 1;
            temp_segments[segment_size].next = NULL;
            temp_segments[segment_size].ref_name = NULL;

            // INFO
            token = strtok_r(NULL, "\t", &saveptr); // process rest 
            while (token) {
                // process SN:Z:chr22	SO:i:10510000	SR:i:0
                char *info_saveptr;
                char *tag = strtok_r(token, ":", &info_saveptr);  // get TAG
                char *type = strtok_r(NULL, ":", &info_saveptr);  // get TYPE
                char *value = strtok_r(NULL, ":", &info_saveptr); // get VALUE

                if (tag && type && value) {
                    if (strcmp(tag, "SN") == 0 && strcmp(type, "Z") == 0) {
                        // temp_segments[segment_size].seq_name = strdup(value);
                    } else if (strcmp(tag, "SO") == 0 && strcmp(type, "i") == 0) {
                        temp_segments[segment_size].start = atoi(value);
                        temp_segments[segment_size].end = atoi(value) + strlen(temp_segments[segment_size].seq);
                    } else if (strcmp(tag, "SR") == 0 && strcmp(type, "i") == 0) {
                        temp_segments[segment_size].rank = atoi(value);
                    } 
                }

                token = strtok_r(NULL, "\t", &saveptr); // move to the next field
            }
            
            khint_t k; int absent;
            k = map32_put(h, temp_segments[segment_size].id, &absent);
            kh_val(h, k) = segment_size; // set value as index in segments
            segment_size++;
        } else if (line[0] == 'L') {
            uint64_t id1, id2;
            char strand1, strand2;
            size_t overlap = 0;
            sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM", &id1, &strand1, &id2, &strand2, &overlap);

            if (overlap) {
                fprintf(stderr, "[ERROR] Overlaps are not zero, cannot make conversion.\n");
                exit(EXIT_FAILURE);
            }
        } else if (line[0] == 'P') {
            char *token = strtok(line, "\t");
            token = strtok(NULL, "\t"); // skip path ID

            if (path_size == path_capacity) {
                path_capacity *= 2;
                char **temp1 = (char **)realloc(paths, path_capacity * sizeof(char *));
                uint64_t *temp2 = (uint64_t *)realloc(path_lens, path_capacity * sizeof(uint64_t));
                if (temp1 == NULL || temp2 == NULL) {
                    fprintf(stderr, "[ERROR] Realloc failed.\n");
                    free(line); free_segments(&temp_segments, segment_size); map32_destroy(h);
                    exit(EXIT_FAILURE);
                }
                paths = temp1;
                path_lens = temp2;
            }
            paths[path_size] = strdup(token);
            
            token = strtok(NULL, "\t"); // process segment IDs
            if (!token) {
                fprintf(stderr, "[ERROR] Missing segment list in path %s\n", paths[path_size-1]);
                path_size++;
                continue;
            }

            char *segment_token = strtok(token, ",");
            segment *prev_segment = NULL;
            int ref_pos = 0;

            while (segment_token) {
                size_t len = strlen(segment_token);
                if (segment_token[len - 1] == '+') segment_token[len - 1] = '\0'; // remove the trailing '+'

                uint64_t seg_id = strtoull(segment_token, NULL, 10);
                khint_t k = map32_get(h, seg_id);
                segment *current_segment = temp_segments + kh_val(h, k);
                current_segment->ref_name = paths[path_size];

                if (prev_segment) {
                    prev_segment->rank = 0;
                    prev_segment->next = current_segment;
                }
                if (current_segment->seq != NULL) {
                    path_lens[path_size] += strlen(current_segment->seq);
                    current_segment->start = ref_pos;
                    ref_pos += strlen(current_segment->seq);
                    current_segment->end = ref_pos; // for now, not needed.
                }
                prev_segment = current_segment; // update previous
                segment_token = strtok(NULL, ","); // next segment
            }
            if (prev_segment) {
                prev_segment->rank = 0;
            }

            path_size++;
        }
    }
    io_close(file, &line);

    *segments = temp_segments;
    *size = segment_size;

    write_sam_hdr(sam_out, paths, path_size, path_lens);

    *paths_ptr = paths;
    *path_size_ptr = path_size;
    *path_lens_ptr = path_lens;

    return 1;
}

int gaf2sam(int argc, char* argv[]) {
    if (argc < 5) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [OUTPUT file]\n");
        return -1;
    }

    if (!ends_with(argv[2], ".gfa") && !ends_with(argv[2], ".rgfa")) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [OUTPUT file]\n");
        return -1;
    }

    if (!ends_with(argv[3], ".gaf")) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [OUTPUT file]\n");
        return -1;
    }

    if (!ends_with(argv[4], ".sam")) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [OUTPUT file]\n");
        return -1;
    }

    FILE *sam = fopen(argv[4], "w");
    if (sam == NULL) {
        fprintf(stderr, "[ERROR] Couldn't open output file.\n");
        exit(EXIT_FAILURE);
    }

    int isGFA = ends_with(argv[2], ".gfa");


    map32_t *h = map32_init();
    segment *segments;
    int size = 0;

    int path_size = 0;
    char **paths = NULL;
    uint64_t *path_lens;

    if (gaf2sam_parse_gfa(argv[2], &segments, &size, h, sam, &paths, &path_size, &path_lens)) {
        printf("[INFO] Processed %s\n", isGFA ? "GFA" : "rGFA");
        if (gaf2sam_parse_gaf(argv[3], segments, h, sam)) {
            printf("[INFO] Processed GAF\n");
        }
    }

    free_segments(&segments, size);
    map32_destroy(h);
    fclose(sam);

    if (paths != NULL) {
        for (int i=0; i<path_size; i++) {
            free(paths[i]);
        }
        free(paths);
        free(path_lens);
    }

    return 0;
}