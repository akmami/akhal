#include "gaf2sam.h"

// KHASHL_MAP_INIT(SCOPE, HType, prefix, khkey_t, kh_val_t, __hash_fn, __hash_eq)
KHASHL_MAP_INIT(KH_LOCAL, map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, strmap_t, strmap, const char*, char*, kh_hash_str, kh_eq_str)

void gaf2sam_parse_gaf(const char* file_path, segment *segments, map32_t *h1, strmap_t *h2, int read_count, char **rg_headers, int rg_headers_size, int simplify, FILE *out_sam) {

    int invalid_no_read = 0, invalid_segment = 0, invalid_mixed_strand = 0, invalid_rank = 0, invalid_no_qual = 0, invalid_other = 0;

    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;

    int fwd_aln_count = 0, rc_aln_count = 0, inv_aln_count = 0;
    int rank_0_aln_count = 0, rank_1_aln_count = 0, rank_both_aln_count = 0, rank_inv_aln_count = 0;

    while ((line_len = io_read(file, &line, &line_cap))) {

        alignment aln;
        memset(&aln, 0, sizeof(alignment));

        char *saveptr, *token;
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
                    if (strncmp(token, "NM:i:", 5) == 0) aln.xdi = atoi(token + 5);
                    else if (strncmp(token, "AS:f:", 5) == 0) aln.score = atof(token + 5);
                    else if (strncmp(token, "dv:f:", 5) == 0) aln.divergence = atof(token + 5);
                    else if (strncmp(token, "id:f:", 5) == 0) aln.identity = atof(token + 5);
                    else if (strncmp(token, "cg:Z:", 5) == 0) aln.cigar = strdup(token + 5);
                    break;
            }
            field++;
            token = strtok_r(NULL, "\t", &saveptr);
        }

        // Discard reads that has no quality
        if (aln.qual == 0) {
            invalid_no_qual++;
            if (aln.readName) free(aln.readName);
            if (aln.path) free(aln.path);
            if (aln.cigar) free(aln.cigar);
            continue;
        }
        
        // Get sequence from hash table that stores reads from fastq/fasta
        khint_t k = strmap_get(h2, aln.readName);
        if (k == kh_end(h2)) {
            fprintf(stderr, "[ERROR] Read %s does not exists\n", aln.readName);
            invalid_no_read++;
            if (aln.readName) free(aln.readName);
            if (aln.path) free(aln.path);
            if (aln.cigar) free(aln.cigar);
            continue;
        }
        aln.read = kh_val(h2, k);
        
        // Get nodes and store in array. It is needed in case if it is mapped as reverse complement
        char *ref_name = NULL; // to which chromosome it aligned
        uint64_t *nodes = (uint64_t *)malloc(sizeof(uint64_t) * MAX_CIGAR);
        int node_size = 0, path_index, fwd_strand_count = 0, rev_strand_count = 0;
        const char *path = aln.path; uint64_t id = 0; char strand = aln.cigar[0];
        int is_alt = 0, is_ref = 0, segment_not_found = 0;

        while ((path_index = next_node(path, &id, &strand))) {
            path += path_index;
            k = map32_get(h1, id);
            if (k == kh_end(h1)) {
                fprintf(stderr, "[ERROR] Segment %lu in read %s not found.\n", id, aln.readName); 
                invalid_segment++;
                if (aln.readName) free(aln.readName);
                if (aln.path) free(aln.path);
                if (aln.cigar) free(aln.cigar);
                free(nodes);
                segment_not_found = 1;
                break;
            }
            nodes[node_size++] = id;
            // stat 1
            if (strand == '>') fwd_strand_count++;
            else rev_strand_count++;
            // stat 2
            segment *temp = segments + kh_val(h1, k);
            if (temp->rank == 0) is_ref = 1;
            else if (temp->rank == 1) is_alt = 1;
            // Set ref_name
            if (temp->ref_name != NULL) ref_name = temp->ref_name;
        }
        // Continue if the wile loop terminated because of missing segment
        if (segment_not_found) continue;
        // Validation of valid strands
        if (fwd_strand_count && rev_strand_count) {
            fprintf(stderr, "[ERROR] Read %s aligned in mixed strands. %d %d\n", aln.readName, fwd_strand_count, rev_strand_count);
            invalid_mixed_strand++;
            if (aln.readName) free(aln.readName);
            if (aln.path) free(aln.path);
            if (aln.cigar) free(aln.cigar);
            free(nodes);
            continue;
        }
        // Stats of Segment Ranks
        if (is_ref && !is_alt) rank_0_aln_count++;
        else if (!is_ref && is_alt) rank_1_aln_count++;
        else if (is_alt && is_ref) rank_both_aln_count++;
        else {
            rank_inv_aln_count++;
            fprintf(stderr, "[ERROR] Read %s contains invalid/none ranks.\n", aln.readName);
            invalid_rank++;
            if (aln.readName) free(aln.readName);
            if (aln.path) free(aln.path);
            if (aln.cigar) free(aln.cigar);
            free(nodes);
            continue;
        }
        // Stats of Read alignment Strand
        if (fwd_strand_count) { fwd_aln_count++; aln.strand = '+'; }
        else { rc_aln_count++; aln.strand = '-'; }
        // Discard the reads that are mapped in alt paths only
        if (!is_ref) {
            if (aln.readName) free(aln.readName);
            if (aln.path) free(aln.path);
            if (aln.cigar) free(aln.cigar);
            free(nodes);
            continue;
        }

        // If read mapped as rc, reverse the order of the nodes
        if (rev_strand_count) {
            int i = 0;
            while (i < node_size / 2) {
                uint64_t temp = nodes[i];
                nodes[i] = nodes[node_size-i-1];
                nodes[node_size-i-1] = temp;
                i++;
            }
            int pathPreLen = aln.pathStart;
            int pathPostLen = aln.pathLen - aln.pathEnd;
            int pathLen = aln.pathLen - pathPreLen - pathPostLen;
            if (pathLen + pathPreLen + pathPostLen != aln.pathLen) {
                fprintf(stderr, "[ERROR] Read %s (path) has indices mismatch.\n", aln.readName);
                invalid_other++;
                if (aln.readName) free(aln.readName);
                if (aln.path) free(aln.path);
                if (aln.cigar) free(aln.cigar);
                free(nodes);
            }
            aln.pathStart = pathPostLen;
            aln.pathEnd = pathPostLen + pathLen;

            int readPreLen = aln.readStart;
            int readPostLen = aln.readLen - aln.readEnd;
            int readLen = aln.readLen - readPreLen - readPostLen;
            if (readLen + readPreLen + readPostLen != aln.readLen) {
                fprintf(stderr, "[ERROR] Read %s (read) has indices mismatch.\n", aln.readName);
                invalid_other++;
                if (aln.readName) free(aln.readName);
                if (aln.path) free(aln.path);
                if (aln.cigar) free(aln.cigar);
                free(nodes);
            }
            aln.readStart = readPostLen;
            aln.readEnd = readPostLen + readLen;

            aln.read = strdup(aln.read);
            reverseComplement(aln.read, strlen(aln.read));
        }

        // Initialize the segments
        int node_index = 0;
        k = map32_get(h1, nodes[node_index]);
        segment *seg = segments + kh_val(h1, k);
        segment *temp_prev_seg = (seg->rank == 0 ? seg : NULL);
        int p_length = strlen(seg->seq)-aln.pathStart;

        // Initialize the cigar strings
        char ops[MAX_CIGAR]; int op_index = 0;

        // Soft clip for the prefix part of the read
        for (int i = 0; i < aln.readStart; i++) ops[op_index++] = CIGAR_SOFT_CLIP;

        char cigar_ops[MAX_CIGAR];
        int n_ops = parse_cigar(aln.cigar, cigar_ops, MAX_CIGAR, rev_strand_count);
        if (n_ops < 0) {
            fprintf(stderr, "[ERROR] Unable to parse CIGAR %s in read %s\n", aln.cigar, aln.readName);
            invalid_other++;
            if (aln.readName) free(aln.readName);
            if (aln.path) free(aln.path);
            if (aln.cigar) free(aln.cigar);
            free(nodes);
            if (rev_strand_count) free(aln.read);
            continue;
        }

        int cigar_op_index = 0;
        int reference_start = seg->start+aln.pathStart+1;
        int is_reference_start_set = 0;

        while (cigar_op_index < n_ops) {
            char op = cigar_ops[cigar_op_index++];

            if (!is_reference_start_set && seg->rank == 0 && CIGAR_REF(op)) {
                if (op == CIGAR_SEQUENCE_MATCH) is_reference_start_set = 1;
                else reference_start++;
            }

            if (seg->rank == 0) ops[op_index++] = op;
            else if (CIGAR_QUE(op)) ops[op_index++] = CIGAR_INSERTION;

            if (CIGAR_REF(op)) {
                p_length--;

                if (p_length) continue;

                node_index++;
                if (node_index == node_size) break;
                khint_t k = map32_get(h1, nodes[node_index]);
                seg = segments + kh_val(h1, k);
                p_length = strlen(seg->seq);
                if (seg->rank == 0) {
                    if(!is_reference_start_set) reference_start = seg->start;

                    if (temp_prev_seg != NULL) {
                        for (int index = temp_prev_seg->end; index < seg->start; index++)
                            ops[op_index++] = CIGAR_DELETION;
                    }

                    temp_prev_seg = seg;
                }
            }
        }

        // Soft clip for the end part of the read
        for (int i = aln.readEnd; i < aln.readLen; i++) ops[op_index++] = CIGAR_SOFT_CLIP;

        // Validate cigar
        int cigar_query_count = 0;
        for (int i=0; i<op_index; i++) {
            if (CIGAR_QUE(ops[i])) cigar_query_count++;
        }
        if (cigar_query_count != aln.readLen) fprintf(stderr, "[ERROR] CIGAR (query) mismatch in length. %d %d\n", cigar_query_count, aln.readLen);

        write_sam_record(out_sam, &aln, ops, op_index, ref_name, reference_start, rg_headers, rg_headers_size, simplify);
        
        // Cleanup
        if (aln.readName) free(aln.readName);
        if (aln.path) free(aln.path);
        if (aln.cigar) free(aln.cigar);
        free(nodes);
        if (rev_strand_count) free(aln.read);
    }    
    fclose(file);
    free(line);

    printf("[INFO] Stats:\n");
    printf("[INFO] # of reads: %d\n", read_count);
    printf("[INFO] # of reads aligned in forward strands: %d\n", fwd_aln_count);
    printf("[INFO] # of reads aligned in reverse strands: %d\n", rc_aln_count);
    printf("[INFO] # of reads aligned provided in invalid format: %d\n", inv_aln_count);

    printf("[INFO] # of alignments that only mapped to segments with rank 0: %d\n", rank_0_aln_count);
    printf("[INFO] # of alignments that only mapped to segments with rank 1: %d\n", rank_1_aln_count);
    printf("[INFO] # of alignments that mapped to segments with ranks 0 and 1: %d\n", rank_both_aln_count);
    printf("[INFO] # of alignments that mapped to segments with ranks 1<: %d\n", rank_inv_aln_count);

    printf("[INFO] Stats: read/segment/strand/rank/qual/other: %d/%d/%d/%d/%d/%d\n", invalid_no_read, invalid_segment, invalid_mixed_strand, invalid_rank, invalid_no_qual, invalid_other);
}

int gaf2sam_parse_gfa(const char* file_path, segment **segments, int *size, map32_t *h1, FILE *sam_out, char ***paths_ptr, int *path_size_ptr, uint64_t **path_lens_ptr, char **rg_headers, int rg_headers_size) {

    size_t line_cap = MAX_LINE;
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
        free(line); free_segments(&temp_segments, segment_size); map32_destroy(h1);
        exit(EXIT_FAILURE);
    }

    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == 'S') {
            if (segment_size == segment_capacity) {
                segment_capacity = segment_capacity * 2;
                segment *temp = realloc(temp_segments, segment_capacity * sizeof(segment));
                if (temp == NULL) {
                    fprintf(stderr, "[ERROR] Realloc failed.\n");
                    free(line); free_segments(&temp_segments, segment_size); map32_destroy(h1);
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
            k = map32_put(h1, temp_segments[segment_size].id, &absent);
            kh_val(h1, k) = segment_size; // set value as index in segments
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
                    free(line); free_segments(&temp_segments, segment_size); map32_destroy(h1);
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
                khint_t k = map32_get(h1, seg_id);
                segment *current_segment = temp_segments + kh_val(h1, k);
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
                } else {
                    fprintf(stderr, "[ERROR] Segment seq is null\n");
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

    write_sam_hdr(sam_out, paths, path_size, path_lens, rg_headers, rg_headers_size);

    *paths_ptr = paths;
    *path_size_ptr = path_size;
    *path_lens_ptr = path_lens;

    return 1;
}

int gaf2sam_parse_fa(const char *file_path, strmap_t *h2, int *read_count, char ***rg_headers) {
    
    int rg_headers_size = 0, rg_headers_capacity = 128;
    char **temp_rg_headers = (char **)malloc(sizeof(char *) * rg_headers_capacity);
    size_t line_cap = MAX_LINE;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;
    int count = 0;
    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == '>') {
            int absent;
            khint_t k = strmap_put(h2, line + 1, &absent);
            if (absent) {
                kh_key(h2, k) = strdup(line + 1);
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
                line_len = io_read(file, &line, &line_cap);
                if (line_len) {
                    kh_val(h2, k) = strdup(line);
                } else {
                    kh_val(h2, k) = strdup("");
                }
            }
            count++;          
        }
    }

    *rg_headers = temp_rg_headers;

    io_close(file, &line);

    *read_count = count;
    return rg_headers_size;
}

int gaf2sam(int argc, char* argv[]) {
    if (argc < 6) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [FASTA FILE] [OUTPUT file]\n");
        return -1;
    }

    if (!ends_with(argv[2], ".gfa") && !ends_with(argv[2], ".rgfa")) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [FASTA FILE] [OUTPUT file]\n");
        return -1;
    }

    if (!ends_with(argv[3], ".gaf")) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [FASTA FILE] [OUTPUT file]\n");
        return -1;
    }

    if (!ends_with(argv[4], ".fa") && !ends_with(argv[4], ".fasta")) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [FASTA FILE] [OUTPUT file]\n");
        return -1;
    }

    if (!ends_with(argv[5], ".sam")) {
        printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [FASTA FILE] [OUTPUT file]\n");
        return -1;
    }

    FILE *sam = fopen(argv[5], "w");
    if (sam == NULL) {
        fprintf(stderr, "[ERROR] Couldn't open output file %s.\n", argv[5]);
        exit(EXIT_FAILURE);
    }

    int simplify = 0;

    if (argc == 7) {
        if (strcmp(argv[6], "--simple") == 0) {
            simplify = 1;
        } else {
            printf("[ERROR] Usage: ./akhal gam2sam [r/GFA file] [GAF file] [FASTA FILE] [OUTPUT file] [--simple]\n");
            return -1;
        }
    }

    int isGFA = ends_with(argv[2], ".gfa");

    map32_t *h1 = map32_init();
    segment *segments;
    int size = 0;

    strmap_t *h2 = strmap_init();

    int path_size = 0;
    char **paths = NULL;
    uint64_t *path_lens;
    int read_count;

    char **rg_headers;
    int rg_headers_size;

    rg_headers_size = gaf2sam_parse_fa(argv[4], h2, &read_count, &rg_headers);
    printf("[INFO] Processed FA\n");
    if (gaf2sam_parse_gfa(argv[2], &segments, &size, h1, sam, &paths, &path_size, &path_lens, rg_headers, rg_headers_size)) {
        printf("[INFO] Processed %s\n", isGFA ? "GFA" : "rGFA");
        gaf2sam_parse_gaf(argv[3], segments, h1, h2, read_count, rg_headers, rg_headers_size, simplify, sam);
        printf("[INFO] Processed GAF.\n");
    }

    free_segments(&segments, size);
    // cleanup h1
    map32_destroy(h1);
    // cleanup h2
    khint_t k;
    kh_foreach(h2, k) {
        free((char*)kh_val(h2, k));
        free((char*)kh_key(h2, k));
    }
    strmap_destroy(h2);
    fclose(sam);

    if (paths != NULL) {
        for (int i=0; i<path_size; i++) {
            free(paths[i]);
        }
        free(paths);
        free(path_lens);
    }

    for (int i = 0; i < rg_headers_size; i++) free(rg_headers[i]);
    free(rg_headers);

    return 0;
}