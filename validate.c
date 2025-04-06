#include "validate.h"

KHASHL_MAP_INIT(KH_LOCAL, val_map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

int validate_parse_gfa(const char* file_path, segment **segments, int *size, int rGFA, val_map32_t *h) {

    size_t line_cap = 1048576;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;
    
    int segment_size = 0, segment_capacity = 1000000;
    segment *temp_segments = (segment *)malloc(segment_capacity * sizeof(segment));

    int main_segment_count = 0;
    int main_segment_count_rgfa = 0;

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
            if (strlen(token) == 0) {
                fprintf(stderr, "[ERROR] (S) Sequence %lu is empty\n", temp_segments[segment_size].id);
                temp_segments[segment_size].seq = NULL;
            } else {
                temp_segments[segment_size].seq = malloc(strlen(token) + 1);
                if (!temp_segments[segment_size].seq) {
                    fprintf(stderr, "[ERROR] Memory allocation failed\n");
                    exit(EXIT_FAILURE);
                }
                strcpy(temp_segments[segment_size].seq, token);
            }
            
            temp_segments[segment_size].rank = 1; // assume rank is 1 for all segments
            temp_segments[segment_size].next = NULL;

            // INFO
            token = strtok_r(NULL, "\t", &saveptr); // process rest if any
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
                        // temp_segments[segment_size].rank = atoi(value); // rank will be set based on Path
                        if (atoi(value) == 0) main_segment_count_rgfa++;
                    } 
                }

                token = strtok_r(NULL, "\t", &saveptr); // move to the next field
            }
            
            khint_t k; int absent;
            k = map32_put(h, temp_segments[segment_size].id, &absent);
            if (absent) kh_val(h, k) = segment_size; // set value as index in segments
            else fprintf(stderr, "[ERROR] (S) Duplicate id for sequence %lu\n", temp_segments[segment_size].id);
            segment_size++;
        } else if (line[0] == 'L') {
            uint64_t id1, id2;
            char strand1, strand2;
            size_t overlap = 0;
            sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM", &id1, &strand1, &id2, &strand2, &overlap);
            
            khint_t k = map32_get(h, id1); // query the hash table
            if (k < kh_end(h)) {
                uint32_t seg1_index = kh_val(h, k);
                k = map32_get(h, id2); // query the hash table
                if (k < kh_end(h)) {
                    uint32_t seg2_index = kh_val(h, k);
                    if (0 < overlap) {
                        if (overlap < strlen(temp_segments[seg1_index].seq) && overlap < strlen(temp_segments[seg2_index].seq)) {
                            if (strncmp(temp_segments[seg1_index].seq, temp_segments[seg2_index].seq + (strlen(temp_segments[seg2_index].seq) - overlap), overlap) != 0) {
                                fprintf(stderr, "[ERROR] (L) Overlap mismatch %lu -> %lu with length %lu\n", id1, id2, overlap);
                            }
                        } else {
                            fprintf(stderr, "[ERROR] (L) Overlap mismatch %lu -> %lu with length %lu\n", id1, id2, overlap);
                        }
                    }
                } 
                else fprintf(stderr, "[ERROR] (L-2) Segment %lu not found\n", id2);
            }
            else fprintf(stderr, "[ERROR] (L-1) Segment %lu not found\n", id1);
        } else if (line[0] == 'P') {
            char *path_id;
            char *save_ptr;
            char *token = strtok_r(line, "\t", &save_ptr);
            
            token = strtok_r(NULL, "\t", &save_ptr); // skip path ID
            path_id = strdup(token);
            
            token = strtok_r(NULL, "\t", &save_ptr); // process segment IDs
            if (!token) {
                fprintf(stderr, "[ERROR] (P) Missing segment list in path %s\n", path_id);
                continue;
            }

            char *segment_ptr;
            char *segment_token = strtok_r(token, ",", &segment_ptr);
            segment *prev_segment = NULL;

            while (segment_token) {
                main_segment_count++;
                size_t len = strlen(segment_token);
                if (segment_token[len - 1] == '+') segment_token[len - 1] = '\0'; // remove the trailing '+'

                uint64_t seg_id = strtoull(segment_token, NULL, 10);
                khint_t k = map32_get(h, seg_id);
                if (k == kh_end(h)) fprintf(stderr, "[ERROR] (P) Segment %lu not found in path %s\n", seg_id, path_id);
                else {
                    segment *current_segment = temp_segments + kh_val(h, k);
                    current_segment->rank = 0;
                    if (prev_segment) {
                        prev_segment->next = current_segment;     
                    }               
                    prev_segment = current_segment; // update previous
                }

                segment_token = strtok_r(NULL, ",", &segment_ptr); // next segment
            }
            if (prev_segment)
                prev_segment->next = NULL;

            free(path_id);
        }
    }

    rewind(file);

    // Now validate the path and links (whether there are links are missing)
    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == 'L') {
            uint64_t id1, id2;
            char strand1, strand2;
            size_t overlap = 0;
            sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM", &id1, &strand1, &id2, &strand2, &overlap);
            
            khint_t k = map32_get(h, id1); // query the hash table
            if (k < kh_end(h)) {
                uint32_t seg1_index = kh_val(h, k);
                if (temp_segments[seg1_index].rank == 0) {
                    k = map32_get(h, id2); // query the hash table
                    if (k < kh_end(h)) {
                        uint32_t seg2_index = kh_val(h, k);
                        // If segment present in Path and not the last one
                        if (temp_segments[seg2_index].rank == 0) {
                            if (temp_segments[seg1_index].next != NULL) { 
                                if (temp_segments[seg1_index].next->id == temp_segments[seg2_index].id) 
                                    temp_segments[seg1_index].next = NULL;
                            }
                        }
                    } 
                }
            }
        }
    }
    io_close(file, &line);

    for (int i=0; i<segment_size; i++) {
        if (temp_segments[i].rank == 0 && temp_segments[i].next != NULL)
            fprintf(stderr, "[ERROR] (V) No link found %lu -> %lu which is presented in path\n", temp_segments[i].id, temp_segments[i].next->id);
    }

    *segments = temp_segments;
    *size = segment_size;

    if (rGFA && main_segment_count != main_segment_count_rgfa)
        fprintf(stderr, "[ERROR] (P) Number of rank-0 segments mismatch in segment info and path\n");

    return 1;
}

int validate(int argc, char* argv[]) {

    if (argc < 3) {
        fprintf(stderr, "[ERROR] Usage: ./akhal parse [r/GFA]\n");
        exit(EXIT_FAILURE);
    }

    val_map32_t *h = map32_init();
    segment *segments;
    int size = 0;
    int rGFA = ends_with(argv[2], ".rgfa");

    if (validate_parse_gfa(argv[2], &segments, &size, rGFA, h)) {
        printf("[INFO] Parsed %s successfully\n", rGFA ? "rGFA" : "GFA");
    }

    free_segments(&segments, size);
    map32_destroy(h);

    return 0;
}