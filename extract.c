#include "validate.h"

KHASHL_MAP_INIT(KH_LOCAL, ext_map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

int extract_parse_gfa(const char* file_path, const char *fa_out, segment **segments, int *size, ext_map32_t *h) {

    size_t line_cap = 1048576;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;

    FILE *out = fopen(fa_out, "w");
    if (out == NULL) {
        fprintf(stderr, "[ERROR] Output file couldn't be opened.\n");
        exit(EXIT_FAILURE);
    }
    
    int segment_size = 0, segment_capacity = 1000000;
    segment *temp_segments = (segment *)malloc(segment_capacity * sizeof(segment));

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
                fprintf(stderr, "Error: Sequence %lu is empty\n", temp_segments[segment_size].id);
                temp_segments[segment_size].seq = NULL;
            } else {
                temp_segments[segment_size].seq = malloc(strlen(token) + 1);
                if (!temp_segments[segment_size].seq) {
                    fprintf(stderr, "[ERROR] Memory allocation failed\n");
                    exit(EXIT_FAILURE);
                }
                strcpy(temp_segments[segment_size].seq, token);
            }
            
            temp_segments[segment_size].next = NULL;
            
            khint_t k; int absent;
            k = map32_put(h, temp_segments[segment_size].id, &absent);
            if (absent) kh_val(h, k) = segment_size; // set value as index in segments
            else fprintf(stderr, "[ERROR] Duplicate id for sequence %lu.\n", temp_segments[segment_size].id);
            segment_size++;
        } else if (line[0] == 'P') {
            char *path_id;
            char *save_ptr;
            char *token = strtok_r(line, "\t", &save_ptr);
            
            token = strtok_r(NULL, "\t", &save_ptr); // skip path ID
            path_id = strdup(token);

            fprintf(out, ">%s\n", path_id);
            
            token = strtok_r(NULL, "\t", &save_ptr); // process segment IDs
            if (!token) {
                fprintf(stderr, "[ERROR] Missing segment list in path %s\n", path_id);
                continue;
            }

            char *segment_ptr;
            char *segment_token = strtok_r(token, ",", &segment_ptr);

            int chars_in_line = 0;

            while (segment_token) {
                size_t len = strlen(segment_token);
                if (segment_token[len - 1] == '+') segment_token[len - 1] = '\0'; // remove the trailing '+'

                uint64_t seg_id = strtoull(segment_token, NULL, 10);
                khint_t k = map32_get(h, seg_id);
                if (k == kh_end(h)) fprintf(stderr, "[ERROR] Segment %lu not found in path %s\n", seg_id, path_id);
                
                segment *current_segment = temp_segments + kh_val(h, k);
                if (current_segment->seq != NULL) {
                    const char *seq = current_segment->seq;
                    size_t seq_len = strlen(seq);
                    size_t printed = 0;

                    while (printed < seq_len) {
                        size_t space_left = FASTA_WRAP_SIZE - chars_in_line;
                        size_t to_print = (seq_len - printed > space_left) ? space_left : (seq_len - printed);

                        fwrite(seq + printed, 1, to_print, out);
                        printed += to_print;
                        chars_in_line += to_print;

                        if (chars_in_line >= FASTA_WRAP_SIZE) {
                            fputc('\n', out);
                            chars_in_line = 0;
                        }
                    }
                }
                segment_token = strtok_r(NULL, ",", &segment_ptr); // next segment
            }

            if (chars_in_line)
                fputc('\n', out);
            
            chars_in_line = 0;
            free(path_id);
        }
    }
    io_close(file, &line);
    fclose(out);

    *segments = temp_segments;
    *size = segment_size;

    return 1;
}

int extract(int argc, char* argv[]) {

    if (argc < 5) {
        fprintf(stderr, "[ERROR] Usage: ./akhal extract fa [r/GFA] [Output FA]\n");
        exit(EXIT_FAILURE);
    }

    if (!ends_with(argv[3], ".gfa") && !ends_with(argv[3], ".rgfa")) {
        fprintf(stderr, "[ERROR] Usage: ./akhal extract fa [r/GFA] [Output FA]\n");
        exit(EXIT_FAILURE);
    }

    if (!ends_with(argv[4], ".fa") && !ends_with(argv[4], ".fasta")) {
        fprintf(stderr, "[ERROR] Usage: ./akhal extract fa [r/GFA] [Output FA]\n");
        exit(EXIT_FAILURE);
    }

    ext_map32_t *h = map32_init();
    segment *segments;
    int size = 0;

    if (strcmp(argv[2], "fa") == 0) {
        extract_parse_gfa(argv[3], argv[4], &segments, &size, h);
        printf("[INFO] Extracted reference successfully\n");
    }

    free_segments(&segments, size);
    map32_destroy(h);

    return 0;
}