#include "analyze.h"

KHASHL_MAP_INIT(KH_LOCAL, stat_map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

int analize_parse_gfa(const char* file_path, segment **segments, int *size, stat_map32_t *h) {

    size_t line_cap = 1048576;
    char *line = NULL;
    FILE *file = io_open(file_path, &line, line_cap);
    int line_len = 0;
    
    int segment_size = 0, segment_capacity = 1000000;
    segment *temp_segments = (segment *)malloc(segment_capacity * sizeof(segment));
    size_t *segment_sizes = (size_t *)malloc(segment_capacity*sizeof(size_t));

    int overlap_size = 0, overlap_capacity = 1000000;
    size_t *overlap_sizes = (size_t *)malloc(overlap_capacity*sizeof(size_t));

    int main_segment_count = 0;

    while ((line_len = io_read(file, &line, &line_cap))) {
        if (line[0] == 'S') {
            if (segment_size == segment_capacity) {
                segment_capacity = segment_capacity * 2;
                segment *temp1 = (segment *)realloc(temp_segments, segment_capacity * sizeof(segment));
                if (temp1 == NULL) {
                    fprintf(stderr, "[ERROR] Realloc of segments array failed.\n");
                    free(line); free_segments(&temp_segments, segment_size); map32_destroy(h);
                    exit(EXIT_FAILURE);
                }
                temp_segments = temp1;

                size_t *temp2 = (size_t *)realloc(segment_sizes, segment_capacity * sizeof(size_t));
                if (temp2 == NULL) {
                    fprintf(stderr, "[ERROR] Realloc of segment sizes array failed.\n");
                    exit(EXIT_FAILURE);
                }
                segment_sizes = temp2;
            }
            
            char *saveptr;
            char *token = strtok_r(line, "\t", &saveptr); // Skip 'S'
            
            // ID
            token = strtok_r(NULL, "\t", &saveptr); // ID
            temp_segments[segment_size].id = strtoull(token, NULL, 10);

            // Sequence
            token = strtok_r(NULL, "\t", &saveptr); // Sequence string
            segment_sizes[segment_size] = strlen(token);
            temp_segments[segment_size].seq = strdup(token);
            
            temp_segments[segment_size].out_degree = 0;
            temp_segments[segment_size].in_degree = 0;
            
            khint_t k; int absent;
            k = map32_put(h, temp_segments[segment_size].id, &absent);
            kh_val(h, k) = segment_size; // set value as index in segments
            segment_size++;
        } else if (line[0] == 'L') {
            uint64_t id1, id2;
            char strand1, strand2;
            size_t overlap = 0;
            sscanf(line, "L\t%lu\t%c\t%lu\t%c\t%luM", &id1, &strand1, &id2, &strand2, &overlap);
            
            if (overlap_capacity == overlap_size) {
                overlap_capacity = overlap_capacity * 2;
                size_t *temp = (size_t *)realloc(overlap_sizes, overlap_capacity * sizeof(size_t));
                if (temp == NULL) {
                    fprintf(stderr, "[ERROR] Realloc of overlap array failed.\n");
                    exit(EXIT_FAILURE);
                }
                overlap_sizes = temp;
            }
            overlap_sizes[overlap_size++] = overlap;
            
            khint_t k = map32_get(h, id1); // query the hash table
            uint32_t seg1_index = kh_val(h, k);
            k = map32_get(h, id2); // query the hash table
            uint32_t seg2_index = kh_val(h, k);

            temp_segments[seg1_index].out_degree++;
            temp_segments[seg2_index].in_degree++;
        } else if (line[0] == 'P') {
            char *save_ptr;
            char *token = strtok_r(line, "\t", &save_ptr);
            token = strtok_r(NULL, "\t", &save_ptr); // skip path ID            
            token = strtok_r(NULL, "\t", &save_ptr); // process segment IDs
            char *segment_ptr;
            char *segment_token = strtok_r(token, ",", &segment_ptr);

            while (segment_token) {
                main_segment_count++;
                segment_token = strtok_r(NULL, ",", &segment_ptr); // next segment
            }
        }
    }
    io_close(file, &line);


    double segment_mean = calculate_mean(segment_sizes, segment_size);
    double segment_variance = calculate_variance(segment_sizes, segment_size, segment_mean);
    double segment_std_dev = calculate_std_dev(segment_variance);

    size_t min_segment_len = segment_sizes[0];

    for (int i = 1; i < segment_size; i++) {
        min_segment_len = min_segment_len < strlen(temp_segments[i].seq) ? min_segment_len : strlen(temp_segments[i].seq);
    }

    size_t max_segment_len = segment_sizes[0];
    for (int i = 1; i < segment_size; i++) {
        max_segment_len = max_segment_len > strlen(temp_segments[i].seq) ? max_segment_len : strlen(temp_segments[i].seq);
    }

    double overlap_mean = calculate_mean(overlap_sizes, overlap_size);
    double overlap_variance = calculate_variance(overlap_sizes, overlap_size, overlap_mean);
    double overlap_std_dev = calculate_std_dev(overlap_variance);

    int min_in_degree, max_in_degree;
    find_in_degrees(temp_segments, segment_size, &min_in_degree, &max_in_degree);
    int min_out_degree, max_out_degree;
    find_out_degrees(temp_segments, segment_size, &min_out_degree, &max_out_degree);

    printf("Segment count: %d\n", segment_size);
    printf("Rank 0 segment count: %d\n", main_segment_count);
    printf("Rank 0< segment count: %d\n", segment_size-main_segment_count);
    printf("Segment avg length: %f\n", segment_mean);
    printf("Segment std length: %f\n", segment_std_dev);
    printf("Segment min. length %lu\n", min_segment_len);
    printf("Segment max. length %lu\n", max_segment_len);
    printf("Link count: %d\n", overlap_size);
    printf("Link overlapping avg length: %f\n", overlap_mean);
    printf("Link overlapping std length: %f\n", overlap_std_dev);
    printf("Minimum in degree: %d\n", min_in_degree);
    printf("Maximum in degree: %d\n", max_in_degree);
    printf("Minimum out degree: %d\n", min_out_degree);
    printf("Maximum out degree: %d\n", max_out_degree);

    free(overlap_sizes);
    free(segment_sizes);

    *segments = temp_segments;
    *size = segment_size;

    return 1;
}

int analyze(int argc, char* argv[]) {

    if (argc < 3) {
        fprintf(stderr, "[ERROR] Usage: ./akhal analize [r/GFA]\n");
        exit(EXIT_FAILURE);
    }

    stat_map32_t *h = map32_init();
    segment *segments;
    int size = 0;

    analize_parse_gfa(argv[2], &segments, &size, h);

    free_segments(&segments, size);
    map32_destroy(h);

    return 0;
}