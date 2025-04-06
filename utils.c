#include "utils.h"

void free_segments(segment **segments, int segment_size) {
    segment *temp = *segments;
    for (int i = 0; i < segment_size; i++) {
        free(temp[i].seq);
    }
    free(*segments);
}

void print_segments(segment *segments, int segment_count) {
    for (int i = 0; i < segment_count; i++) {
        printf("Segment ID: %lu, String: %s, Next: %lu\n", 
               segments[i].id, segments[i].seq, 
               segments[i].next != NULL ? segments[i].next->id : 0);
    }
}

void print_ref(segment *segments, int segment_count) {
    for (int i = 0; i < segment_count; i++) {
        if (segments[i].next != NULL) {
            printf("Segment ID: %lu, String: %s, Next: %lu\n", 
                segments[i].id, segments[i].seq, segments[i].next->id);
        }
    }
}

void write_sam_hdr(FILE *file, char **paths, int path_size, uint64_t *lens) {

    fprintf(file, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");

    for (int i=0; i<path_size; i++) {
        fprintf(file, "@SQ\tSN:%s\tLN:%lu\n", paths[i], lens[i]);
    }

    fprintf(file, "@PG\tID:LCPan\tPN:LCPan\tNV:1.0\n");
    fprintf(file, "@PG\tID:GraphAligner\tPN:GraphAligner\n");

    // fprintf(file, "@SQ\tSN:%s\tLN:%d\tAH:%s:%d-%d\n", id, len, chr, start, end); // for large ins or alt
}

double calculate_mean(size_t *arr, size_t s) {
    if (s == 0) return 0.0;

    size_t sum = 0;
    for (size_t i = 0; i < s; i++) {
        sum += arr[i];
    }
    return (double)sum / s;
}

double calculate_variance(size_t *arr, size_t s, double mean) {
    if (s == 0) return 0.0;

    double sum_sq_diff = 0.0;
    for (size_t i = 0; i < s; i++) {
        double diff = (double)arr[i] - mean;
        sum_sq_diff += diff * diff;
    }
    return sum_sq_diff / s;
}

double calculate_std_dev(double variance) {
    return sqrt(variance);
}

void find_in_degrees(segment *segments, int segment_count, int *min_degree, int *max_degree) {
    *min_degree = -1;
    *max_degree = -1;
    if (segment_count == 0) {
        return;
    }
    
    int min = INT32_MAX;
    int max = 0;
    for (int i = 1; i < segment_count; i++) {
        if (segments[i].in_degree)
            min = min < segments[i].in_degree ? min : segments[i].in_degree;
        if (segments[i].in_degree)
            max = max > segments[i].in_degree ? max : segments[i].in_degree;
    }

    *min_degree = min;
    *max_degree = max;
}

void find_out_degrees(segment *segments, int segment_count, int *min_degree, int *max_degree) {
    *min_degree = -1;
    *max_degree = -1;
    if (segment_count == 0) {
        return;
    }
    
    int min = INT32_MAX;
    int max = 0;
    for (int i = 0; i < segment_count; i++) {
        if (segments[i].out_degree)
            min = min < segments[i].out_degree ? min : segments[i].out_degree;
        if (segments[i].out_degree)
            max = max > segments[i].out_degree ? max : segments[i].out_degree;
    }

    *min_degree = min;
    *max_degree = max;
}