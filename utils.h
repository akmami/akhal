#ifndef __UTILS_H__
#define __UTILS_h__

#include <math.h>
#include <string.h>
#include "struct_def.h"

#define ends_with(str, suffix) ((strlen(str) < strlen(suffix)) ? 0 : (strcmp(str+strlen(str)-strlen(suffix), suffix) == 0))

void free_segments(segment **segments, int segment_size);
void print_segments(segment *segments, int segment_count);
void print_ref(segment *segments, int segment_count);
void write_sam_hdr(FILE *file, char **paths, int path_size, uint64_t *lens);
double calculate_mean(size_t *arr, size_t s);
double calculate_variance(size_t *arr, size_t s, double mean);
double calculate_std_dev(double variance);
void find_in_degrees(segment *segments, int segment_count, int *min_degree, int *max_degree);
void find_out_degrees(segment *segments, int segment_count, int *min_degree, int *max_degree);  

#endif