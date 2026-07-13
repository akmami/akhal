#ifndef __UTILS_H__
#define __UTILS_h__

#include "struct_def.h"
#include <math.h>
#include <string.h>
#include <ctype.h>


#define ends_with(str, suffix) ((strlen(str) < strlen(suffix)) ? 0 : (strcmp(str+strlen(str)-strlen(suffix), suffix) == 0))

#define complement(base) ((base & 0xDF) == 'A' ? 'T' : \
                          (base & 0xDF) == 'T' ? 'A' : \
                          (base & 0xDF) == 'C' ? 'G' : \
                          (base & 0xDF) == 'G' ? 'C' : 'N')

#define abs(a) (a < 0 ? -a : a)

void free_segments(segment **segments, int segment_size);
void print_segments(segment *segments, int segment_count);
void free_ref_seq(struct ref_seq *seqs);
void print_ref(segment *segments, int segment_count);
void write_sam_hdr(FILE *file, char **paths, int path_size, uint64_t *lens);
double calculate_mean(size_t *arr, size_t s);
double calculate_variance(size_t *arr, size_t s, double mean);
double calculate_std_dev(double variance);
void find_in_degrees(segment *segments, int segment_count, int *min_degree, int *max_degree);
void find_out_degrees(segment *segments, int segment_count, int *min_degree, int *max_degree);  
int parse_cigar(char *cigar, char *ops, int max_ops, int rev);
void reverseComplement(char *sequence, int length);
char *parse_rg_prefix(const char *header);
int next_node(const char *path, uint64_t *id, char *strand);

#endif