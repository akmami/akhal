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

void free_ref_seq(struct ref_seq *seqs) {
	if (seqs->size) {
		for (int i=0; i<seqs->size; i++) {
			free(seqs->chrs[i].seq_name);
			free(seqs->chrs[i].seq);
		}
		free(seqs->chrs);
		seqs->size = 0;
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

int chrom_rank(const char *chr) {
    if (strncmp(chr, "chr", 3) == 0) {
        chr += 3;
    }

    if (isdigit(chr[0])) {
        int n = atoi(chr);
        if (n >= 1 && n <= 22) return n;
    }

    if (strcmp(chr, "X") == 0) return 23;
    if (strcmp(chr, "Y") == 0) return 24;
    if (strcmp(chr, "M") == 0 || strcmp(chr, "MT") == 0) return 25;

    return 1000;  // non-canonical contigs go last
}

void write_sam_hdr(FILE *file, char **paths, int path_size, uint64_t *lens) {

    fprintf(file, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");

    int *is_printed = calloc(path_size, sizeof(int));

    for (int printed = 0; printed < path_size; printed++) {

        int min_i = -1, min_rank = 0;

        for (int i = 0; i < path_size; i++) {
            if (is_printed[i]) continue;

            int r = chrom_rank(paths[i]);

            if (min_i == -1 || r < min_rank || (r == min_rank && strcmp(paths[i], paths[min_i]) < 0)) {
                min_i = i;
                min_rank = r;
            }
        }

        if (min_i == -1)
            break;

        fprintf(file, "@SQ\tSN:%s\tLN:%lu\n", paths[min_i], lens[min_i]);

        is_printed[min_i] = 1;
    }

    free(is_printed);

    fprintf(file, "@PG\tID:LCPan\tPN:LCPan\tNV:1.0\n");
    fprintf(file, "@PG\tID:GraphAligner\tPN:GraphAligner\n");
    fprintf(file, "@RG\tID:%s.%d\tPL:%s\tPU:%s\tSM:%s\n", "akhal", 0, "UNKNOWN", "UNKNOWN", "UNKNOWN");
    
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

int parse_cigar(char *cigar, char *ops, int max_ops, int rev) {
    int i = 0, j = 0, num = 0;
    while (cigar[i]) {
        if ('0' <= cigar[i] && cigar[i] <= '9') {
            num = num * 10 + (cigar[i] - '0');
        } else {
            if (j + num >= max_ops) return -1; // out of range
            for (int k = 0; k < num; k++) {
                ops[j++] = cigar[i];
            }
            num = 0;
        }
        i++;
    }

    if (j > max_ops) return -1;
    
    if (rev) {
        int k=0;
        while (k < j/2) {
            char temp = ops[k];
            ops[k] = ops[j-k-1];
            ops[j-k-1] = temp;
            k++;
        }
    }
    return j;  // number of ops
}

void reverseComplement(char *sequence, int length) {
    if (length == 0) {
        sequence[0] = '\0';
        return;
    }
    int i = 0;
    while (i < length / 2) {
        char temp1 = sequence[i];
        char temp2 = sequence[length-i-1];
        sequence[i] = complement(temp2);
        sequence[length-i-1] = complement(temp1);
        i++;
    }

    if (length % 2 == 1) {
        sequence[length / 2] = complement(sequence[length / 2]);
    }
}

char *parse_rg_prefix(const char *header) {
    if (header[0] == '@' || header[0] == '>') header++;

    const char *dot = strchr(header, '.');
    const char *slash = strchr(header, '/');
    const char *uline = strchr(header, '_');
    const char *end = NULL;
    if (dot && slash) end = dot < slash ? dot : slash;
    else if (dot) end = dot;
    else if (slash) end = slash;
    else if (uline) end = uline;
    else end = header + strlen(header); // no delimiter found, use full string

    size_t len = end - header;
    char *prefix = (char *)malloc(len + 1);
    if (!prefix) return NULL;

    strncpy(prefix, header, len);
    prefix[len] = '\0';
    return prefix;
}

int next_node(const char *path, uint64_t *id, char *strand) {
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