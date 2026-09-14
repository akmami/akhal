#include "akhal/util.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// case-insensitive complement of one DNA base; unknown maps to 'N'
char ak_complement(char base) {
    switch (base & 0xDF) {   // uppercase
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 'N';
    }
}

void ak_revcomp(char *seq, size_t len) {
    if (!seq || len == 0) return;

    size_t i = 0, j = len - 1;
    while (i < j) {
        char a = seq[i], b = seq[j];
        seq[i] = ak_complement(b);
        seq[j] = ak_complement(a);
        i++;
        j--;
    }
    if (i == j) {
        seq[i] = ak_complement(seq[i]);
    }
}

int ak_ends_with(const char *str, const char *suffix) {
    size_t ls = strlen(str), lf = strlen(suffix);
    if (lf > ls) return 0;
    return strcmp(str + ls - lf, suffix) == 0;
}

// whole-string strtol: rejects "", trailing junk, whitespace-only and overflow
int ak_str2int(const char *str, int *out) {
    if (!str || !*str) return 0;

    errno = 0;
    char *end;
    long v = strtol(str, &end, 10);
    if (*end || end == str) return 0;
    if (errno == ERANGE || v < INT_MIN || v > INT_MAX) return 0;

    *out = (int)v;
    return 1;
}

// Welford's online update: the mean and the sum of squared deviations are
// carried along together, so one pass gives both without a second walk and
// without the cancellation a naive sum-of-squares suffers on large values
void ak_dist_add(ak_dist_t *d, double x) {
    if (d->n == 0) {
        d->min = d->max = x;
    } else if (x < d->min) {
        d->min = x;
    } else if (x > d->max) {
        d->max = x;
    }
    d->n++;
    double delta = x - d->mean;
    d->mean += delta / (double)d->n;
    d->m2   += delta * (x - d->mean);
}

// population variance: divides by n, not n - 1
double ak_dist_variance(const ak_dist_t *d) {
    return d->n ? d->m2 / (double)d->n : 0.0;
}

double ak_dist_sd(const ak_dist_t *d) {
    return sqrt(ak_dist_variance(d));
}
