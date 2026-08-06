#include "akhal/util.h"

#include <math.h>
#include <string.h>

/** Case-insensitive complement of one DNA base; unknown maps to 'N'. */
char ak_complement(char base) {
	switch (base & 0xDF) {   // uppercase
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
		default:  return 'N';
	}
}

/**
 * Reverse-complement a sequence in place.
 * @param seq Sequence buffer to transform.
 * @param len Number of bases in seq.
 */
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
	if (i == j) seq[i] = ak_complement(seq[i]);
}

/** @return 1 if str ends with suffix, else 0. */
int ak_ends_with(const char *str, const char *suffix) {
	size_t ls = strlen(str), lf = strlen(suffix);
	if (lf > ls) return 0;
	return strcmp(str + ls - lf, suffix) == 0;
}

/**
 * Arithmetic mean of an array of values.
 * @param a Array of values.
 * @param n Number of values; 0 is safe.
 * @return The mean, or 0.0 when n is 0.
 */
double ak_mean(const size_t *a, size_t n) {
	if (n == 0) return 0.0;
	double sum = 0.0;
	for (size_t i = 0; i < n; i++) sum += (double)a[i];
	return sum / (double)n;
}

/**
 * Population variance of an array of values.
 * @param a Array of values.
 * @param n Number of values; 0 is safe.
 * @param mean Precomputed mean (see ak_mean).
 * @return The variance, or 0.0 when n is 0.
 */
double ak_variance(const size_t *a, size_t n, double mean) {
	if (n == 0) return 0.0;
	double acc = 0.0;
	for (size_t i = 0; i < n; i++) {
		double d = (double)a[i] - mean;
		acc += d * d;
	}
	return acc / (double)n;
}

/** @return Square root of the given variance. */
double ak_stddev(double variance) {
	return sqrt(variance);
}
