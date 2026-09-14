#ifndef AKHAL_UTIL_H
#define AKHAL_UTIL_H

#include <stddef.h>
#include <stdint.h>

/**
 * Small, dependency-free helpers shared across the library and CLI
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Case-insensitive complement of a single DNA base
 * @param base A nucleotide character in any case
 * @return The complementary base, or 'N' if unrecognized
 */
char ak_complement(char base);

/**
 * Reverse-complement a sequence in place
 * @param seq Sequence buffer to transform
 * @param len Number of bases in seq (excludes any NUL terminator)
 */
void ak_revcomp(char *seq, size_t len);

/**
 * Test whether a string ends with a given suffix
 * @param str String to inspect
 * @param suffix Suffix to look for
 * @return 1 if str ends with suffix, otherwise 0
 */
int ak_ends_with(const char *str, const char *suffix);

/**
 * Parse a string as a base-10 int, requiring the whole string to be consumed
 * @param str String to parse; NULL or empty is rejected
 * @param out Set to the parsed value on success; untouched on failure
 * @return 1 if str is a valid int, otherwise 0
 */
int ak_str2int(const char *str, int *out);

/**
 * Running summary of a distribution: count, mean, population variance and
 * extremes, accumulated one value at a time with Welford's update so nothing
 * has to be collected into an array first. Zero-initialize and feed values
 * with ak_dist_add(); every query is safe on an empty accumulator.
 */
typedef struct {
    int64_t n;           // values seen
    double  mean;        // running mean
    double  m2;          // sum of squared deviations from the running mean
    double  min, max;    // extremes; meaningless while n == 0
} ak_dist_t;

/**
 * Add one value to a running distribution
 * @param d Accumulator, zero-initialized before the first call
 * @param x The value
 */
void ak_dist_add(ak_dist_t *d, double x);

/**
 * Population variance (divides by n, not n - 1) of the values seen so far
 * @param d Accumulator
 * @return The variance, or 0.0 while fewer than one value has been added
 */
double ak_dist_variance(const ak_dist_t *d);

/**
 * Population standard deviation of the values seen so far
 * @param d Accumulator
 * @return sqrt of ak_dist_variance(), so 0.0 on an empty accumulator
 */
double ak_dist_sd(const ak_dist_t *d);

#ifdef __cplusplus
}
#endif

#endif
