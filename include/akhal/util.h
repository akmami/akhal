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
 * Arithmetic mean of an array of values
 * @param a Array of values
 * @param n Number of values; 0 is safe
 * @return The mean, or 0.0 when n is 0
 */
double ak_mean(const size_t *a, size_t n);

/**
 * Population variance of an array of values
 * @param a Array of values
 * @param n Number of values; 0 is safe
 * @param mean Precomputed mean of the array (see ak_mean)
 * @return The variance, or 0.0 when n is 0
 */
double ak_variance(const size_t *a, size_t n, double mean);

/**
 * Standard deviation derived from a variance
 * @param variance A variance value (see ak_variance)
 * @return The square root of variance
 */
double ak_stddev(double variance);

#ifdef __cplusplus
}
#endif

#endif
