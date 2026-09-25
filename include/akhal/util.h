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

// Selections

/**
 * What ak_select() and ak_select_check() made of a selection. Anything but
 * AK_SELECT_OK is a verdict rather than a failure: the selection cannot stand,
 * and `what` points at the text that broke it
 */
enum {
    AK_SELECT_OK = 0,
    AK_SELECT_NONE,      // there are no candidates to pick from
    AK_SELECT_EMPTY,     // a list entry is empty ("a,,b", a trailing comma)
    AK_SELECT_TWICE,     // a name is listed more than once
    AK_SELECT_MISSING,   // a listed name is not among the candidates
    AK_SELECT_SHARED     // under "all", two candidates carry the same name
};

/**
 * Check a selection's own shape, before there is anything to select from:
 * that no list entry is empty and none is given twice. ak_select() applies
 * the same rules, so this is for failing fast - a command can reject a bad
 * list before spending the time to read its input
 * @param spec NULL, "all", or a comma-separated list of names
 * @param what Set to the offending text when the verdict is not AK_SELECT_OK
 * @param what_len Its length; the text is not NUL-terminated
 * @return AK_SELECT_OK, AK_SELECT_EMPTY or AK_SELECT_TWICE
 */
int ak_select_check(const char *spec, const char **what, int *what_len);

/**
 * Pick names out of a set of candidates, the way a --ref style option names
 * them: NULL picks the first candidate, "all" picks every candidate in order,
 * and anything else is a comma-separated list of names, blanks around each
 * ignored, picked in the order listed.
 *
 * Every listed name must be a candidate - the first one carrying it is taken
 * - and none may be listed twice. Under "all" no two candidates may share a
 * name: for a graph's paths that is a path written as fragments, which a name
 * would only ever select the first piece of. The candidates are any array of
 * names, such as a graph's path names (g->path) or a FASTA's record names.
 * `spec` is not modified
 * @param cand Candidate names
 * @param n_cand How many there are
 * @param spec Selection, as above
 * @param pick Receives indices into cand, in selection order; room for n_cand
 * @param n_pick Set to how many were picked; 0 unless the verdict is AK_SELECT_OK
 * @param what Set to the offending text when the verdict is not AK_SELECT_OK
 * @param what_len Its length; the text is not NUL-terminated
 * @return AK_SELECT_OK, or the rule the selection broke
 */
int ak_select(const char *const *cand, int32_t n_cand, const char *spec, int32_t *pick, int32_t *n_pick, const char **what, int *what_len);

/**
 * Word a selection verdict for a message: "no P line named 'chrZ'"
 * @param buf Buffer to write into
 * @param size Its size
 * @param verdict What ak_select() or ak_select_check() returned
 * @param what The offending text they pointed at
 * @param what_len Its length
 * @param noun What the candidates are, singular: "P line", "record"
 * @return buf
 */
char *ak_select_msg(char *buf, size_t size, int verdict, const char *what, int what_len, const char *noun);

/**
 * Number formatting with thousands separators, for readable output. Each
 * writes into a buffer the caller owns and returns it, so the call can sit
 * inline in a printf() argument list. AK_NUM_LEN is enough for any value.
 */
#define AK_NUM_LEN 48

/**
 * Format an unsigned integer with commas: 543020649 -> "543,020,649"
 * @param buf Buffer of at least AK_NUM_LEN bytes
 * @param v Value
 * @return buf
 */
char *ak_format_u64(char *buf, uint64_t v);

/**
 * Format a signed integer with commas: -1234567 -> "-1,234,567"
 * @param buf Buffer of at least AK_NUM_LEN bytes
 * @param v Value
 * @return buf
 */
char *ak_format_i64(char *buf, int64_t v);

/**
 * Format a real with commas in the integer part and a fixed number of
 * decimals: 12345.678 with prec 2 -> "12,345.68". NaN and infinities are
 * written as printf() would
 * @param buf Buffer of at least AK_NUM_LEN bytes
 * @param v Value
 * @param prec Digits after the point, 0..20
 * @return buf
 */
char *ak_format_f64(char *buf, double v, int prec);

/**
 * Format a byte count the way a person reads it: 1536 -> "1.5 KB",
 * 2199023255552 -> "2.0 TB". Powers of 1024, one decimal above the first unit
 * @param buf Buffer of at least AK_NUM_LEN bytes
 * @param bytes Value
 * @return buf
 */
char *ak_format_bytes(char *buf, uint64_t bytes);

// Process measurements, for code that reports what it cost

/**
 * Seconds from a monotonic clock, for timing a span. The origin is arbitrary,
 * so only differences mean anything
 * @return The reading, or 0.0 where no monotonic clock is available
 */
double ak_realtime(void);

/**
 * Peak resident set size of this process so far, in bytes - the high-water
 * mark, not what is held right now, so it never falls
 * @return The peak, or 0 where it cannot be read
 */
size_t ak_peak_rss(void);

/**
 * Resident set size of this process right now, in bytes. Unlike the peak this
 * can fall - but only as far as the allocator lets it, since a freed block it
 * keeps for reuse stays resident, so this measures what the process holds
 * rather than what any one structure in it costs
 * @return The size, or 0 where it cannot be read
 */
size_t ak_rss(void);

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
