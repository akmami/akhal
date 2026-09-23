#include "akhal/util.h"

#include <errno.h>
#include <sys/resource.h>
#include <unistd.h>
#if defined(__APPLE__)
#include <mach/mach.h>
#endif
#include <time.h>
#include <stdio.h>
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

// commas every three digits, written right to left from a small scratch of
// plain digits; the sign, if any, goes on first
static char *format_digits(char *buf, const char *digits, size_t n, int neg) {
    char *o = buf;
    if (neg) *o++ = '-';
    for (size_t i = 0; i < n; i++) {
        if (i && (n - i) % 3 == 0) *o++ = ',';
        *o++ = digits[i];
    }
    *o = '\0';
    return buf;
}

char *ak_format_u64(char *buf, uint64_t v) {
    char d[24];
    int n = snprintf(d, sizeof d, "%llu", (unsigned long long)v);
    return format_digits(buf, d, (size_t)n, 0);
}

char *ak_format_i64(char *buf, int64_t v) {
    // negate as unsigned so INT64_MIN does not overflow
    uint64_t u = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    char d[24];
    int n = snprintf(d, sizeof d, "%llu", (unsigned long long)u);
    return format_digits(buf, d, (size_t)n, v < 0);
}

char *ak_format_f64(char *buf, double v, int prec) {
    if (prec < 0) prec = 0;
    if (prec > 20) prec = 20;
    if (isnan(v) || isinf(v)) {
        snprintf(buf, AK_NUM_LEN, "%f", v);
        return buf;
    }
    // let printf do the rounding, then re-emit the integer part with commas
    // and copy the fraction through untouched
    char s[64];
    snprintf(s, sizeof s, "%.*f", prec, v);
    const char *p = s;
    int neg = (*p == '-');
    if (neg) p++;
    const char *dot = strchr(p, '.');
    size_t n = dot ? (size_t)(dot - p) : strlen(p);
    format_digits(buf, p, n, neg);
    if (dot) strcat(buf, dot);
    return buf;
}

char *ak_format_bytes(char *buf, uint64_t bytes) {
    static const char *unit[] = { "B", "KB", "MB", "GB", "TB", "PB" };
    double v = (double)bytes;
    int u = 0;
    while (v >= 1024.0 && u < 5) {
        v /= 1024.0;
        u++;
    }
    // whole bytes below the first step up; a fraction of one is noise
    if (u == 0) snprintf(buf, AK_NUM_LEN, "%llu B", (unsigned long long)bytes);
    else        snprintf(buf, AK_NUM_LEN, "%.3f %s", v, unit[u]);
    return buf;
}

double ak_realtime(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) return 0.0;
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

size_t ak_peak_rss(void) {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0;
    // ru_maxrss is kilobytes on Linux and bytes on the BSDs, macOS included
#ifdef __APPLE__
    return (size_t)ru.ru_maxrss;
#else
    return (size_t)ru.ru_maxrss * 1024;
#endif
}

size_t ak_rss(void) {
#if defined(__APPLE__)
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t n = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &n) != KERN_SUCCESS) return 0;
    return (size_t)info.resident_size;
#elif defined(__linux__)
    // field 2 of statm is the resident page count
    FILE *f = fopen("/proc/self/statm", "r");
    if (!f) return 0;
    long total = 0, res = 0;
    int ok = (fscanf(f, "%ld %ld", &total, &res) == 2);
    fclose(f);
    long page = sysconf(_SC_PAGESIZE);
    return (ok && res > 0 && page > 0) ? (size_t)res * (size_t)page : 0;
#else
    return 0;
#endif
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
