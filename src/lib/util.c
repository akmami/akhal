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

// selections

static int is_blank(char c) {
    return c == ' ' || c == '\t';
}

// the next field of a comma-separated list, trimmed of the blanks around it
static int next_field(const char **cur, const char **s, int *len) {
    const char *p = *cur;
    if (!p) return 0;
    const char *comma = strchr(p, ',');
    const char *e = comma ? comma : p + strlen(p);
    while (p < e && is_blank(*p)) p++;
    while (e > p && is_blank(e[-1])) e--;
    *s = p;
    *len = (int)(e - p);
    *cur = comma ? comma + 1 : NULL;
    return 1;
}

// the first candidate named exactly by the text s[0 .. len), or -1
static int32_t find_name(const char *const *cand, int32_t n_cand, const char *s, int len) {
    for (int32_t k = 0; k < n_cand; k++)
        if ((int)strlen(cand[k]) == len && !memcmp(cand[k], s, (size_t)len)) return k;
    return -1;
}

// order pick[] by the names it points at: a shellsort
static void sort_by_name(const char *const *cand, int32_t *pick, int32_t n) {
    for (int32_t gap = n / 2; gap > 0; gap /= 2) {
        for (int32_t i = gap; i < n; i++) {
            int32_t v = pick[i], j = i;
            for (; j >= gap && strcmp(cand[pick[j - gap]], cand[v]) > 0; j -= gap) pick[j] = pick[j - gap];
            pick[j] = v;
        }
    }
}

int ak_select_check(const char *spec, const char **what, int *what_len) {
    if (!spec || !strcmp(spec, "all")) return AK_SELECT_OK;

    const char *cur = spec, *s;
    int len;
    while (next_field(&cur, &s, &len)) {
        if (len == 0) {
            *what = s;
            *what_len = 0;
            return AK_SELECT_EMPTY;
        }
        // against every field before this one
        const char *prev = spec, *t;
        int tlen;
        while (next_field(&prev, &t, &tlen) && t < s) {
            if (tlen == len && !memcmp(t, s, (size_t)len)) {
                *what = s;
                *what_len = len;
                return AK_SELECT_TWICE;
            }
        }
    }
    return AK_SELECT_OK;
}

int ak_select(const char *const *cand, int32_t n_cand, const char *spec, int32_t *pick, int32_t *n_pick, const char **what, int *what_len) {
    *n_pick = 0;
    int verdict = ak_select_check(spec, what, what_len);
    if (verdict != AK_SELECT_OK) return verdict;
    if (n_cand <= 0) {
        *what = "";
        *what_len = 0;
        return AK_SELECT_NONE;
    }

    if (!spec) {
        pick[0] = 0;
        *n_pick = 1;
        return AK_SELECT_OK;
    }

    if (!strcmp(spec, "all")) {
        // shared names show up side by side once the picks are ordered by name; the picks are then put back in file order
        for (int32_t k = 0; k < n_cand; k++) pick[k] = k;
        sort_by_name(cand, pick, n_cand);
        for (int32_t i = 1; i < n_cand && verdict == AK_SELECT_OK; i++) {
            if (!strcmp(cand[pick[i - 1]], cand[pick[i]])) {
                *what = cand[pick[i]];
                *what_len = (int)strlen(*what);
                verdict = AK_SELECT_SHARED;
            }
        }
        if (verdict != AK_SELECT_OK) return verdict;
        for (int32_t k = 0; k < n_cand; k++) pick[k] = k;
        *n_pick = n_cand;
        return AK_SELECT_OK;
    }

    // a list of distinct names each picks a different candidate, so pick[] never needs more room than there are candidates
    const char *cur = spec, *s;
    int len;
    while (next_field(&cur, &s, &len)) {
        int32_t k = find_name(cand, n_cand, s, len);
        if (k < 0) {
            *n_pick = 0;
            *what = s;
            *what_len = len;
            return AK_SELECT_MISSING;
        }
        pick[(*n_pick)++] = k;
    }
    return AK_SELECT_OK;
}

char *ak_select_msg(char *buf, size_t size, int verdict, const char *what, int what_len, const char *noun) {
    switch (verdict) {
    case AK_SELECT_NONE:
        snprintf(buf, size, "there is no %s to choose from", noun);
        break;
    case AK_SELECT_EMPTY:
        snprintf(buf, size, "the list has an empty name");
        break;
    case AK_SELECT_TWICE:
        snprintf(buf, size, "'%.*s' is listed more than once", what_len, what);
        break;
    case AK_SELECT_MISSING:
        snprintf(buf, size, "no %s named '%.*s'", noun, what_len, what);
        break;
    case AK_SELECT_SHARED:
        snprintf(buf, size, "'%.*s' names more than one %s", what_len, what, noun);
        break;
    default:
        snprintf(buf, size, "the selection stands");
        break;
    }
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
