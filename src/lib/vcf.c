#include "akhal/vcf.h"
#include "akhal/io.h"
#include "akhal/kstr.h"
#include "akhal/error.h"

#include <stdlib.h>
#include <string.h>

struct vcf_reader {
    ak_file  *f;
    kstring_t line;
    long      lineno;
};

// zero a record; see akhal/vcf.h
void vcf_rec_init(vcf_rec_t *r) {
    memset(r, 0, sizeof(*r));
}

// free owned strings and reset; see akhal/vcf.h
void vcf_rec_clear(vcf_rec_t *r) {
    free(r->chrom);
    free(r->id);
    free(r->ref);
    free(r->alt);
    memset(r, 0, sizeof(*r));
}

// open a VCF for streaming; see akhal/vcf.h
vcf_reader_t *vcf_open(const char *fn) {
    ak_file *f = ak_open(fn);
    if (!f) return NULL;
    vcf_reader_t *r = (vcf_reader_t *)calloc(1, sizeof(*r));
    if (!r) {
        ak_close(f);
        ak_log(AK_LOG_ERROR, "vcf", "out of memory");
        return NULL;
    }
    r->f = f;
    return r;
}

// strdup, except a lone "." yields NULL rather than a copy
static int dup_field(const char *tok, char **out) {
    if (tok[0] == '.' && tok[1] == '\0') {
        *out = NULL;
        return AK_OK;
    }
    *out = strdup(tok);
    return *out ? AK_OK : AK_ENOMEM;
}

// parse the next data record; see akhal/vcf.h
int vcf_read1(vcf_reader_t *r, vcf_rec_t *rec) {
    long len;
    while ((len = ak_getline(r->f, &r->line)) >= 0) {
        r->lineno++;
        if (len == 0 || r->line.s[0] == '#') continue;   // header / empty

        char *save;
        char *chrom = strtok_r(r->line.s, "\t", &save);
        char *pos   = strtok_r(NULL, "\t", &save);
        char *id    = strtok_r(NULL, "\t", &save);
        char *ref   = strtok_r(NULL, "\t", &save);
        char *alt   = strtok_r(NULL, "\t", &save);
        if (!chrom || !pos || !id || !ref || !alt) {
            ak_log(AK_LOG_WARN, "vcf", "line %ld: fewer than 5 fields, skipped", r->lineno);
            continue;
        }

        vcf_rec_clear(rec);
        rec->pos = strtoll(pos, NULL, 10);
        if (rec->pos <= 0) {
            ak_log(AK_LOG_WARN, "vcf", "line %ld: bad POS '%s', skipped", r->lineno, pos);
            continue;
        }
        rec->chrom = strdup(chrom);
        if (!rec->chrom) return AK_ENOMEM;
        if (dup_field(id,  &rec->id)  != AK_OK) return AK_ENOMEM;
        rec->ref = strdup(ref);
        if (!rec->ref) return AK_ENOMEM;
        if (dup_field(alt, &rec->alt) != AK_OK) return AK_ENOMEM;
        return 1;
    }
    return 0;
}

// close a reader; see akhal/vcf.h
void vcf_close(vcf_reader_t *r) {
    if (!r) return;
    ak_close(r->f);
    ks_free(&r->line);
    free(r);
}
