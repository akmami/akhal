#ifndef AKHAL_VCF_H
#define AKHAL_VCF_H

#include <stdint.h>

/**
 * Minimal streaming VCF reader.
 *
 * Only the fields the toolkit needs are parsed (CHROM, POS, ID, REF, ALT);
 * QUAL, FILTER, INFO and the sample columns are skipped. The reader follows
 * the same streaming shape as the GAF reader: vcf_open() + vcf_read1() pull
 * one record at a time into a caller-owned, reusable vcf_rec_t, so memory
 * stays flat for arbitrarily large variant files.
 *
 * Header lines (starting with '#') are skipped. Malformed data lines are
 * logged and skipped; only allocation failures abort a read.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char   *chrom;   // owned: reference/contig name
    int64_t pos;     // 1-based position of the first REF base
    char   *id;      // owned: variant id (e.g. an rsID), NULL when '.'
    char   *ref;     // owned: reference allele
    char   *alt;     // owned: alternate allele(s), comma-separated; NULL when '.'
} vcf_rec_t;

/**
 * Zero-initialize a record. Does not allocate
 * @param r Record to initialize
 */
void vcf_rec_init(vcf_rec_t *r);

/**
 * Free any owned strings and reset the record for reuse
 * @param r Record to clear
 */
void vcf_rec_clear(vcf_rec_t *r);

typedef struct vcf_reader vcf_reader_t;

/**
 * Open a VCF file for streaming
 * @param fn Path to the .vcf file
 * @return A reader, or NULL on failure (logged)
 */
vcf_reader_t *vcf_open(const char *fn);

/**
 * Parse the next data record. The record is cleared first, so the same one
 * can be reused across the whole file. Header and malformed lines are
 * skipped transparently
 * @param r Open reader
 * @param rec Destination record (reused each call)
 * @return 1 on success, 0 at end of file, or a negative AK_E* code on error
 */
int vcf_read1(vcf_reader_t *r, vcf_rec_t *rec);

/**
 * Close a reader. Safe with NULL. Does not touch caller-owned records
 * @param r Reader to close
 */
void vcf_close(vcf_reader_t *r);

#ifdef __cplusplus
}
#endif

#endif
