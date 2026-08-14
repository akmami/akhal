#ifndef AKHAL_SAM_H
#define AKHAL_SAM_H

#include <stdio.h>
#include <stdint.h>

/**
 * SAM output helpers and CIGAR utilities.
 *
 * This module owns the SAM-side conventions used by the toolkit: the CIGAR
 * operation alphabet, expansion of a CIGAR string into a per-base op array,
 * and writing a SAM header and alignment records. Commands assemble the
 * fields and hand them here.
 */

#ifdef __cplusplus
extern "C" {
#endif

// Upper bound on a per-base CIGAR op array / rendered CIGAR string.
#define SAM_MAX_CIGAR 1048576

// CIGAR operations                       consumes query / reference
#define CIGAR_ALIGNMENT_MATCH   'M'   // yes  yes
#define CIGAR_INSERTION         'I'   // yes  no
#define CIGAR_DELETION          'D'   // no   yes
#define CIGAR_SKIPPED           'N'   // no   yes
#define CIGAR_SOFT_CLIP         'S'   // yes  no
#define CIGAR_HARD_CLIP         'H'   // no   no
#define CIGAR_PADDING           'P'   // no   no
#define CIGAR_SEQUENCE_MATCH    '='   // yes  yes
#define CIGAR_SEQUENCE_MISMATCH 'X'   // yes  yes

// 1 if op consumes the query sequence.
#define CIGAR_QUE(x) ((x)==CIGAR_ALIGNMENT_MATCH || (x)==CIGAR_INSERTION || \
                      (x)==CIGAR_SOFT_CLIP || (x)==CIGAR_SEQUENCE_MATCH || \
                      (x)==CIGAR_SEQUENCE_MISMATCH)
// 1 if op consumes the reference.
#define CIGAR_REF(x) ((x)==CIGAR_ALIGNMENT_MATCH || (x)==CIGAR_DELETION || \
                      (x)==CIGAR_SKIPPED || (x)==CIGAR_SEQUENCE_MATCH || \
                      (x)==CIGAR_SEQUENCE_MISMATCH)

// SAM FLAG bits used by the toolkit.
#define SAM_FMUNMAP        0x4
#define SAM_FREVERSE       0x10
#define SAM_FSECONDARY     0x100
#define SAM_FSUPPLEMENTARY 0x800

/**
 * Expand a CIGAR string ("3=1I2X") into a per-base op array ("===I=XX...")
 * @param cigar CIGAR string to expand
 * @param ops Destination buffer for per-base ops
 * @param max_ops Capacity of ops
 * @param rev If non-zero, the resulting op array is reversed
 * @return The number of ops written, or -1 on overflow
 */
int sam_cigar_expand(const char *cigar, char *ops, int max_ops, int rev);

/**
 * Write a SAM header: @HD, one @SQ per reference sequence (in chromosome
 * order), a @PG line, and a default @RG
 * @param out Destination stream
 * @param names Reference sequence names (parallel to lens), length n
 * @param n Number of reference sequences
 * @param lens Reference sequence lengths (parallel to names)
 * @param pg Program name for the @PG / @RG id
 */
void sam_write_header(FILE *out, char **names, int n, const uint64_t *lens, const char *pg);

/**
 * One alignment record to emit. `ops` is a per-base op array of length n_ops
 * and may be mutated (leading/trailing insert/mismatch runs become soft clips);
 * if simplify is set, '='/'X' are rendered as 'M'
 */
typedef struct {
    const char *qname;
    int         flag;
    const char *rname;
    int         pos;      // 1-based reference position
    int         mapq;
    char       *ops;
    int         n_ops;
    int         simplify;
    const char *seq;
    int         nm;       // NM:i
    double      as;       // AS:f
    double      dv;       // dv:f
    double      id;       // id:f
    const char *rg;       // read-group id for RG:Z
} sam_rec_t;

/**
 * Render and write a single SAM alignment line
 * @param out Destination stream
 * @param r Record to emit (its ops array may be mutated in place)
 */
void sam_write_record(FILE *out, sam_rec_t *r);

/**
 * Read-group prefix of a read name: the leading component before the first
 * '.', '/', or '_'. Accepts a bare name or a full SAM/FASTA line (a leading
 * '@' or '>' is skipped)
 * @param name Read name or line to derive the prefix from
 * @return A newly allocated prefix string (caller frees), or NULL on OOM
 */
char *sam_rg_prefix(const char *name);

#ifdef __cplusplus
}
#endif

#endif
