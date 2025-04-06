#ifndef __STRUCT_DEF__
#define __STRUCT_DEF__

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

//                                Consumes     query   reference   Op
#define CIGAR_ALIGNMENT_MATCH       'M'    //  yes     yes         M
#define CIGAR_INSERTION             'I'    //  yes     no          I
#define CIGAR_DELETION              'D'    //  no      yes         D
#define CIGAR_SKIPPED               'N'    //  no      yes         N
#define CIGAR_SOFT_CLIP             'S'    //  yes     no          S
#define CIGAR_HARD_CLIP             'H'    //  no      no          H
#define CIGAR_PADDING               'P'    //  no      no          P
#define CIGAR_SEQUENCE_MATCH        '='    //  yes     yes         =
#define CIGAR_SEQUENCE_MISMATCH     'X'    //  yes     yes         X

#define CIGAR_QUE(x) (x==CIGAR_ALIGNMENT_MATCH || x==CIGAR_INSERTION || x==CIGAR_SOFT_CLIP || x==CIGAR_SEQUENCE_MATCH || x==CIGAR_SEQUENCE_MISMATCH )
#define CIGAR_REF(x) (x==CIGAR_ALIGNMENT_MATCH || x==CIGAR_DELETION || x==CIGAR_SKIPPED || x==CIGAR_SEQUENCE_MATCH || x==CIGAR_SEQUENCE_MISMATCH )

#define __FLAG_MULTIPLE_SEGMENTS          0x1
#define __FLAG_SECONDARY_ALIGNMENT        0x100
#define __FLAG_SUPPLEMENTARY_ALIGNMENT    0x800

#define FASTA_WRAP_SIZE 80

typedef struct alignment {
    char *readName;
    int readLen;
    int readStart; // The start position of the query (read) in the alignment.
    int readEnd; // The end position of the query (read) in the alignment.
    char strand; // + or -
    char *path; // node names in < or > directions
    int pathLen;
    int pathStart; // The start position of the alignment on the reference node path.
    int pathEnd; // The end position of the alignment on the reference node path.
    int matches; // The number of matching bases between the query and the reference.
    int blockLen; // The length of the aligned block (number of bases in the alignment block).
    int qual; // The mapping quality score, 255 not avaialble 
    int xdi; // NM:i:22 The number of mismatches, deletions, and insertions in the alignment (represented as an integer value).
    float score; // AS:f:14239.3 (optional) The alignment score as a floating-point number, if provided (this field is optional).
    float divergence; // dv:f:0.00153685 The divergence score, which quantifies how different the alignment is by calculating a ratio of mismatches, deletions, and insertions relative to the total alignment length.
    float identity; // id:f:0.998463 The identity score, which is the ratio of matches to the total alignment length (matches, mismatches, deletions, and insertions).
    char *cigar; //  cg:Z:3=1I (optional)
} alignment;

typedef struct segment {
    uint64_t id;
    char *seq;
    char *ref_name;
    int start;
    int end;
    int rank;
    struct segment *next;
    int out_degree;
    int in_degree;
} segment;

#endif