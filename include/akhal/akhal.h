#ifndef AKHAL_H
#define AKHAL_H

/**
 * libakhal umbrella header.
 *
 * akhal is split into a reusable library (libakhal, under src/lib/) and a
 * thin command-line front end (src/cli/), mirroring the htslib + samtools
 * separation. Include this single header to pull in the whole public API,
 * or include the individual headers you need.
 *
 *   akhal/error.h   diagnostics + return-code conventions
 *   akhal/kstr.h    growable string buffer
 *   akhal/io.h      line-oriented input handle
 *   akhal/util.h    sequence / statistics helpers
 *   akhal/gfa.h     (r)GFA graph model, reader and traversal
 *   akhal/gaf.h     GAF alignment records, streaming + batch readers
 *   akhal/fasta.h   FASTA loaded into a name-indexed store
 *   akhal/sam.h     SAM output + CIGAR utilities
 */

#include "akhal/version.h"
#include "akhal/error.h"
#include "akhal/kstr.h"
#include "akhal/io.h"
#include "akhal/util.h"
#include "akhal/gfa.h"
#include "akhal/gaf.h"
#include "akhal/fasta.h"
#include "akhal/sam.h"

#endif  // AKHAL_H
