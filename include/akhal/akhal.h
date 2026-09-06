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
 *   akhal/vcf.h     streaming VCF reader
 *   akhal/vg.h      vg native (.vg protobuf) graph reader
 *   akhal/annot.h   node annotation store, builders and file format
 *   akhal/call.h    backbone labelling and variant discovery on a graph
 *   akhal/diff.h    structural comparison of two graphs, ids aside
 *   akhal/rgfa.h    stable-sequence labelling: a GFA turned into an rGFA
 *   akhal/compact.h folding non-branching runs of segments into single nodes
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
#include "akhal/vcf.h"
#include "akhal/vg.h"
#include "akhal/annot.h"
#include "akhal/call.h"
#include "akhal/diff.h"
#include "akhal/rgfa.h"
#include "akhal/compact.h"

#endif
