#ifndef AKHAL_CLI_H
#define AKHAL_CLI_H

/**
 * Subcommand entry points.
 *
 * Each command has the signature `int cmd(int argc, char **argv)` where
 * argv[1] is the command name. A command returns 0 on success or non-zero on
 * failure; main() forwards that as the process exit status. Every command is
 * built on libakhal.
 */

/**
 * `stats` - print summary statistics for an r/GFA graph
 * @param argc Argument count
 * @param argv Arguments; argv[2] is the graph path
 * @return 0 on success, non-zero on failure
 */
int cmd_stats(int argc, char **argv);

/**
 * `extract` - pull information out of a graph. The target is argv[2]:
 * "fa" writes each path's sequence as FASTA, "path" stitches a fragmented
 * path back together first, and "vcf" reports the graph's variations over
 * its reference backbone
 * @param argc Argument count
 * @param argv Arguments; the target, the graph path, an output path, then any
 *             target-specific arguments
 * @return 0 on success, non-zero on failure
 */
int cmd_extract(int argc, char **argv);

/**
 * `compare` - compare two graphs that need not agree on segment ids: segments
 * are matched by sequence content, links through the labelling that produces,
 * and paths by the bases they spell
 * @param argc Argument count
 * @param argv Arguments; the two graph paths, then optional --verbose
 * @return 0 when the graphs match, 1 when they differ, 2 when the comparison
 *         could not be made
 */
int cmd_compare(int argc, char **argv);

/**
 * `parse` - validate an r/GFA graph
 * @param argc Argument count
 * @param argv Arguments; argv[2] is the graph path
 * @return 0 if valid, non-zero if issues were found or on error
 */
int cmd_parse(int argc, char **argv);

/**
 * `sort` - topologically sort a graph and renumber nodes 1..N
 * @param argc Argument count
 * @param argv Arguments; input .gfa and optional output .gfa
 * @return 0 on success, non-zero on failure
 */
int cmd_sort(int argc, char **argv);

/**
 * `rank` - rewrite a graph's SR ranks against its backbone, consolidating
 * fragmented paths or replacing them with a traced reference
 * @param argc Argument count
 * @param argv Arguments; input .gfa, optional output .gfa, optional
 *             --fasta/--ref
 * @return 0 on success, non-zero on failure
 */
int cmd_rank(int argc, char **argv);

/**
 * `vg2gfa` - convert vg's native .vg format to GFA
 * @param argc Argument count
 * @param argv Arguments; input .vg and optional output .gfa
 * @return 0 on success, non-zero on failure
 */
int cmd_vg2gfa(int argc, char **argv);

/**
 * `gaf2sam` - convert GAF alignments to SAM
 * @param argc Argument count
 * @param argv Arguments; graph, GAF, reads FASTA, output SAM, optional --simple
 * @return 0 on success, non-zero on failure
 */
int cmd_gaf2sam(int argc, char **argv);

/**
 * `sampoke` - validate SAM CIGAR/positions against a reference
 * @param argc Argument count
 * @param argv Arguments; reference FASTA, input SAM, optional output SAM
 * @return 0 on success, non-zero on failure
 */
int cmd_sampoke(int argc, char **argv);

/**
 * `annotate` - trace node origins (backbone / VCF variants / FASTA walks)
 * and write a binary .annot file
 * @param argc Argument count
 * @param argv Arguments; graph, output .annot, optional --vcf/--fasta/--ref
 * @return 0 on success, non-zero on failure
 */
int cmd_annotate(int argc, char **argv);

/**
 * `annotget` - look up node annotations in a .annot file
 * @param argc Argument count
 * @param argv Arguments; the .annot file, then node ids (none = dump all)
 * @return 0 on success, non-zero on failure
 */
int cmd_annotget(int argc, char **argv);

#endif
