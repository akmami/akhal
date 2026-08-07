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
 * `stats` - print summary statistics for an r/GFA graph.
 * @param argc Argument count.
 * @param argv Arguments; argv[2] is the graph path.
 * @return 0 on success, non-zero on failure.
 */
int cmd_stats(int argc, char **argv);

/**
 * `extract` - write each path's sequence as FASTA.
 * @param argc Argument count.
 * @param argv Arguments; expects "fa", the graph path, and an output path.
 * @return 0 on success, non-zero on failure.
 */
int cmd_extract(int argc, char **argv);

/**
 * `parse` - validate an r/GFA graph.
 * @param argc Argument count.
 * @param argv Arguments; argv[2] is the graph path.
 * @return 0 if valid, non-zero if issues were found or on error.
 */
int cmd_parse(int argc, char **argv);

/**
 * `sort` - topologically sort a graph and renumber nodes 1..N.
 * @param argc Argument count.
 * @param argv Arguments; input .gfa and optional output .gfa.
 * @return 0 on success, non-zero on failure.
 */
int cmd_sort(int argc, char **argv);

/**
 * `vg2gfa` - convert vg's native .vg format to GFA.
 * @param argc Argument count.
 * @param argv Arguments; input .vg and optional output .gfa.
 * @return 0 on success, non-zero on failure.
 */
int cmd_vg2gfa(int argc, char **argv);

/**
 * `gaf2sam` - convert GAF alignments to SAM.
 * @param argc Argument count.
 * @param argv Arguments; graph, GAF, reads FASTA, output SAM, optional --simple.
 * @return 0 on success, non-zero on failure.
 */
int cmd_gaf2sam(int argc, char **argv);

/**
 * `sampoke` - validate SAM CIGAR/positions against a reference.
 * @param argc Argument count.
 * @param argv Arguments; reference FASTA, input SAM, optional output SAM.
 * @return 0 on success, non-zero on failure.
 */
int cmd_sampoke(int argc, char **argv);

#endif  // AKHAL_CLI_H
