#include "akhal/version.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <string.h>

typedef int (*cmd_fn)(int, char **);

typedef struct {
    const char *name;
    cmd_fn      fn;
    const char *group;
    const char *desc;
} command_t;

static const command_t commands[] = {
    { "stats",    cmd_stats,    "graph", "summary statistics for an r/GFA graph"      },
    { "parse",    cmd_parse,    "graph", "validate an r/GFA graph"                    },
    { "extract",  cmd_extract,  "graph", "extract FASTA, paths, or a VCF"             },
    { "compare",  cmd_compare,  "graph", "compare two graphs"                         },
    { "compact",  cmd_compact,  "graph", "fold non-branching runs into one node"      },
    { "sort",     cmd_sort,     "graph", "topologically sort and renumber a graph"    },
    { "rank",     cmd_rank,     "graph", "rewrite SR ranks against the backbone"      },
    { "vg2gfa",   cmd_vg2gfa,   "graph", "convert vg to GFA"                          },
    { "gfa2rgfa", cmd_gfa2rgfa, "graph", "label a GFA as rGFA (SN/SO/SR tags)"        },
    { "gfa2dot",  cmd_gfa2dot,  "graph", "write a graph as Graphviz, to look at it"   },
    { "gaf2sam",  cmd_gaf2sam,  "align", "convert GAF alignments to SAM"              },
    { "sampoke",  cmd_sampoke,  "align", "validate SAM CIGAR/positions against a ref" },
    { "annotate", cmd_annotate, "annot", "trace node origins (VCF/FASTA) to a .annot file" },
    { "annotget", cmd_annotget, "annot", "look up node annotations in a .annot file" },
    { NULL, NULL, NULL, NULL }
};

// print the top-level usage / command list
static void usage(FILE *out) {
    fprintf(out, "akhal %s - assembly graph analysis toolkit\n\n", AKHAL_VERSION_STR);
    fprintf(out, "Usage: akhal <command> [options]\n\n");

    const char *group = NULL;
    for (const command_t *c = commands; c->name; c++) {
        if (!group || strcmp(group, c->group) != 0) {
            group = c->group;
            fprintf(out, "%s:\n", group);
        }
        fprintf(out, "  %-10s %s\n", c->name, c->desc);
    }
    fprintf(out, "\nRun `akhal <command>` with no arguments for command-specific help.\n");
}

// dispatch argv[1] to a subcommand
int main(int argc, char **argv) {
    if (argc < 2) {
        usage(stderr);
        return 1;
    }

    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        usage(stdout);
        return 0;
    }
    if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--version")) {
        printf("akhal %s\n", AKHAL_VERSION_STR);
        return 0;
    }

    for (const command_t *c = commands; c->name; c++) {
        if (!strcmp(argv[1], c->name)) return c->fn(argc, argv);
    }

    ak_log(AK_LOG_ERROR, NULL, "unknown command: %s", argv[1]);
    usage(stderr);
    return 1;
}
