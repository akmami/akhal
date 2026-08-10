/**
 * akhal command-line front end.
 *
 * This file is deliberately thin: it maps a subcommand name to a function and
 * hands off. All real work lives in the per-command files and in libakhal. To
 * add a command, implement `int cmd_x(int, char**)`, declare it in cli.h, and
 * add one row to the table below.
 */

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
    { "stats",   cmd_stats,   "graph", "summary statistics for an r/GFA graph"      },
    { "parse",   cmd_parse,   "graph", "validate an r/GFA graph"                    },
    { "extract", cmd_extract, "graph", "extract the reference FASTA from an r/GFA"  },
    { "sort",    cmd_sort,    "graph", "topologically sort and renumber a graph"    },
    { "vg2gfa",  cmd_vg2gfa,  "graph", "convert vg (protobuf .vg) to GFA"           },
    { "gaf2sam", cmd_gaf2sam, "align", "convert GAF alignments to SAM"              },
    { "sampoke", cmd_sampoke, "align", "validate SAM CIGAR/positions against a ref" },
    { NULL, NULL, NULL, NULL }
};

/** 
 * @brief Print the top-level usage / command list. 
 */
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

/** 
 * @brief Dispatch argv[1] to a subcommand. 
 * @return the command's exit status. 
 */
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
        if (!strcmp(argv[1], c->name))
            return c->fn(argc, argv);
    }

    ak_log(AK_LOG_ERROR, NULL, "unknown command: %s", argv[1]);
    usage(stderr);
    return 1;
}
