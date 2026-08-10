/**
 * `akhal parse <r/GFA>` - validate a graph.
 *
 * gfa_read(GFA_VALIDATE) already checks link overlap consistency and
 * referential integrity while reading (issues are logged). This command adds
 * the structural checks that need the assembled model:
 *   - every pair of consecutive segments in a path is joined by a link;
 *   - for rGFA, the number of rank-0 (SR:i:0) segments matches the number of
 *     segment occurrences laid out on paths.
 */

#include "akhal/gfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>

/** 
 * @brief `parse` entry point; see cli.h. 
 */
int cmd_parse(int argc, char **argv) {
    if (argc < 3) {
        ak_log(AK_LOG_ERROR, NULL, "usage: akhal parse <r/GFA>");
        return 1;
    }
    const char *fn = argv[2];
    if (!ak_ends_with(fn, ".gfa") && !ak_ends_with(fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
        return 1;
    }
    int is_rgfa = ak_ends_with(fn, ".rgfa");

    gfa_t *g = gfa_read(fn, GFA_LINKS | GFA_PATHS | GFA_VALIDATE);
    if (!g) return 1;

    long issues = 0;

    // Consecutive path segments must be connected by a link.
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        const uint32_t *segs;
        int n = gfa_path_segs(g, k, &segs);
        for (int i = 1; i < n; i++) {
            uint32_t a = segs[i - 1], b = segs[i];
            if (a == GFA_NIL || b == GFA_NIL) continue;
            if (!gfa_has_arc(g, (int32_t)a, (int32_t)b)) {
                ak_log(AK_LOG_WARN, "parse", "no link %lu -> %lu present in path %s", (unsigned long)gfa_seg_at(g, (int32_t)a)->id, (unsigned long)gfa_seg_at(g, (int32_t)b)->id, gfa_path_name(g, k));
                issues++;
            }
        }
    }

    // rGFA: rank-0 segment count vs path occurrence count.
    if (is_rgfa) {
        uint64_t n_ref = 0;
        for (int32_t i = 0; i < gfa_n_seg(g); i++)
            if (gfa_seg_at(g, i)->rank == 0) n_ref++;
        if (n_ref != g->n_path_seg) {
            ak_log(AK_LOG_WARN, "parse", "rank-0 segment count (%lu) != path segment count (%lu)", (unsigned long)n_ref, (unsigned long)g->n_path_seg);
            issues++;
        }
    }

    gfa_destroy(g);

    if (issues == 0) {
        printf("[INFO] Parsed %s successfully\n", is_rgfa ? "rGFA" : "GFA");
        return 0;
    }
    ak_log(AK_LOG_WARN, "parse", "%ld issue(s) found (overlap warnings, if any, logged above)", issues);
    return 1;
}
