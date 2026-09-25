#include "akhal/gfa.h"
#include "akhal/util.h"
#include "akhal/error.h"
#include "cli.h"

#include "khashl.h"

#include <stdio.h>
#include <string.h>


#define CHECK_LINKS    0x1   // L lines naming segments no S line defines
#define CHECK_OVERLAPS 0x2   // the bases either side of an overlap agree (loads every sequence; implies CHECK_LINKS)
#define CHECK_PATHS    0x4   // P steps name defined segments and consecutive steps are joined by a link; reused path names are pointed out
#define CHECK_RANKS    0x8   // rGFA only: rank-0 segments match the path steps
#define CHECK_ALL      (CHECK_LINKS | CHECK_OVERLAPS | CHECK_PATHS | CHECK_RANKS)
#define CHECK_BASIC    (CHECK_LINKS | CHECK_PATHS)

// path name -> index of the first path that carried it
KHASHL_MAP_INIT(KH_LOCAL, pathmap_t, pathmap, const char *, int32_t, kh_hash_str, kh_eq_str)

static void usage(void) {
    ak_log(AK_LOG_ERROR, NULL, "usage: akhal parse <r/GFA> [--basic] [--no-links] [--no-overlaps] [--no-paths] [--no-ranks] [--verbose] [--footprint]");
}

// P lines sharing a name
static void check_path_names(const gfa_t *g) {
    pathmap_t *h = pathmap_init();
    if (!h) return;
    long dup = 0;
    for (int32_t k = 0; k < gfa_n_path(g); k++) {
        int absent;
        khint_t it = pathmap_put(h, gfa_path_name(g, k), &absent);
        if (absent) kh_val(h, it) = k;
        else dup++;
    }
    if (dup > 0) {
        ak_log(AK_LOG_WARN, "parse", "%ld of %d P lines reuse a name already given to another path (%d distinct name(s)); fine for a path written as fragments, otherwise a problem", dup, gfa_n_path(g), (int)kh_size(h));
    }
    pathmap_destroy(h);
}

// `parse` entry point; see cli.h
int cmd_parse(int argc, char **argv) {
    const char *fn = NULL;
    int checks = CHECK_ALL;
    int report = cli_take_report(&argc, argv, 2);

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--basic")) {
            checks = CHECK_BASIC;
        } else if (!strcmp(argv[i], "--no-links")) {
            checks &= ~CHECK_LINKS;
        } else if (!strcmp(argv[i], "--no-overlaps")) {
            checks &= ~CHECK_OVERLAPS;
        } else if (!strcmp(argv[i], "--no-paths")) {
            checks &= ~CHECK_PATHS;
        } else if (!strcmp(argv[i], "--no-ranks")) {
            checks &= ~CHECK_RANKS;
        } else if (argv[i][0] == '-') {
            ak_log(AK_LOG_ERROR, NULL, "unknown option: %s", argv[i]);
            usage();
            return 1;
        } else if (!fn) {
            fn = argv[i];
        } else {
            usage();
            return 1;
        }
    }
    if (!fn) {
        usage();
        return 1;
    }
    if (!ak_ends_with(fn, ".gfa") && !ak_ends_with(fn, ".rgfa")) {
        ak_log(AK_LOG_ERROR, NULL, "expected a .gfa or .rgfa file: %s", fn);
        return 1;
    }
    int is_rgfa = ak_ends_with(fn, ".rgfa");
    if (!is_rgfa) checks &= ~CHECK_RANKS;   // nothing to compare against
    if (checks & CHECK_OVERLAPS) checks |= CHECK_LINKS;

    // read only what the requested checks need: 
    //  1) the link and overlap checks run inside gfa_read() under GFA_VALIDATE, the overlap one only once the sequences are loaded; 
    //  2) the path checks need the steps (in path) and the adjacency;
    //  3) the rank check needs the steps (in path)
    int flags = GFA_SEGS | report;
    if (checks & CHECK_LINKS)    flags |= GFA_LINKS | GFA_VALIDATE;
    if (checks & CHECK_OVERLAPS) flags |= GFA_LINKS | GFA_VALIDATE | GFA_SEQ;
    if (checks & CHECK_PATHS)    flags |= GFA_PATHS | GFA_LINKS | GFA_ARCS;
    if (checks & CHECK_RANKS)    flags |= GFA_PATHS;

    gfa_t *g = gfa_read(fn, flags);
    if (!g) return 1;

    long issues = 0;

    // consecutive path steps must be joined by a link, on the strand the path walks them: "2-,1-" is a valid walk of "L 1 + 2 +"
    if (checks & CHECK_PATHS) {
        for (int32_t k = 0; k < gfa_n_path(g); k++) {
            const uint32_t *segs;
            int n = gfa_path_segs(g, k, &segs);
            const char *ori = g->path_ori + g->path_off[k];
            for (int i = 1; i < n; i++) {
                uint32_t a = segs[i - 1], b = segs[i];
                if (a == GFA_NIL || b == GFA_NIL) continue;
                if (!gfa_has_link(g, (int32_t)a, ori[i - 1], (int32_t)b, ori[i])) {
                    ak_log(AK_LOG_WARN, "parse", "no link %llu%c -> %llu%c present in path %s",
                           (unsigned long long)gfa_seg_at(g, (int32_t)a)->id, ori[i - 1],
                           (unsigned long long)gfa_seg_at(g, (int32_t)b)->id, ori[i], gfa_path_name(g, k));
                    issues++;
                }
            }
        }
        check_path_names(g);
    }

    // rGFA: rank-0 segment count vs path occurrence count
    if (checks & CHECK_RANKS) {
        uint64_t n_ref = 0;
        for (int32_t i = 0; i < gfa_n_seg(g); i++) {
            if (gfa_seg_at(g, i)->rank == 0) {
                n_ref++;
            }
        }
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
    ak_log(AK_LOG_WARN, "parse", "%ld issue(s) found (unknown-id and overlap warnings, if any, logged above)", issues);
    return 1;
}
