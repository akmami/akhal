#include "akhal/annot.h"
#include "akhal/error.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>

// `annotget` entry point; see cli.h
int cmd_annotget(int argc, char **argv) {
    if (argc < 3) {
        ak_log(AK_LOG_ERROR, NULL, "usage: akhal annotget <file.annot> [node-id ...]");
        return 1;
    }

    annot_t *an = annot_read(argv[2]);
    if (!an) return 1;

    if (argc == 3) {   // no ids: dump everything
        for (int64_t i = 0; i < annot_n(an); i++) {
            const annot_rec_t *r = annot_at(an, i);
            const char *info = annot_info_at(an, i);
            printf("%llu\t%s%s%s\n", (unsigned long long)r->id,
                   annot_kind_name(r->kind),
                   info ? "\t" : "", info ? info : "");
        }
    } else {
        for (int i = 3; i < argc; i++) {
            char *end;
            uint64_t id = strtoull(argv[i], &end, 10);
            if (end == argv[i] || *end) {
                ak_log(AK_LOG_WARN, NULL, "not a node id: %s", argv[i]);
                continue;
            }
            const char *info;
            int kind = annot_get(an, id, &info);
            printf("%llu\t%s%s%s\n", (unsigned long long)id,
                   annot_kind_name(kind),
                   info ? "\t" : "", info ? info : "");
        }
    }

    annot_destroy(an);
    return 0;
}
