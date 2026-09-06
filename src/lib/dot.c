#include "akhal/dot.h"
#include "akhal/error.h"

// shades for the ranks: the reference backbone, then everything hanging off it
#define DOT_FILL_REF "#cfe3f7"
#define DOT_FILL_ALT "#f7e2cf"

// a label is DOT-quoted, so the two characters that would end it early are
// escaped. Sequences never carry them; a stray byte in a malformed file might
static void put_escaped(FILE *out, const char *s, uint32_t len) {
    for (uint32_t i = 0; i < len; i++) {
        if (s[i] == '"' || s[i] == '\\') fputc('\\', out);
        fputc(s[i], out);
    }
}

// the id, and either the bases or how many of them there are
static void put_label(FILE *out, const gfa_seg_t *s, int flags) {
    fprintf(out, "%llu", (unsigned long long)s->id);

    if (flags & DOT_IDS) return;
    if (s->seq && s->len && s->len <= DOT_SEQ_MAX) {
        fputs("\\n", out);
        put_escaped(out, s->seq, s->len);
    } else if (s->len) {
        fprintf(out, "\\n%u bp", s->len);
    }
}

// one box per segment, shaded by rank where the file ranks itself
static void put_segs(const gfa_t *g, FILE *out, int flags) {
    for (int32_t i = 0; i < gfa_n_seg(g); i++) {
        const gfa_seg_t *s = gfa_seg_at(g, i);

        fprintf(out, "  %llu [label=\"", (unsigned long long)s->id);
        put_label(out, s, flags);
        fputc('"', out);

        if (g->has_sr && s->rank >= 0) {
            fprintf(out, ", style=\"rounded,filled\", fillcolor=\"%s\"", s->rank == 0 ? DOT_FILL_REF : DOT_FILL_ALT);
        }
        fputs("];\n", out);
    }
}

// one arrow per link, labelled only where it has something to say
static void put_links(const gfa_t *g, FILE *out) {
    for (int32_t i = 0; i < gfa_n_link(g); i++) {
        const gfa_link_t *e = gfa_link_at(g, i);
        int flip = e->from_orient != '+' || e->to_orient != '+';

        fprintf(out, "  %llu -> %llu",
                (unsigned long long)gfa_seg_at(g, (int32_t)e->v)->id,
                (unsigned long long)gfa_seg_at(g, (int32_t)e->w)->id);

        if (flip || e->overlap) {
            fputs(" [label=\"", out);
            if (flip) {
                fprintf(out, "%c/%c", e->from_orient, e->to_orient);
            }
            if (e->overlap) {
                fprintf(out, "%s%uM", flip ? " " : "", e->overlap);
            }
            fputs("\"]", out);
        }
        fputs(";\n", out);
    }
}

// write a graph as a digraph; see akhal/dot.h
int dot_write(const gfa_t *g, FILE *out, int flags) {
    if (gfa_n_seg(g) > DOT_BIG) {
        ak_log(AK_LOG_WARN, "dot", "%d nodes is a lot to lay out; graphviz may take a very long time over it", gfa_n_seg(g));
    }

    fputs("digraph gfa {\n", out);
    fputs("  rankdir=LR;\n", out);
    fputs("  node [shape=box, style=rounded, fontname=\"monospace\"];\n", out);
    fputs("  edge [fontname=\"monospace\", fontsize=9];\n\n", out);

    put_segs(g, out, flags);
    if (gfa_n_link(g)) {
        fputc('\n', out);
        put_links(g, out);
    }

    fputs("}\n", out);

    if (ferror(out)) {
        ak_log(AK_LOG_ERROR, "dot", "write failed");
        return AK_EIO;
    }
    return AK_OK;
}
