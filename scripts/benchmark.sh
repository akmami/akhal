#!/usr/bin/env bash
#
# Usage:
#   ./benchmark.sh [options]
#
#   --config FILE     read settings from FILE (default: ./benchmark.conf, if it exists)
#   --data DIR        directory holding the inputs (default: data)
#   --out DIR         where results, logs and scratch go (default: out)
#   --threads N       threads for the tools that take them (default: 1)
#   --repeats N       runs per measurement, median reported (default: 1)
#   --timeout SECS    give up on a run after this long (default: 86400, one day; 0 = never)
#   --only REGEX      run only the tasks whose name matches
#   --tools LIST      which tools to measure, comma or space separated
#   --install         fetch the missing competitors into <out>/bin first
#   --keep            keep the scratch files instead of deleting them at the end
#   --dry-run         print what would be run, measure nothing
#   -h, --help        this text
#
# The defaults assume the layout described in the README:
#
#   data/human_v38.gfa    the graph
#   data/human_v38.vg     the same graph in vg's native format
#   data/human_v38.gaf    alignments against it
#   data/human_v38.2.gaf  a second, different set of alignments
#   data/human_v38.fa     the reads those alignments came from (gaf2sam only)
#

# the helpers every benchmark script shares: timing, rows, the summary
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/utils.sh"

# defaults

AKHAL="./akhal"
DATA_DIR="data"
GFA=""
VG_FILE=""
GAF=""
GAF_B=""
READS=""
REF_CSV=""
REF_P=""
REF_ALL=""
GAFTOOLS_REF=""

OUTDIR="out"
THREADS=1
REPEATS=1
TIMEOUT=86400

TOOLS="all"

CONFIG=""
DO_INSTALL=0
KEEP=0
DRY_RUN=0
ONLY=""

# argument parsing

usage() {
    # the header comment block is the help text; print it up to the first line that is not a comment
    awk 'NR>1 && !/^#/{exit} NR>1{sub(/^# ?/, ""); print}' "$0"
    exit "${1:-0}"
}

# the config is read before the options, so that the command line overrides it
for ((i = 1; i < $#; i++)); do
    j=$((i + 1))
    [ "${!i}" = "--config" ] && CONFIG=${!j}
done

# the config file is optional: without one the defaults above stand
if [ -z "$CONFIG" ] && [ -f "./benchmark.conf" ]; then
    CONFIG="./benchmark.conf"
fi
if [ -n "$CONFIG" ]; then
    if [ ! -f "$CONFIG" ]; then
        echo "config file not found: $CONFIG" >&2
        exit 2
    fi
    # shellcheck disable=SC1090
    . "$CONFIG"
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --config|--data|--out|--threads|--repeats|--timeout|--only|--tools)
            [ $# -ge 2 ] || { echo "$1 needs a value" >&2; usage 2; } ;;
    esac
    case "$1" in
        --config)  shift 2 ;;
        --data)    DATA_DIR=$2; shift 2 ;;
        --out)     OUTDIR=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --repeats) REPEATS=$2; shift 2 ;;
        --timeout) TIMEOUT=$2; shift 2 ;;
        --only)    ONLY=$2; shift 2 ;;
        --tools)   TOOLS=$2; shift 2 ;;
        --install) DO_INSTALL=1; shift ;;
        --keep)    KEEP=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage 0 ;;
        *) echo "unknown option: $1" >&2; usage 2 ;;
    esac
done

if [ -n "$CONFIG" ]; then
    echo "config: $CONFIG"
else
    echo "config: none, using defaults"
fi

# which tools are in play
ALL_TOOLS="akhal gfatools odgi vg gaftools"
case "$(printf '%s' "$TOOLS" | tr 'A-Z,' 'a-z ')" in
    *all*) TOOLS="$ALL_TOOLS" ;;
    *)     TOOLS="$(printf '%s' "$TOOLS" | tr 'A-Z,' 'a-z ' | tr -s ' ')" ;;
esac
for t in $TOOLS; do
    case " $ALL_TOOLS " in
        *" $t "*) ;;
        *) echo "unknown tool: $t (choose from $ALL_TOOLS, or all)" >&2; exit 2 ;;
    esac
done
if [ -z "$(printf '%s' "$TOOLS" | tr -d ' ')" ]; then
    echo "no tools selected: --tools takes one or more of $ALL_TOOLS, or all" >&2
    exit 2
fi

# anything the config did not set falls back to the documented layout
: "${GFA:=$DATA_DIR/human_v38.gfa}"
: "${VG_FILE:=$DATA_DIR/human_v38.vg}"
: "${GAF:=$DATA_DIR/human_v38.gaf}"
: "${GAF_B:=$DATA_DIR/human_v38.2.gaf}"
: "${READS:=$DATA_DIR/human_v38.fa}"

WORK="$OUTDIR/work"
LOGS="$OUTDIR/logs"
BINDIR="$OUTDIR/bin"
RESULTS="$OUTDIR/results.tsv"
ENVFILE="$OUTDIR/env.txt"

OG="$WORK/graph.og"
VGP="$WORK/graph.packed.vg"
SORTED="$WORK/akhal.sorted.gfa"
ANNOT="$WORK/akhal.annot"
REUSED=" $OG $VGP $SORTED $ANNOT "

STARTED=$(date +%s)
mkdir -p "$WORK" "$LOGS" "$BINDIR"
PATH="$(cd "$BINDIR" && pwd):$PATH"
export PATH

# how to measure

find_timer
find_timeout

# is it actually runnable here? akhal is a path out of the config rather than a name on PATH, so it answers differently from the rest
present() {
    case "$1" in
        akhal) [ -x "$AKHAL" ] || have "$AKHAL" ;;
        *)     have "$1" ;;
    esac
}

# installing the competitors

install_tools() {
    local d
    d="$(cd "$BINDIR" && pwd)"

    if want gfatools && ! have gfatools; then
        echo "installing gfatools ..."
        ( set -e
          tmp=$(mktemp -d)
          git clone --depth 1 https://github.com/lh3/gfatools "$tmp/gfatools" >/dev/null 2>&1
          make -C "$tmp/gfatools" -j "$THREADS" >/dev/null 2>&1
          cp "$tmp/gfatools/gfatools" "$d/"
          rm -rf "$tmp" ) || echo "  gfatools: build failed (needs git, make, gcc, zlib headers)"
    fi

    if want vg && ! have vg; then
        echo "installing vg (static release binary) ..."
        # vg publishes one static linux x86_64 binary per release, named `vg`
        if curl -fsSL https://github.com/vgteam/vg/releases/latest/download/vg -o "$d/vg"; then
            chmod +x "$d/vg"
        else
            rm -f "$d/vg"
            echo "  vg: download failed (the release binary is linux x86_64 only)"
        fi
    fi

    if want gaftools && ! have gaftools; then
        echo "installing gaftools (pip, into $OUTDIR/venv) ..."
        ( set -e
          python3 -m venv "$OUTDIR/venv" >/dev/null 2>&1
          "$OUTDIR/venv/bin/pip" install --quiet --upgrade pip >/dev/null 2>&1
          "$OUTDIR/venv/bin/pip" install --quiet gaftools >/dev/null 2>&1
          ln -sf "$(cd "$OUTDIR/venv/bin" && pwd)/gaftools" "$d/gaftools" ) \
          || echo "  gaftools: pip install failed (needs python >= 3.9)"
    fi

    if want odgi && ! have odgi; then
        # odgi is a cmake build with heavy dependencies, so the package manager
        # is the only sane route here
        if have mamba; then
            echo "installing odgi (mamba) ..."
            mamba install -y -c conda-forge -c bioconda odgi >/dev/null 2>&1 || echo "  odgi: mamba install failed"
        elif have conda; then
            echo "installing odgi (conda) ..."
            conda install -y -c conda-forge -c bioconda odgi >/dev/null 2>&1 || echo "  odgi: conda install failed"
        else
            echo "  odgi: not installed and no conda/mamba found - conda install -c bioconda odgi"
        fi
    fi
    echo
}

[ "$DO_INSTALL" = 1 ] && install_tools

# what is here

AVAIL="" ABSENT=""
for t in $TOOLS; do
    if present "$t"; then AVAIL="$AVAIL $t"; else ABSENT="$ABSENT $t"; fi
done
want akhal && echo "akhal:  $AKHAL"
echo "tools: ${AVAIL:- none of the requested tools are installed}"
[ -n "$ABSENT" ] && echo "        asked for but missing:$ABSENT (their rows are marked skipped)"
echo

# only fatal when akhal is one of the tools being measured
if want akhal && ! present akhal; then
    echo "akhal not found at '$AKHAL' - build it with make, or set AKHAL in the config" >&2
    exit 2
fi

# A missing inputs are not fatal: the tasks that need it are skipped
missing=""
for f in "$GFA" "$VG_FILE" "$GAF"; do
    [ -f "$f" ] || missing="$missing $f"
done
[ -n "$missing" ] && echo "note: missing input(s):$missing (their tasks will be skipped)"

# the backbones for vcf and gfa2rgfa, as the config set them: REF_CSV and REF_P (the same names as vg's -p list) fall back to
# the first P line's name, which is what akhal itself defaults to, and REF_ALL to REF_CSV
if [ -z "$REF_CSV" ] && [ -f "$GFA" ]; then
    REF_CSV=$(awk '$1=="P"{print $2; exit}' "$GFA" 2>/dev/null)
    [ -n "$REF_CSV" ] && REF_P=" -p '$REF_CSV'"
fi
: "${REF_ALL:=$REF_CSV}"
if [ -n "$REF_CSV" ] && [ -z "$REF_P" ]; then
    echo "REF_CSV is set but REF_P is not: build REF_P from REF_CSV in the config (see benchmark.conf.example)" >&2
    exit 2
fi
[ -n "$REF_CSV" ] && echo "backbones: $REF_CSV (akhal: $REF_ALL; vg:$REF_P)"

ref_note() {  # a --ref value, short enough for the results table
    case "$1" in
        all) echo "backbone: every path" ;;
        *,*) echo "backbone: ${1%%,*} and $(printf '%s' "$1" | tr -cd , | wc -c | tr -d ' ') more" ;;
        *)   echo "backbone: $1" ;;
    esac
}

# a second GAF is what makes `compare gaf` interesting; without one the comparison still runs, against a copy of the first
GAF_B_NOTE=""
if [ ! -f "$GAF_B" ] && [ -f "$GAF" ]; then
    GAF_B="$GAF"
    GAF_B_NOTE="second GAF missing: compared against itself"
    echo "note: $GAF_B_NOTE"
fi
echo

# the machine, and what every tool answers to --version, so a results file can be read a year from now
{
    echo "date:    $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "host:    $(uname -a)"
    if [ -r /proc/cpuinfo ]; then
        echo "cpu:     $(awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo)"
        echo "cores:   $(grep -c ^processor /proc/cpuinfo)"
        echo "memory:  $(awk '/MemTotal/{printf "%.1f GB\n", $2/1048576}' /proc/meminfo)"
    elif have sysctl; then
        echo "cpu:     $(sysctl -n machdep.cpu.brand_string 2>/dev/null)"
        echo "cores:   $(sysctl -n hw.ncpu 2>/dev/null)"
        echo "memory:  $(sysctl -n hw.memsize 2>/dev/null | awk '{printf "%.1f GB\n", $1/1073741824}')"
    fi
    echo "threads: $THREADS"
    echo "repeats: $REPEATS"
    echo "timeout: ${TIMEOUT}s"
    echo "tools:   $TOOLS"
    echo
    want akhal    && echo "akhal:    $("$AKHAL" --version 2>&1 | head -1)"
    want gfatools && have gfatools && echo "gfatools: $(gfatools version 2>&1 | head -1)"
    want odgi     && have odgi     && echo "odgi:     $(odgi version 2>&1 | head -1)"
    want vg       && have vg       && echo "vg:       $(vg version 2>&1 | head -1)"
    want gaftools && have gaftools && echo "gaftools: $(gaftools --version 2>&1 | head -1)"
    echo
    echo "inputs:"
    for f in "$GFA" "$VG_FILE" "$GAF" "$GAF_B" "$READS"; do
        if [ -f "$f" ]; then
            echo "  $f  $(wc -c < "$f" | tr -d ' ') bytes"
        fi
    done
} > "$ENVFILE" 2>/dev/null

results_init

# preparation
#
# odgi and vg both prefer their own on-disk format. Neither conversion is
# hidden: it is the price of every odgi/vg row that follows, and it is measured
# like everything else

# the only rows here belong to odgi and vg, so with neither asked for the
# section does not exist at all. A conversion an earlier run left behind is used as it is, as long as the GFA has not
# changed since; one that fails is removed, so a later run builds it again rather than picking up a broken file
if want odgi || want vg; then
    heading "== prep: the formats the other tools want ==" "prep"
    if [ -f "$GFA" ]; then
        if fresh "$OG" "$GFA"; then
            skip "prep" "odgi" "reusing $OG from an earlier run"
        else
            try "prep" "odgi" "odgi build -g '$GFA' -o '$OG' -t $THREADS" "GFA -> .og, needed by every odgi row" "$OG" || done_with 0 "$OG"
        fi
        if fresh "$VGP" "$GFA"; then
            skip "prep" "vg" "reusing $VGP from an earlier run"
        else
            try "prep" "vg" "vg convert -g -p '$GFA' > '$VGP'" "GFA -> packed graph, needed by every vg row" "$VGP" || done_with 0 "$VGP"
        fi
    else
        skip "prep" "all" "no GFA at $GFA"
    fi
fi

# odgi and vg read a GFA directly too, just slower; if the conversion failed,
# fall back to that rather than dropping the tool from the whole benchmark
[ -s "$OG" ]  || { done_with 0 "$OG";  OG="$GFA"; }
[ -s "$VGP" ] || { done_with 0 "$VGP"; VGP="$GFA"; }

# 1. stats

heading "== stats: count nodes, edges and sequence ==" "stats"
if [ -f "$GFA" ]; then
    measure "stats" "akhal"    "$AKHAL stats '$GFA'"
    try "stats" "gfatools" "gfatools stat '$GFA'"
    try "stats" "odgi"     "odgi stats -i '$OG' -S -t $THREADS" "on the prebuilt .og"
    try "stats" "vg"       "vg stats -z -l '$VGP'" "on the converted graph"
else
    skip "stats" "all" "no GFA at $GFA"
fi

# 2. validate

heading "== validate: is the graph well formed ==" "validate"
if [ -f "$GFA" ]; then
    measure "validate" "akhal" "$AKHAL parse '$GFA'"
    skip "validate" "gfatools" "no equivalent subcommand"
    try "validate" "odgi" "odgi validate -i '$OG' -t $THREADS" "on the prebuilt .og"
    try "validate" "vg"   "vg validate '$VGP'" "on the converted graph"
else
    skip "validate" "all" "no GFA at $GFA"
fi

# 3. sort

heading "== sort: topological order, ids renumbered ==" "sort"
if [ -f "$GFA" ]; then
    measure "sort" "akhal" "$AKHAL sort '$GFA' '$SORTED'" "" "$SORTED"
    skip "sort" "gfatools" "no equivalent subcommand"
    try "sort" "odgi" "odgi sort -i '$OG' -o '$WORK/odgi.sorted.og' -z -t $THREADS" "depth-first topological sort" "$WORK/odgi.sorted.og"
    try "sort" "vg"   "vg ids -s '$VGP' > '$WORK/vg.sorted.vg'" "generalized topological order" "$WORK/vg.sorted.vg"
else
    skip "sort" "all" "no GFA at $GFA"
fi

# 4. compact / unchop

heading "== compact: fold non-branching runs into one node ==" "compact"
if [ -f "$GFA" ]; then
    measure "compact" "akhal" "$AKHAL compact '$GFA' '$WORK/akhal.compact.gfa'" "" "$WORK/akhal.compact.gfa"
    skip "compact" "gfatools" "asm -u builds unitigs, which is a different operation"
    try "compact" "odgi" "odgi unchop -i '$OG' -o '$WORK/odgi.unchop.og' -t $THREADS" "unchop" "$WORK/odgi.unchop.og"
    try "compact" "vg"   "vg mod -u '$VGP' > '$WORK/vg.unchop.vg'" "mod -u" "$WORK/vg.unchop.vg"
else
    skip "compact" "all" "no GFA at $GFA"
fi

# 5. paths as FASTA

heading "== gfa2fa: write every path as FASTA ==" "gfa2fa"
if [ -f "$GFA" ]; then
    measure "gfa2fa" "akhal"    "$AKHAL extract fa '$GFA' '$WORK/akhal.paths.fa'" "one record per P line" "$WORK/akhal.paths.fa"
    try "gfa2fa" "gfatools" "gfatools gfa2fa -s '$GFA' > '$WORK/gfatools.paths.fa'" "-s: stable sequences" "$WORK/gfatools.paths.fa"
    try "gfa2fa" "odgi"     "odgi paths -i '$OG' -f -t $THREADS > '$WORK/odgi.paths.fa'" "on the prebuilt .og" "$WORK/odgi.paths.fa"
    try "gfa2fa" "vg"       "vg paths -x '$VGP' -F > '$WORK/vg.paths.fa'" "on the converted graph" "$WORK/vg.paths.fa"
else
    skip "gfa2fa" "all" "no GFA at $GFA"
fi

# 6. variants as VCF

heading "== vcf: variation off the reference backbone ==" "vcf"
if [ -f "$GFA" ] && [ -n "$REF_CSV" ]; then
    measure "vcf" "akhal" "$AKHAL extract vcf '$GFA' '$WORK/akhal.vcf' --ref '$REF_ALL'" "$(ref_note "$REF_ALL")" "$WORK/akhal.vcf"
    skip "vcf" "gfatools" "no equivalent subcommand"
    skip "vcf" "odgi" "no equivalent subcommand"
    try "vcf" "vg" "vg deconstruct $REF_P -t $THREADS '$VGP' > '$WORK/vg.vcf'" "$(ref_note "$REF_CSV")" "$WORK/vg.vcf"
else
    skip "vcf" "all" "no GFA, or no P line to use as the backbone"
fi

# 7. vg -> GFA

heading "== vg2gfa: vg's native format to GFA ==" "vg2gfa"
if [ -f "$VG_FILE" ]; then
    measure "vg2gfa" "akhal" "$AKHAL vg2gfa '$VG_FILE' '$WORK/akhal.fromvg.gfa'" "" "$WORK/akhal.fromvg.gfa"
    try "vg2gfa" "odgi" "odgi view -i '$OG' -g > '$WORK/odgi.view.gfa'" "og -> GFA; odgi cannot read vg's protobuf, so this is the nearest operation" "$WORK/odgi.view.gfa"
    try "vg2gfa" "vg" "vg convert -f '$VG_FILE' > '$WORK/vg.fromvg.gfa'" "" "$WORK/vg.fromvg.gfa"
else
    skip "vg2gfa" "all" "no .vg at $VG_FILE"
fi
done_with 1 "$OG"   # odgi's last row

# 8. GFA -> rGFA

heading "== gfa2rgfa: label a GFA with SN/SO/SR ==" "gfa2rgfa"
if [ -f "$GFA" ]; then
    if [ -n "$REF_CSV" ]; then
        measure "gfa2rgfa" "akhal" "$AKHAL gfa2rgfa '$GFA' '$WORK/akhal.rgfa' --ref '$REF_CSV'" "$(ref_note "$REF_CSV")" "$WORK/akhal.rgfa"
    else
        measure "gfa2rgfa" "akhal" "$AKHAL gfa2rgfa '$GFA' '$WORK/akhal.rgfa'" "" "$WORK/akhal.rgfa"
    fi
    # gaftools names a genome, not paths: the sample its W lines carry, else whatever the H lines say
    if [ -n "$GAFTOOLS_REF" ]; then
        try "gfa2rgfa" "gaftools" "gaftools gfa2rgfa '$GFA' --reference-name '$GAFTOOLS_REF' --output '$WORK/gaftools.rgfa'" "backbone: genome $GAFTOOLS_REF" "$WORK/gaftools.rgfa"
    else
        try "gfa2rgfa" "gaftools" "gaftools gfa2rgfa '$GFA' --output '$WORK/gaftools.rgfa'" "backbone: genome named in the H lines" "$WORK/gaftools.rgfa"
    fi
else
    skip "gfa2rgfa" "all" "no GFA at $GFA"
fi

# 9. GFA -> dot

heading "== gfa2dot: the graph as Graphviz ==" "gfa2dot"
if [ -f "$GFA" ]; then
    measure "gfa2dot" "akhal" "$AKHAL gfa2dot '$GFA' '$WORK/akhal.dot'" "" "$WORK/akhal.dot"
    try "gfa2dot" "vg" "vg view -d '$VGP' > '$WORK/vg.dot'" "on the converted graph" "$WORK/vg.dot"
else
    skip "gfa2dot" "all" "no GFA at $GFA"
fi

# 10. GAF -> SAM

heading "== gaf2sam: alignments as SAM ==" "gaf2sam"
if [ -f "$GFA" ] && [ -f "$GAF" ] && [ -f "$READS" ]; then
    measure "gaf2sam" "akhal" "$AKHAL gaf2sam '$GFA' '$GAF' '$READS' '$WORK/akhal.sam'" "" "$WORK/akhal.sam"
    try "gaf2sam" "vg" "vg surject -x '$VGP' -G -s -t $THREADS '$GAF' > '$WORK/vg.sam'" "surject onto reference paths - related, not identical" "$WORK/vg.sam"
else
    skip "gaf2sam" "all" "needs the graph, the GAF and the reads FASTA ($READS)"
fi
done_with 1 "$VGP"  # vg's last row

# 11. sorting a GAF
#
# akhal sorts a GAF inside `compare gaf` rather than as a command of its own,
# so the honest comparison is the whole comparison against gaftools' sort -
# with gaftools' prerequisite (BO/NO tags, from order_gfa) timed separately

heading "== gafsort: putting a GAF in order ==" "gafsort"
if [ -f "$GAF" ]; then
    measure "gafsort" "akhal" "$AKHAL compare gaf '$GAF' '$GAF'" "sorts both files, then compares them" "" "0 1"
    if want gaftools && have gaftools; then
        ORDERED="$WORK/ordered"
        PREP_RGFA="$WORK/gaftools.prep.rgfa"
        if grep -m1 '^S' "$GFA" 2>/dev/null | grep -q 'SN:Z:'; then
            PREP_CMD="gaftools order_gfa --outdir '$ORDERED' '$GFA'"
            PREP_NOTE="order_gfa: adds the BO/NO tags gaftools sort requires"
        else
            PREP_CMD="gaftools gfa2rgfa '$GFA'${GAFTOOLS_REF:+ --reference-name '$GAFTOOLS_REF'} --output '$PREP_RGFA' && gaftools order_gfa --outdir '$ORDERED' '$PREP_RGFA'"
            PREP_NOTE="gfa2rgfa + order_gfa: the rGFA with BO/NO tags gaftools sort requires"
        fi
        try "gafsort" "gaftools-prep" "$PREP_CMD" "$PREP_NOTE" ""
        ORDERED_GFA=$(ls "$ORDERED"/*.gfa 2>/dev/null | head -1)
        if [ -n "${ORDERED_GFA:-}" ]; then
            try "gafsort" "gaftools" "gaftools sort '$GAF' '$ORDERED_GFA' --outgaf '$WORK/gaftools.sorted.gaf'" "needs the ordered GFA above" "$WORK/gaftools.sorted.gaf"
        else
            skip "gafsort" "gaftools" "order_gfa produced no GFA to sort against"
        fi
        done_with "$KEEP" "$ORDERED" "$PREP_RGFA"
        try "gafstat" "gaftools" "gaftools stat '$GAF' -o '$WORK/gaftools.stat.txt'" "GAF parsing reference point; akhal has no gaf stats command" "$WORK/gaftools.stat.txt"
    else
        skip "gafsort" "gaftools" "gaftools is not installed"
    fi
else
    skip "gafsort" "all" "no GAF at $GAF"
fi

# 12. akhal on its own
#
# Nothing else here compares two graphs or two alignment sets, so these are
# timed rather than raced. The graph comparison doubles as a correctness check:
# a graph and its own sort output must come out identical, and the exit status
# says whether they did

# nothing else here compares two graphs or two alignment sets, so with akhal
# left out there is no section at all
if want akhal; then
    heading "== akhal only: no equivalent in the other tools ==" "compare|rank|annotate"

    # the comparison needs the sorted graph, which --only may have skipped past;
    # make it here rather than dropping the task, untimed since the sort has
    # already been measured on its own
    if [ ! -s "$SORTED" ] && [ -f "$GFA" ] && [ "$DRY_RUN" = 0 ]; then
        "$AKHAL" sort "$GFA" "$SORTED" >/dev/null 2>&1
    fi

    if [ -f "$GFA" ] && [ -s "$SORTED" ]; then
        measure "compare" "akhal" "$AKHAL compare gfa '$GFA' '$SORTED'" "the graph against its own sort output: must be identical" "" "0 1"
        st=$(awk -F'\t' '$1=="compare" && $2=="akhal"{print $4}' "$RESULTS" | tail -1)
        case "$st" in
            "") ;;   # --only filtered the row out, so there is nothing to judge
            0) echo "             correctness: PASS - sorting renumbered every node and the graph still compares identical" ;;
            1) echo "             correctness: FAIL - the sorted graph differs from the original, see $LOGS/compare.akhal.log" ;;
            *) echo "             correctness: could not be established (exit $st)" ;;
        esac
        done_with "$KEEP" "$SORTED"
    else
        skip "compare" "akhal" "needs the GFA and a successful sort"
    fi

    if [ -f "$GAF" ] && [ -f "$GAF_B" ]; then
        measure "comparegaf" "akhal" "$AKHAL compare gaf '$GAF' '$GAF_B'" "${GAF_B_NOTE:-two different alignment sets}" "" "0 1"
        st=$(awk -F'\t' '$1=="comparegaf" && $2=="akhal"{print $4}' "$RESULTS" | tail -1)
        case "$st" in
            0) echo "             the two GAF files place every read the same way" ;;
            1) echo "             the two GAF files differ - the counts are in $LOGS/comparegaf.akhal.log" ;;
        esac
    else
        skip "comparegaf" "akhal" "needs two GAF files"
    fi

    if [ -f "$GFA" ]; then
        # without --fasta, rank takes every P line as backbone; --ref would only name a record in that FASTA
        measure "rank" "akhal" "$AKHAL rank '$GFA' '$WORK/akhal.ranked.gfa'" "backbone: every path" "$WORK/akhal.ranked.gfa"
        measure "annotate" "akhal" "$AKHAL annotate '$GFA' '$ANNOT'" "" "$ANNOT"
        [ -s "$ANNOT" ] && measure "annotget" "akhal" "$AKHAL annotget '$ANNOT' > /dev/null" "dump every node's annotation" ""
        done_with "$KEEP" "$ANNOT"
    fi
fi

# the summary

if [ "$DRY_RUN" = 1 ]; then
    echo
    echo "dry run: nothing was measured"
    exit 0
fi

echo
echo "== summary =="
echo
summary "akhal"

printf 'total: %d s of benchmarking\n\n' "$(( $(date +%s) - STARTED ))"
echo "results: $RESULTS"
echo "logs:    $LOGS/ (*.log is each run's output, *.time the timer's own report)"
echo "machine: $ENVFILE"

clean_work
