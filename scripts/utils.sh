#!/usr/bin/env bash
#
# Shared helpers for the benchmark scripts: find a timer, run a command under it, record a row, summarize the rows.
# Source it, don't run it:
#
#   . "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/utils.sh"
#
# Nothing here runs on sourcing; every helper reads its settings from these globals at call time, so they can be set
# before or after the source:
#
#   TOOLS      the tools asked for, space separated (want)
#   ONLY       a regex; only the tasks it matches run (skip, measure, heading); empty for all
#   DRY_RUN    1 to print each command instead of running it (measure)
#   KEEP       1 to keep outputs, 0 to delete them once measured (the default for measure's and done_with's keep)
#   REPEATS    runs per measurement; the median is reported (measure)
#   TIMEOUT    seconds a run may take; 0 for no limit (find_timeout)
#   WORK       scratch directory: the timer's report goes here, and done_with deletes nothing outside it
#   LOGS       directory for <task>.<tool>.log (the run's output) and <task>.<tool>.time (the timer's report)
#   RESULTS    the results table (results_init, row, summary)
#   KEPT       set by done_with: the files it was told to keep, which clean_work leaves in place
#   REUSED     outputs a later row reads back, space separated with a space at each end; measure leaves them for
#              done_with to delete once their last reader is done. Empty or unset, every output goes once measured
#
# present <tool> says whether a tool can run here, and defaults to "is it on PATH". A script whose tool is a path
# rather than a name redefines it after sourcing.

# tools

have() {
    command -v "$1" >/dev/null 2>&1
}

present() {
    have "$1"
}

# was this tool asked for? 
want() {
    case "$1" in all) return 0 ;; esac
    case " $TOOLS " in *" ${1%-prep} "*) return 0 ;; esac
    return 1
}

# is this tool's subcommand actually there? versions differ, and a missing one should read as "unsupported", not as a failed benchmark
has_sub() {
    local tool=$1 sub=$2
    have "$tool" || return 1
    "$tool" --help 2>&1 | grep -qw -- "$sub" && return 0
    "$tool" 2>&1 | grep -qw -- "$sub"
}

# timing

# GNU time reports peak RSS, BSD time reports it differently, and a shell builtin reports none at all -
# so the mode is settled once, here. Sets TIME_BIN and TIME_MODE (gnu, bsd or none) and says which it found
find_timer() {
    TIME_BIN=""
    TIME_MODE="none"
    if command -v gtime >/dev/null 2>&1 && gtime -v true 2>/dev/null >/dev/null; then
        TIME_BIN="$(command -v gtime)"; TIME_MODE="gnu"
    elif [ -x /usr/bin/time ] && /usr/bin/time -v true 2>/dev/null >/dev/null; then
        TIME_BIN="/usr/bin/time"; TIME_MODE="gnu"
    elif [ -x /usr/bin/time ] && /usr/bin/time -l true 2>/dev/null >/dev/null; then
        TIME_BIN="/usr/bin/time"; TIME_MODE="bsd"
    fi

    case "$TIME_MODE" in
        gnu) echo "timer:  $TIME_BIN -v (wall clock and peak RSS)" ;;
        bsd) echo "timer:  $TIME_BIN -l (wall clock and peak RSS)" ;;
        none) echo "timer:  shell only - install GNU time for peak RSS (apt install time / brew install gnu-time)" ;;
    esac
}

# a run that outlives the budget is killed - the whole process group, so a pipeline goes with it - and its row says
# "timeout". GNU timeout is coreutils on linux and gtimeout from brew's coreutils on a mac. Sets TIMEOUT_BIN
find_timeout() {
    TIMEOUT_BIN=""
    if [ "${TIMEOUT:-0}" -gt 0 ] 2>/dev/null; then
        if have timeout; then TIMEOUT_BIN="timeout"
        elif have gtimeout; then TIMEOUT_BIN="gtimeout"
        fi
        if [ -n "$TIMEOUT_BIN" ]; then
            echo "budget: $TIMEOUT s per run ($TIMEOUT_BIN)"
        else
            echo "budget: none - no timeout binary found (coreutils); runs are not limited"
        fi
    fi
}

# results

# the table's header; row() writes the lines under it, in the same order
results_init() {
    printf 'task\ttool\tstatus\texit\twall_min\tmax_rss_gb\tout_bytes\tnote\tcommand\n' > "$RESULTS"
}

# median of the numbers on stdin
median() {
    sort -g | awk '{v[NR]=$1} END{ if(NR==0){print "NA"} else if(NR%2){printf "%.3f\n", v[(NR+1)/2]} else {printf "%.3f\n", (v[NR/2]+v[NR/2+1])/2} }'
}

# GNU time writes h:mm:ss or m:ss.ss; seconds is what a table wants
to_seconds() {
    awk -F: '{ if (NF==3) printf "%.3f\n", $1*3600+$2*60+$3; else if (NF==2) printf "%.3f\n", $1*60+$2; else printf "%.3f\n", $1 }'
}

row() {  # task tool status exit wall rss bytes note command
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "$RESULTS"
}

skip() {  # task tool reason
    want "$2" || return 0
    [ -n "$ONLY" ] && ! printf '%s' "$1" | grep -Eq "$ONLY" && return 0
    row "$1" "$2" "skipped" "NA" "NA" "NA" "NA" "$3" ""
    printf '  %-10s %-9s %s\n' "$1" "$2" "skipped: $3"
}

# seconds with a fraction where the platform offers one
now() {
    local t
    t=$(date +%s.%N 2>/dev/null)
    case "$t" in *N*|"") date +%s ;; *) printf '%s' "$t" ;; esac
}

# measure <task> <tool> <command> [note] [output file] [accepted exit codes] [keep]
#
# Returns 0 when the row is ok, or when it was filtered out or only printed; non-zero when the run failed or timed out.
#
# `accepted exit codes` defaults to 0. 
# The compare commands answer 1 for "the two files differ", which is a result rather than a failure, so those rows pass "0 1" and the real code is kept in its own column
measure() {
    local task=$1 tool=$2 cmd=$3 note=${4:-} out=${5:-} accept=${6:-0} keep=${7:-$KEEP}

    want "$tool" || return 0
    if [ -n "$ONLY" ] && ! printf '%s' "$task" | grep -Eq "$ONLY"; then
        return 0
    fi
    if [ "$DRY_RUN" = 1 ]; then
        printf '  %-10s %-9s %s\n' "$task" "$tool" "$cmd"
        return 0
    fi

    local log="$LOGS/$task.$tool.log"
    local timelog="$LOGS/$task.$tool.time"
    local tf="$WORK/.time.$$"
    local walls="" rsss="" status="ok" rc=0
    : > "$log"
    : > "$timelog"

    # what actually runs: the command under the budget, if there is one.
    # -k gives a stubborn process a minute after TERM before it gets KILL
    local -a run
    if [ -n "$TIMEOUT_BIN" ]; then
        run=("$TIMEOUT_BIN" -k 60 "$TIMEOUT" bash -c "$cmd")
    else
        run=(bash -c "$cmd")
    fi

    local i w r
    for i in $(seq 1 "$REPEATS"); do
        printf '== run %s of %s: %s\n' "$i" "$REPEATS" "$cmd" >> "$timelog"
        w="" r="NA"
        case "$TIME_MODE" in
            gnu)
                "$TIME_BIN" -v -o "$tf" "${run[@]}" >>"$log" 2>>"$log"
                rc=$?
                cat "$tf" >> "$timelog"
                w=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$tf" | to_seconds | awk '{printf "%.4f\n", $1/60}')
                r=$(awk '/Maximum resident set size/{printf "%.4f\n", $NF/1024/1024}' "$tf")
                ;;
            bsd)
                "$TIME_BIN" -l "${run[@]}" >>"$log" 2>"$tf"
                rc=$?
                cat "$tf" >> "$log"
                cat "$tf" >> "$timelog"
                w=$(awk '/^ *real/{print $1}' "$tf" | tail -1 | awk '{printf "%.4f\n", $1/60}')
                r=$(awk '/maximum resident set size/{printf "%.4f\n", $1/1073741824}' "$tf")
                ;;
            *)
                local t0 t1
                t0=$(now)
                "${run[@]}" >>"$log" 2>>"$log"
                rc=$?
                t1=$(now)
                w=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.4f\n", (b-a)/60}')
                printf 'no timer installed; wall clock taken from the shell: %s min\n' "$w" >> "$timelog"
                ;;
        esac
        printf 'exit status: %s\n\n' "$rc" >> "$timelog"
        walls="$walls$w"
        rsss="$rsss$r"
        # 124 is timeout's own code for a run it had to stop (137 if it
        # needed KILL); whatever the run wrote by then is not a result
        if [ -n "$TIMEOUT_BIN" ] && { [ "$rc" = 124 ] || [ "$rc" = 137 ]; }; then
            status="timeout"
            printf 'killed after %s s\n\n' "$TIMEOUT" >> "$timelog"
            [ -n "$out" ] && rm -f "$out"
            break
        fi
        if ! printf '%s' " $accept " | grep -q " $rc "; then
            status="failed"
            break
        fi
    done
    rm -f "$tf"

    local wall rss bytes
    wall=$(printf '%s' "$walls" | grep -v '^$' | median)
    rss=$(printf '%s' "$rsss" | grep -v '^$' | grep -v NA | median)
    [ -z "$rss" ] && rss="NA"
    if [ -n "$out" ] && [ -f "$out" ]; then
        bytes=$(wc -c < "$out" | tr -d ' ')
    else
        bytes="NA"
    fi
    case "$REUSED" in *" $out "*) ;; *) done_with "$keep" "$out" ;; esac

    row "$task" "$tool" "$status" "$rc" "$wall" "$rss" "$bytes" "$note" "$cmd"
    if [ "$status" = "timeout" ]; then
        local timeout_min
        timeout_min=$(awk -v t="$TIMEOUT" 'BEGIN{printf "%.2f", t/60}')
        printf '  %-10s %-9s %8s   %8s GB  timeout: killed after %s m\n' "$task" "$tool" ">${timeout_min}m" "$rss" "$timeout_min"
    else
        printf '  %-10s %-9s %8sm %8sGB  %s\n' "$task" "$tool" "$wall" "$rss" "$status(exit $rc)"
    fi
    [ "$status" != "ok" ] && printf '             see %s\n' "$log"
    # 0 when the run succeeded (or was only printed), so a caller can act on a failure: measure ... || ...
    [ "$status" = "ok" ]
}

# done_with <keep> <file>... - scratch that nothing further reads, deleted unless <keep> is 1.
# A file kept this way is remembered in KEPT, so clean_work leaves it too. Only files under $WORK are touched -
# an input stood in for a failed conversion (OG and VGP fall back to the GFA) is neither deleted nor listed
done_with() {
    local keep=$1 f
    shift
    for f in "$@"; do
        case "$f" in "$WORK"/*) ;; *) continue ;; esac
        if [ "$keep" = 1 ]; then
            KEPT="${KEPT:- }$f "
        else
            rm -rf "$f"
        fi
    done
}

# the end of a run: whatever is left in scratch goes, except the files done_with was told to keep; with KEEP=1 it all stays
clean_work() {
    if [ "$KEEP" = 1 ]; then
        echo "scratch: $WORK/"
        return 0
    fi
    local f
    for f in "$WORK"/* "$WORK"/.[!.]*; do
        [ -e "$f" ] || continue
        case "${KEPT:-}" in *" $f "*) ;; *) rm -rf "$f" ;; esac
    done
    rmdir "$WORK" 2>/dev/null || echo "kept:   ${KEPT% }"
}

# fresh <output> <input> - the output is there, non-empty and newer than the input it was made from, so building it
# again would give the same file. A later edit of the input makes it stale
fresh() {
    [ -s "$1" ] && [ "$1" -nt "$2" ]
}

# measure a competitor's row
try() {  # task tool cmd [note] [out] [accept] [keep]
    want "$2" || return 0
    if ! present "${2%-prep}"; then
        skip "$1" "$2" "not installed"
        return 0
    fi
    measure "$@"
}

heading() {
    [ -n "$ONLY" ] && ! printf '%s' "$2" | grep -Eq "$ONLY" && return 0
    printf '\n%s\n' "$1"
}

# summary <baseline tool> - the results table as a report, each tool's wall time as a multiple of the baseline's
summary() {
    local base_tool=$1
    awk -F'\t' -v base_tool="$base_tool" '
    NR == 1 { next }
    {
        task=$1; tool=$2; status=$3; wall=$5; rss=$6; note=$8
        if (!(task in seen)) { order[++n]=task; seen[task]=1 }
        key=task SUBSEP tool
        st[key]=status; w[key]=wall; r[key]=rss; nt[key]=note
        tools[task]=tools[task] " " tool
        if (tool==base_tool && status=="ok") base[task]=wall
    }
    END {
        fmt = "%-11s %-13s %10s %11s  %-9s %s\n"
        printf fmt, "task", "tool", "wall (m)", "peak (GB)", "vs " base_tool, "note"
        printf fmt, "-----------", "-------------", "----------", "-----------", "---------", "----"
        for (i=1; i<=n; i++) {
            task=order[i]
            c=split(tools[task], tl, " ")
            for (j=1; j<=c; j++) {
                tool=tl[j]; if (tool=="") continue
                key=task SUBSEP tool
                # a prep row is what a tool needs before it can start, not a
                # competitor for the same work, so it gets no ratio
                ratio="-"
                if (tool!=base_tool && tool !~ /-prep$/ && st[key]=="ok" && (task in base) && base[task]+0 > 0)
                    ratio=sprintf("%.2fx", w[key]/base[task])
                note=nt[key]
                if (length(note) > 52) note=substr(note, 1, 49) "..."
                if (st[key]=="skipped") {
                    printf fmt, task, tool, "-", "-", "skipped", note
                } else if (st[key]!="ok") {
                    printf fmt, task, tool, w[key], r[key], st[key], note
                } else {
                    printf fmt, task, tool, w[key], r[key], ratio, note
                }
            }
            printf "\n"
        }
    }' "$RESULTS"

    # /usr/bin/time counts in hundredths of a second, so anything this quick is being reported by the clock rather than measured by it
    if awk -F'\t' 'NR>1 && $3=="ok" && $5*60 < 0.02 {found=1} END{exit !found}' "$RESULTS"; then
        echo "note: some rows finished under 0.02 s, which is the timer's own resolution -"
        echo "      those numbers say \"too fast to measure here\", not what they literally read."
        echo
    fi
}
