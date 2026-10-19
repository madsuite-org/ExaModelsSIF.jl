#!/bin/bash
# Parallel restart-on-chunk / resume-on-crash driver for test/sweep.jl.
# Usage: sweep_driver.sh [NSHARDS]   (default 8)
# Env: SIF_DIR or MASTSIF names the collection; SWEEP_ONLY=<file> restricts to
# names; JULIA overrides the julia executable.
cd "$(dirname "$0")/.." || exit 1
N=${1:-8}
JULIA="${JULIA:-julia}"
onesweep() {
    local shard=$1
    for i in $(seq 1 60); do
        SWEEP_SHARD=$shard SWEEP_NSHARDS=$N \
            "$JULIA" --project=. --heap-size-hint=6G test/sweep.jl \
            >> "sweep_progress.$shard.log" 2>&1
        rc=$?
        echo "driver[$shard]: pass $i rc=$rc" >> "sweep_progress.$shard.log"
        [ "$rc" -eq 42 ] && continue
        [ "$rc" -eq 0 ] && return 0
        sleep 5
    done
    echo "driver[$shard]: retry budget exhausted" >> "sweep_progress.$shard.log"
    return 1
}
pids=()
for s in $(seq 0 $((N-1))); do
    onesweep "$s" &
    pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
cat sweep_results.*.tsv > sweep_results.tsv 2>/dev/null
echo "merged $(wc -l < sweep_results.tsv) results; fail=$fail"
exit $fail
