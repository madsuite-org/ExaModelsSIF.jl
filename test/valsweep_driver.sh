#!/bin/bash
# Parallel driver for test/valsweep.jl (same scheme as sweep_driver.sh).
cd "$(dirname "$0")/.." || exit 1
N=${1:-12}
JULIA="$HOME/.juliaup/bin/julia"
onesweep() {
    local shard=$1
    for i in $(seq 1 80); do
        SWEEP_SHARD=$shard SWEEP_NSHARDS=$N \
            "$JULIA" --project=. --heap-size-hint=6G test/valsweep.jl \
            >> "val_progress.$shard.log" 2>&1
        rc=$?
        echo "valdriver[$shard]: pass $i rc=$rc" >> "val_progress.$shard.log"
        [ "$rc" -eq 42 ] && continue
        [ "$rc" -eq 0 ] && return 0
        sleep 5
    done
    return 1
}
pids=()
for s in $(seq 0 $((N-1))); do
    onesweep "$s" &
    pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "valsweep merged $(cat val_results.val*.tsv | wc -l) results; fail=$fail"
exit $fail
