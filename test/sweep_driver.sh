#!/bin/bash
# Restart-on-chunk / resume-on-crash driver for test/sweep.jl.
cd "$(dirname "$0")/.." || exit 1
for i in $(seq 1 40); do
    julia --project=. --heap-size-hint=16G test/sweep.jl >> sweep_progress.log 2>&1
    rc=$?
    echo "driver: pass $i exited rc=$rc" >> sweep_progress.log
    [ "$rc" -eq 42 ] && continue          # chunk boundary: restart julia
    [ "$rc" -eq 0 ] && exit 0             # sweep done
    sleep 5                               # crash: resume (culprit recorded)
done
echo "driver: retry budget exhausted" >> sweep_progress.log
exit 1
