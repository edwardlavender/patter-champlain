#!/usr/bin/env bash
set -euo pipefail

NROW=210
N_CPU_JULIA=20
LOGDIR="./data/output/sim/main/logs"
mkdir -p "$LOGDIR"

export JULIA_NUM_THREADS=1

seq 1 "$NROW" |
  xargs -P "$N_CPU_JULIA" -I {} \
    bash -c '
      LOGFILE="'"$LOGDIR"'/log-$1.log"
      julia "./Julia/001-run-algorithms.jl" "$1" >"$LOGFILE" 2>&1 \
        || echo "Julia failed on row $1" >>"$LOGFILE"
    ' _ {}