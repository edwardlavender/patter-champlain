#!/usr/bin/env bash
set -euo pipefail

ANALYSIS="sim"
NROW=4
LOGDIR="./data/output/$ANALYSIS/main/logs"
NCPU=2
mkdir -p "$LOGDIR"

export JULIA_NUM_THREADS=1

seq 1 "$NROW" |
  xargs -P "$NCPU" -I {} \
    bash -c '
      LOGFILE="'"$LOGDIR"'/log-$1.log"
      julia "./Julia/001-run-algorithms.jl" "$1" >"$LOGFILE" 2>&1 \
        || echo "Julia failed on row $1" >>"$LOGFILE"
    ' _ {}