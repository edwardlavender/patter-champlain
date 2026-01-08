#!/usr/bin/env bash
set -euo pipefail

SCRIPT="./Julia/002-run-algorithms.jl"
# SCRIPT="./Julia/001-run-filter.jl" 
ANALYSIS="sim"
NROW=4
NCPU=2
LOGDIR="./data/output/$ANALYSIS/main/logs"
mkdir -p "$LOGDIR"

export JULIA_NUM_THREADS=1

seq 1 "$NROW" |
  xargs -P "$NCPU" -I {} \
    bash -c '
      LOGFILE="'"$LOGDIR"'/log-$1.log"
      julia "'"$SCRIPT"'" "$1" >"$LOGFILE" 2>&1 \
        || echo "Julia failed on row $1" >>"$LOGFILE"
    ' _ {}
