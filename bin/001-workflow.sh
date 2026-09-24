#!/usr/bin/env bash
set -euo pipefail

# SCRIPT="./Julia/002-run-filter.jl" 
SCRIPT="./Julia/003-run-algorithms.jl"
ANALYSIS="validation"
NROW=665
NCPU=40
LOGDIR="./data/output/$ANALYSIS/main/logs"
mkdir -p "$LOGDIR"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OPENBLAS_DEFAULT_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

seq 1 "$NROW" |
  xargs -P "$NCPU" -I {} \
    bash -c '
      LOGFILE="'"$LOGDIR"'/log-$1.log"
      julia "'"$SCRIPT"'" "$1" >"$LOGFILE" 2>&1 \
        || echo "Julia failed on row $1" >>"$LOGFILE"
    ' _ {}
