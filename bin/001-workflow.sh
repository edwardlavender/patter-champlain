#!/usr/bin/env bash
set -euo pipefail

# Define global settings 
NROW=210
NROW_PER_BATCH=100
N_CPU_JULIA=20
N_CPU_R=5

export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Loop over batches: 1:100, 101:200, etc.
for (( BATCH_START = 1; BATCH_START <= NROW; BATCH_START += NROW_PER_BATCH )); do
  BATCH_END=$(( BATCH_START + NROW_PER_BATCH - 1 ))
  (( BATCH_END > NROW )) && BATCH_END=$NROW

  echo "Rows ${BATCH_START}-${BATCH_END}"

  # --------------------------------
  # A) Run Julia workflow I (parallel)
  # --------------------------------
  seq "$BATCH_START" "$BATCH_END" |
    xargs -n 1 -P "$N_CPU_JULIA" -I {} \
      bash -c 'julia ./Julia/001-run-algorithms.jl "$1" || echo "Julia I error (row $1)" >&2' _ {}

  # --------------------------------
  # B) Run Julia workflow II (per row, serial)
  # --------------------------------
  # Use 1 CPU as this script is more memory intensive
  for (( i = BATCH_START; i <= BATCH_END; i++ )); do
    julia ./Julia/002-collate-states.jl "$i" || \
      echo "Julia II error (row $i)" >&2
  done

  # --------------------------------
  # C) Run R workflow (per row, limited parallelism)
  # --------------------------------
  # Each iteration requires > 1.80 GB of memory
  seq "$BATCH_START" "$BATCH_END" |
    xargs -n 1 -P "$N_CPU_R" -I {} \
      bash -c 'Rscript ./R/018-run-patter.R "$1" || echo "R error (row $1)" >&2' _ {}

done
