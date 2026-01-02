#!/usr/bin/env bash
set -euo pipefail
trap 'kill -- -$$' EXIT

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
  for (( i = BATCH_START; i <= BATCH_END; i++ )); do
    (
      julia ./Julia/001-run-algorithms.jl "$i"
    ) || echo "Julia I error (row $i)" >&2 &

    (( $(jobs -r | wc -l) >= N_CPU_JULIA )) && wait -n
  done
  wait

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
  for (( i = BATCH_START; i <= BATCH_END; i++ )); do
    (
      Rscript ./R/018-run-patter.R "$i"
    ) || echo "R error (row $i)" >&2 &

    (( $(jobs -r | wc -l) >= N_CPU_R )) && wait -n
  done
  wait
done
