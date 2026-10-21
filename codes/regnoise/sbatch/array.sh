#!/bin/bash
# Job-array alternative to node.sh: one single-CPU task per chunk.
#
#   NET=TS MODE=Additive sbatch --array=1-200 sbatch/array.sh
#
# NCHUNKS must equal the array size. Validate afterwards with scripts/validate.jl.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=2-00:00:00
#SBATCH --partition=ctbp
#SBATCH --job-name=regnoise
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
: "${NET:?set NET}" "${MODE:?set MODE}"
CODE="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${CONFIG:-$CODE/config/production.toml}"
NCHUNKS="${SLURM_ARRAY_TASK_COUNT:?run as an array job}"
EXTRA=()
[[ -n "${DATA_DIR:-}" ]] && EXTRA+=(--data-dir "$DATA_DIR")
[[ -n "${RESULTS_ROOT:-}" ]] && EXTRA+=(--results-root "$RESULTS_ROOT")
module load julia/1.10.4
cd "$CODE"
julia --project="$CODE" --threads=1 scripts/run.jl --config "$CONFIG" "${EXTRA[@]}" \
      --network "$NET" --mode "$MODE" --chunk "$SLURM_ARRAY_TASK_ID/$NCHUNKS"
