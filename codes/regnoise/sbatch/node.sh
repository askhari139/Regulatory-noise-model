#!/bin/bash
# Run one or more (network, noise mode) combinations on a single node. Each combination
# gets NCHUNKS single-threaded Julia processes (default: cores / combinations), each
# handling every NCHUNKS-th selected ParamID. Selection is made first if it is missing.
# run.jl skips finished ParamIDs, so resubmitting the same command resumes.
#
#   NETS=TT,TT4 MODES=Additive,Multiplicative CONFIG=config/pilot.toml \
#     DATA_DIR=... RESULTS_ROOT=... sbatch --cpus-per-task=128 sbatch/node.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=300GB
#SBATCH --time=7-00:00:00
#SBATCH --partition=ctbp
#SBATCH --job-name=regnoise
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${NETS:?set NETS, e.g. TT,TT4}" "${MODES:?set MODES, e.g. Additive}"
CODE="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${CONFIG:-$CODE/config/production.toml}"
EXTRA=()
[[ -n "${DATA_DIR:-}" ]] && EXTRA+=(--data-dir "$DATA_DIR")
[[ -n "${RESULTS_ROOT:-}" ]] && EXTRA+=(--results-root "$RESULTS_ROOT")
IFS=, read -ra NETL <<< "$NETS"; IFS=, read -ra MODEL <<< "$MODES"
NCOMBO=$(( ${#NETL[@]} * ${#MODEL[@]} ))
NCHUNKS="${NCHUNKS:-$(( SLURM_CPUS_PER_TASK / NCOMBO ))}"
(( NCHUNKS >= 1 )) || { echo "fewer cores than combinations"; exit 1; }

module load julia/1.10.4
export JULIA_NUM_THREADS=1
cd "$CODE"
echo "$(date) start on $(hostname): nets=$NETS modes=$MODES chunks/combo=$NCHUNKS commit=$(git rev-parse HEAD)"
julia --project="$CODE" scripts/select.jl --config "$CONFIG" "${EXTRA[@]}" --network "$NETS" --skip-existing

TASKS=$(mktemp)
for net in "${NETL[@]}"; do for mode in "${MODEL[@]}"; do
    mkdir -p "logs/$(basename "$CONFIG" .toml)/$net/$mode"
    for k in $(seq 1 "$NCHUNKS"); do echo "$net $mode $k"; done
done; done > "$TASKS"

run_chunk() {
    local net=$1 mode=$2 k=$3
    local log="logs/$(basename "$CONFIG" .toml)/$net/$mode/chunk_$k.log"
    julia --project="$CODE" --threads=1 scripts/run.jl --config "$CONFIG" "${EXTRA[@]}" \
        --network "$net" --mode "$mode" --chunk "$k/$NCHUNKS" > "$log" 2>&1 \
        || echo "$(date) FAILED $net $mode chunk $k (see $log)"
}
export -f run_chunk; export CODE CONFIG NCHUNKS; export EXTRA_STR="${EXTRA[*]}"
# EXTRA is an array; re-split inside the subshell
xargs -P "$SLURM_CPUS_PER_TASK" -L 1 bash -c 'EXTRA=($EXTRA_STR); run_chunk "$@"' _ < "$TASKS"
rm -f "$TASKS"

echo "$(date) all chunks finished; validating"
status=0
for net in "${NETL[@]}"; do for mode in "${MODEL[@]}"; do
    julia --project="$CODE" scripts/validate.jl --config "$CONFIG" "${EXTRA[@]}" --network "$net" --mode "$mode" || status=1
done; done
exit $status
