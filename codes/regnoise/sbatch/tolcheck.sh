#!/bin/bash
# Tolerance check: rerun the first PER_CLASS selected ParamIDs of each class (per network)
# with config/tolcheck.toml (reltol 1e-6). The pilot selection must already exist; the
# tolcheck config has the same selection settings, so it selects the same ParamIDs.
#
#   NETS=TT,TT4 MODES=Additive,Multiplicative DATA_DIR=... RESULTS_ROOT=... sbatch sbatch/tolcheck.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100GB
#SBATCH --time=1-00:00:00
#SBATCH --partition=ctbp
#SBATCH --job-name=regnoise-tol
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${NETS:?}" "${MODES:?}" "${DATA_DIR:?}" "${RESULTS_ROOT:?}"
PER_CLASS="${PER_CLASS:-4}"
CODE="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="$CODE/config/tolcheck.toml"
module load julia/1.10.4
cd "$CODE"
julia --project="$CODE" scripts/select.jl --config "$CONFIG" --data-dir "$DATA_DIR" --results-root "$RESULTS_ROOT" --network "$NETS" --skip-existing
RUN_ID=$(grep '^run_id' "$CONFIG" | cut -d'"' -f2)
mkdir -p logs/tolcheck
TASKS=$(mktemp)
IFS=, read -ra NETL <<< "$NETS"; IFS=, read -ra MODEL <<< "$MODES"
for net in "${NETL[@]}"; do
    ids=$(awk -F, -v n="$PER_CLASS" 'NR>1 && c[$2]++ < n {print $1}' "$RESULTS_ROOT/$RUN_ID/selection/$net.csv")
    for mode in "${MODEL[@]}"; do for id in $ids; do echo "$net $mode $id"; done; done
done > "$TASKS"
xargs -P "$SLURM_CPUS_PER_TASK" -L 1 bash -c \
  'julia --project="'"$CODE"'" --threads=1 scripts/run.jl --config "'"$CONFIG"'" --data-dir "'"$DATA_DIR"'" --results-root "'"$RESULTS_ROOT"'" --network $0 --mode $1 --params $2 > logs/tolcheck/$0_$1_$2.log 2>&1 || echo "FAILED $0 $1 $2"' < "$TASKS"
rm -f "$TASKS"
echo "$(date) tolerance check finished"
