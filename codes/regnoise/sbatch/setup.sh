#!/bin/bash
# One-time setup on the cluster: install the environment (Julia 1.10.4) and run the tests.
#   DATA_DIR=/scratch/a.hari/regnoise/inputs sbatch sbatch/setup.sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --partition=ctbp
#SBATCH --job-name=regnoise-setup
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
: "${DATA_DIR:?set DATA_DIR to the folder holding the RACIPE inputs}"
module load julia/1.10.4
CODE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$CODE"
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'
REGNOISE_DATA="$DATA_DIR" julia --project=. test/runtests.jl
