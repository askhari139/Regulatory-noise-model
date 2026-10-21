# regnoise

Stochastic regulatory-strength (λ) simulations on RACIPE networks for the CSBJ revision.
**SPEC.md is the specification**; `config/production.toml` holds every number.

## Layout

| Path | Purpose |
|---|---|
| `config/production.toml` | the only config used for results (hash stamped into every output) |
| `config/test.toml` | tiny sizes for the end-to-end test; never for results |
| `src/` | the `RegNoise` module |
| `scripts/select.jl` | choose parameter sets per network (once per run) |
| `scripts/run.jl` | all conditions for a chunk of ParamIDs (resumable) |
| `scripts/validate.jl` | check outputs against the config; exit 1 on any problem |
| `scripts/collect.jl` | validate, then merge into `collected/*.csv.gz` |
| `sbatch/` | cluster setup and job scripts |
| `test/` | unit + end-to-end tests (`legacy_ode/` holds reference ODE files) |

## Cluster workflow (Explorer)

```bash
# 0. code: a clean, committed checkout (run.jl refuses a dirty tree)
git clone https://github.com/askhari139/Regulatory-noise-model.git /scratch/a.hari/regnoise_code
cd /scratch/a.hari/regnoise_code/codes/regnoise

# 1. inputs: the 7 networks' RACIPE files (checked against config [inputs.sha256] at every start)
mkdir -p /scratch/a.hari/regnoise/inputs   # copy TS*, TSSA*, TT*, TTSA*, TS4*, TT4*, TT4SA* .prs/.topo/_parameters.dat/_solution.dat

# 2. environment + tests
DATA_DIR=/scratch/a.hari/regnoise/inputs bash sbatch/setup.sh

# 3. select parameter sets (writes <results_root>/<run_id>/selection/)
export DATA_DIR=/scratch/a.hari/regnoise/inputs RESULTS_ROOT=/scratch/a.hari/regnoise/results
module load julia/1.10.4
julia --project=. scripts/select.jl --data-dir $DATA_DIR --results-root $RESULTS_ROOT

# 4. run, one node per (network, mode); resubmit the same line to resume
NET=TS MODE=Additive sbatch --cpus-per-task=128 sbatch/node.sh

# 5. collect (validates first)
julia --project=. scripts/collect.jl --data-dir $DATA_DIR --results-root $RESULTS_ROOT --network TS --mode Additive
```

## Local tests

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
REGNOISE_DATA=/path/to/RACIPEdata/data julia --project=. test/runtests.jl
```
