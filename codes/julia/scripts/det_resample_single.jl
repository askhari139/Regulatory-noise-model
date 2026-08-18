"""
Worker script: called by each SLURM array task.
Runs the resampling analysis for a single parameter ID and saves a CSV.

Usage (direct):
    julia scripts/det_resample_single.jl \\
        network=TS \\
        param-id=42 \\
        sigma-values=0.0,0.01,0.05,0.1 \\
        work-dir=det_resample/TS \\
        num-sims=100 \\
        tspan=0.0,500.0
"""

using DataFrames
using CSV

sourceDir = get(ENV, "RACIPE_SOURCE", pwd()) * "/"
include(sourceDir * "scripts/stochastic_setup.jl")
include(sourceDir * "scripts/parse_args.jl")
include(sourceDir * "src/simulate_racipe_resampled.jl")

# ── arguments ─────────────────────────────────────────────────────────────────

args = parse_named_args(defaults=Dict(
    "network"      => "TS",
    "param-id"     => "0",
    "sigma-values" => "0.0,0.01,0.05,0.1",
    "num-sims"     => "100",
    "tspan"        => "0.0,500.0",
    "work-dir"     => "det_resample"
))

NETWORK      = get_arg(args, "network",      "TS")
PARAM_ID     = get_arg(args, "param-id",     0,   type=Int)
SIGMA_VALUES = get_arg(args, "sigma-values", [0.0, 0.01, 0.05, 0.1],
                       type=Vector{Float64})
NUM_SIMS     = get_arg(args, "num-sims",     100, type=Int)
tspan_vec    = get_arg(args, "tspan",        [0.0, 500.0], type=Vector{Float64})
TSPAN        = (tspan_vec[1], tspan_vec[2])
WORK_DIR     = get_arg(args, "work-dir",     joinpath("det_resample", NETWORK))

PARAM_ID == 0 && error("Must specify param-id=XXX")
RESAMPLE = get_arg(args, "resample", true, type=Bool)

println("=" ^ 60)
println("WORKER: network=$NETWORK  param-id=$PARAM_ID")
println("=" ^ 60)

# ── load data ─────────────────────────────────────────────────────────────────

cd("data")
PRS_FILE        = NETWORK * ".prs"
PARAMETERS_FILE = NETWORK * "_parameters.dat"
TOPO_FILE       = NETWORK * ".topo"
ODE_FILE        = sourceDir * "src/" * NETWORK * "_ode.jl"

prs         = read_prs(PRS_FILE)
params_data = read_parameters(PARAMETERS_FILE, prs)

if !isfile(ODE_FILE)
    generate_ode_function(TOPO_FILE, PRS_FILE, ODE_FILE)
end
include(ODE_FILE)
cd("..")

# ── run ───────────────────────────────────────────────────────────────────────

t_start = time()

df, pDf = analyze_single_param(
    ode_system!,
    PARAM_ID,
    params_data,
    prs,
    NODE_NAMES;
    sigma_values = SIGMA_VALUES,
    num_sims     = NUM_SIMS,
    tspan        = TSPAN,
    # seed         = 42,
    resample      = RESAMPLE
)

elapsed = round(time() - t_start, digits=1)
println("Done in $(elapsed)s  ($(nrow(df)) rows)")

# ── save ──────────────────────────────────────────────────────────────────────

results_dir = joinpath(WORK_DIR, "results")
mkpath(results_dir)
out_file = joinpath(results_dir, "sol_$(PARAM_ID).csv")
CSV.write(out_file, df)
out_file = joinpath(results_dir, "param_$(PARAM_ID).csv")
CSV.write(out_file, pDf)
println("Saved: $out_file")