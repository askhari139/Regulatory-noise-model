"""
5_plot_results.jl

For each selected parameter set, run short simulations at the lowest and
highest non-zero sigma levels and plot expression trajectories.

Optionally saves the raw Float32 trajectory matrices to JLD2 (binary) so
they can be loaded in R (via RJulia / JLD2.jl) or analysed further in Julia.

Usage:
    julia --project=. scripts/5_plot_results.jl \\
        work-dir=Additive/TS

    # Select specific parameter IDs instead of auto-selecting top/bottom 10
    julia --project=. scripts/5_plot_results.jl \\
        work-dir=Additive/TS \\
        param-ids=12,45,78,200

    # Change number of params shown, ICs per plot, and traj length
    julia --project=. scripts/5_plot_results.jl \\
        work-dir=Additive/TS \\
        top-n=5 \\
        num-ics=8 \\
        tmax=500.0

    # Also save raw Float32 trajectories to JLD2
    julia --project=. scripts/5_plot_results.jl \\
        work-dir=Additive/TS \\
        save-trajectories=true

    # Save trajectories only (skip plot generation)
    julia --project=. scripts/5_plot_results.jl \\
        work-dir=Additive/TS \\
        save-trajectories=true \\
        plots-only=false
"""

using CSV
using DataFrames
using Statistics
using JLD2
using Plots
using Random

function has_slurm()
    try
        run(pipeline(`which sbatch`, stdout=devnull, stderr=devnull))
        return true
    catch
        return false
    end
end
# sourceDir = get(ENV, "RACIPE_SOURCE", pwd())
if has_slurm()
    sourceDir = "/scratch/a.hari/RACIPEdata"
else
    sourceDir = dirname(@__DIR__)  # portable: resolves to this codes/julia checkout
end
scriptsDir = joinpath(sourceDir, "scripts")
include(joipath(scriptsDir, "stochastic_setup.jl"))
include(joipath(scriptsDir, "parse_args.jl"))
include(joipath(scriptsDir, "config_utils.jl"))
dataFolder = joinpath(sourceDir, "data")
dataUniform = joinpath(sourceDir, "dataUniform")
ENV["GKSwstype"] = "nul"
gr(show=false)

# ============================================================
# ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "work-dir"           => "work",
    "param-ids"          => "",        # comma-separated; empty = auto-select
    "top-n"              => "10",      # params to show at each sigma extreme
    "num-ics"            => "6",       # trajectories per param per sigma
    "tmax"               => "500.0",
    "saveat"             => "1.0",
    "seed"               => "42",
    "save-trajectories"  => "false",   # write Float32 JLD2 files
    "plots-only"         => "true",    # set false to skip PNG generation
    "sigma-low"          => "",        # override lowest non-zero sigma
    "sigma-high"         => "",        # override highest sigma
    "use-uniform"        => "false",
))

WORK_DIR           = get_arg(args, "work-dir",          "work")
PARAM_IDS_ARG      = get_arg(args, "param-ids",         "")
TOP_N              = get_arg(args, "top-n",             10,     type=Int)
NUM_ICS            = get_arg(args, "num-ics",           6,      type=Int)
TMAX               = get_arg(args, "tmax",              500.0,  type=Float64)
SAVEAT             = get_arg(args, "saveat",            1.0,    type=Float64)
SEED               = get_arg(args, "seed",              42,     type=Int)
SAVE_TRAJECTORIES  = get_arg(args, "save-trajectories", false,  type=Bool)
PLOTS_ONLY         = get_arg(args, "plots-only",        true,   type=Bool)   # when false, skip PNGs
SIGMA_LOW_ARG      = get_arg(args, "sigma-low",         "")
SIGMA_HIGH_ARG     = get_arg(args, "sigma-high",        "")
USE_UNIFORM        = get_arg(args, "use-uniform", false, type=Bool)
if USE_UNIFORM
    dataFolder = dataUniform
end
Random.seed!(SEED)

println("=" ^ 70)
println("TRAJECTORY PLOTS")
println("=" ^ 70)
println("Work dir          : $WORK_DIR")
println("Top-N per extreme : $TOP_N")
println("ICs per param     : $NUM_ICS")
println("tmax              : $TMAX")
println("Save trajectories : $SAVE_TRAJECTORIES")
println("Generate plots    : $PLOTS_ONLY")
println("=" ^ 70)

# ============================================================
# LOAD CONFIGURATION
# ============================================================

param_file = joinpath(WORK_DIR, "selected_parameters.jld2")
job_file   = joinpath(WORK_DIR, "job_config.jld2")

(isfile(param_file) && isfile(job_file)) ||
    error("Config files not found in $WORK_DIR — run steps 1 and 2 first")

param_data = load_parameter_config(param_file)
job_config = load_job_config(job_file)

NETWORK           = param_data["network"]
NODE_NAMES        = param_data["node_names"]
lambda_indices    = param_data["lambda_indices"]
max_lambdas       = param_data["max_lambdas"]
min_lambdas       = param_data["min_lambdas"]
racipe_thresholds = param_data["racipe_thresholds"]
NOISE_MODE        = job_config["noise_mode"]
NOISE_LEVELS      = job_config["noise_levels"]
PARAM_TYPE        = Symbol(job_config["param_type"])

n_nodes = length(NODE_NAMES)

# ============================================================
# LOAD RACIPE DATA + ODE
# ============================================================
original_dir = pwd()
cd(dataFolder)
PRS_FILE        = NETWORK * ".prs"
PARAMETERS_FILE = NETWORK * "_parameters.dat"
SOLUTIONS_FILE  = NETWORK * "_solution.dat"
ODE_FILE        = sourceDir * "src/" * NETWORK * "_ode.jl"

prs          = read_prs(PRS_FILE)
params_data  = read_parameters(PARAMETERS_FILE, prs)
sol_data     = read_solutions(SOLUTIONS_FILE, prs)
include(ODE_FILE)
cd(original_dir)

# ============================================================
# DETERMINE SIGMA PAIR
# ============================================================

nonzero_sigmas = sort(filter(x -> x > 0.0, NOISE_LEVELS))
isempty(nonzero_sigmas) && error("No non-zero sigma values in job config")

sigma_low  = isempty(SIGMA_LOW_ARG)  ? nonzero_sigmas[1]   : parse(Float64, SIGMA_LOW_ARG)
sigma_high = isempty(SIGMA_HIGH_ARG) ? nonzero_sigmas[end] : parse(Float64, SIGMA_HIGH_ARG)

println("Sigma low  : $sigma_low")
println("Sigma high : $sigma_high")

# ============================================================
# SELECT PARAMETER IDs
# ============================================================

if !isempty(PARAM_IDS_ARG)
    # Explicit list from CLI
    target_ids = parse.(Int, split(PARAM_IDS_ARG, ","))
    println("Using supplied param-ids: $target_ids")
else
    # Auto-select: load all_parameters_results.csv, rank by mean MRT spread
    results_file = joinpath(WORK_DIR, "results", "all_parameters_results.csv")
    if !isfile(results_file)
        @warn "all_parameters_results.csv not found — using all param_ids from config"
        target_ids = param_data["param_ids"]
    else
        df = CSV.read(results_file, DataFrame)

        # Score each param: max per-state MRT spread across noise levels
        spread_df = combine(
            groupby(df, [:ParamID, :State]),
            :MRT => (x -> maximum(x) - minimum(x)) => :MRT_spread
        )
        param_scores = combine(
            groupby(spread_df, :ParamID),
            :MRT_spread => maximum => :max_spread
        )
        sort!(param_scores, :max_spread, rev=true)

        n_avail    = nrow(param_scores)
        n_each     = min(TOP_N, n_avail)
        top_ids    = param_scores.ParamID[1:n_each]
        bottom_ids = param_scores.ParamID[max(1, n_avail - n_each + 1):n_avail]
        # Remove overlap
        bottom_ids = setdiff(bottom_ids, top_ids)
        target_ids = unique(vcat(top_ids, bottom_ids))

        println("Auto-selected $(length(top_ids)) high-spread + $(length(bottom_ids)) low-spread params")
    end
end

# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

fig_dir   = joinpath(WORK_DIR, "figures", "trajectories")
traj_dir  = joinpath(WORK_DIR, "trajectories")
PLOTS_ONLY && mkpath(fig_dir)
SAVE_TRAJECTORIES && mkpath(traj_dir)

# ============================================================
# NODE COLOURS (cycle for up to 8 nodes)
# ============================================================

NODE_COLORS = [:steelblue, :firebrick, :seagreen, :darkorange,
               :mediumpurple, :saddlebrown, :teal, :hotpink]

# ============================================================
# DISCRETE-STATE HELPERS
# ============================================================

all_states_enum = [Tuple(reverse(digits(i, base=2, pad=n_nodes)))
                   for i in 0:(2^n_nodes - 1)]
state_labels    = ["(" * join(s, ",") * ")" for s in all_states_enum]
state_index_map = Dict(all_states_enum[i] => i - 1 for i in eachindex(all_states_enum))

@inline function disc_idx(u, thr)
    get(state_index_map, Tuple(u[j] > thr[j] ? 1 : 0 for j in eachindex(thr)), 0)
end

# ============================================================
# MAIN LOOP
# ============================================================

println("\nProcessing $(length(target_ids)) parameter sets …")

for param_id in target_ids
    param_row = params_data.data[params_data.data.ParamID .== param_id, :]
    nrow(param_row) > 0 || (println("  ⚠ ParamID $param_id not found — skip"); continue)

    param_vector = [param_row[1, Symbol(nm)] for nm in prs.names]

    # Build u0List from RACIPE steady states; fall back to random
    u0List = Vector{Vector{Float64}}()
    for row in eachrow(sol_data.data[sol_data.data.ParamID .== param_id, :])
        push!(u0List, [2.0^row[Symbol(nd)] for nd in NODE_NAMES])
    end
    if isempty(u0List)
        upper = [param_row[1, Symbol("Prod_of_$nd")] / param_row[1, Symbol("Deg_of_$nd")]
                 for nd in NODE_NAMES]
        u0List = [rand(n_nodes) .* upper for _ in 1:NUM_ICS]
    end
    # Cycle to NUM_ICS
    u0List = u0List[mod1.(1:NUM_ICS, length(u0List))]

    tspan = (0.0, TMAX)

    # ── run at both sigmas ────────────────────────────────────────────────────
    traj_store = Dict{Float64, NamedTuple}()   # sigma => (times, expr, disc)

    for sigma in (sigma_low, sigma_high)
        all_times   = Float64[]
        all_expr    = [Float64[] for _ in NODE_NAMES]
        all_disc    = Int[]
        all_sim_ids = Int[]

        for (i, u0) in enumerate(u0List)
            sol = simulate_with_noise(
                ode_system!, u0, tspan, copy(param_vector), lambda_indices;
                sigma      = sigma,
                saveat     = SAVEAT,
                max_lambdas = max_lambdas,
                min_lambdas = min_lambdas,
                noise_mode  = NOISE_MODE
            )

            t_vec   = collect(sol.t)
            ex_mat  = Array(sol)   # nodes × timepoints
            d_idx   = [disc_idx(ex_mat[:, k], racipe_thresholds) for k in axes(ex_mat, 2)]

            append!(all_times,   t_vec)
            append!(all_disc,    d_idx)
            append!(all_sim_ids, fill(i, length(t_vec)))
            for j in 1:n_nodes
                append!(all_expr[j], ex_mat[j, :])
            end
        end

        traj_store[sigma] = (
            times   = all_times,
            expr    = all_expr,       # Vector{Vector{Float64}}, one per node
            disc    = all_disc,
            sim_ids = all_sim_ids,
        )
    end

    # ── save trajectories to JLD2 (Float32 for space) ───────────────────────
    if SAVE_TRAJECTORIES
        out_jld = joinpath(traj_dir, "param_$(param_id).jld2")
        jldsave(out_jld;
            param_id     = param_id,
            node_names   = NODE_NAMES,
            sigma_low    = sigma_low,
            sigma_high   = sigma_high,
            thresholds   = Float32.(racipe_thresholds),
            state_labels = state_labels,
            # low sigma
            low_times    = Float32.(traj_store[sigma_low].times),
            low_sim_ids  = traj_store[sigma_low].sim_ids,
            low_disc     = traj_store[sigma_low].disc,
            low_expr     = [Float32.(v) for v in traj_store[sigma_low].expr],
            # high sigma
            high_times   = Float32.(traj_store[sigma_high].times),
            high_sim_ids = traj_store[sigma_high].sim_ids,
            high_disc    = traj_store[sigma_high].disc,
            high_expr    = [Float32.(v) for v in traj_store[sigma_high].expr],
        )
        println("  ✓ Trajectories saved → $out_jld")
    end

    # ── generate plots ────────────────────────────────────────────────────────
    if !PLOTS_ONLY
        println("  (plots skipped for param $param_id)")
        continue
    end

    # One figure per sigma level
    for sigma in (sigma_low, sigma_high)
        td = traj_store[sigma]
        sim_ids_unique = unique(td.sim_ids)

        # Layout: n_nodes expression rows + 1 discrete-state row
        n_rows = n_nodes + 1
        fig = plot(layout=(n_rows, 1),
                   size=(900, 220 * n_rows),
                   left_margin=8Plots.mm,
                   bottom_margin=4Plots.mm)

        # Expression panels
        for nd in 1:n_nodes
            for sim_id in sim_ids_unique
                mask = td.sim_ids .== sim_id
                plot!(fig[nd],
                      td.times[mask], td.expr[nd][mask];
                      color     = NODE_COLORS[mod1(nd, length(NODE_COLORS))],
                      alpha     = 0.45,
                      linewidth = 0.9,
                      label     = "")
            end
            hline!(fig[nd], [racipe_thresholds[nd]];
                   color=NODE_COLORS[mod1(nd, length(NODE_COLORS))],
                   linestyle=:dash, linewidth=1.2, alpha=0.8, label="")
            ylabel!(fig[nd], NODE_NAMES[nd])
        end

        # Discrete-state panel
        n_states = length(all_states_enum)
        for sim_id in sim_ids_unique
            mask = td.sim_ids .== sim_id
            plot!(fig[n_rows],
                  td.times[mask], td.disc[mask];
                  color=:steelblue, alpha=0.4, linewidth=0.9, label="")
        end
        plot!(fig[n_rows];
              yticks  = (0:(n_states - 1), state_labels),
              ylims   = (-0.5, n_states - 0.5),
              xlabel  = "Time",
              ylabel  = "State")

        sigma_str = replace(string(sigma), "." => "p")
        plot!(fig;
              plot_title = "ParamID $param_id | $(NOISE_MODE) | σ=$sigma | $(length(sim_ids_unique)) ICs",
              plot_titlefontsize = 10)

        fname = joinpath(fig_dir, "param$(param_id)_sigma$(sigma_str).png")
        savefig(fig, fname)
        println("  ✓ $(basename(fname))")
    end
end

println("\n" * "=" ^ 70)
println("✅ DONE")
println("=" ^ 70)
PLOTS_ONLY       && println("Figures  → $(joinpath(WORK_DIR, "figures", "trajectories"))/")
SAVE_TRAJECTORIES && println("JLD2     → $(joinpath(WORK_DIR, "trajectories"))/")
println()
println("To load a trajectory file in Julia:")
println("  using JLD2")
println("  d = load(\"$(joinpath(WORK_DIR, "trajectories"))/param_<ID>.jld2\")")
println("  # keys: low_expr, low_times, low_sim_ids, low_disc,")
println("  #       high_expr, high_times, high_sim_ids, high_disc,")
println("  #       node_names, thresholds, state_labels")