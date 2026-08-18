println("Running the script now.")
using JSON
using DataFrames
using CSV
using Statistics
using JLD2
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
# dataFolder = joinpath(sourceDir, "data")
# dataUniform = joinpath(sourceDir, "dataUniform")
scriptsDir = joinpath(sourceDir, "scripts")

y1 = @elapsed include(joinpath(scriptsDir, "stochastic_setup.jl"))
println("Loading libraries took $(round(y1, digits=2)) seconds.")
include(joinpath(scriptsDir, "parse_args.jl"))
include(joinpath(scriptsDir, "config_utils.jl"))
include(joinpath(scriptsDir, "lambda_sampler.jl"))

# ============================================================
# PARSE ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "param-id"   => "0",
    "work-dir"   => ".",
    "network"    => "TS",
    "noise-mode" => "Additive",
    "det-iters"  => "10"
    # "use-uniform" => "false",
))

PARAM_ID   = get_arg(args, "param-id",   0,         type=Int)
NOISE_MODE = get_arg(args, "noise-mode", "Additive")
NETWORK    = get_arg(args, "network",    "TS")
USE_UNIFORM= get_arg(args, "use-uniform", false, type=Bool)
# FIX #11: actually use the work-dir argument instead of hard-coding it
WORK_DIR = get_arg(args, "work-dir", joinpath(NOISE_MODE, NETWORK))
DET_ITERS = get_arg(args, "det-iters", 10, type = Int)
if WORK_DIR == "."
    WORK_DIR = joinpath(NOISE_MODE, NETWORK)
end
println(WORK_DIR)
# if USE_UNIFORM
#     dataFolder = dataUniform
# end
PARAM_ID == 0 && error("Must specify param-id=XXX")

println("=" ^ 70)
println("ANALYZING PARAMETER $PARAM_ID")
println("=" ^ 70)

# ============================================================
# LOAD CONFIGURATION
# ============================================================

param_file = joinpath(WORK_DIR, "selected_parameters.jld2")
job_file   = joinpath(WORK_DIR, "job_config.jld2")

(isfile(param_file) && isfile(job_file)) ||
    error("Configuration files not found in $WORK_DIR\nRun steps 1 and 2 first")

param_data = load_parameter_config(param_file)
job_config = load_job_config(job_file)

NETWORK           = param_data["network"]
NODE_NAMES        = param_data["node_names"]
lambda_indices         = param_data["lambda_indices"]
max_lambdas            = param_data["max_lambdas"]
min_lambdas            = param_data["min_lambdas"]
racipe_thresholds = param_data["racipe_thresholds"]
dataFolder = param_data["data_folder"]
NOISE_MODE   = job_config["noise_mode"]
NUM_SIMS     = job_config["num_sims"]
NOISE_LEVELS = job_config["noise_levels"]
PARAM_TYPE   = Symbol(job_config["param_type"])
DT_vals           = job_config["dt"]
println(DT_vals)

param_types = param_data["param_types"]

# FIX #10: JSON deserialises integer keys as strings — always look up as string
# param_type_key = string(PARAM_ID)
param_type_key = haskey(param_types, string(PARAM_ID)) ? string(PARAM_ID) : PARAM_ID
# param_type = param_types[param_type_key]
haskey(param_types, param_type_key) ||
    error("Parameter $PARAM_ID not found in param_types. Available: $(keys(param_types))")
param_type = param_types[param_type_key]

println("  Network       : $NETWORK")
println("  Parameter type: $param_type")
println("  Noise mode    : $NOISE_MODE")
println("=" ^ 70)

# ============================================================
# LOAD DATA
# ============================================================
original_dir = pwd()
cd(dataFolder)

PRS_FILE        = NETWORK * ".prs"
PARAMETERS_FILE = NETWORK * "_parameters.dat"
SOLUTIONS_FILE  = NETWORK * "_solution.dat"
ODE_FILE        = joinpath(sourceDir , "src" , NETWORK * "_ode.jl")

prs             = read_prs(PRS_FILE)
params_data     = read_parameters(PARAMETERS_FILE, prs)
solutions_data  = read_solutions(SOLUTIONS_FILE, prs)

include(ODE_FILE)

cd(original_dir)
println("Thresholds used to discretize : "* string(racipe_thresholds))
println("RACIPE data sourced from "* dataFolder)
# ============================================================
# GET PARAMETER DATA
# ============================================================

param_row = params_data.data[params_data.data.ParamID .== PARAM_ID, :]
nrow(param_row) > 0 || error("Parameter $PARAM_ID not found in data")

param_vector = [param_row[1, Symbol(name)] for name in prs.names]
upper_limits = [param_row[1, Symbol("Prod_of_"*node)]/param_row[1, Symbol("Deg_of_"*node)] for node in NODE_NAMES]
A = rand(100, length(NODE_NAMES))
u0List = [A[i, :].*upper_limits for i in 1:size(A, 1)]

function replace_lambdas(param_vector, lambda_indices, max_lambdas,
                         noise_level, lambda_inh, lambda_act, noise_mode)
    p = copy(param_vector)
    for l in eachindex(lambda_indices)
        idx = lambda_indices[l]
        cache = max_lambdas[l] <= 1.0 ? lambda_inh : lambda_act
        p[idx] = sample_lambda(cache, noise_mode, noise_level, p[idx])
    end
    return p
end

function replace_lambdas_uniform(param_vector, lambda_indices, max_lambdas)
    p = copy(param_vector)
    r = rand(length(lambda_indices))        # never mutate the base vector
    for l in eachindex(lambda_indices)
        idx = lambda_indices[l]
        lambda = p[idx]
        rnd = r[l]
        if max_lambdas[l] <= 1.0 
            lambda = rnd * (0.99) + 0.01
        else
            lambda = rnd * 99.1 + 1.0
        end
        p[idx] = lambda
    end
    return p
end

# ============================================================
# RUN ANALYSIS
# ============================================================

println("\nRunning deterministic simulations...")
start_time = time()

TSPAN = (0.0, 1000.0)
results_all = DataFrame[]
transitions_all = DataFrame[]
stats_all = DataFrame[]
lambdas_all = DataFrame[]
try
    for iter in 1:DET_ITERS
        Random.seed!(rand(1:10000))  # different seed per iteration for lambda sampling
        all_results = Dict{Tuple{Float64, Float64}, Vector{StochasticResult}}()
        replaced_lambdas = Dict{Float64, Vector{Float64}}()
        y1 = 0.0
        for dt in DT_vals
            jld = joinpath(dataFolder, "lambda_hist_"*NOISE_MODE*"_dt$(dt).jld2")
            println(jld)
            if (isfile(jld))
                @load jld lambda_inh lambda_act
            else
                error("Lambda distribution not found.")
            end
        
            for noise in NOISE_LEVELS
                # println("Running $(PARAM_ID) DT $(dt) at Sigma $(noise) for the $(iter) time")
                # global y1
                p_noise = replace_lambdas(param_vector, lambda_indices, max_lambdas, noise, lambda_inh, lambda_act, NOISE_MODE)
                if noise == 0.0
                    p_noise = copy(param_vector)
                end
                if noise == maximum(NOISE_LEVELS)
                    println(p_noise[lambda_indices])
                end
                replaced_lambdas[noise] = p_noise[lambda_indices]
                y = 0.0
                results = run_multiple_stochastic_simulations(
                    ode_system!, p_noise, prs.names, lambda_indices, NODE_NAMES,
                    num_sims=100, tspan=TSPAN, sigma=0.0,
                    u0List = u0List, racipe_thresholds = racipe_thresholds,
                    max_lambdas=max_lambdas, min_lambdas=min_lambdas, noise_mode = NOISE_MODE,
                    saveat=1.0, dt = dt
                )
                all_results[(noise, dt)] = results
                y1 = y1 +y
            end
        end
        summary = create_summary_dataframe(all_results, NODE_NAMES)

        println("Simulations took $(round(y1, digits=2)) seconds")
        
        println("aggregate_mrt took $(round(y1, digits=2)) seconds")

        # ---- Results DataFrame ----
        # WITH THIS
        results_df = DataFrame(
            ParamID      = Int[],
            ParamType    = String[],
            NoiseLevel   = Float64[],
            DT           = Float64[],
            State        = String[],
            MRT          = Float64[],
            RelativeMRT  = Float64[],
            MeanSwitches = Float64[],
            StdSwitches  = Float64[]
        )
        lambdas_df = DataFrame(
            ParamID    = Int[],
            ParamType  = String[],
            NoiseLevel = Float64[],
            DT         = Float64[]
        )
        for (l, idx) in enumerate(lambda_indices)
            lambdas_df[!, Symbol("lambda_$(prs.names[idx])")] .= Float64[]
        end

        transitions_df = DataFrame(
            ParamID    = Int[],
            ParamType  = String[],
            NoiseLevel = Float64[],
            DT         = Float64[],
            Transition = String[],
            Count      = Float64[]
        )

        # WITH THIS
        stats_df = DataFrame(
            ParamID    = Int[],
            ParamType  = String[],
            NoiseLevel = Float64[],
            DT         = Float64[],
            Node       = String[],
            Mean       = Float64[],
            Count      = Int[]
        )
        for dt in DT_vals
            deterministic_mrt = aggregate_mrt(all_results[(0.0, dt)])
            for sigma in NOISE_LEVELS
                results   = all_results[(sigma,dt)]
                mean_mrt  = aggregate_mrt(results)
                mean_counts = aggregate_transitions(results)
                switches  = [r.switches for r in results]

                for (state, mrt) in mean_mrt
                    det_mrt = get(deterministic_mrt, state, 0.0)
                    rel_mrt = 2.0^(mrt - det_mrt)
                    
                    # WITH THIS
                    push!(results_df, (
                        ParamID      = PARAM_ID,
                        ParamType    = param_type,
                        NoiseLevel   = sigma,
                        DT           = dt,
                        State        = string(state),
                        MRT          = mrt,
                        RelativeMRT  = rel_mrt,
                        MeanSwitches = mean(switches),
                        StdSwitches  = std(switches)
                    ))
                end

                # ADD THIS after the for (state, mrt) loop
                lambda_row = Dict{Symbol, Any}(
                    :ParamID    => PARAM_ID,
                    :ParamType  => param_type,
                    :NoiseLevel => sigma,
                    :DT         => dt
                )
                for (l, idx) in enumerate(lambda_indices)
                    lambda_row[Symbol("lambda_$(prs.names[idx])")] = replaced_lambdas[sigma][l]
                end
                push!(lambdas_df, lambda_row)

                for (trans, count) in mean_counts
                    push!(transitions_df, (
                        ParamID    = PARAM_ID,
                        ParamType  = param_type,
                        NoiseLevel = sigma,
                        DT         = dt,
                        Transition = string(trans),
                        Count      = count
                    ))
                end

                # WITH THIS
                raw_stats = DataFrame(Node = String[], Mean = Float64[])
                for r in results
                    for (j, node) in enumerate(NODE_NAMES)
                        push!(raw_stats, (
                            Node = node,
                            Mean = round(r.trajectory_stats["Mean"][j], digits=2)
                        ))
                    end
                end
                grouped = combine(groupby(raw_stats, [:Node, :Mean]), nrow => :Count)
                for row in eachrow(grouped)
                    push!(stats_df, (
                        ParamID    = PARAM_ID,
                        ParamType  = param_type,
                        NoiseLevel = sigma,
                        DT         = dt,
                        Node       = row.Node,
                        Mean       = row.Mean,
                        Count      = row.Count
                    ))
                end
            end
        end
        # WITH THIS
        push!(results_all, results_df)
        push!(transitions_all, transitions_df)
        push!(stats_all, stats_df)
        push!(lambdas_all, lambdas_df)
    end
    output_dir  = joinpath(WORK_DIR, "results_det")
    mkpath(output_dir)
    output_file = joinpath(output_dir, "param_$(PARAM_ID).csv")
    # WITH THIS
    combined_results = vcat(results_all...)
    # Build full grid of all combinations
    all_states    = unique(combined_results.State)
    all_keys      = unique(combined_results[:, [:ParamID, :ParamType, :NoiseLevel, :DT]])
    full_grid     = crossjoin(all_keys, DataFrame(State = all_states))

    # Left join and fill missing with 0
    combined_results = leftjoin(full_grid, combined_results, on = [:ParamID, :ParamType, :NoiseLevel, :DT, :State])
    for col in [:MRT, :MeanSwitches, :StdSwitches]
        combined_results[!, col] = coalesce.(combined_results[!, col], 0.0)
    end

    # Sum over iterations then divide by DET_ITERS
    combined_results = combine(
        groupby(combined_results, [:ParamID, :ParamType, :NoiseLevel, :DT, :State]),
        :MRT          => sum => :MRT,
        :MeanSwitches => sum => :MeanSwitches,
        :StdSwitches  => sum => :StdSwitches
    )
    combined_results.MRT          ./= DET_ITERS
    combined_results.MeanSwitches ./= DET_ITERS
    combined_results.StdSwitches  ./= DET_ITERS
    CSV.write(output_file, combined_results)

    combined_trans = vcat(transitions_all...)

    # Build full grid of all combinations
    all_transitions = unique(combined_trans.Transition)
    all_keys_trans  = unique(combined_trans[:, [:ParamID, :ParamType, :NoiseLevel, :DT]])
    full_grid_trans = crossjoin(all_keys_trans, DataFrame(Transition = all_transitions))

    # Left join and fill missing with 0
    combined_trans = leftjoin(full_grid_trans, combined_trans, on = [:ParamID, :ParamType, :NoiseLevel, :DT, :Transition])
    combined_trans.Count = coalesce.(combined_trans.Count, 0.0)

    # Sum over iterations then divide by DET_ITERS
    combined_trans = combine(
        groupby(combined_trans, [:ParamID, :ParamType, :NoiseLevel, :DT, :Transition]),
        :Count => sum => :Count
    )
    combined_trans.Count ./= DET_ITERS
    CSV.write(replace(output_file, ".csv" => "_transitions.csv"), combined_trans)

    combined_stats = vcat(stats_all...)
    combined_stats = combine(
        groupby(combined_stats, [:ParamID, :ParamType, :NoiseLevel, :DT, :Node, :Mean]),
        :Count => sum => :Count
    )
    CSV.write(replace(output_file, ".csv" => "_stats.csv"), combined_stats)

    CSV.write(replace(output_file, ".csv" => "_lambdas.csv"), vcat(lambdas_all...))

    elapsed = time() - start_time
    println("\n" * "=" ^ 70)
    println("✓ ANALYSIS COMPLETE")
    println("=" ^ 70)
    println("Time   : $(round(elapsed, digits=1)) seconds")
    println("Results: $output_file")
catch e
    println("\n" * "=" ^ 70)
    println("✗ ANALYSIS FAILED")
    println("=" ^ 70)
    println("Error: $e")
    println(stacktrace(catch_backtrace()))
    exit(1)
end