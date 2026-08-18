using JSON
using DataFrames
using CSV
using Statistics

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

# ============================================================
# PARSE ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "param-id"   => "0",
    "work-dir"   => ".",
    "network"    => "TS",
    "noise-mode" => "Additive",
    "resample"   => "true",
    # "use-uniform" => "false",
))

PARAM_ID   = get_arg(args, "param-id",   0,         type=Int)
NOISE_MODE = get_arg(args, "noise-mode", "Additive")
NETWORK    = get_arg(args, "network",    "TS")
RESAMPLE   = get_arg(args, "resample",   true,      type=Bool)
# USE_UNIFORM= get_arg(args, "use-uniform", false, type=Bool)

# FIX #11: actually use the work-dir argument instead of hard-coding it
WORK_DIR = get_arg(args, "work-dir", joinpath(NOISE_MODE, NETWORK))
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

# u0List = Vector{Vector{Float64}}()
# for row in eachrow(solutions_data.data[solutions_data.data.ParamID .== PARAM_ID, :])
#     push!(u0List, [2.0^row[Symbol(node)] for node in NODE_NAMES])
# end

upper_limits = [param_row[1, Symbol("Prod_of_"*node)]/param_row[1, Symbol("Deg_of_"*node)] for node in NODE_NAMES]
A = rand(NUM_SIMS, length(NODE_NAMES))
u0List = [A[i, :].*upper_limits for i in 1:size(A, 1)]

# ============================================================
# RUN ANALYSIS
# ============================================================

println("\nRunning stochastic simulations...")
start_time = time()

TSPAN = (0.0, 1000.0)

try
    # if RESAMPLE
    #     for (l, i) in enumerate(lambda_indices)
    #         if min_lambdas[l] < 1
    #             # println("Resampling lambda[$i] from range: $(min_lambdas[i]) - $(max_lambdas[i])")
    #             param_vector[i] = rand() * (max_lambdas[l] - min_lambdas[l]) + min_lambdas[l]
    #         end
    #     end
    # end
    y1 = @elapsed all_results, summary_df = analyze_noise_effects(
        ode_system!,
        param_vector,
        prs.names,
        lambda_indices,
        NODE_NAMES;
        sigma_values         = NOISE_LEVELS,
        num_sims         = NUM_SIMS,
        tspan            = TSPAN,
        u0List           = u0List,
        racipe_thresholds= racipe_thresholds,
        max_lambdas           = max_lambdas,
        min_lambdas           = min_lambdas,
        noise_mode       = NOISE_MODE,
        saveat           = 1.0,
        dt_vals          = DT_vals
    )
    println("Simulations took $(round(y1, digits=2)) seconds")

    
    println("aggregate_mrt took $(round(y1, digits=2)) seconds")

    # ---- Results DataFrame ----
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

    transitions_df = DataFrame(
        ParamID    = Int[],
        ParamType  = String[],
        NoiseLevel = Float64[],
        DT         = Float64[],
        Transition = String[],
        Count      = Float64[]
    )

    stats_df = DataFrame(
        ParamID      = Int[],
        ParamType    = String[],
        NoiseLevel   = Float64[],
        DT           = Float64[],
        Node         = String[],
        Mean         = Float64[],
        StDev        = Float64[],
        Count        = Int[],
        MeanSkewness = Float64[],
        MeanKurtosis = Float64[]
    )

    for dt in DT_vals
        deterministic_mrt = aggregate_mrt(all_results[(0.0, dt)])
        for sigma in NOISE_LEVELS
            results   = all_results[(sigma, dt)]
            mean_mrt  = aggregate_mrt(results)
            mean_counts = aggregate_transitions(results)
            switches  = [r.switches for r in results]

            for (state, mrt) in mean_mrt
                det_mrt = get(deterministic_mrt, state, 0.0)
                # FIX #12: consistent RelativeMRT formula — log2 ratio used throughout
                rel_mrt = 2.0^(mrt - det_mrt)
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
            raw_stats = DataFrame(
                Node     = String[],
                Mean     = Float64[],
                StDev    = Float64[],
                Skewness = Float64[],
                Kurtosis = Float64[]
            )
            for r in results
                for (j, node) in enumerate(NODE_NAMES)
                    push!(raw_stats, (
                        Node     = node,
                        Mean     = round(r.trajectory_stats["Mean"][j],     digits=2),
                        StDev    = round(r.trajectory_stats["StDev"][j],    digits=2),
                        Skewness = r.trajectory_stats["Skewness"][j],
                        Kurtosis = r.trajectory_stats["Kurtosis"][j],
                    ))
                end
            end
            grouped = combine(groupby(raw_stats, [:Node, :Mean, :StDev]),
                nrow         => :Count,
                :Skewness    => mean => :MeanSkewness,
                :Kurtosis    => mean => :MeanKurtosis
            )
            for row in eachrow(grouped)
                push!(stats_df, (
                    ParamID      = PARAM_ID,
                    ParamType    = param_type,
                    NoiseLevel   = sigma,
                    DT           = dt,
                    Node         = row.Node,
                    Mean         = row.Mean,
                    StDev        = row.StDev,
                    Count        = row.Count,
                    MeanSkewness = row.MeanSkewness,
                    MeanKurtosis = row.MeanKurtosis
                ))
            end
        end
    end

    output_dir  = joinpath(WORK_DIR, "results")
    mkpath(output_dir)
    output_file = joinpath(output_dir, "param_$(PARAM_ID).csv")
    CSV.write(output_file, results_df)
    CSV.write(replace(output_file, ".csv" => "_transitions.csv"), transitions_df)
    CSV.write(replace(output_file, ".csv" => "_stats.csv"),       stats_df)

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