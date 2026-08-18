using JLD2
using DataFrames
using CSV
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
dataFolder = joinpath(sourceDir, "data")
dataUniform = joinpath(sourceDir, "dataUniform")
scriptsDir = joinpath(sourceDir, "scripts")
include(joinpath(scriptsDir, "stochastic_setup.jl"))
include(joinpath(scriptsDir, "parse_args.jl"))
include(joinpath(scriptsDir, "config_utils.jl"))
# ============================================================
# PARSE ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "network" => "TS",
    "param-type" => "FoldChange",
    "max-params" => "10",
    "stability-class" => "TTFF",
    "output-dir" => "work",
    "noise-mode" => "Additive",
    "threshold-by" => "Solutions", 
    "param-ids-include" => "1",
    "filter-no-reach" => "false",
    "use-uniform" => "false",
    # "dt" => "0.01",
))

NETWORK = get_arg(args, "network", "TS")
PARAM_TYPE = get_arg(args, "param-type", "FoldChange", type=Symbol)
MAX_PARAMS_PER_TYPE = get_arg(args, "max-params", 10, type=Int)
STABILITY_CLASS = get_arg(args, "stability-class", [true, true, false, false], type=Vector{Bool})
OUTPUT_DIR = get_arg(args, "output-dir", "work")
NOISE_MODE = get_arg(args, "noise-mode", "Additive")
THRESHOLD_BY = get_arg(args, "threshold-by", "Solutions")
# ADDITIONAL_PARAMS = get_arg(args, "param-ids-include", "1")
# OUTPUT_DIR = joinpath(NOISE_MODE, NETWORK, OUTPUT_DIR)
FILTER_NO_REACH = get_arg(args, "filter-no-reach", false; type=Bool)
USE_UNIFORM= get_arg(args, "use-uniform", false, type=Bool)
# DT=get_arg(args, "dt", 0.01, type=Float64)
stability_cats = ["monostable", "bistable", "tristable", "tetrastable"]
param_req = Dict(stability_cats .=> STABILITY_CLASS)

if USE_UNIFORM
    dataFolder = dataUniform
end

println("="^70)
println("PARAMETER SELECTION")
println("="^70)
println("Network: $NETWORK")
println("Parameter type: $PARAM_TYPE")
println("Max parameters per type: $MAX_PARAMS_PER_TYPE")
println("Stability classes: $STABILITY_CLASS")
println("="^70)

# ============================================================
# LOAD RACIPE DATA
# ============================================================
original_dir = pwd()
cd(dataFolder)

PRS_FILE = NETWORK*".prs"
TOPO_FILE = NETWORK*".topo"
PARAMETERS_FILE = NETWORK*"_parameters.dat"
SOLUTIONS_FILE = NETWORK*"_solution.dat"
ODE_FILE = joinpath(sourceDir,"src",NETWORK*"_ode.jl")

# Run RACIPE if needed
if !isfile(PARAMETERS_FILE)
    println("\nRunning RACIPE...")
    run_racipe(TOPO_FILE, num_paras=10000, threads=threads, keep_only_essentials=true)
end

# Generate ODE file
generate_ode_function(TOPO_FILE, PRS_FILE, ODE_FILE)

# Load data
println("\nLoading RACIPE data...")
prs = read_prs(PRS_FILE)
topo = read_topo(TOPO_FILE)
params_data = read_parameters(PARAMETERS_FILE, prs)
solutions_data = read_solutions(SOLUTIONS_FILE, prs)

include(ODE_FILE)
NODE_NAMES = solutions_data.node_names
println("  Nodes: $NODE_NAMES")
println("  Total parameter sets: $(nrow(params_data.data))")

# Get lambda indices
lambda_indices, max_lambdas, min_lambdas = get_lambda_indices_from_topology(prs.names; param_type=PARAM_TYPE)

# ============================================================
# IDENTIFY ATTRACTORS
# ============================================================

println("\n" * "="^70)
println("IDENTIFYING ATTRACTORS")
println("="^70)

# Calculate thresholds

# if THRESHOLD_BY == "Solutions"
racipe_means = get_mean_expression(solutions_data, weight_by_frequency=true)
racipe_thresholds = [first(racipe_means[racipe_means.Node .== node, :MeanExpression]) 
                    for node in NODE_NAMES]
#     df = DataFrame(ParamID = params_data.data.ParamID)
#     for (i, node) in enumerate(NODE_NAMES)
#         df[!, Symbol(node)] .= racipe_thresholds[i]
#     end
#     racipe_thresholds = copy(df)
# else 
#     racipe_thresholds = get_mean_expressions(params_data)
# end

if FILTER_NO_REACH
    n_before = nrow(params_data.data)
    par_keep = Int[]
    for row in eachrow(params_data.data)
        uppers = [row[Symbol("Prod_of_$nd")] / row[Symbol("Deg_of_$nd")]
                  for nd in NODE_NAMES]
        # thresholds_this = collect(racipe_thresholds[
        #     racipe_thresholds.ParamID .== row.ParamID, Symbol.(NODE_NAMES)][1, :])
        if all(uppers .> racipe_thresholds)
            push!(par_keep, row.ParamID)
        end
    end
    keep_set = Set(par_keep)
    params_data = Parameters(
        filter(r -> r.ParamID ∈ keep_set, params_data.data),
        params_data.param_names
    )
    solutions_data = Solutions(
        filter(r -> r.ParamID ∈ keep_set, solutions_data.data),
        solutions_data.node_names
    )
    # racipe_thresholds = filter(r -> r.ParamID ∈ keep_set, racipe_thresholds)
    discrete_df = discretize_racipe_states(solutions_data; thresholds=racipe_thresholds)
    println("  After reachability filter: $(length(par_keep)) / $n_before parameter sets retained")
end
# println("\nThresholds for discretization:")
# for (i, node) in enumerate(NODE_NAMES)
#     println("  $node: $(round(racipe_thresholds[i], digits=3))")
# end

# Identify attractors
attractors = identify_attractors(solutions_data; thresholds=racipe_thresholds)
println("\nFound $(length(attractors)) unique attractor structures")

# Classify by stability
classified = classify_attractors_by_stability(attractors)

param_requirement = Dict(
    cls => requested && haskey(classified, cls)
    for (cls, requested) in param_req
)

println("\nAttractor classification:")
for (class, attractor_list) in collect(classified)
    total_params = sum(length(a[2]["param_ids"]) for a in attractor_list)
    println("  $class: $total_params parameters in $(length(attractor_list)) attractor(s)")
end

# ============================================================
# SAMPLE PARAMETERS
# ============================================================

println("\n" * "="^70)
println("SAMPLING PARAMETERS")
println("="^70)

param_types = Dict{Int, String}()
all_param_ids = Int[]

for stability_class in stability_cats
    if haskey(classified, stability_class) && param_requirement[stability_class]
        sampled = sample_balanced_parameters(
            attractors,
            stability_class,
            MAX_PARAMS_PER_TYPE,
            weight_by_frequency=true,
            prefer_balanced=true
        )
        
        println("\n$stability_class: Sampled $(length(sampled)) parameters")
        
        for param_id in sampled
            param_types[param_id] = stability_class
        end
        append!(all_param_ids, sampled)
    end
end

# append!(all_param_ids, ADDITIONAL_PARAMS)
# for param in ADDITIONAL_PARAMS
#     param_types[param] = "Additional"
# end
sort!(all_param_ids)

println("\nTotal parameters selected: $(length(all_param_ids)) from folder $(dataFolder)")

# ============================================================
# SAVE CONFIGURATION
# ============================================================

cd(original_dir)
mkpath(OUTPUT_DIR)

# Save parameter list
param_data = Dict(
    "network" => NETWORK,
    "param_ids" => all_param_ids,
    "param_types" => param_types,
    "node_names" => NODE_NAMES,
    "lambda_indices" => lambda_indices,
    "max_lambdas" => max_lambdas,
    "min_lambdas" => min_lambdas,
    "racipe_thresholds" => racipe_thresholds,
    "param_type" => string(PARAM_TYPE),
    "timestamp" => string(time()), 
    "data_folder" => dataFolder
)
println(param_data)

save_parameter_config(joinpath(OUTPUT_DIR, "selected_parameters.jld2"), param_data)

println("\n" * "="^70)
println("✓ PARAMETER SELECTION COMPLETE")
println("="^70)
println("Saved to: $(OUTPUT_DIR)/selected_parameters.jld2")
println("\nNext step: julia scripts/2_submit_jobs.jl")