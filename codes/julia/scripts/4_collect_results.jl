using CSV
using DataFrames
using JSON

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
scriptDir = joinpath(sourceDir, "scripts")
include(joinpath(scriptDir, "parse_args.jl"))
include(joinpath(scriptDir, "config_utils.jl"))

# ============================================================
# PARSE ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "work-dir" => "work",
    "output-dir" => "results",
    "network" => "TS",
    "noise-mode" => "Additive",
    "deterministic" => "false"
))
# NOISE_MODE = get_arg(args, "noise-mode", "Additive")
# NETWORK = get_arg(args, "network", "TS")
WORK_DIR = get_arg(args, "work-dir", "work")
println(WORK_DIR)
# WORK_DIR = joinpath(NOISE_MODE, NETWORK, WORKDIR)


# OUTPUT_DIR = joinpath(NOISE_MODE, NETWORK, OUTPUT_DIR)
DETERMINISTIC = get_arg(args, "deterministic",   false, type=Bool)
if DETERMINISTIC
    OUTPUT_DIR = joinpath(WORK_DIR, "results_det")
else
    OUTPUT_DIR = joinpath(WORK_DIR, "results")
end
println(OUTPUT_DIR)
println("="^70)
println("COLLECTING RESULTS")
println("="^70)

# ============================================================
# LOAD CONFIGURATION
# ============================================================

param_file = joinpath(WORK_DIR, "selected_parameters.jld2")
if !isfile(param_file)
    error("Parameter file not found: $param_file")
end

param_data = load_parameter_config(param_file)
all_param_ids = param_data["param_ids"]
println("Expected parameters: $(length(all_param_ids))")

# ============================================================
# COLLECT RESULTS
# ============================================================

results_dir = OUTPUT_DIR

if !isdir(results_dir)
    error("Results directory not found: $results_dir\nRun analysis first")
end

result_files = [joinpath(results_dir, "param_$(id).csv") for id in all_param_ids]
transition_files = [joinpath(results_dir, "param_$(id)_transitions.csv") for id in all_param_ids]
stats_files = [joinpath(results_dir, "param_$(id)_stats.csv") for id in all_param_ids]
if DETERMINISTIC
    lambda_files = [joinpath(results_dir, "param_$(id)_lambdas.csv") for id in all_param_ids]
end
found_files = filter(isfile, result_files)
found_files_trans = filter(isfile, transition_files)
found_files_stats = filter(isfile, stats_files)

println("Found results: $(length(found_files)) / $(length(all_param_ids))")
println("Found transitions: $(length(found_files_trans)) / $(length(all_param_ids))")
println("Found stats: $(length(found_files_stats)) / $(length(all_param_ids))")

if length(found_files) == 0
    error("No result files found!")
end

if length(found_files) < length(all_param_ids)
    missing = setdiff(all_param_ids, [parse(Int, match(r"param_(\d+)", f).captures[1]) for f in found_files])
    @warn "Missing results for parameters: $missing"
end

# Load and combine
println("\nCombining results...")
# all_results = DataFrame()
# transition_results = DataFrame()
# stats_results = DataFrame()
# for file in found_files
#     df = CSV.read(file, DataFrame)
#     if any(ismissing, df.ParamID)
#         @warn "Bad file: $file"
#         continue
#     end
#     append!(all_results, df, promote = true)
# end

# for file in found_files_trans
#     df = CSV.read(file, DataFrame)
#     append!(transition_results, df, promote = true)
# end

# for file in found_files_stats
#     df = CSV.read(file, DataFrame)
#     append!(stats_results, df, promote = true)
# end

# println("Total rows: $(nrow(all_results))")

function stich_files(file_vector, output_name)
    open(output_name, "w") do output_io
        # Handle the very first file (keep the header)
        if !isempty(file_vector) && isfile(file_vector[1])
            open(file_vector[1], "r") do input_io
                write(output_io, input_io) # Blazing fast raw byte copy
            end
            write(output_io, '\n')
        end

        # Handle the thousands of remaining files (drop the header)
        for i in 2:length(file_vector)
            file = file_vector[i]
            isfile(file) || continue
            
            open(file, "r") do input_io
                if !eof(input_io)
                    readline(input_io) # Reads and discards the header line instantly
                    write(output_io, input_io) # Streams the remaining bytes directly to disk
                end
            end
            write(output_io, '\n')
        end
    end
    for f in file_vector
        rm(f)
    end
end

# ============================================================
# SAVE COMBINED RESULTS
# ============================================================

mkpath(OUTPUT_DIR)
output_file = joinpath(OUTPUT_DIR, "all_parameters_results.csv")
stich_files(found_files, output_file)
stich_files(found_files_trans, replace(output_file, "_results.csv" => "_transitions.csv"))
stich_files(found_files_stats, replace(output_file, "_results.csv" => "_stats.csv"))
if DETERMINISTIC
    stich_files(lambda_files, replace(output_file, "_results.csv" => "_lambdas.csv"))
end


# CSV.write(output_file, all_results)
# CSV.write(replace(output_file, "_results.csv" => "_transitions.csv"), transition_results)
# CSV.write(replace(output_file, "_results.csv" => "_stats.csv"),       stats_results)
println("\n" * "="^70)
println("✓ RESULTS COLLECTED")
println("="^70)
println("Saved to: $output_file")
println("Deleting individual files ...")

println("\nNext step: julia scripts/5_create_plots.jl")