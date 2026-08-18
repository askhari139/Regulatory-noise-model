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
# dataFolder = joinpath(sourceDir, "data")
# dataUniform = joinpath(sourceDir, "dataUniform")
scriptsDir = joinpath(sourceDir, "scripts")
include(joinpath(scriptsDir, "parse_args.jl"))
include(joinpath(scriptsDir, "parallel_run_utils.jl"))
include(joinpath(scriptsDir, "config_utils.jl"))
scriptSubmit = joinpath(scriptsDir, "1_select_parameters.jl")
# ============================================================
# PARSE ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "work-dir" => "work",
    "noise-mode" => "Additive",
    "num-sims" => "25",
    "noise-levels" => "0.0,0.001,0.005,0.01,0.05,0.1",
    "run-mode" => "parallel",
    "num-threads" => "4",
    "log-dir" => "logs", 
    "network" => "TS", 
    "deterministic" => "false",
    "resample" => "true",
    # "use-uniform" => "false",
    "partition" => "ctbp",
    "det-iters" => "10",
    "dt" => "0.01",
))


NOISE_MODE = get_arg(args, "noise-mode", "Additive")
NUM_SIMS = get_arg(args, "num-sims", 25, type=Int)
NOISE_LEVELS = get_arg(args, "noise-levels", [0.0, 0.001, 0.005, 0.01, 0.05, 0.1], type=Vector{Float64})
RUN_MODE = get_arg(args, "run-mode", "parallel")
NUM_THREADS = get_arg(args, "num-threads", 4, type=Int)
LOG_DIR = get_arg(args, "log-dir", "logs")
NETWORK = get_arg(args, "network", "TS")
WORK_DIR = get_arg(args, "work-dir", "work")
Deterministic = get_arg(args, "deterministic", false, type=Bool)
RESAMPLE = get_arg(args, "resample", true, type=Bool)
# USE_UNIFORM= get_arg(args, "use-uniform", false, type=Bool)
PARTITION = get_arg(args, "partition", "ctbp")
DET_ITERS = get_arg(args, "det-iters", 10, type = Int)
DT = get_arg(args, "dt", [0.01], type = Vector{Float64})
# WORK_DIR = joinpath(NOISE_MODE, NETWORK, WORK_DIR)
# ============================================================
# LOAD PARAMETERS
# ============================================================
# dataFolder = get_arg(args, "data-folder")
param_file = joinpath(WORK_DIR, "selected_parameters.jld2")

if !isfile(param_file)
    error("Parameter file not found: $param_file\nRun: julia $scriptSelect first")
end

println("="^70)
println("JOB SUBMISSION")
println("="^70)
println("Loading parameters from: $param_file")

param_data = load_parameter_config(param_file)

all_param_ids = param_data["param_ids"]
param_types = param_data["param_types"]
st_cats = unique(values(param_types))
st_cats_all = ["monostable", "bistable", "tristable", "tetrastable"]
stability_class = [x in st_cats_all for x in st_cats]
NETWORK = param_data["network"]
PARAM_TYPE = Symbol(param_data["param_type"])
dataFolder = param_data["data_folder"]

println("  Network: $NETWORK")
println("  Parameters to analyze: $(length(all_param_ids))")
println("  Noise mode: $NOISE_MODE")
println("  Simulations per level: $NUM_SIMS")
println("="^70)

# ============================================================
# SAVE JOB CONFIGURATION
# ============================================================

job_config = Dict(
    "network" => NETWORK,
    "noise_mode" => NOISE_MODE,
    "param_type" => string(PARAM_TYPE),
    "num_sims" => NUM_SIMS,
    "noise_levels" => NOISE_LEVELS,
    "param_ids" => all_param_ids,
    "param_types" => param_types,
    "timestamp" => string(time()),
    "min_lambdas" => param_data["min_lambdas"],
    "max_lambdas" => param_data["max_lambdas"],
    "dt" => DT
)

mkpath(WORK_DIR)
save_job_config(joinpath(WORK_DIR, "job_config.jld2"), job_config)
# ============================================================
# SUBMIT JOBS
# ============================================================
script = joinpath(scriptsDir, "3_analyze_single_parameter.jl")
if Deterministic
    script = joinpath(scriptsDir, "3_analyze_single_parameter_det.jl")
end
if RUN_MODE == "serial"
    println("\n⚠️  Serial mode - running locally one by one")
    println("This will take approximately $(length(all_param_ids) * 40) seconds")
    println("\nConsider using: run-mode=parallel")
    
    # Run serially
    for param_id in all_param_ids
        cmd = ```
            julia --project=. 
            $script 
            param-id=$param_id
            work-dir=$WORK_DIR
            ```
        println("\nProcessing parameter $param_id...")
        run(cmd)
    end
    
elseif RUN_MODE == "parallel" || RUN_MODE == "shell"
    println("\nSubmitting parallel jobs...")
    function has_slurm()
        try
            run(pipeline(`which sbatch`, stdout=devnull, stderr=devnull))
            return true
        catch
            return false
        end
    end
    # run-mode=shell forces shell-based (subprocess) parallelism even on a
    # SLURM-enabled system -- for when you're already inside a single
    # sbatch allocation with many cores (e.g. sim_shell.script) and want to
    # fan the per-parameter work out across THIS job's cores directly,
    # instead of has_slurm()==true triggering submit_slurm_job_array()'s
    # one-cpu-per-parameter job array (a separate sbatch submission per
    # parameter, which schedules poorly at this parameter count and adds
    # per-job queueing/startup overhead on top of Julia's own per-process
    # startup cost). run-mode=parallel keeps the old auto-detect behavior.
    HAS_SLURM = RUN_MODE == "shell" ? false : has_slurm()

    # For parallel mode, we need stability class info
    # Extract from param_types
    # stability_class = [true, true, false, false]  # Default

    if HAS_SLURM
        println("Using SLURM scheduler")
        mkpath(LOG_DIR)
        submit_slurm_job_array(
            all_param_ids, NOISE_MODE, NETWORK, PARAM_TYPE, NUM_SIMS,
            NOISE_LEVELS, stability_class, LOG_DIR;
            sourceDir = sourceDir, deterministic = Deterministic, resample = RESAMPLE,
            # use_uniform = USE_UNIFORM,
            det_iters = DET_ITERS,
            cpus=1, mem="8G", time="2:00:00",
            partition = PARTITION
        )
    else
        println("Using shell parallelization")
        NUM_THREADS = min(NUM_THREADS, Sys.CPU_THREADS)
        mkpath(LOG_DIR)
        run_parallel_shell(
            all_param_ids, NOISE_MODE, NETWORK, PARAM_TYPE, NUM_SIMS,
            NOISE_LEVELS, stability_class;
            log_dir=LOG_DIR, max_parallel=NUM_THREADS,
            deterministic = Deterministic, resample = RESAMPLE, det_iters = DET_ITERS,
            source_dir = sourceDir
        )
    end
    println("Submitted jobs to use $script file, because Deterministic is set to $Deterministic")
else
    error("Unknown run-mode: $RUN_MODE. Use 'serial', 'parallel', or 'shell'")
end

println("\n" * "="^70)
println("✓ JOB SUBMISSION COMPLETE")
println("="^70)

if RUN_MODE == "parallel"
    if HAS_SLURM
        println("Monitor with: squeue -u \$USER")
        println("Check logs: $LOG_DIR/param_*.out")
    else
        println("Check logs: $LOG_DIR/param_*.log")
    end
    println("\nAfter jobs complete, run: julia scripts/4_collect_results.jl")
else
    println("\nNext step: julia scripts/4_collect_results.jl")
end