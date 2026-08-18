using JSON
using Dates

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
include(joinpath(scriptsDir , "parse_args.jl"))
script0 = joinpath(scriptsDir, "0_run_all.jl")
scriptExisting = joinpath(scriptsDir, "0_run_existing.jl")
scriptSelect = joinpath(scriptsDir, "1_select_parameters.jl")
scriptSubmit = joinpath(scriptsDir, "2_submit_jobs.jl")
scriptAnalyze = joinpath(scriptsDir, "3_analyze_single_parameter.jl")
scriptsAnalyzeDet = joinpath(scriptsDir, "3_analyze_single_parameter_det.jl")
scriptsMonitor = joinpath(scriptsDir, "monitor_and_collect.jl")
scriptsCollect = joinpath(scriptsDir, "4_collect_results.jl")
scriptsPlot = joinpath(scriptsDir, "5_plot_results.jl")

# ============================================================
# PARSE ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "networks"         => "TS",
    "noise-modes"      => "Additive",
    "param-type"       => "FoldChange",
    "max-params"       => "10",
    "stability-class"  => "TTFF",
    "num-sims"         => "25",
    "noise-levels"     => "0.0,0.001,0.005,0.01,0.05,0.1",
    "saveat"           => "10.0",
    "run-mode"         => "parallel",
    "num-threads"      => "4",
    "poll-interval"    => "900",   # seconds; passed to monitor_and_collect.jl
    "max-wait"         => "86400",
    "collect-only"     => "false",
    "status"           => "false",
    "deterministic"    => "false",
    "cleanup"          => "true",
    "det-iters"        => "100",
    "DT"               => "0.01"
))

NETWORKS          = get_arg(args, "networks",        ["TS"],   type=Vector{String})
NOISE_MODES       = get_arg(args, "noise-modes",     ["Additive"], type=Vector{String})
PARAM_TYPE        = get_arg(args, "param-type",      "FoldChange")
MAX_PARAMS        = get_arg(args, "max-params",      10,  type=Int)
STABILITY_CLASS   = get_arg(args, "stability-class", [true, true, false, false], type=Vector{Bool})
NUM_SIMS          = get_arg(args, "num-sims",        25,  type=Int)
NOISE_LEVELS      = get_arg(args, "noise-levels",    [0.0, 0.001, 0.005, 0.01, 0.05, 0.1], type=Vector{Float64})
SAVEAT            = get_arg(args, "saveat",          10.0, type=Float64)
RUN_MODE          = get_arg(args, "run-mode",        "parallel")
NUM_THREADS       = get_arg(args, "num-threads",     4,   type=Int)
POLL_INTERVAL     = get_arg(args, "poll-interval",   900, type=Int)
MAX_WAIT          = get_arg(args, "max-wait",        86400, type=Int)
COLLECT_ONLY      = get_arg(args, "collect-only",    false, type=Bool)
STATUS_CHECK      = get_arg(args, "status",          false, type=Bool)
Deterministic     = get_arg(args, "deterministic",   false, type=Bool)
CLEANUP           = get_arg(args, "cleanup", true, type=Bool)
DET_ITERS         = get_arg(args, "det-iters", 100, type = Int)
DT                = get_arg(args, "DT", 0.01, type = Vector{Float64})
combinations = [(net, mode) for net in NETWORKS for mode in NOISE_MODES]
work_dirs    = [joinpath(mode, net) for (net, mode) in combinations]

println("=" ^ 70)
println("MASTER ORCHESTRATION SCRIPT")
println("=" ^ 70)
println("Networks:    $(join(NETWORKS,    ", "))")
println("Noise modes: $(join(NOISE_MODES, ", "))")
println("Combinations: $(length(combinations))")
println("=" ^ 70)

# if Deterministic
#     for (NOISE_MODE,dt) in [(x, y) for x in NOISE_MODES for y in DT]
#         jld = joinpath(dataFolder,
#                        "lambda_hist_$(NOISE_MODE)_dt$(dt).jld2")
#         if !isfile(jld)
#             lambda_inh = build_lambda_cache(NOISE_MODE, NOISE_LEVELS,
#                                             0.01, 1.0; dt=dt)
#             lambda_act = build_lambda_cache(NOISE_MODE, NOISE_LEVELS,
#                                             1.0, 100.0; dt=dt)
#             @save jld lambda_inh lambda_act
#         end
#     end
# end

# ============================================================
# STEP 2 — SUBMIT JOBS
# ============================================================

println("\n" * "=" ^ 70)
println("STEP 2 / 2 — SUBMITTING JOBS - STOCHASTIC")
println("=" ^ 70)
failed_step2 = String[]

for (work_dir, (network, noise_mode)) in zip(work_dirs, combinations)
    log_dir = joinpath(work_dir, "logs")
    println("\n--- $noise_mode / $network ---")
    if CLEANUP
        if Deterministic
            res_dir = joinpath(work_dir, "results_det")
        else
            res_dir = joinpath(work_dir, "results")
        end
        if ispath(res_dir)
            rm(res_dir, recursive=true, force=true)
        end
    end

    cmd = `julia --project=. $scriptSubmit \
            work-dir=$work_dir \
            noise-mode=$noise_mode \
            network=$network \
            num-sims=$NUM_SIMS \
            noise-levels=$(join(NOISE_LEVELS, ",")) \
            run-mode=$RUN_MODE \
            deterministic=$Deterministic \
            num-threads=$NUM_THREADS \
            log-dir=$log_dir \
            saveat=$SAVEAT \
            det-iters=$DET_ITERS \
            dt=$(join(DT, ","))`
    try
        run(cmd)
        println("  ✓ Jobs submitted")
    catch e
        @error "Job submission failed for $noise_mode/$network: $e"
        push!(failed_step2, work_dir)
    end
end

if !isempty(failed_step2)
    msg = "❌ Job submission failed for: $(join(failed_step2, ", ")). Aborting monitor."
    @error msg
    send_discord_message(JSON.json(Dict("content" => msg)))
    exit(1)
end

# ============================================================
# HAND OFF TO MONITOR
# ============================================================

println("\n" * "=" ^ 70)
println("HANDING OFF TO MONITOR")
println("=" ^ 70)

work_dirs_arg = join(work_dirs, ",")

monitor_cmd = `julia --project=. $scriptsMonitor \
                work-dirs=$work_dirs_arg \
                poll-interval=$POLL_INTERVAL \
                deterministic=$Deterministic \
                max-wait=$MAX_WAIT`

println("Starting monitor (poll every $(POLL_INTERVAL)s, max wait $(MAX_WAIT)s)...")
println("Command: $monitor_cmd")

if RUN_MODE == "parallel" && occursin("slurm", lowercase(get(ENV, "SCHEDULER", "")))
    # FIX #18: redirect stdout/stderr to a log file so output isn't lost when detached
    log_path = "monitor_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log"
    monitor_proc = run(pipeline(monitor_cmd,
                                stdout=log_path,
                                stderr=log_path);
                    wait=false)   # detach without blocking
    println("✓ Monitor detached — output → $log_path  (PID: $(getpid()))")
else
    # Local / shell-parallel: run inline (jobs already blocking or finished)
    run(monitor_cmd)
end

println("\n" * "=" ^ 70)
println("✅ 0_run_all.jl COMPLETE")
println("=" ^ 70)