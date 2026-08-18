"""
Master orchestration script

Usage:
    julia scripts/0_run_all.jl \\
        networks=TS,EMT \\
        noise-modes=Additive,OU \\
        run-mode=parallel

    # Check status manually
    julia scripts/0_run_all.jl ... status=true

    # Collect without re-running
    julia scripts/0_run_all.jl ... collect-only=true
"""

using JSON
using Dates
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

include(joinpath(scriptsDir, "parse_args.jl"))
include(joinpath(scriptsDir,"discord_notifier.jl"))
include(joinpath(scriptsDir, "lambda_sampler.jl"))

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
    "resample"         => "false",
    "filter-no-reach"  => "false",
    "run-plots"        => "false",
    "use-uniform"      => "false",
    "partition"        => "ctbp",
    "det-iters"        => "10",
    "dt"               => "0.01",
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
RESAMPLE          = get_arg(args, "resample",        true,  type=Bool)
FILTER_NO_REACH   = get_arg(args, "filter-no-reach", false, type = Bool)
RUN_PLOTS         = get_arg(args, "run-plots",        false, type = Bool)
USE_UNIFORM       = get_arg(args, "use-uniform",     false, type = Bool)
PARTITION         = get_arg(args, "partition", "ctbp")
DET_ITERS         = get_arg(args, "det-iters", 10, type = Int)
DT                = get_arg(args, "dt", [0.01], type = Vector{Float64})

combinations = [(net, mode) for net in NETWORKS for mode in NOISE_MODES]
work_dirs    = [joinpath(mode, net) for (net, mode) in combinations]
if USE_UNIFORM
    dataFolder = dataUniform
end
println("=" ^ 70)
println("MASTER ORCHESTRATION SCRIPT")
println("=" ^ 70)
println("Networks:    $(join(NETWORKS,    ", "))")
println("Noise modes: $(join(NOISE_MODES, ", "))")
println("Combinations: $(length(combinations))")
println("=" ^ 70)

function check_move(src, tgt)
    if ispath(src) mv(src, tgt, force = true) end 
end

function get_lambda_dist(noise_mode::String, noise_level::Float64, min_lambda::Float64, max_lambda::Float64)
    lambda_init = collect(1:100)*1.0
    if min_lambda < 1
        lambda_init = lambda_init./100
    end
    # lambda_init_key = round.(lambda_init, digits = 0)
    lambda_trajectories = Dict{Float64, Vector{Float64}}()
    for i in eachindex(lambda_init)
        lambda = lambda_init[i]
        key = copy(lambda_init[i])
        lt = zeros(10000)
        trajectories = Float64[]
        for iter in 1:100
            noise = noise_level*randn(10000)
            lt[1] = lambda
            for t in 2:10000
                if noise_mode == "Additive"
                    ln = lt[t-1] + noise[t]
                elseif noise_mode == "Multiplicative"
                    ln = lt[t-1]*(1 + noise[t])
                elseif noise_mode == "Fluctuating"
                    ln = lambda + noise[t]
                elseif noise_mode == "Jumping"
                    ln = lt[t-1] + (noise[t] > 0 ? noise_level : -1*noise_level)
                elseif noise_mode == "Extreme"
                    ln = lt[t-1] + (noise[t] > 0 ? 1 : -1)
                end
                lt[t] = max(min(ln, max_lambda), min_lambda)
            end
            trajectories = vcat(trajectories, lt[5000:end])
        end
        lambda_trajectories[key] = trajectories
    end
    return(lambda_trajectories)
end


# ============================================================
# STATUS CHECK MODE  (manual invocation only)
# ============================================================

if STATUS_CHECK
    include(sourceDir * "/scripts/config_utils.jl")

    println("\n" * "=" ^ 70)
    println("STATUS CHECK")
    println("=" ^ 70)
    println("\n" * rpad("Network", 12) * rpad("Noise Mode", 18) * rpad("Results", 14) * "Status")
    println("-" ^ 60)

    for (work_dir, (network, noise_mode)) in zip(work_dirs, combinations)
        param_file = joinpath(work_dir, "selected_parameters.jld2")
        results_dir = joinpath(work_dir, "results")

        if !isfile(param_file)
            status = "❌ Not started"
            results_str = "-"
        else
            param_data = load_parameter_config(param_file)
            expected = length(param_data["param_ids"])

            n_done = isdir(results_dir) ?
                length(filter(
                    f -> startswith(f, "param_") && endswith(f, ".csv") && !occursin("transitions", f),
                    readdir(results_dir)
                )) : 0

            results_str = "$n_done/$expected"
            has_combined = isfile(joinpath(results_dir, "all_parameters_results.csv"))

            status = if n_done == 0
                "⏳ Waiting"
            elseif n_done < expected
                "⏳ In progress"
            elseif !has_combined
                "⚠️  Ready to collect"
            else
                "✅ Complete"
            end
        end

        println(rpad(network, 12) * rpad(noise_mode, 18) * rpad(results_str, 14) * status)
    end

    println("\n" * "=" ^ 70)
    exit(0)
end

# ============================================================
# COLLECT-ONLY MODE
# ============================================================

if COLLECT_ONLY
    println("\n" * "=" ^ 70)
    println("COLLECT & PLOT MODE")
    println("=" ^ 70)

    work_dirs_arg = join(work_dirs, ",")
    collect_cmd   = `julia --project=. $scriptsMonitor \
                        work-dirs=$work_dirs_arg \
                        poll-interval=0 \
                        max-wait=0 \
                        deterministic=$Deterministic`
    # poll-interval=0 / max-wait=0 → the monitor will see jobs already done
    # and go straight to collection without sleeping
    run(collect_cmd)
    exit(0)
end

if Deterministic
    for (NOISE_MODE,dt) in [(x, y) for x in NOISE_MODES for y in DT]
        jld = joinpath(dataFolder,
                       "lambda_hist_$(NOISE_MODE)_dt$(dt).jld2")
        if !isfile(jld)
            lambda_inh = build_lambda_cache(NOISE_MODE, NOISE_LEVELS,
                                            0.01, 1.0; dt=dt)
            lambda_act = build_lambda_cache(NOISE_MODE, NOISE_LEVELS,
                                            1.0, 100.0; dt=dt)
            @save jld lambda_inh lambda_act
        end
    end
end

# ============================================================
# FULL PIPELINE: STEP 1 — SELECT PARAMETERS
# ============================================================

println("\n" * "=" ^ 70)
println("STEP 1 / 2 — SELECTING PARAMETERS")
println("=" ^ 70)

notify_pipeline_start(NETWORKS, NOISE_MODES)

failed_step1 = String[]

for (work_dir, (network, noise_mode)) in zip(work_dirs, combinations)
    println("\n--- $noise_mode / $network ---")
    if ispath(work_dir)
        foreach(x -> rm(x; recursive=true, force=true), readdir(work_dir, join=true))
    else
        mkpath(work_dir)
    end

    cmd = `julia --project=. $scriptSelect \
            network=$network \
            noise-mode=$noise_mode \
            param-type=$PARAM_TYPE \
            max-params=$MAX_PARAMS \
            stability-class=$(join([x ? "T" : "F" for x in STABILITY_CLASS], "")) \
            output-dir=$work_dir \
            filter-no-reach=$FILTER_NO_REACH \
            use-uniform=$USE_UNIFORM`
    try
        run(cmd)
        println("  ✓ Parameters selected")
    catch e
        @error "Parameter selection failed for $noise_mode/$network: $e"
        push!(failed_step1, work_dir)
    end
end

if !isempty(failed_step1)
    # Notify and abort — no point submitting jobs with missing configs
    msg = "❌ Parameter selection failed for: $(join(failed_step1, ", ")). Aborting."
    @error msg
    send_discord_message(JSON.json(Dict("content" => msg)))
    exit(1)
end

# ============================================================
# STEP 2 — SUBMIT JOBS
# ============================================================

println("\n" * "=" ^ 70)
println("STEP 2 / 2 — SUBMITTING JOBS")
println("=" ^ 70)
failed_step2 = String[]

for (work_dir, (network, noise_mode)) in zip(work_dirs, combinations)
    log_dir = joinpath(work_dir, "logs")
    println("\n--- $noise_mode / $network ---")

    cmd = `julia --project=. $scriptSubmit \
            work-dir=$work_dir \
            noise-mode=$noise_mode \
            network=$network \
            num-sims=$NUM_SIMS \
            noise-levels=$(join(NOISE_LEVELS, ",")) \
            run-mode=$RUN_MODE \
            num-threads=$NUM_THREADS \
            log-dir=$log_dir \
            saveat=$SAVEAT \
            deterministic=$Deterministic \
            resample=$RESAMPLE \
            partition=$PARTITION \
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
                max-wait=$MAX_WAIT \
                run-plots=$RUN_PLOTS \
                deterministic=$Deterministic`

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

if Deterministic
    # for (work_dir, (network, noise_mode)) in zip(work_dirs, combinations)
    #     original_res = joinpath(work_dir, "results")
    #     original_figs = joinpath(work_dir, "figures")
    #     check_move(original_figs, original_figs*"_det")
    #     check_move(original_res, original_res*"_det")
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

        cmd = `julia --project=. $scriptSubmit \
                work-dir=$work_dir \
                noise-mode=$noise_mode \
                network=$network \
                num-sims=$NUM_SIMS \
                noise-levels=$(join(NOISE_LEVELS, ",")) \
                run-mode=$RUN_MODE \
                num-threads=$NUM_THREADS \
                log-dir=$log_dir \
                saveat=$SAVEAT \
                partition=$PARTITION \
                deterministic=false \
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
                    max-wait=$MAX_WAIT \
                    run-plots=$RUN_PLOTS \
                    deterministic=false`

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

end