"""
monitor_and_collect.jl

Polls job progress across all work directories, sends Discord updates every
15 minutes, then triggers result collection once all simulations are complete.

Usage (called automatically by 0_run_all.jl, or standalone):
    julia scripts/monitor_and_collect.jl \\
        work-dirs=Additive/TS,OU/TS \\
        poll-interval=900 \\
        max-wait=86400
"""

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
include(joinpath(scriptsDir, "parse_args.jl"))
include(joinpath(scriptsDir, "discord_notifier.jl"))  # defines fmt_duration, discord_embed, send_discord_message
include(joinpath(scriptsDir, "config_utils.jl"))
script0 = joinpath(scriptsDir, "0_run_all.jl")
scriptExisting = joinpath(scriptsDir, "0_run_existing.jl")
scriptSelect = joinpath(scriptsDir, "1_select_parameters.jl")
scriptSubmit = joinpath(scriptsDir, "2_submit_jobs.jl")
scriptAnalyze = joinpath(scriptsDir, "3_analyze_single_parameter.jl")
scriptsAnalyzeDet = joinpath(scriptsDir, "3_analyze_single_parameter_det.jl")
scriptsMonitor = joinpath(scriptsDir, "monitor_and_collect.jl")
scriptsCollect = joinpath(scriptsDir, "4_collect_results.jl")
scriptsPlot = joinpath(scriptsDir, "5_plot_results.jl")
# FIX #19: fmt_duration and discord_embed are defined in discord_notifier.jl — do NOT redefine here

# ============================================================
# PARSE ARGUMENTS
# ============================================================

args = parse_named_args(defaults=Dict(
    "work-dirs"      => "work",
    "poll-interval"  => "900",    # seconds between polls (default: 15 min)
    "max-wait"       => "86400",  # bail-out after 24 h; 0 = no timeout (collect-only)
    "collect-script" => scriptsCollect,
    "plot-script"    => scriptsPlot,
    "run-plots"      => "false",
    "deterministic"  => "false",
))

WORK_DIRS      = get_arg(args, "work-dirs",      ["work"],   type=Vector{String})
POLL_INTERVAL  = get_arg(args, "poll-interval",  900,        type=Int)
MAX_WAIT       = get_arg(args, "max-wait",        86400,      type=Int)   # 0 = unlimited
COLLECT_SCRIPT = get_arg(args, "collect-script", scriptsCollect)
PLOT_SCRIPT    = get_arg(args, "plot-script",    scriptsPlot)
RUN_PLOTS      = get_arg(args, "run-plots",       false; type=Bool)
DETERMINISTIC  = get_arg(args, "deterministic", false; type=Bool)

# ============================================================
# HELPERS
# ============================================================

"""Count result files present in a work directory."""
function count_results(work_dir::String; DETERMINISTIC::Bool=DETERMINISTIC)
    if DETERMINISTIC
        results_dir = joinpath(work_dir, "results_det")
    else
        results_dir = joinpath(work_dir, "results")
    end
    isdir(results_dir) || return 0
    length(filter(
        f -> startswith(f, "param_") && endswith(f, ".csv") && occursin("transitions", f),
        readdir(results_dir)
    ))
end

"""Load expected parameter count from the saved config."""
function expected_count(work_dir::String)
    param_file = joinpath(work_dir, "selected_parameters.jld2")
    isfile(param_file) || return 0
    length(load_parameter_config(param_file)["param_ids"])
end

"""
Build a status snapshot across all work directories.
Returns (total_done, total_expected, per_dir_status)
"""
function snapshot(work_dirs::Vector{String})
    total_done     = 0
    total_expected = 0
    per_dir        = Dict{String, NamedTuple}()

    for wd in work_dirs
        done = count_results(wd)
        exp  = expected_count(wd)
        total_done     += done
        total_expected += exp
        per_dir[wd] = (done=done, expected=exp, complete=(done >= exp > 0))
    end

    return total_done, total_expected, per_dir
end

# ============================================================
# DISCORD PAYLOAD BUILDERS
# (send_discord_embed / fmt_duration come from discord_notifier.jl)
# ============================================================

function notify_monitor_start(work_dirs, total_expected, poll_interval)
    fields = Dict{String,Any}[
        Dict{String,Any}("name" => wd, "value" => "$(expected_count(wd)) parameters", "inline" => true)
        for wd in work_dirs
    ]
    send_discord_embed(
        "🚀 Monitor Started",
        "Watching **$(length(work_dirs))** directory/directories. " *
        "Total parameters: **$total_expected**. " *
        "Polling every **$(fmt_duration(poll_interval))**.",
        color=0x57F287,
        fields=fields
    )
end

function notify_status_poll(work_dirs, total_done, total_expected,
                            elapsed_secs, per_dir, poll_num)
    pct    = total_expected > 0 ? round(100 * total_done / total_expected, digits=1) : 0.0
    filled = round(Int, 20 * total_done / max(total_expected, 1))
    bar    = "█"^filled * "░"^(20 - filled)

    fields = Dict{String,Any}[
        Dict{String,Any}(
            "name"   => basename(wd),
            "value"  => "$(per_dir[wd].done)/$(per_dir[wd].expected) " *
                        (per_dir[wd].complete ? "✅" : "⏳"),
            "inline" => true
        ) for wd in work_dirs
    ]
    push!(fields, Dict{String,Any}("name" => "Elapsed", "value" => fmt_duration(elapsed_secs), "inline" => true))

    send_discord_embed(
        "⏱️ Status Update #$poll_num",
        "`[$bar]` **$pct%** ($total_done / $total_expected)",
        color=0xFEE75C,
        fields=fields
    )
end

function notify_collection_starting(work_dirs, total_expected, elapsed_secs)
    send_discord_embed(
        "📦 Collecting Results",
        "All **$total_expected** simulations complete in **$(fmt_duration(elapsed_secs))**. " *
        "Starting collection and plotting for **$(length(work_dirs))** directory/directories...",
        color=0xEB459E
    )
end

function notify_collection_done(work_dirs, elapsed_secs, failed_dirs)
    n_ok = length(work_dirs) - length(failed_dirs)
    desc = "Collection finished in **$(fmt_duration(elapsed_secs))**. " *
           "**$n_ok/$(length(work_dirs))** directories OK."
    !isempty(failed_dirs) && (desc *= "\n\n⚠️ **Failed:** " * join(failed_dirs, ", "))
    send_discord_embed(
        isempty(failed_dirs) ? "✅ Pipeline Complete" : "⚠️ Pipeline Complete (with errors)",
        desc,
        color=isempty(failed_dirs) ? 0x57F287 : 0xED4245
    )
end

function notify_timeout(elapsed_secs, total_done, total_expected)
    send_discord_embed(
        "❌ Monitor Timed Out",
        "Gave up after **$(fmt_duration(elapsed_secs))**. " *
        "Only **$total_done/$total_expected** results present.",
        color=0xED4245
    )
end

# ============================================================
# COLLECT + PLOT FOR ONE WORK DIRECTORY
# ============================================================

"""
Run collection then plotting for a single work directory.
Returns true on success.
"""
function collect_and_plot(work_dir::String, collect_script::String, plot_script::String; 
    DETERMINISTIC::Bool=DETERMINISTIC)
    ok = true

    println("  [collect] $work_dir")
    try
        run(`julia --project=. $collect_script work-dir=$work_dir deterministic=$DETERMINISTIC`)
    catch e
        @error "Collection failed for $work_dir: $e"
        ok = false
    end

    ok || return false
    if RUN_PLOTS
        println("  [plot]    $work_dir")
        try
            run(`julia --project=. $plot_script work-dir=$work_dir`)
        catch e
            @error "Plotting failed for $work_dir: $e"
            ok = false
        end
    else
        println("RUN_PLOTS is set to false. Not generating plots for $work_dir")
    end
    return ok
end
using Printf

function check_stoch_jobs()
    output = read(`squeue --user a.hari`, String)
    lines = split(output, '\n')
    # Skip header line, check only data lines
    has_stoch = any(lines[2:end]) do line
        parts = split(line)
        length(parts) >= 3 && startswith(parts[3], "stoch_")
    end
    return !has_stoch
end

# ============================================================
# MAIN MONITOR LOOP
# ============================================================

# Replace everything from "MAIN MONITOR LOOP" to the end of the file with this:

# ============================================================
# MAIN MONITOR LOOP
# ============================================================

function run_monitor(work_dirs, poll_interval, max_wait, collect_script, plot_script)
    println("=" ^ 70)
    println("MONITOR & COLLECT")
    println("=" ^ 70)
    println("Watching: $(join(work_dirs, ", "))")
    println("Poll interval: $(fmt_duration(poll_interval))")
    println("Max wait: $(fmt_duration(max_wait))")
    println("=" ^ 70)
    flush(stdout)

    start_time = time()
    poll_count = 0

    total_done, total_expected, per_dir = snapshot(work_dirs)
    notify_monitor_start(work_dirs, total_expected, poll_interval)

    all_complete = total_expected > 0 && total_done >= total_expected
    POLL_INTERVAL = copy(poll_interval)
    while !all_complete
        elapsed = time() - start_time

        if max_wait > 0 && elapsed >= max_wait
            notify_timeout(elapsed, total_done, total_expected)
            println("⚠️  Max wait exceeded. Exiting without collection.")
            exit(1)
        end

        sleep_remaining = poll_interval
        while sleep_remaining > 0
            sleep(min(30, sleep_remaining))
            sleep_remaining -= 30
        end

        poll_count += 1
        total_done, total_expected, per_dir = snapshot(work_dirs)
        elapsed = time() - start_time

        println("\n[Poll #$poll_count | $(fmt_duration(elapsed))] $total_done/$total_expected done")
        for wd in work_dirs
            s = per_dir[wd]
            println("  $(rpad(basename(wd), 20)) $(s.done)/$(s.expected) $(s.complete ? "✅" : "⏳")")
        end
        flush(stdout)

        notify_status_poll(work_dirs, total_done, total_expected, elapsed, per_dir, poll_count)

        all_complete = total_expected > 0 && total_done >= total_expected
        if !all_complete
            fraction_done = total_done/total_expected
            time_required = (1-fraction_done)*elapsed/fraction_done
            if time_required - POLL_INTERVAL < 60
                poll_interval = max(time_required, 60)
            end
            # check_stoch_jobs()/squeue is a SLURM-specific safety net (checks
            # for lingering "stoch_"-named jobs in case file-count tracking
            # undercounts) -- squeue doesn't exist off-cluster, so this
            # crashed with ENOENT on macProHome once a run was slow enough to
            # need more than one poll (fast combos dodged it by finishing
            # before the first poll). Only run it where SLURM is actually
            # present; file-count-based all_complete above is authoritative
            # otherwise.
            if has_slurm()
                x = check_stoch_jobs()
                if x
                    output = read(`squeue --user a.hari`, String)
                    println(output)
                    all_complete = copy(x)
                end
            end
        end



    end

    # ============================================================
    # COLLECT & PLOT
    # ============================================================

    elapsed_before_collect = time() - start_time
    notify_collection_starting(work_dirs, total_expected, elapsed_before_collect)
    println("\n" * "=" ^ 70)
    println("ALL SIMULATIONS DONE — COLLECTING RESULTS")
    println("=" ^ 70)
    flush(stdout)

    failed_dirs = String[]
    for wd in work_dirs
        success = collect_and_plot(wd, collect_script, plot_script)
        success || push!(failed_dirs, wd)
    end

    total_elapsed = time() - start_time
    notify_collection_done(work_dirs, total_elapsed, failed_dirs)

    println("\n" * "=" ^ 70)
    isempty(failed_dirs) ? println("✅ PIPELINE COMPLETE") : println("⚠️  PIPELINE COMPLETE (with errors)")
    println("Total time: $(fmt_duration(total_elapsed))")
    println("=" ^ 70)
    flush(stdout)

    isempty(failed_dirs) || exit(1)
end

run_monitor(WORK_DIRS, POLL_INTERVAL, MAX_WAIT, COLLECT_SCRIPT, PLOT_SCRIPT)