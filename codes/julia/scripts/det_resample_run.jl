"""
Driver script: reads RACIPE data for a network, submits SLURM array jobs,
monitors completion, then collects results.

Usage:
    julia scripts/det_resample_run.jl network=TS sigma-values=0.0,0.001,0.005,0.01,0.05,0.1 num-sims=100 tspan=0.0,1000.0 work-dir=det_resample/TS mem=4G time=1:00:00 resample=false poll-interval=300 max-wait=3600
"""

using DataFrames
using CSV

sourceDir = get(ENV, "RACIPE_SOURCE", pwd()) * "/"
include(sourceDir * "scripts/stochastic_setup.jl")
include(sourceDir * "scripts/parse_args.jl")
include(sourceDir * "scripts/discord_notifier.jl")
include(sourceDir * "src/simulate_racipe_resampled.jl")

# ── arguments ─────────────────────────────────────────────────────────────────

args = parse_named_args(defaults=Dict(
    "network"       => "TS",
    "sigma-values"  => "0.0,0.001,0.005,0.01,0.05,0.1",
    "num-sims"      => "100",
    "tspan"         => "0.0,500.0",
    "work-dir"      => "det_resample",
    "mem"           => "4G",
    "time"          => "1:00:00",
    "cpus"          => "1",
    "poll-interval" => "900",
    "max-wait"      => "86400"
))

NETWORK       = get_arg(args, "network",       "TS")
SIGMA_VALUES  = get_arg(args, "sigma-values",  [0.0, 0.001, 0.005, 0.01, 0.05, 0.1], type=Vector{Float64})
NUM_SIMS      = get_arg(args, "num-sims",      100,   type=Int)
tspan_vec     = get_arg(args, "tspan",         [0.0, 500.0], type=Vector{Float64})
TSPAN         = (tspan_vec[1], tspan_vec[2])
WORK_DIR      = get_arg(args, "work-dir",      joinpath("det_resample", NETWORK))
MEM           = get_arg(args, "mem",           "4G")
TIME_LIMIT    = get_arg(args, "time",          "1:00:00")
CPUS          = get_arg(args, "cpus",          1,     type=Int)
POLL_INTERVAL = get_arg(args, "poll-interval", 900,   type=Int)
MAX_WAIT      = get_arg(args, "max-wait",      86400, type=Int)
RESAMPLE       = get_arg(args, "resample",      true,  type=Bool)

println("=" ^ 70)
println("DETERMINISTIC RESAMPLE — SUBMISSION")
println("=" ^ 70)
println("Network      : $NETWORK")
println("Sigma values : $SIGMA_VALUES")
println("Sims per set : $NUM_SIMS")
println("tspan        : $TSPAN")
println("Work dir     : $WORK_DIR")
println("=" ^ 70)

# ── load RACIPE data ──────────────────────────────────────────────────────────

cd("data")
PRS_FILE        = NETWORK * ".prs"
PARAMETERS_FILE = NETWORK * "_parameters.dat"
TOPO_FILE       = NETWORK * ".topo"
ODE_FILE        = sourceDir * "src/" * NETWORK * "_ode.jl"

prs         = read_prs(PRS_FILE)
params_data = read_parameters(PARAMETERS_FILE, prs)
generate_ode_function(TOPO_FILE, PRS_FILE, ODE_FILE)
include(ODE_FILE)
cd("..")

param_ids       = params_data.data.ParamID
total_expected  = length(param_ids)
println("Parameter sets loaded: $total_expected")

# ── submit ────────────────────────────────────────────────────────────────────

submit_resampled_slurm(
    NETWORK, param_ids, SIGMA_VALUES, WORK_DIR;
    num_sims   = NUM_SIMS,
    tspan      = TSPAN,
    mem        = MEM,
    time_limit = TIME_LIMIT,
    cpus       = CPUS,
    resample    = RESAMPLE
)

# ── helpers ───────────────────────────────────────────────────────────────────

fmt_duration(s) = s < 60 ? "$(round(Int, s))s" :
                  s < 3600 ? "$(round(s/60, digits=1))m" :
                  "$(round(s/3600, digits=2))h"

function count_done(work_dir, param_ids)
    results_dir = joinpath(work_dir, "results")
    return count(id -> isfile(joinpath(results_dir, "param_$(id).csv")), param_ids)
end

# ── monitor ───────────────────────────────────────────────────────────────────

function monitor_jobs(param_ids, work_dir; poll_interval=900, max_wait=86400)
    start_time = time()
    poll_num   = 0

    send_discord_embed(
        "🚀 Monitor Started",
        "Total parameters: **$total_expected**. " *
        "Polling every **$(fmt_duration(poll_interval))**.",
        color=0x57F287
    )

    while true
        elapsed    = time() - start_time
        total_done = count_done(work_dir, param_ids)
        pct        = round(Int, 100 * total_done / total_expected)
        bar        = repeat("█", pct ÷ 5) * repeat("░", 20 - pct ÷ 5)
        poll_num  += 1

        if total_done >= total_expected
            send_discord_embed(
                "✅ All Simulations Complete",
                "All **$total_expected** results present in **$(fmt_duration(elapsed))**.",
                color=0x57F287
            )
            println("✓ All jobs complete.")
            break
        end

        if elapsed >= max_wait
            send_discord_embed(
                "❌ Monitor Timed Out",
                "Gave up after **$(fmt_duration(elapsed))**. " *
                "Only **$total_done/$total_expected** results present.",
                color=0xED4245
            )
            println("⚠️  Max wait exceeded. Exiting.")
            exit(1)
        end

        send_discord_embed(
            "⏱️ Status Update #$poll_num",
            "`[$bar]` **$pct%** ($total_done / $total_expected) — elapsed $(fmt_duration(elapsed))",
            color=0xFEE75C
        )
        println("[$bar] $pct% ($total_done/$total_expected) — sleeping $(fmt_duration(poll_interval))")
        sleep(poll_interval)
    end
end

monitor_jobs(param_ids, WORK_DIR; poll_interval=POLL_INTERVAL, max_wait=MAX_WAIT)

# ── collect ───────────────────────────────────────────────────────────────────

function collect_results(work_dir, param_ids)
    results_dir = joinpath(work_dir, "results")
    all_results = DataFrame()

    for id in param_ids
        file = joinpath(results_dir, "param_$(id).csv")
        if isfile(file)
            append!(all_results, CSV.read(file, DataFrame), promote=true)
        else
            @warn "Missing result file: $file"
        end
    end

    CSV.write(joinpath(work_dir, "all_results.csv"), all_results)
    println("✓ Saved: $(joinpath(work_dir, "all_results.csv"))  ($(nrow(all_results)) rows)")

    # Sigma=0 subset in RACIPE-compatible tab-delimited format
    noise0 = filter(row -> row.NoiseLevel == 0.0, all_results)
    CSV.write(joinpath(work_dir, "noise0_solutions.dat"), noise0;
              delim='\t', writeheader=false)
    println("✓ Saved: $(joinpath(work_dir, "noise0_solutions.dat"))  ($(nrow(noise0)) rows)")
end

collect_results(WORK_DIR, param_ids)