# Configuration loading. The TOML file is the single source of truth (SPEC §0.2).
# Every key is required and unknown keys are rejected, so a typo can never fall
# back to a silent default.

const SCHEMA = Dict(
    "run"        => Dict("run_id" => String, "seed_root" => Integer,
                         "networks" => Vector, "noise_modes" => Vector),
    "paths"      => Dict("data_dir" => String, "results_root" => String),
    "selection"  => Dict("n_per_class" => Integer, "classes" => Vector),
    "model"      => Dict("lambda_inh_bounds" => Vector, "lambda_act_bounds" => Vector),
    "noise"      => Dict("sigmas" => Vector, "dts" => Vector),
    "sim"        => Dict("t_end" => Real, "t_burn" => Real, "saveat" => Real,
                         "reltol" => Real, "abstol" => Real, "n_ics" => Integer,
                         "n_traj" => Integer, "n_blocks" => Integer,
                         "threshold_multipliers" => Vector),
    "controls"   => Dict("modes" => Vector, "n_det_draws" => Integer,
                         "convergence_ns" => Vector, "n_bootstrap" => Integer,
                         "ctrl_sigmas" => Vector, "hold_times" => Vector,
                         "pool_dt" => Real, "acf_lags" => Vector,
                         "lambda_quantiles" => Vector),
    "inputs"     => Dict("sha256" => Dict),
)

const NOISE_MODES = ("Additive", "Multiplicative", "Stationary")
const STABILITY_CLASSES = ("monostable", "bistable", "tristable", "tetrastable")

struct Config
    raw::Dict{String,Any}
    hash::String
    path::String
    # resolved, typed fields
    run_id::String
    seed_root::Int
    networks::Vector{String}
    noise_modes::Vector{String}
    data_dir::String
    results_root::String
    n_per_class::Int
    classes::Vector{String}
    inh_bounds::Tuple{Float64,Float64}
    act_bounds::Tuple{Float64,Float64}
    sigmas::Vector{Float64}
    dts::Vector{Float64}
    t_end::Float64
    t_burn::Float64
    saveat::Float64
    reltol::Float64
    abstol::Float64
    n_ics::Int
    n_traj::Int
    n_blocks::Int
    threshold_multipliers::Vector{Float64}
    ctrl_modes::Vector{String}
    n_det_draws::Int
    convergence_ns::Vector{Int}
    n_bootstrap::Int
    ctrl_sigmas::Vector{Float64}
    hold_times::Vector{Float64}
    pool_dt::Float64
    acf_lags::Vector{Int}
    lambda_quantiles::Vector{Float64}
    input_sha256::Dict{String,String}
end

function check_schema(raw::Dict)
    for (sec, keys_) in SCHEMA
        haskey(raw, sec) || error("config: missing section [$sec]")
        for (k, T) in keys_
            haskey(raw[sec], k) || error("config: missing key [$sec].$k")
            raw[sec][k] isa T || error("config: [$sec].$k must be $T, got $(typeof(raw[sec][k]))")
        end
        extra = setdiff(keys(raw[sec]), keys(keys_))
        isempty(extra) || error("config: unknown key(s) in [$sec]: $(join(extra, ", "))")
    end
    extra = setdiff(keys(raw), keys(SCHEMA))
    isempty(extra) || error("config: unknown section(s): $(join(extra, ", "))")
end

"Hash of the config with [paths] removed, so moving data does not change the hash."
function config_hash(raw::Dict)
    hashed = Dict(k => v for (k, v) in raw if k != "paths")
    io = IOBuffer()
    TOML.print(io, hashed; sorted=true)
    return bytes2hex(sha256(take!(io)))
end

function load_config(path::AbstractString; data_dir=nothing, results_root=nothing)
    raw = TOML.parsefile(path)
    check_schema(raw)
    r, s, c = raw["run"], raw["sim"], raw["controls"]

    modes = String.(r["noise_modes"])
    for m in vcat(modes, String.(c["modes"]))
        m in NOISE_MODES || error("config: unknown noise mode '$m' (allowed: $(join(NOISE_MODES, ", ")))")
    end
    classes = String.(raw["selection"]["classes"])
    for cl in classes
        cl in STABILITY_CLASSES || error("config: unknown stability class '$cl'")
    end
    sig = Float64.(raw["noise"]["sigmas"])
    issorted(sig) && sig[1] == 0.0 || error("config: sigmas must be sorted and start at 0.0")
    for cs in c["ctrl_sigmas"]
        cs in sig || error("config: ctrl_sigma $cs not in sigmas")
    end
    c["pool_dt"] in raw["noise"]["dts"] || error("config: pool_dt must be one of noise.dts")
    maximum(c["convergence_ns"]) == c["n_det_draws"] ||
        error("config: the largest convergence_n must equal n_det_draws")
    0 < s["t_burn"] < s["t_end"] || error("config: need 0 < t_burn < t_end")
    s["n_traj"] <= s["n_ics"] || error("config: n_traj must be <= n_ics (trajectory r uses IC r)")
    window = s["t_end"] - s["t_burn"]
    isapprox(rem(window, s["n_blocks"] * s["saveat"]), 0; atol=1e-9) ||
        error("config: analysis window must split into n_blocks whole multiples of saveat")

    inh = Tuple(Float64.(raw["model"]["lambda_inh_bounds"]))
    act = Tuple(Float64.(raw["model"]["lambda_act_bounds"]))
    (0 < inh[1] < inh[2] <= 1) || error("config: bad lambda_inh_bounds")
    (1 <= act[1] < act[2]) || error("config: bad lambda_act_bounds")

    p = raw["paths"]
    base = dirname(abspath(path))
    resolve(x) = isabspath(x) ? x : normpath(joinpath(base, "..", x))
    ddir = data_dir === nothing ? resolve(p["data_dir"]) : abspath(data_dir)
    rroot = results_root === nothing ? resolve(p["results_root"]) : abspath(results_root)

    return Config(raw, config_hash(raw), abspath(path),
        r["run_id"], r["seed_root"], String.(r["networks"]), modes,
        ddir, rroot,
        raw["selection"]["n_per_class"], classes, inh, act,
        sig, Float64.(raw["noise"]["dts"]),
        s["t_end"], s["t_burn"], s["saveat"], s["reltol"], s["abstol"],
        s["n_ics"], s["n_traj"], s["n_blocks"], Float64.(s["threshold_multipliers"]),
        String.(c["modes"]), c["n_det_draws"], Int.(c["convergence_ns"]), c["n_bootstrap"],
        Float64.(c["ctrl_sigmas"]), Float64.(c["hold_times"]), c["pool_dt"],
        Int.(c["acf_lags"]), Float64.(c["lambda_quantiles"]),
        Dict{String,String}(k => String(v) for (k, v) in raw["inputs"]["sha256"]))
end

"Verify the RACIPE input files for `net` against the checksums in the config."
function verify_inputs(cfg::Config, net::AbstractString)
    for f in input_files(net)
        haskey(cfg.input_sha256, f) || error("config: no sha256 recorded for $f")
        path = joinpath(cfg.data_dir, f)
        isfile(path) || error("input file missing: $path")
        got = open(io -> bytes2hex(sha256(io)), path)
        got == cfg.input_sha256[f] || error("checksum mismatch for $path (expected $(cfg.input_sha256[f]), got $got)")
    end
end

input_files(net) = [net * ".prs", net * ".topo", net * "_parameters.dat", net * "_solution.dat"]

run_dir(cfg::Config) = joinpath(cfg.results_root, cfg.run_id)
