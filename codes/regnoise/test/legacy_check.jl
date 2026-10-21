# Regression check against the old pipeline (SPEC §12.6).
# Runs C0 + C1 with the old λ bounds for ParamIDs that the old final/ run also used, and
# compares occupancy with the old MRT. Differences should look like Monte-Carlo noise.
#
#   julia --project test/legacy_check.jl --data-dir D --old-results R --network TS --mode Additive --n 5
# where R is the old results root (e.g. .../RACIPEdata/final).
using RegNoise, CSV, DataFrames, Statistics, Printf
include(joinpath(@__DIR__, "..", "scripts", "cli.jl"))
o = parse_cli(ARGS, Dict{String,Any}("config" => joinpath(PKG, "config", "legacy_check.toml"),
    "data-dir" => nothing, "results-root" => mktempdir(), "old-results" => nothing,
    "network" => "TS", "mode" => "Additive", "n" => "5"))
o["old-results"] === nothing && error("--old-results is required")
cfg = cfg_from(o)
netname, mode = o["network"], o["mode"]
select_parameters(cfg, netname)
net = load_network(cfg, netname)
sel, θ, _ = load_selection(cfg, netname)
params = read_parameters(cfg, net)

old = CSV.read(joinpath(o["old-results"], mode, netname, "results", "all_parameters_results.csv"), DataFrame)
old = filter(r -> r.DT == 0.01, old)
old.State = [replace(s, r"[^01]" => "") for s in string.(old.State)]
shared = intersect(Set(sel.ParamID), Set(old.ParamID))
ids = Int[]
for cls in unique(sel.StabilityClass)       # balance the classes
    pool = sort([p for p in shared if sel.StabilityClass[findfirst(==(p), sel.ParamID)] == cls])
    append!(ids, first(pool, cld(parse(Int, o["n"]), 2)))
end

rows = DataFrame()
for pid in ids
    t = @elapsed u = run_unit(cfg, net, mode, pid, params[pid], θ)
    occ = filter(r -> r.Block == 0 && r.ThresholdMult == 1.0 && r.Condition in ("C0", "C1"), DataFrame(u[:occupancy]))
    for r in eachrow(occ)
        σ = r.Condition == "C0" ? 0.0 : r.Sigma
        m = filter(x -> x.ParamID == pid && x.NoiseLevel == σ && x.State == r.State, old)
        oldv = nrow(m) == 0 ? 0.0 : m.MRT[1]     # old files omitted never-visited states
        push!(rows, (ParamID=pid, Sigma=σ, State=r.State, New=r.Mean, Old=oldv, SE=r.SE))
    end
    println("ParamID $pid done in $(round(t, digits=1)) s")
end
rows.Diff = rows.New .- rows.Old
rows.Z = [s > 0 ? d / (sqrt(2) * s) : (abs(d) < 1e-9 ? 0.0 : Inf) for (d, s) in zip(rows.Diff, rows.SE)]
println("\nper-σ summary (|z| uses √2·SE_new as the SE of the difference):")
for g in groupby(rows, :Sigma)
    @printf("σ=%-6g  max|diff|=%.3f  mean|diff|=%.4f  frac|z|>3 = %.2f\n", g.Sigma[1],
            maximum(abs.(g.Diff)), mean(abs.(g.Diff)), mean(abs.(g.Z) .> 3))
end
CSV.write(joinpath(cfg.results_root, "legacy_check_$(netname)_$(mode).csv"), rows)
println("details: ", joinpath(cfg.results_root, "legacy_check_$(netname)_$(mode).csv"))
