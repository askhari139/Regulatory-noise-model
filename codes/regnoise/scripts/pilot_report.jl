# Summaries of a collected run, used to check the pilot (SPEC §13).
# julia --project scripts/pilot_report.jl --collected DIR [--tolcheck DIR] [--out FILE.md]
# DIR holds <NET>_<MODE>_<table>.csv.gz from collect.jl. Prevalence-weighted means use
# the class prevalences stored in the selection (ClassPrevalence column).
using CSV, DataFrames, CodecZlib, Statistics, Printf
include(joinpath(@__DIR__, "cli.jl"))
o = parse_cli(ARGS, Dict{String,Any}("collected" => nothing, "tolcheck" => nothing,
    "selection" => nothing, "out" => nothing))
o["collected"] === nothing && error("--collected is required")

readgz(path) = CSV.read(GzipDecompressorStream(open(path)), DataFrame;
                        types=Dict(:State => String, :From => String, :To => String), validate=false)
files = readdir(o["collected"])
combos = unique([(m[1], m[2]) for f in files
                 for m in (match(r"^(\w+?)_(Additive|Multiplicative|Stationary)_occupancy\.csv\.gz$", f),)
                 if m !== nothing])
out = IOBuffer()
say(args...) = (println(out, args...); println(args...))

allhigh(n) = "H"^n

for (net, mode) in sort(combos)
    occ = readgz(joinpath(o["collected"], "$(net)_$(mode)_occupancy.csv.gz"))
    nn = length(first(occ.State))
    AH = allhigh(nn)
    whole = filter(r -> r.Block == 0 && r.ThresholdMult == 1.0, occ)
    say("\n## $net — $mode   ($(length(unique(occ.ParamID))) parameter sets)\n")

    # 1. C2 convergence
    conv = readgz(joinpath(o["collected"], "$(net)_$(mode)_convergence.csv.gz"))
    nmax = maximum(conv.NDraws)
    ref = filter(r -> r.NDraws == nmax, conv)
    refd = Dict((r.ParamID, r.Sigma, r.DT, r.State) => r.Estimate for r in eachrow(ref))
    say("### C2 convergence: |estimate(n) − estimate($nmax)| over all parameter sets, σ, dt, states\n")
    say("| n draws | median | 95th pct | max |\n|---|---|---|---|")
    for n in sort(unique(conv.NDraws))
        n == nmax && continue
        d = [abs(r.Estimate - refd[(r.ParamID, r.Sigma, r.DT, r.State)]) for r in eachrow(conv) if r.NDraws == n]
        say(@sprintf("| %d | %.4f | %.4f | %.4f |", n, median(d), quantile(d, 0.95), maximum(d)))
    end
    say(@sprintf("\nBootstrap SE at n = %d: median %.4f, max %.4f\n", nmax, median(ref.BootSE), maximum(ref.BootSE)))

    # 2. Decomposition for the all-high state at the control σ values (pool dt)
    say("### All-high ($AH) occupancy, mean over parameter sets (unweighted)\n")
    c3 = filter(r -> r.Condition == "C3" && r.State == AH, whole)
    holds = sort(unique(c3.THold))
    ctrlσ = sort(unique(c3.Sigma))
    pooldt = isempty(c3) ? 0.01 : first(c3.DT)
    say("| σ | C0 | C4 mean-λ | C2 frozen | " * join(["C3 T=$h" for h in holds], " | ") * " | C1 stochastic |")
    say("|" * repeat("---|", 4 + length(holds) + 1))
    m(cond; σ=nothing, th=nothing) = begin
        d = filter(r -> r.Condition == cond && r.State == AH &&
                        (σ === nothing || r.Sigma == σ) && (cond in ("C0",) || r.DT == pooldt || cond == "C0") &&
                        (th === nothing || r.THold == th), whole)
        isempty(d) ? NaN : mean(d.Mean)
    end
    for σ in ctrlσ
        row = [m("C0"), m("C4"; σ), m("C2"; σ), [m("C3"; σ, th=h) for h in holds]..., m("C1"; σ)]
        say(@sprintf("| %g | ", σ) * join([@sprintf("%.3f", v) for v in row], " | ") * " |")
    end

    # 3. C1 vs C2 across all σ at each dt (the old Fig. 4D comparison)
    say("\n### All-high occupancy: C1 (stochastic) vs C2 (frozen, resampled) by σ and dt\n")
    say("| dt | σ | C2 | C1 | C1 − C2 |\n|---|---|---|---|---|")
    for dt in sort(unique(filter(r -> r.Condition == "C1", whole).DT)), σ in sort(unique(filter(r -> r.Condition == "C1", whole).Sigma))
        a = filter(r -> r.Condition == "C1" && r.State == AH && r.Sigma == σ && r.DT == dt, whole)
        b = filter(r -> r.Condition == "C2" && r.State == AH && r.Sigma == σ && r.DT == dt, whole)
        (isempty(a) || isempty(b)) && continue
        say(@sprintf("| %g | %g | %.3f | %.3f | %+.3f |", dt, σ, mean(b.Mean), mean(a.Mean), mean(a.Mean) - mean(b.Mean)))
    end

    # 4. Block convergence of C1 (pool dt)
    say("\n### C1 all-high occupancy by 100-unit block (dt = $pooldt)\n")
    blk = filter(r -> r.Condition == "C1" && r.State == AH && r.Block > 0 && r.ThresholdMult == 1.0 && r.DT == pooldt, occ)
    nb = maximum(blk.Block)
    say("| σ | " * join(["block $b" for b in 1:nb], " | ") * " | max drift |")
    say("|" * repeat("---|", nb + 2))
    for σ in sort(unique(blk.Sigma))
        v = [mean(filter(r -> r.Sigma == σ && r.Block == b, blk).Mean) for b in 1:nb]
        say(@sprintf("| %g | ", σ) * join([@sprintf("%.3f", x) for x in v], " | ") * @sprintf(" | %.3f |", maximum(v) - minimum(v)))
    end
end

# 5. Tolerance check
if o["tolcheck"] !== nothing
    say("\n## Tolerance check (reltol 1e-6 vs 1e-4), same ParamIDs\n")
    say("| network | mode | cells compared | median |Δ| | max |Δ| | frac |Δ| > 2·SE |\n|---|---|---|---|---|---|")
    for (net, mode) in sort(combos)
        f = joinpath(o["tolcheck"], "$(net)_$(mode)_occupancy.csv.gz")
        isfile(f) || continue
        t = filter(r -> r.Block == 0 && r.ThresholdMult == 1.0, readgz(f))
        p = filter(r -> r.Block == 0 && r.ThresholdMult == 1.0 && r.ParamID in Set(t.ParamID),
                   readgz(joinpath(o["collected"], "$(net)_$(mode)_occupancy.csv.gz")))
        j = innerjoin(select(t, :ParamID, :Condition, :Sigma, :DT, :THold, :State, :Mean => :Tight, :SE => :SETight),
                      select(p, :ParamID, :Condition, :Sigma, :DT, :THold, :State, :Mean => :Loose, :SE => :SELoose),
                      on=[:ParamID, :Condition, :Sigma, :DT, :THold, :State])
        d = abs.(j.Tight .- j.Loose)
        se = sqrt.(coalesce.(j.SETight, 0.0) .^ 2 .+ coalesce.(j.SELoose, 0.0) .^ 2)
        frac = mean((d .> 2 .* se) .& (d .> 1e-9))
        say(@sprintf("| %s | %s | %d | %.4f | %.4f | %.3f |", net, mode, nrow(j), median(d), maximum(d), frac))
    end
end

o["out"] === nothing || write(o["out"], take!(out))
