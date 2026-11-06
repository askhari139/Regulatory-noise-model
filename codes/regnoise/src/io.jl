# Provenance, raw unit files, validation (SPEC §10, §11).

const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))

function git_info()
    try
        commit = readchomp(`git -C $PKG_ROOT rev-parse HEAD`)
        dirty = !isempty(readchomp(`git -C $PKG_ROOT status --porcelain -- .`))
        return commit, dirty
    catch
        return "unknown", true
    end
end

function provenance()
    commit, dirty = git_info()
    mf = joinpath(PKG_ROOT, "Manifest.toml")
    return Dict("git_commit" => commit, "git_dirty" => dirty,
                "julia_version" => string(VERSION),
                "manifest_sha256" => isfile(mf) ? open(io -> bytes2hex(sha256(io)), mf) : "none",
                "hostname" => gethostname(),
                "timestamp" => string(now()))
end

"Refuse to produce results from uncommitted code unless explicitly allowed."
function check_clean_tree(allow_dirty::Bool)
    commit, dirty = git_info()
    if dirty && !allow_dirty
        error("code tree at $PKG_ROOT has uncommitted changes (commit $commit). Commit first, or pass --allow-dirty for tests.")
    end
end

raw_path(cfg::Config, net, mode, pid) = joinpath(run_dir(cfg), "raw", net, mode, "param_$(pid).jld2")

"Write one unit atomically (tmp file + rename), so an existing file is always complete."
function write_unit(cfg::Config, net::Network, mode, pid, class, θ, rows, wall)
    path = raw_path(cfg, net.name, mode, pid)
    mkpath(dirname(path))
    meta = merge(provenance(), Dict(
        "config_hash" => cfg.hash, "run_id" => cfg.run_id, "network" => net.name,
        "noise_mode" => String(mode), "param_id" => pid, "stability_class" => String(class),
        "nodes" => net.nodes, "edges" => [e.name for e in net.edges],
        "thresholds" => θ, "wall_seconds" => wall))
    tables = Dict(String(t) => (isempty(rows[t]) ? DataFrame() : pool_strings!(DataFrame(rows[t]))) for t in keys(rows))
    tmp = path * ".tmp"
    jldsave(tmp; meta=meta, tables=tables)
    mv(tmp, path; force=true)
    return path
end

read_unit(path) = jldopen(f -> (f["meta"], f["tables"]), path, "r")

"Store repeated string columns (Condition, State, Edge, ...) as pooled integer codes."
function pool_strings!(df::DataFrame)
    for c in names(df)
        eltype(df[!, c]) <: AbstractString && (df[!, c] = PooledArray(String.(df[!, c]); compress=true))
    end
    return df
end

"Status of an existing unit file: :missing, :done, or error on a config mismatch."
function unit_status(cfg::Config, net, mode, pid)
    path = raw_path(cfg, net, mode, pid)
    isfile(path) || return :missing
    meta = jldopen(f -> f["meta"], path, "r")
    meta["config_hash"] == cfg.hash ||
        error("$path was produced with config $(meta["config_hash"]), current config is $(cfg.hash). " *
              "Use a new run_id or delete the old results explicitly.")
    return :done
end

# ---------------------------------------------------------------------------
# Validation of one unit (SPEC §11)
# ---------------------------------------------------------------------------

function validate_unit(cfg::Config, path, mode, lo, hi)
    problems = String[]
    meta, T = read_unit(path)
    meta["config_hash"] == cfg.hash || push!(problems, "config hash mismatch")
    occ = T["occupancy"]

    have = Set(zip(occ.Condition, occ.Sigma, occ.DT, occ.THold, occ.NReplicates))
    for c in expected_cells(cfg, mode)
        c in have || push!(problems, "missing or wrong-size cell $(c)")
    end
    length(have) == length(expected_cells(cfg, mode)) ||
        push!(problems, "unexpected cells: $(length(have)) present, $(length(expected_cells(cfg, mode))) expected")

    must_be_finite = Dict("occupancy" => ["Mean"], "cells" => ["MeanSwitches"],
                          "expression" => ["Mean"], "lambda" => ["Mean", "Min", "Max"],
                          "convergence" => ["Estimate"])
    for (t, cols) in must_be_finite
        haskey(T, t) || (push!(problems, "missing table $t"); continue)
        nrow(T[t]) == 0 && continue
        for c in cols
            any(isnan, T[t][!, c]) && push!(problems, "NaN in $t.$c")
        end
    end
    any(x -> x < -1e-12 || x > 1 + 1e-12, occ.Mean) && push!(problems, "occupancy outside [0,1]")

    g = combine(groupby(occ, [:Condition, :Sigma, :DT, :THold, :ThresholdMult, :Block]), :Mean => sum => :Total)
    bad = filter(r -> abs(r.Total - 1) > 1e-9, g)
    nrow(bad) > 0 && push!(problems, "occupancies do not sum to 1 in $(nrow(bad)) groups")

    # per-replicate tables: first `replicates_stored` replicates × nodes per cell;
    # the windowed fold-change only for the time-varying conditions C1 and C3
    for (tname, conds) in (("replicates", ("C0", "C1", "C2", "C3", "C4", "C5")), ("lfc", ("C1", "C3")))
        haskey(T, tname) || (push!(problems, "missing table $tname"); continue)
        t = T[tname]
        cnt = nrow(t) == 0 ? Dict() :
            Dict((r.Condition, r.Sigma, r.DT, r.THold) => r.n for r in eachrow(combine(groupby(t, [:Condition, :Sigma, :DT, :THold]), nrow => :n)))
        for c in expected_cells(cfg, mode)
            c[1] in conds || continue
            want = min(c[5], cfg.replicates_stored)
            get(cnt, c[1:4], 0) == want || push!(problems, "$tname: cell $(c[1:4]) has $(get(cnt, c[1:4], 0)) rows, expected $want")
        end
    end

    n5 = count(c -> c[1] == "C5", expected_cells(cfg, mode))
    if n5 > 0
        nb = haskey(T, "basins") ? nrow(T["basins"]) : 0
        nb == n5 * cfg.n_c5_draws || push!(problems, "basins: $nb rows, expected $(n5 * cfg.n_c5_draws)")
    end

    lam = T["lambda"]
    for r in eachrow(lam)
        k = findfirst(==(r.Edge), meta["edges"])
        (r.Min >= lo[k] - 1e-12 && r.Max <= hi[k] + 1e-12) ||
            push!(problems, "λ outside bounds on edge $(r.Edge) at σ=$(r.Sigma), dt=$(r.DT)")
    end
    return problems
end
