# Observables (SPEC §9). Replicates are accumulated into a CellAccumulator, then turned
# into long-format rows.

"Accumulates per-replicate observables for one condition cell."
mutable struct CellAccumulator
    nstates::Int
    n_nodes::Int
    occ::Vector{Vector{Float64}}             # per replicate: occupancy per state (m = 1, whole window)
    block_occ::Vector{Matrix{Float64}}       # per replicate: nstates × n_blocks
    mult_occ::Vector{Matrix{Float64}}        # per replicate: nstates × n_multipliers
    dwell_sum::Vector{Float64}               # completed visits only
    dwell_n::Vector{Int}
    trans::Matrix{Int}                       # from × to, summed over replicates
    switches::Vector{Int}
    node_mean::Vector{Vector{Float64}}
    node_sd::Vector{Vector{Float64}}
    store_lfc::Bool                          # windowed fold-change kept for this cell (C1, C3)
    lfc::Vector{Matrix{UInt16}}              # per stored replicate: n_nodes × (lags × thresholds)
end

CellAccumulator(n_nodes::Int; store_lfc::Bool=false) = (ns = 1 << n_nodes;
    CellAccumulator(ns, n_nodes, Vector{Float64}[], Matrix{Float64}[], Matrix{Float64}[],
                    zeros(ns), zeros(Int, ns), zeros(Int, ns, ns), Int[],
                    Vector{Float64}[], Vector{Float64}[], store_lfc, Matrix{UInt16}[]))

"Fraction stored as UInt16 (resolution 1.5e-5)."
encode_frac(f) = round(UInt16, clamp(f, 0, 1) * typemax(UInt16))
decode_frac(u) = u / typemax(UInt16)

"""
Windowed log fold-change of each node: for each lag τ and log2 threshold c, the fraction of
time points t in the window with |log2 x(t+τ) − log2 x(t)| ≥ c. Columns ordered lag-major.
"""
function lfc_fractions(cfg::Config, x::AbstractMatrix)
    lx = log2.(max.(x, 1e-12))
    nt = size(lx, 2)
    out = zeros(UInt16, size(lx, 1), length(cfg.lfc_lags) * length(cfg.lfc_thresholds))
    for (a, τ) in enumerate(cfg.lfc_lags)
        L = round(Int, τ / cfg.saveat)
        for i in axes(lx, 1)
            d = abs.(view(lx, i, 1+L:nt) .- view(lx, i, 1:nt-L))
            for (b, c) in enumerate(cfg.lfc_thresholds)
                out[i, (a - 1) * length(cfg.lfc_thresholds) + b] = encode_frac(mean(d .>= c))
            end
        end
    end
    return out
end

nreplicates(a::CellAccumulator) = length(a.occ)

"Block index (1-based) of each saved time in the analysis window; the last block is closed."
function block_index(cfg::Config, t::Vector{Float64})
    len = (cfg.t_end - cfg.t_burn) / cfg.n_blocks
    return [min(cfg.n_blocks, 1 + floor(Int, (ti - cfg.t_burn) / len + 1e-9)) for ti in t]
end

function add_replicate!(a::CellAccumulator, cfg::Config, tr::Trajectory, θ::Vector{Float64})
    nt = size(tr.x, 2)
    codes = [state_code(view(tr.x, :, j), θ) for j in 1:nt]

    occ = zeros(a.nstates)
    for c in codes
        occ[c+1] += 1
    end
    push!(a.occ, occ ./ nt)

    bi = block_index(cfg, tr.t)
    bocc = zeros(a.nstates, cfg.n_blocks)
    bn = zeros(Int, cfg.n_blocks)
    for (c, b) in zip(codes, bi)
        bocc[c+1, b] += 1
        bn[b] += 1
    end
    push!(a.block_occ, bocc ./ bn')

    mocc = zeros(a.nstates, length(cfg.threshold_multipliers))
    for (k, m) in enumerate(cfg.threshold_multipliers), j in 1:nt
        mocc[state_code(view(tr.x, :, j), θ, m)+1, k] += 1
    end
    push!(a.mult_occ, mocc ./ nt)

    # runs of identical states; first and last runs are censored and not counted as dwell
    run_start = 1
    nsw = 0
    for j in 2:nt
        if codes[j] != codes[j-1]
            a.trans[codes[j-1]+1, codes[j]+1] += 1
            nsw += 1
            if run_start > 1
                a.dwell_sum[codes[j-1]+1] += (j - run_start) * cfg.saveat
                a.dwell_n[codes[j-1]+1] += 1
            end
            run_start = j
        end
    end
    push!(a.switches, nsw)
    push!(a.node_mean, vec(mean(tr.x; dims=2)))
    push!(a.node_sd, vec(std(tr.x; dims=2)))
    a.store_lfc && length(a.lfc) < cfg.replicates_stored && push!(a.lfc, lfc_fractions(cfg, tr.x))
    return a
end

# ---------------------------------------------------------------------------
# Rows. Every table shares the key columns below.
# ---------------------------------------------------------------------------

cellkey(cond, sigma, dt, thold, nrep) =
    (Condition=cond, Sigma=Float64(sigma), DT=Float64(dt), THold=Float64(thold), NReplicates=nrep)

se(x) = length(x) > 1 ? std(x) / sqrt(length(x)) : NaN

function occupancy_rows!(rows, key, a::CellAccumulator, cfg::Config)
    n = a.n_nodes
    O = reduce(hcat, a.occ)                       # nstates × nrep
    for s in 1:a.nstates
        v = O[s, :]
        push!(rows, (; key..., State=state_label(s - 1, n), ThresholdMult=1.0, Block=0,
                     Mean=mean(v), SD=std(v), SE=se(v)))
        for b in 1:cfg.n_blocks
            vb = [B[s, b] for B in a.block_occ]
            push!(rows, (; key..., State=state_label(s - 1, n), ThresholdMult=1.0, Block=b,
                         Mean=mean(vb), SD=std(vb), SE=se(vb)))
        end
        for (k, m) in enumerate(cfg.threshold_multipliers)
            vm = [M[s, k] for M in a.mult_occ]
            push!(rows, (; key..., State=state_label(s - 1, n), ThresholdMult=m, Block=0,
                         Mean=mean(vm), SD=std(vm), SE=se(vm)))
        end
    end
end

function transition_rows!(rows, key, a::CellAccumulator)
    for i in 1:a.nstates, j in 1:a.nstates
        a.trans[i, j] > 0 &&
            push!(rows, (; key..., From=state_label(i - 1, a.n_nodes),
                         To=state_label(j - 1, a.n_nodes), Count=a.trans[i, j]))
    end
end

function dwell_rows!(rows, key, a::CellAccumulator)
    for s in 1:a.nstates
        push!(rows, (; key..., State=state_label(s - 1, a.n_nodes),
                     MeanDwell=a.dwell_n[s] > 0 ? a.dwell_sum[s] / a.dwell_n[s] : NaN,
                     NVisits=a.dwell_n[s]))
    end
end

function cell_rows!(rows, key, a::CellAccumulator)
    push!(rows, (; key..., MeanSwitches=mean(a.switches), SDSwitches=std(a.switches)))
end

"""
Per-replicate time-averaged mean and SD of every node (columns Mean_<node>, SD_<node>), one row
per replicate, for the first `replicates_stored` replicates.
"""
function replicate_rows!(rows, key, a::CellAccumulator, net::Network, cfg::Config)
    names = Tuple(Symbol(s, "_", n) for n in net.nodes for s in ("Mean", "SD"))
    for r in 1:min(nreplicates(a), cfg.replicates_stored)
        vals = Tuple(Float32(s == 1 ? a.node_mean[r][i] : a.node_sd[r][i]) for i in eachindex(net.nodes) for s in 1:2)
        push!(rows, (; key..., Replicate=Int32(r), NamedTuple{names}(vals)...))
    end
end

lfc_colnames(cfg::Config) = [Symbol("Frac_lag", L == round(L) ? Int(L) : L, "_", Int(2^c), "fold")
                             for L in cfg.lfc_lags for c in cfg.lfc_thresholds]

"""
Per-replicate windowed fold-change fractions, one row per replicate, columns <node>_Frac_lag<τ>_<k>fold
(stored as UInt16; divide by 65535).
"""
function lfc_rows!(rows, key, a::CellAccumulator, net::Network, cfg::Config)
    a.store_lfc || return
    cols = lfc_colnames(cfg)
    names = Tuple(Symbol(n, "_", c) for n in net.nodes for c in cols)
    for (r, M) in enumerate(a.lfc)
        push!(rows, (; key..., Replicate=Int32(r), NamedTuple{names}(Tuple(vec(permutedims(M))))...))
    end
end

function expression_rows!(rows, key, a::CellAccumulator, net::Network)
    M = reduce(hcat, a.node_mean); S = reduce(hcat, a.node_sd)
    for (i, node) in enumerate(net.nodes)
        push!(rows, (; key..., Node=node, Mean=mean(M[i, :]), MeanSD=mean(S[i, :])))
    end
end

# ---------------------------------------------------------------------------
# λ pool (SPEC §7, §9.4)
# ---------------------------------------------------------------------------

"Autocorrelation of x at `lag` (NaN when x is constant or too short)."
function acf(x::AbstractVector, lag::Int)
    n = length(x)
    lag < n || return NaN
    μ = mean(x); d = x .- μ
    v = sum(abs2, d)
    v > 0 || return NaN
    return sum(d[1:n-lag] .* d[1+lag:n]) / v
end

function lambda_rows!(rows, key, trajs::Vector{Trajectory}, net::Network, lo, hi, cfg::Config)
    pool = reduce(hcat, (tr.lam for tr in trajs))
    for (k, e) in enumerate(net.edges)
        v = pool[k, :]
        q = quantile(v, cfg.lambda_quantiles)
        acfs = [mean_or_nan(filter(!isnan, [acf(tr.lam[k, :], L) for tr in trajs]))
                for L in cfg.acf_lags]
        qcols = NamedTuple{Tuple(Symbol("Q", lpad(round(Int, 100p), 2, '0')) for p in cfg.lambda_quantiles)}(q)
        acols = NamedTuple{Tuple(Symbol("ACF_lag", L) for L in cfg.acf_lags)}(acfs)
        push!(rows, (; key..., Edge=e.name, EdgeType=e.act ? "act" : "inh",
                     Mean=mean(v), SD=std(v), Min=minimum(v), Max=maximum(v),
                     FracAtMin=mean(v .== lo[k]), FracAtMax=mean(v .== hi[k]),
                     qcols..., acols...))
    end
    return pool
end

mean_or_nan(x) = isempty(x) ? NaN : mean(x)

# ---------------------------------------------------------------------------
# C2 convergence (SPEC §9.5)
# ---------------------------------------------------------------------------

function convergence_rows!(rows, key, a::CellAccumulator, cfg::Config, rng)
    O = reduce(hcat, a.occ)                       # nstates × ndraws (in draw order)
    nd = size(O, 2)
    nd == cfg.n_det_draws || error("convergence: have $nd draws, expected $(cfg.n_det_draws)")
    boot = [mean(O[:, rand(rng, 1:nd, nd)]; dims=2) for _ in 1:cfg.n_bootstrap]
    B = reduce(hcat, boot)
    for s in 1:a.nstates
        for n in cfg.convergence_ns
            push!(rows, (; key..., State=state_label(s - 1, a.n_nodes), NDraws=n,
                         Estimate=mean(O[s, 1:n]),
                         BootSE=n == nd ? std(B[s, :]) : NaN))
        end
    end
end
