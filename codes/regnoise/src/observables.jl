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
end

CellAccumulator(n_nodes::Int) = (ns = 1 << n_nodes;
    CellAccumulator(ns, n_nodes, Vector{Float64}[], Matrix{Float64}[], Matrix{Float64}[],
                    zeros(ns), zeros(Int, ns), zeros(Int, ns, ns), Int[],
                    Vector{Float64}[], Vector{Float64}[]))

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
