# Seeds, λ update rules and single-trajectory simulation (SPEC §3, §8, §10).

"""
Derive a reproducible 64-bit seed from labelled parts. Uses SHA-256 (stable across Julia
versions, unlike `hash`). Floats are written with `repr` so 0.1 and 0.10000000000000002
never collide.
"""
function derive_seed(parts...)
    s = join((p isa AbstractFloat ? repr(Float64(p)) : string(p) for p in parts), "|")
    return reinterpret(UInt64, sha256(s)[1:8])[1]
end

rng_for(parts...) = Xoshiro(derive_seed(parts...))

# ---------------------------------------------------------------------------
# λ updaters: callables used as PeriodicCallback affects
# ---------------------------------------------------------------------------

"Noise process of SPEC §3 (C1)."
struct ProcessUpdater
    mode::Symbol              # :Additive, :Multiplicative, :Stationary, :LogAdditive, :Geometric
    sigma::Float64
    lam_idx::Vector{Int}
    lo::Vector{Float64}
    hi::Vector{Float64}
    sigma_eff::Vector{Float64}   # Additive/Stationary step scale per edge
    target::Vector{Float64}      # Multiplicative target λ* per edge
    lambda0::Vector{Float64}     # initial (RACIPE) λ, for Stationary
    rng::Xoshiro
end

function ProcessUpdater(mode::AbstractString, sigma, net::Network, lo, hi, p0, rng)
    acts = [e.act for e in net.edges]
    m = Symbol(mode)
    m in (:Additive, :Multiplicative, :Stationary, :LogAdditive, :Geometric) || error("unknown noise mode $mode")
    # Additive/Stationary: step in λ units (σ for inhibition, σ·λ_max for activation).
    # LogAdditive: step in log10 λ, scaled by the edge's log10 range (2 decades for both types).
    seff = m === :LogAdditive ?
        [sigma * (log10(h) - log10(l)) for (l, h) in zip(lo, hi)] :
        [a ? sigma * h : sigma for (a, h) in zip(acts, hi)]
    return ProcessUpdater(m, sigma, net.lam_idx, lo, hi, seff,
        [a ? h : l for (a, l, h) in zip(acts, lo, hi)],
        p0[net.lam_idx], rng)
end

"One update of a single λ (exposed for unit tests)."
@inline function update_lambda(mode::Symbol, λ, λ0, sigma, sigma_eff, target, lo, hi, ξ)
    new = if mode === :Additive
        λ + sigma_eff * ξ
    elseif mode === :Multiplicative
        λ + abs(target - λ) * sigma * ξ
    elseif mode === :LogAdditive
        λ * exp10(sigma_eff * ξ)          # log10 λ ← log10 λ + σ_eff ξ; clamping λ = clamping log λ
    elseif mode === :Geometric
        # The fold-change term that full regulation removes (λ for inhibition, 1/λ for
        # activation) is multiplied by (1 + σξ). Its target 0 lies outside the bounds, so λ is
        # clamped at the full-strength bound (not reflected); the bound is not absorbing, since a
        # later factor > 1 moves λ back. log λ drifts by -σ²/2 per update, so λ piles up near
        # that bound (stationary density ∝ 1/μ², independent of σ; SPEC §3). A factor ≤ 0 sends
        # the edge to full strength. `target` is λ_max for activation edges (> 1), λ_min otherwise.
        f = 1 + sigma * ξ
        target > 1 ? (f <= 0 ? hi : λ / f) : λ * f
    else # :Stationary
        λ0 + sigma_eff * ξ
    end
    return clamp(new, lo, hi)
end

function (u::ProcessUpdater)(integrator)
    p = integrator.p
    @inbounds for k in eachindex(u.lam_idx)
        i = u.lam_idx[k]
        p[i] = update_lambda(u.mode, p[i], u.lambda0[k], u.sigma, u.sigma_eff[k],
                             u.target[k], u.lo[k], u.hi[k], randn(u.rng))
    end
    u_modified!(integrator, false)
end

"Hold-resample control (C3): every T_hold, λ ← a whole vector drawn from the pool."
struct HoldUpdater
    lam_idx::Vector{Int}
    pool::Matrix{Float64}     # n_edges × n_pool
    rng::Xoshiro
end

function (u::HoldUpdater)(integrator)
    j = rand(u.rng, 1:size(u.pool, 2))
    @inbounds for k in eachindex(u.lam_idx)
        integrator.p[u.lam_idx[k]] = u.pool[k, j]
    end
    u_modified!(integrator, false)
end

# ---------------------------------------------------------------------------
# One trajectory
# ---------------------------------------------------------------------------

struct Trajectory
    t::Vector{Float64}              # saved times (analysis window)
    x::Matrix{Float64}              # n_nodes × n_t
    lam::Union{Nothing,Matrix{Float64}}   # n_edges × n_t, only when recorded
end

window_times(cfg::Config) = collect(cfg.t_burn:cfg.saveat:cfg.t_end)

nonnegative_violated(u, p, t) = any(x -> x < 0, u)

"""
Integrate one trajectory. `updater === nothing` gives a deterministic run with λ = p0;
otherwise `updater` is applied every `interval` time units (first at t = interval).
For C3 the caller sets the initial pool draw in `p0` before calling.
"""
function simulate(net::Network, cfg::Config, p0::Vector{Float64}, u0::Vector{Float64};
                  updater=nothing, interval::Float64=0.0, record_lambda::Bool=false,
                  reltol=cfg.reltol, abstol=cfg.abstol)
    p = copy(p0)
    ts = window_times(cfg)
    cbs = Any[]
    # The recorder must come before the updater: callbacks in a CallbackSet run in order
    # at the end of each step, and the recorder fills in every save time inside that step,
    # so it has to see λ before an update that happens at the step's end.
    saved = nothing
    if record_lambda
        saved = SavedValues(Float64, Vector{Float64})
        push!(cbs, SavingCallback((u, t, integ) -> integ.p[net.lam_idx], saved; saveat=ts))
    end
    if updater !== nothing
        interval > 0 || error("simulate: interval must be > 0 with an updater")
        push!(cbs, PeriodicCallback(updater, interval; save_positions=(false, false)))
    end
    prob = ODEProblem(RHS(net), copy(u0), (0.0, cfg.t_end), p)
    # Expression is non-negative analytically. With long λ hold intervals the adaptive solver can
    # step a near-zero node below 0, where the Hill terms change sign and the solution diverges.
    # isoutofdomain rejects such steps and retries with a smaller one; accepted steps are unchanged.
    sol = solve(prob, Tsit5(); saveat=ts, reltol=reltol, abstol=abstol,
                callback=isempty(cbs) ? nothing : CallbackSet(cbs...),
                isoutofdomain=nonnegative_violated, maxiters=10^8)
    SciMLBase.successful_retcode(sol) || error("ODE solve failed: retcode $(sol.retcode)")
    length(sol.t) == length(ts) || error("ODE solve saved $(length(sol.t)) points, expected $(length(ts))")
    lam = nothing
    if record_lambda
        length(saved.t) == length(ts) || error("λ recorder saved $(length(saved.t)) points, expected $(length(ts))")
        lam = reduce(hcat, saved.saveval)
    end
    return Trajectory(sol.t, Array(sol), lam)
end

"Initial conditions: x_i(0) = U(0,1)·g_i/k_i, one stream per (network, ParamID) (SPEC §8)."
function initial_conditions(cfg::Config, net::Network, pid::Int, p0::Vector{Float64})
    rng = rng_for(cfg.seed_root, net.name, pid, "ic")
    ub = [p0[g] / p0[k] for (g, k) in zip(net.g_idx, net.k_idx)]
    return [rand(rng, nnodes(net)) .* ub for _ in 1:cfg.n_ics]
end
