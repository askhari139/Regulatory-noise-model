# All conditions for one (network, noise mode, ParamID) unit (SPEC §6, §7).

const TABLES = (:occupancy, :transitions, :dwell, :cells, :expression, :lambda, :convergence,
                :replicates, :lfc)

"Run C0–C4 for one parameter set and noise mode. Returns Dict(table => Vector of rows)."
function run_unit(cfg::Config, net::Network, mode::AbstractString, pid::Int,
                  p0::Vector{Float64}, θ::Vector{Float64})
    rows = Dict(t => NamedTuple[] for t in TABLES)
    lo, hi = lambda_bounds(cfg, net)
    ics = initial_conditions(cfg, net, pid, p0)
    seed(parts...) = rng_for(cfg.seed_root, net.name, mode, pid, parts...)

    function emit!(key, acc)
        occupancy_rows!(rows[:occupancy], key, acc, cfg)
        transition_rows!(rows[:transitions], key, acc)
        dwell_rows!(rows[:dwell], key, acc)
        cell_rows!(rows[:cells], key, acc)
        expression_rows!(rows[:expression], key, acc, net)
        replicate_rows!(rows[:replicates], key, acc, net, cfg)
        lfc_rows!(rows[:lfc], key, acc, net, cfg)
    end

    # C0: deterministic, original λ
    acc = CellAccumulator(nnodes(net))
    for r in 1:cfg.n_ics
        add_replicate!(acc, cfg, simulate(net, cfg, p0, ics[r]), θ)
    end
    emit!(cellkey("C0", 0.0, 0.0, 0.0, nreplicates(acc)), acc)

    run_controls = mode in cfg.ctrl_modes
    for dt in cfg.dts, σ in cfg.sigmas
        σ == 0.0 && continue

        # C1: stochastic; record λ in the analysis window to build the pool
        acc = CellAccumulator(nnodes(net); store_lfc=true)
        trajs = Trajectory[]
        for r in 1:cfg.n_traj
            upd = ProcessUpdater(mode, σ, net, lo, hi, p0, seed("C1", σ, dt, r))
            tr = simulate(net, cfg, p0, ics[r]; updater=upd, interval=dt, record_lambda=true)
            add_replicate!(acc, cfg, tr, θ)
            push!(trajs, tr)
        end
        key1 = cellkey("C1", σ, dt, 0.0, nreplicates(acc))
        emit!(key1, acc)
        pool = lambda_rows!(rows[:lambda], key1, trajs, net, lo, hi, cfg)
        trajs = nothing
        run_controls || continue

        # C2: λ drawn once per draw from the pool and held fixed; one IC per draw
        acc = CellAccumulator(nnodes(net))
        rng = seed("C2", σ, dt)
        for d in 1:cfg.n_det_draws
            p = copy(p0)
            p[net.lam_idx] .= pool[:, rand(rng, 1:size(pool, 2))]
            add_replicate!(acc, cfg, simulate(net, cfg, p, ics[mod1(d, cfg.n_ics)]), θ)
        end
        key2 = cellkey("C2", σ, dt, 0.0, nreplicates(acc))
        emit!(key2, acc)
        convergence_rows!(rows[:convergence], key2, acc, cfg, seed("C2boot", σ, dt))

        (dt == cfg.pool_dt && σ in cfg.ctrl_sigmas) || continue

        # C3: λ redrawn i.i.d. from the pool every T_hold (first draw at t = 0)
        for th in cfg.hold_times
            acc = CellAccumulator(nnodes(net); store_lfc=true)
            for r in 1:cfg.n_traj
                rng = seed("C3", σ, th, r)
                p = copy(p0)
                p[net.lam_idx] .= pool[:, rand(rng, 1:size(pool, 2))]
                upd = HoldUpdater(net.lam_idx, pool, rng)
                add_replicate!(acc, cfg, simulate(net, cfg, p, ics[r]; updater=upd, interval=th), θ)
            end
            emit!(cellkey("C3", σ, dt, th, nreplicates(acc)), acc)
        end

        # C4: every λ fixed at its pool mean
        acc = CellAccumulator(nnodes(net))
        p = copy(p0)
        p[net.lam_idx] .= vec(mean(pool; dims=2))
        for r in 1:cfg.n_ics
            add_replicate!(acc, cfg, simulate(net, cfg, p, ics[r]), θ)
        end
        emit!(cellkey("C4", σ, dt, 0.0, nreplicates(acc)), acc)
    end
    return rows
end

"The condition cells a complete unit must contain (used by validation)."
function expected_cells(cfg::Config, mode::AbstractString)
    cells = [("C0", 0.0, 0.0, 0.0, cfg.n_ics)]
    for dt in cfg.dts, σ in cfg.sigmas
        σ == 0.0 && continue
        push!(cells, ("C1", σ, dt, 0.0, cfg.n_traj))
        mode in cfg.ctrl_modes || continue
        push!(cells, ("C2", σ, dt, 0.0, cfg.n_det_draws))
        if dt == cfg.pool_dt && σ in cfg.ctrl_sigmas
            append!(cells, [("C3", σ, dt, th, cfg.n_traj) for th in cfg.hold_times])
            push!(cells, ("C4", σ, dt, 0.0, cfg.n_ics))
        end
    end
    return cells
end
