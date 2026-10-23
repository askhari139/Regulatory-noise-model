# Run with: REGNOISE_DATA=/path/to/RACIPE/data julia --project test/runtests.jl
using RegNoise, Test, Random, DataFrames, Statistics
using OrdinaryDiffEq
const R = RegNoise

haskey(ENV, "REGNOISE_DATA") || error("set REGNOISE_DATA to the folder holding the RACIPE input files")
const TESTCFG = joinpath(@__DIR__, "..", "config", "test.toml")
# Legacy generated ODE files copied from /scratch/a.hari/RACIPEdata/src (reference fixtures)
const LEGACY_SRC = joinpath(@__DIR__, "legacy_ode")

cfg_tmp() = load_config(TESTCFG; data_dir=ENV["REGNOISE_DATA"], results_root=mktempdir())

@testset "RegNoise" begin

@testset "config" begin
    cfg = cfg_tmp()
    @test cfg.inh_bounds == (0.01, 1.0) && cfg.act_bounds == (1.0, 100.0)
    prod = load_config(joinpath(@__DIR__, "..", "config", "production.toml"))
    @test prod.n_det_draws == 1000 && prod.n_traj == 100 && prod.n_per_class == 1000
    @test prod.ctrl_modes == ["Additive", "Multiplicative", "Stationary"]
    # [paths] does not enter the hash
    a = load_config(TESTCFG; data_dir="/x", results_root="/y")
    b = load_config(TESTCFG; data_dir="/z", results_root="/w")
    @test a.hash == b.hash
end

@testset "λ update rules" begin
    for mode in (:Additive, :Multiplicative, :Stationary, :LogAdditive), ξ in (-2.3, 0.0, 1.7)
        # σ = 0 is the identity whenever λ = λ0
        @test R.update_lambda(mode, 0.3, 0.3, 0.0, 0.0, 0.01, 0.01, 1.0, ξ) == 0.3
    end
    # clipping at both bounds
    @test R.update_lambda(:Additive, 0.5, 0.5, 1.0, 1.0, 0.01, 0.01, 1.0, 10.0) == 1.0
    @test R.update_lambda(:Additive, 0.5, 0.5, 1.0, 1.0, 0.01, 0.01, 1.0, -10.0) == 0.01
    # multiplicative step is proportional to the distance to the target
    @test R.update_lambda(:Multiplicative, 0.21, 0.21, 0.1, 0.0, 0.01, 0.01, 1.0, 1.0) ≈ 0.21 + 0.2 * 0.1
    @test R.update_lambda(:Multiplicative, 20.0, 20.0, 0.1, 0.0, 100.0, 1.0, 100.0, 1.0) ≈ 20.0 + 80.0 * 0.1
    # log-additive: step in log10 λ, clamped to the bounds
    @test R.update_lambda(:LogAdditive, 0.1, 0.1, 0.1, 0.2, 0.01, 0.01, 1.0, 1.0) ≈ 0.1 * 10^0.2
    @test R.update_lambda(:LogAdditive, 0.5, 0.5, 1.0, 2.0, 0.01, 0.01, 1.0, 3.0) == 1.0
    @test R.update_lambda(:LogAdditive, 0.5, 0.5, 1.0, 2.0, 0.01, 0.01, 1.0, -3.0) == 0.01
    # stationary always restarts from λ0
    @test R.update_lambda(:Stationary, 0.9, 0.3, 0.1, 0.1, 0.01, 0.01, 1.0, 1.0) ≈ 0.4
    # activation edges use σ·λ_max for Additive
    cfg = cfg_tmp()
    net = load_network(cfg, "TSSA")
    lo, hi = lambda_bounds(cfg, net)
    upd = R.ProcessUpdater("Additive", 0.1, net, lo, hi, zeros(length(net.param_names)), Xoshiro(1))
    @test all(upd.sigma_eff[k] == (e.act ? 10.0 : 0.1) for (k, e) in enumerate(net.edges))
    @test_throws ErrorException R.ProcessUpdater("Lognormal", 0.1, net, lo, hi, zeros(16), Xoshiro(1))
    # log-additive step = σ × log10 range (2 decades) for both edge types
    updl = R.ProcessUpdater("LogAdditive", 0.1, net, lo, hi, zeros(length(net.param_names)), Xoshiro(1))
    @test all(updl.sigma_eff .≈ 0.2)
end

@testset "RHS matches legacy generated ODE ($n)" for n in ("TS", "TSSA", "TT", "TTSA", "TS4", "TT4", "TT4SA")
    cfg = cfg_tmp()
    net = load_network(cfg, n)
    params = read_parameters(cfg, net)
    m = Module()
    Base.include(m, joinpath(LEGACY_SRC, "$(n)_ode.jl"))
    rng = Xoshiro(7)
    for pid in first(sort(collect(keys(params))), 25)
        p = params[pid]
        u = 20 .* rand(rng, length(net.nodes))
        du1 = zeros(length(u)); du2 = zeros(length(u))
        R.RHS(net)(du1, u, p, 0.0)
        Base.invokelatest(m.ode_system!, du2, u, p, 0.0)
        @test du1 ≈ du2 rtol = 1e-12
    end
end

@testset "seeds" begin
    @test R.derive_seed(1, "TS", 0.1) == R.derive_seed(1, "TS", 0.1)
    @test R.derive_seed(1, "TS", 0.1, 1) != R.derive_seed(1, "TS", 0.1, 2)
    @test R.derive_seed(1, "TS", 0.1) != R.derive_seed(1, "TS", 0.10000000000000002)
end

@testset "hold-resample switches only at multiples of T_hold" begin
    cfg = cfg_tmp()
    net = load_network(cfg, "TS")
    pid = first(sort(collect(keys(read_parameters(cfg, net)))))
    p0 = read_parameters(cfg, net)[pid]
    pool = [0.02 0.5 0.9; 0.03 0.6 0.95]           # 2 edges × 3 vectors
    rng = Xoshiro(3)
    p = copy(p0); p[net.lam_idx] .= pool[:, 1]
    tr = R.simulate(net, cfg, p, [1.0, 1.0]; updater=R.HoldUpdater(net.lam_idx, pool, rng),
                    interval=10.0, record_lambda=true)
    cols = [tr.lam[:, j] for j in axes(tr.lam, 2)]
    @test all(any(c ≈ pool[:, k] for k in 1:3) for c in cols)
    # within each open interval (10k, 10k+10) λ is constant
    for k in 5:9
        js = findall(t -> 10k < t < 10k + 10, tr.t)
        @test all(cols[j] == cols[js[1]] for j in js)
    end
end

@testset "observables on a hand-built trajectory" begin
    cfg = cfg_tmp()                                   # window 50:1:100 → 51 samples, blocks of 10
    t = collect(50.0:1.0:100.0)
    θ = [1.0, 1.0]
    # A high for t < 70, then both high until 90, then B high
    x = reduce(hcat, [ti < 70 ? [2.0, 0.5] : ti < 90 ? [2.0, 2.0] : [0.5, 2.0] for ti in t])
    acc = R.CellAccumulator(2)
    R.add_replicate!(acc, cfg, R.Trajectory(t, x, nothing), θ)
    occ = acc.occ[1]
    # codes: A high only = 1 ("HL"), both high = 3 ("HH"), B high only = 2 ("LH")
    @test occ[1+1] ≈ 20 / 51 && occ[3+1] ≈ 20 / 51 && occ[2+1] ≈ 11 / 51
    @test acc.trans[1+1, 3+1] == 1 && acc.trans[3+1, 2+1] == 1 && sum(acc.trans) == 2
    # only the middle run is a completed visit
    @test acc.dwell_n == [0, 0, 0, 1] && acc.dwell_sum[4] ≈ 20.0
    @test acc.block_occ[1][1+1, 1] == 1.0 && acc.block_occ[1][2+1, 5] ≈ 11 / 11
    @test R.state_label(1, 2) == "HL" && R.state_label(2, 2) == "LH" && R.state_label(3, 2) == "HH"
end

@testset "windowed log fold-change" begin
    cfg = cfg_tmp()                                    # lags 1, 3, 10, 30; thresholds 2- and 4-fold
    nt = 51
    x = vcat([isodd(j) ? 1.0 : 4.0 for j in 1:nt]', fill(3.0, 1, nt))   # A alternates 4-fold, B constant
    F = R.decode_frac.(R.lfc_fractions(cfg, x))
    # columns lag-major: (1,2f) (1,4f) (3,2f) (3,4f) (10,2f) (10,4f) (30,2f) (30,4f)
    @test F[1, :] ≈ [1, 1, 1, 1, 0, 0, 0, 0]          # odd lags see the 4-fold jump, even lags do not
    @test all(F[2, :] .== 0)
    @test R.lfc_colnames(cfg)[1:2] == [:Frac_lag1_2fold, :Frac_lag1_4fold]
end

@testset "selection" begin
    cfg = cfg_tmp()
    s1 = select_parameters(cfg, "TS")
    s2 = select_parameters(cfg, "TS")
    @test s1 == s2
    @test Set(s1.StabilityClass) ⊆ Set(cfg.classes)
    @test all(s1.NSelected .== min.(cfg.n_per_class, s1.ClassSize))
end

@testset "end to end: run, write, validate, determinism" begin
    cfg = cfg_tmp()
    for netname in ("TS", "TSSA")
        select_parameters(cfg, netname)
        net = load_network(cfg, netname)
        sel, θ, _ = load_selection(cfg, netname)
        params = read_parameters(cfg, net)
        lo, hi = lambda_bounds(cfg, net)
        pid = sel.ParamID[1]
        for mode in cfg.noise_modes
            rows = run_unit(cfg, net, mode, pid, params[pid], θ)
            path = write_unit(cfg, net, mode, pid, sel.StabilityClass[1], θ, rows, 0.0)
            @test isempty(validate_unit(cfg, path, mode, lo, hi))
            meta, T = R.read_unit(path)
            @test meta["stability_class"] isa String && meta["noise_mode"] isa String
            @test all(eltype(c) <: Union{Number,String} for t in values(T) for c in eachcol(t))
            @test nrow(T["replicates"]) == sum(min(c[5], cfg.replicates_stored) for c in expected_cells(cfg, mode))
            @test "Mean_A" in names(T["replicates"]) && "A_Frac_lag1_2fold" in names(T["lfc"])
            @test unit_status(cfg, netname, mode, pid) == :done
        end
        rows1 = run_unit(cfg, net, "Additive", pid, params[pid], θ)
        rows2 = run_unit(cfg, net, "Additive", pid, params[pid], θ)
        @test DataFrame(rows1[:occupancy]) == DataFrame(rows2[:occupancy])
        @test DataFrame(rows1[:lambda]) == DataFrame(rows2[:lambda])
        # C2 draws only pool values: its convergence rows exist at every n
        conv = DataFrame(rows1[:convergence])
        @test Set(conv.NDraws) == Set(cfg.convergence_ns)
    end
    # a unit produced under another config is refused
    other = load_config(TESTCFG; data_dir=ENV["REGNOISE_DATA"], results_root=cfg.results_root)
    @test other.hash == cfg.hash
end

end
