# regnoise — simulation specification (CSBJ revision rerun)

Status: **APPROVED for implementation** (decisions in §15 confirmed 2026-10-18).
Owner: Kishore Hari. Last updated: 2026-10-18.

This document specifies a clean rewrite of the RACIPE stochastic-λ pipeline
(`/scratch/a.hari/RACIPEdata/{scripts,src}`) that regenerates every ODE result in the
paper plus the new controls requested by the CSBJ reviewers. It replaces, not patches,
the old pipeline.

---

## 0. Design rules (non-negotiable)

1. **Defaults are the production values.** Every default in code equals the value used
   for the paper. A production run needs no overrides except *which* network / mode /
   parameter chunk to run. Anything that is not a production value (debug sizes, short
   tspan) lives only in `config/test.toml`.
2. **One config file is the source of truth.** `config/production.toml` holds every
   number in this document. Code reads it; nothing numeric is hard-coded elsewhere.
   The config's SHA-256 (`config_hash`) is written into every output file.
3. **Fail loudly.** Unknown noise mode, missing input, wrong vector length, solver
   failure, λ outside bounds, or a config hash mismatch → `error()`, never a silent
   fallback or a `@warn` + continue.
4. **No dead code paths.** Only the modes and conditions in this spec exist. The old
   `Lognormal`, `Jumping`, `Extreme`, `MultiplicativeInvLambda`, `MultiplicativeUnbiased`
   modes are not ported.
5. **Every output is self-describing.** Each file carries network, mode, condition,
   σ, dt, T_hold, ParamID, stability class, counts (n_traj / n_draws), seed root,
   config hash and git commit. You never infer settings from a folder name.
6. **No absolute paths in code.** Paths come from the config (`paths.*`) or CLI.
7. **Validation gate.** Results are not "done" until `validate.jl` passes (§11). The
   figure scripts read only validated, collected outputs.
8. **Paired design.** Within a parameter set, all conditions reuse the same initial
   conditions and the same parameter set, so condition differences are paired.

---

## 1. Scope

**In scope:** the RACIPE ODE networks used in the paper.

| Network | Topology | Figures |
|---|---|---|
| TS | A⊣B, B⊣A | 1–6 |
| TSSA | TS + A→A, B→B | 1–6 |
| TT | 3-node mutual inhibition | 1–6 |
| TTSA | TT + self-activation | 1–6 |
| TS4 | two toggles A⊣B, C⊣D with A⊣D, D⊣A, B⊣C, C⊣B | 5 |
| TT4 | 4-node mutual inhibition | 5, 6 |
| TT4SA | TT4 + self-activation | 5, 6 |

Topology and RACIPE files: `data/<NET>.topo`, `<NET>.prs`, `<NET>_parameters.dat`,
`<NET>_solution.dat` (existing 10,000-set RACIPE ensembles; **not regenerated**).

Noise modes: `Additive`, `Multiplicative`, `Stationary` (old code name `Fluctuating`).

**Out of scope (separate spec):** the Boolean EMT / gonadal simulations (Fig. 7) and
the sign-preservation control requested by R1.2 / R3.3.

---

## 2. Model

For node *i* with inputs *j*:

  dx_i/dt = g_i · Π_j H^S(x_j; x0_ij, n_ij, λ_ij) − k_i · x_i

  H⁻(x) = x0ⁿ / (x0ⁿ + xⁿ), H⁺ = 1 − H⁻
  inhibition (λ ≤ 1): H^S = H⁻ + λ H⁺
  activation (λ > 1): H^S = (H⁻ + λ H⁺) / λ

Parameters g, k, x0, n, λ come from `<NET>_parameters.dat` (λ_inh stored as a fraction in
[0.01, 1]; λ_act in [1, 100]). The ODE right-hand side is generated from the `.topo` and
`.prs` files at load time (port `ODEgenerator.jl`); the generated function is unit-tested
against the existing `TS_ode.jl` / `TT4SA_ode.jl` (§12).

---

## 3. λ noise processes

λ is updated at times t_k = k·dt (k ≥ 1) by a `PeriodicCallback`; between updates λ is
constant and the ODE is integrated deterministically. λ(0) = RACIPE value. Bounds
[λ_min, λ_max] per edge type: [0.01, 1] for inhibition, [1, 100] for activation (the RACIPE
sampling range; §15-D1). After every update λ is **clipped** to the
bounds (not "absorbing": the next update can move it back inside — fixes R3 minor 2).
ξ ~ N(0, 1) independent per edge and per update.

| Mode | Update | σ_eff |
|---|---|---|
| Additive | λ ← clip(λ + σ_eff ξ) | σ for inhibition, σ·λ_max for activation |
| Multiplicative | λ ← clip(λ + \|λ* − λ\| σ ξ), λ* = λ_min (inh), λ_max (act) | raw σ |
| Stationary | λ ← clip(λ(0) + σ_eff ξ) | as Additive |
| LogAdditive (candidate, piloted 2026-10) | log10 λ ← clip(log10 λ + σ_eff ξ) | σ × (log10 λ_max − log10 λ_min) = 2σ decades, both edge types |
| Geometric (candidate, piloted 2026-11) | inh: λ ← clip(λ(1 + σξ)); act: λ ← clip(λ / (1 + σξ)); 1 + σξ ≤ 0 → full-strength bound | raw σ |

LogAdditive is the additive walk on the fold-change scale. RACIPE samples the fold-change
uniformly on [1, 100], so inhibitory λ starts concentrated near 0.01; the linear-scale Additive walk
spreads λ uniformly over [0.01, 1] (stationary mean ≈ 0.5, i.e. at most 2-fold repression), which
the first pilot showed drives most of the additive effect. Whether LogAdditive joins production is
decided after pilot 2.

Geometric restores the old inhibitory rule λ(1 + σξ) and mirrors it for activation on 1/λ (in the
shifted Hill an activating edge leaves its target at 1/λ of maximal production, the counterpart of
λ for inhibition). The term multiplied is driven toward 0, which lies outside the bounds, so λ
reflects off the full-strength bound (λ_min for inh, λ_max for act) rather than being absorbed as
under the target-based Multiplicative rule, where λ* is the bound itself. Stationary law,
independent of σ for σ ≲ 0.1: the full-strength fold change (1/λ for inh, λ for act) is uniform
on [1, 1/λ_min] or [1, λ_max]. Without √Δt scaling the fluctuation rate is set by σ²/Δt.
Production uses an inhibitory floor of 0.001 for this mode (§15-D12). The floor matters: when λ
fluctuates faster than expression, inhibition acts through the time-averaged λ,
E[λ] = λ_min(1 + ln(1/λ_min)), i.e. about 127-fold repression at 0.001 versus 18-fold at 0.01
(pilot 3: at 0.01, TT4/TT4SA lose almost all the single-high enrichment of the old runs).

At σ = 0 every rule is the identity (unit-tested). The manuscript's Eq. 10 must be
corrected to `λ(t+Δt) = λ(t) + |λ* − λ(t)| η` (R3.1 is a typo in the paper; the code
already implements this form).

σ grid: {0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0}.
dt grid (λ update interval): {0.01, 1.0, 10.0}.

---

## 4. Parameter-set selection (`select.jl`, run once per network)

1. Discretize every RACIPE steady state with the thresholds of §5.
2. Stability class of a parameter set = number of **distinct discrete states** among its
   RACIPE steady states (1 = monostable, 2 = bistable, 3 = tristable, 4 = tetrastable).
3. For each class present, draw `n_per_class = 1000` sets without replacement (all sets
   if the class has fewer), uniformly at random (§15-D2). Optional `subsample_of = N`
   (pilots): first make the selection that `n_per_class = N` would make (same RNG stream),
   then draw `n_per_class` of those, so a small pilot runs on a subset of a larger run's sets.
4. Save `selection/<NET>.csv`: ParamID, StabilityClass, selection weight, class size in
   the full ensemble, class prevalence (class size / 10,000), seed.
5. **The same selection is used by every noise mode and every condition.** Selection is
   never re-run implicitly; `run.jl` errors if the selection file is missing.

Classes run: monostable, bistable, tristable (tetrastable only if §15-D3 says so).

---

## 5. Discretization

Threshold for node i: θ_i = 2^( Σ_s f_s log2 x_{i,s} / Σ_s f_s ), the RACIPE
convergence-frequency-weighted geometric mean of node i across all 10,000 parameter sets'
steady states (unchanged from the old code; document as such — R3 minor 3).
A sample is in state b = (1[x_1 > θ_1], …, 1[x_N > θ_N]), written in outputs as letters in
node order (H = high, L = low; e.g. TT all-high = "HHH") so no CSV reader parses it as a number.

Threshold sensitivity (R3.6): every occupancy observable is also computed at
θ × m for m ∈ {0.5, 2.0}; stored alongside m = 1.

---

## 6. Conditions (per network × noise mode × ParamID)

| Code | Condition | λ(t) | Grid | Replicates |
|---|---|---|---|---|
| C0 | Deterministic, original λ | RACIPE λ, fixed | — | 100 ICs |
| C1 | Stochastic | noise process (§3) | σ ∈ grid \ {0} × dt ∈ grid | 100 trajectories |
| C2 | Deterministic resampled | one λ vector per draw, held fixed, drawn from C1's pool (§7) | same (σ, dt) as C1 | 1000 draws × `n_det_ics` ICs each (default 1) |
| C3 | Hold-resample | redrawn i.i.d. from C1's pool every T_hold | σ ∈ S_ctrl × T_hold ∈ {0.01, 0.1, 1, 10, 100}; pool from dt = 0.01 | 100 trajectories |
| C4 | Deterministic, mean λ | each edge fixed at its pool mean | σ ∈ S_ctrl; pool from dt = 0.01 | 100 ICs |

Notes:
- C3 with T_hold = 0.01 is the shuffled-surrogate control (R1.4): same marginal, no memory.
  Large T_hold approaches C2 (frozen λ), but is measured the same way as C1.
- C3 differs from the C1 dt sweep: in C1 a larger dt keeps the per-step size σ (so it
  also changes the diffusion rate and the realized marginal); in C3 only the hold time
  changes (R1.1, R3.2, R3.5.2–3).
- C0 → C4 → C2 → C1 decomposes the effect into mean strength, distribution shape, and
  temporal structure (R3.5.4).
- S_ctrl = {0.01, 0.1, 1.0}. C2–C4 run for all three noise modes (§15-D4, D5).
- C2–C4 depend on C1's λ pool, so a (ParamID, mode) unit runs C1 before C2–C4 in the
  same process; nothing about the pool is written to disk and re-read by another job.

---

## 7. The λ pool (matched marginals)

During every C1 trajectory, record the full λ vector (all edges jointly) at every
**saved time point inside the analysis window** (the value in effect just before any
update scheduled at that same time) (t ∈ [T_burn, T], spacing `saveat`).
Pool(σ, dt) = all recorded vectors over the 100 trajectories (100 × 501 vectors).

- C2 draws whole vectors uniformly from Pool(σ, dt) (preserves cross-edge structure
  and each edge's realized marginal, including any non-stationarity at low σ).
- C3 redraws whole vectors uniformly from Pool(σ, 0.01) every T_hold.
- C4 uses the per-edge mean of Pool(σ, 0.01).
- The old global, parameter-independent histograms (`lambda_hist_*.jld2`, sampled over
  t ∈ [50, 100]) are **not used** anywhere.

Store per (ParamID, σ, dt) the pool summary (§9.4) so marginals can be compared across
conditions and dt (R3.2).

---

## 8. Simulation settings

| Setting | Value |
|---|---|
| tspan | (0, 1000) |
| analysis window | t ∈ [500, 1000] (burn-in = 50%) |
| saveat | 1.0 (uniform samples; occupancy = fraction of samples) |
| solver | Tsit5, reltol 1e-4, abstol 1e-5 (§15-D7) |
| initial conditions | x_i(0) = U(0,1) · g_i / k_i, independent per node |
| IC pairing | IC r is identical across C0, C1 (all σ, dt), C3, C4 for the same ParamID; C2 draw d, IC j uses IC ((d−1)·n_det_ics + j mod 100); replicates are stored draw-major |
| non-negativity | `isoutofdomain`: steps that would make any expression negative are rejected and retried with a smaller step (added 2026-10-24; without it, long λ holds let the solver overshoot a near-zero node below 0, where the Hill terms change sign and self-activating networks diverge) |
| retcode | any non-`Success` retcode → error (logged with seed for replay) |

Blocks for convergence (R1 minor 3, R3.6): occupancy is also reported for the five
100-unit blocks [500,600), …, [900,1000].

---

## 9. Observables

Per trajectory (or draw) then aggregated per condition cell. "Occupancy" replaces the
old term "MRT" everywhere in code and output (R1 minor 1).

9.1 **Occupancy** per state: mean, SD, and SE across replicates; per block (§8); at
threshold multipliers (§5). All 2^N states are written, including zeros (fixes the
old missing-row bug).

9.2 **Dwell times** per state: mean and number of completed visits (true residence times).

9.3 **Transitions**: count matrix of state → state switches (from → to), summed over
replicates, plus n_replicates; total switches mean / SD. Enables in/out rates for the
all-high state (R1.5).

9.4 **λ pool summary** per edge (C1 only): mean, SD, quantiles (1, 5, 25, 50, 75, 95, 99 %),
fraction of samples at λ_min and at λ_max, autocorrelation at lags {1, 5, 10, 50} saved
samples.

9.5 **C2 convergence** (on per-draw occupancy, i.e. averaged over a draw's `n_det_ics` ICs): occupancy estimates from the first n draws for
n ∈ {10, 30, 100, 300, 1000}, plus bootstrap SE (200 resamples) at n = 1000 (R1 minor 2,
R3.5.1).

9.6 **Node expression** per node (analysis window): mean and SD, averaged over replicates.
Note: averaging over replicates blends replicates that sit in different states; use 9.7 for
distributions and joint states.

9.7 **Per-replicate expression** (`replicates` table): for the first `replicates_stored` (100)
replicates of every cell (for C2, the first 100 of the 1,000 draws), the time-averaged mean and the
within-trajectory SD of every node, one row per replicate (columns `Mean_<node>`, `SD_<node>`,
Float32). This replaces the old pipeline's rounded-and-grouped per-node statistics: states and
joint node configurations are kept exactly, and clustering (e.g. round and count) is done in
post-processing.

9.8 **Windowed log fold-change** (`lfc` table, time-varying conditions C1 and C3 only): for each
node, lag τ ∈ {1, 3, 10, 30, 100} and threshold c ∈ {1, 2} (2- and 4-fold), the fraction of time
points t in the analysis window with |log2 x(t+τ) − log2 x(t)| ≥ c, per stored replicate
(columns `<node>_Frac_lag<τ>_<k>fold`, UInt16 = fraction × 65535). The natural lag is the node's
relaxation time 1/k (1–10 time units); storing several lags lets the analysis pick it per node.
Deterministic conditions (C0, C2, C4) are at a fixed point in the window, so the measure is 0 by
construction and is not stored; their per-replicate SD (9.7) shows any oscillation.

---

## 10. RNG and provenance

- Seeds are derived, never global: seed = hash(seed_root, NET, mode, ParamID, condition,
  σ, dt, T_hold, replicate) → `Xoshiro`. Same inputs ⇒ bit-identical outputs on rerun.
- `seed_root` is in the config. The IC stream depends only on (NET, ParamID, replicate),
  which gives the pairing in §8.
- Each output file header/metadata: config_hash, git commit (refuse to run with a dirty
  tree unless `--allow-dirty`), Julia version, `Manifest.toml` hash, hostname, start/end
  time, wall time.
- Package versions are pinned in `Project.toml` / `Manifest.toml` committed with the code
  (R3 minor 8).

---

## 11. Outputs, collection, validation

Layout (one run = one `run_id`, e.g. `2026-10-csbj`):

```
<results_root>/<run_id>/
  config.toml                 # frozen copy of the config used
  selection/<NET>.csv
  raw/<NET>/<MODE>/param_<ID>.jld2    # all conditions for one ParamID
  logs/<NET>/<MODE>/param_<ID>.log
  collected/<NET>_<MODE>_{occupancy,transitions,dwell,lambda,convergence,expression}.csv.gz
  validation/<NET>_<MODE>.txt
```

Long-format columns common to all tables: run_id, config_hash, Network, NoiseMode,
Condition, Sigma, DT, THold, ParamID, StabilityClass, NReplicates, ThresholdMult, Block.

`run.jl` is resumable: it skips a ParamID whose raw file exists **and** has the same
config_hash and is complete; a file with a different hash is an error, not a skip.

`collect.jl` validates first, then merges raw files into one gzipped CSV per table:
`occupancy`, `transitions`, `dwell`, `cells` (switch counts), `expression`, `lambda`,
`convergence`. It refuses to merge mixed config hashes.

`validate.jl` checks, and fails on any violation:
- every selected ParamID present for every expected condition cell;
- NReplicates equals the config value in every row (100 / 1000 — the check that would
  have caught the 10-vs-100 confusion);
- occupancies in [0, 1] and summing to 1 (± 1e-9) per replicate set;
- λ pool within bounds (per-edge min/max of every C1 pool).
- no NaN in occupancy, switches, expression, λ summaries or convergence estimates;
  failed ODE solves raise an error in `run.jl`, so a unit file never exists for them.

Two properties are checked exactly by unit tests instead of statistically on outputs:
σ = 0 leaves λ unchanged in every mode, and C2/C3 draw only whole vectors from the pool.

---

## 12. Tests (`test/runtests.jl`, must pass before any cluster run)

1. Update rules: σ = 0 ⇒ λ unchanged for all modes; clipping at both bounds; Additive
   activation uses σ·λ_max; Multiplicative step ∝ |λ* − λ|.
2. Generated ODE equals the legacy `TS_ode.jl` / `TTSA_ode.jl` RHS on random inputs.
3. Determinism: same seed ⇒ identical C1 output; different replicate ⇒ different.
4. Pool and resampling: C2/C3 draws come only from pool vectors; C3 switches exactly at
   multiples of T_hold.
5. Observables: occupancy / dwell / transitions on hand-built state sequences.
6. Regression: on 5 TS ParamIDs with the legacy settings (bounds, tolerance, histogram
   pool, 10 draws), the new code reproduces the old `final/Additive/TS` occupancies
   within Monte-Carlo error — proves the rewrite changed only what we intended.
7. End-to-end with `config/test.toml` (3 ParamIDs, tspan 100, 5 replicates) through
   select → run → collect → validate.

---

## 13. Execution and cost

CLI: `julia --project run.jl --config config/production.toml --network TS --mode Additive --chunk 7/40`
SLURM: one array job per (network, mode); each task runs a contiguous chunk of
ParamIDs serially, single-threaded (`sbatch/array.sh`). No nested job submission.

Measured cost of one full unit (C0–C4, all σ and dt, production config), single core,
local Apple-silicon machine (2026-10-18): TS ≈ 90–100 s; TT4SA ≈ 270–320 s (one unit per
stability class). Cluster pilot (2026-10-20): TT 4–5 min, TT4 6–8 min per unit (one core).

Estimate for the cluster, assuming a cluster core is ~2× slower than the local one
(the pilot measures this) and interpolating the other networks between TS and TT4SA:

| | CPU-h |
|---|---|
| One noise mode, all 7 networks (~19,400 units) | ~2,200 |
| All three modes | ~6,500 (~50 h on one 128-core node) |
| Raw storage per noise mode (incl. per-replicate tables, measured 2026-10-22: TS 0.86 MB, TT 1.22 MB, TT4SA 1.75 MB per unit) | ~27 GB |
| Raw storage, three modes (four with LogAdditive) | ~80 GB (~110 GB); collected CSV.gz about the same again |

Tolerance check in the pilot: 10 ParamIDs per network rerun at reltol 1e-6 / abstol 1e-8;
occupancy differences must be within Monte-Carlo SE (reported, not a hard gate).

Run order: (1) tests; (2) pilot: TS + TT, Additive, 50 ParamIDs, full conditions →
check validation, timings, and that C2 at n = 1000 has converged; (3) full Additive;
(4) Multiplicative; (5) Stationary.

---

## 14. Reviewer mapping

| Comment | Addressed by |
|---|---|
| R1.1, R3.2 amplitude vs timescale | C3 T_hold ladder (marginal fixed); λ pool ACF/quantiles across dt (§9.4) |
| R1.3 multiplicative drift | C2/C4 for Multiplicative (distribution vs dynamics); §15-D6 |
| R1.4 shuffled surrogate | C3, T_hold = 0.01 |
| R1.5 derepression at low noise | transition matrices → in/out rates of all-high vs σ (§9.3) |
| R1 minor 1 "MRT" | renamed occupancy; dwell times reported (§9.2) |
| R1 minor 2, R3.5.1 draws | C2 with 1000 draws + convergence table (§9.5) |
| R1 minor 3, R3.6 convergence | block occupancies (§8); burn-in sensitivity from blocks |
| R2.2 stochastic model details | §2, §3, §8 become the Methods text |
| R2.3, R3.6 weighting | §4 selection file with class prevalences → reweighted means |
| R3.1 Eq. 10 | §3 (paper typo; formula fixed) |
| R3.5.2–4 matched marginals | C2–C4 drawn from the same pool |
| R3.6 thresholds | threshold multipliers (§5) |
| R3 minor 2 "absorbing" | §3 clipping |
| R3 minor 8 reproducibility | §10 |

---

## 15. Decisions (confirmed 2026-10-18)

| # | Decision | Choice |
|---|---|---|
| D1 | λ bounds | [0.01, 1] inhibition, [1, 100] activation (RACIPE range) |
| D2 | Selection within class | Uniform random without replacement; class prevalences saved for reweighting |
| D3 | Classes | Mono-, bi-, tristable (default; tetrastable not run) |
| D4 | Control σ subset for C3/C4 | {0.01, 0.1, 1.0} |
| D5 | Modes with C2–C4 | All three (Additive, Multiplicative, Stationary) |
| D6 | Extra multiplicative control | None; R1.3 answered by the C0→C4→C2→C1 decomposition (default) |
| D7 | Solver tolerance | Keep reltol 1e-4 / abstol 1e-5. At dt = 0.01 the callback caps step size at 0.01 so tolerance barely matters; elsewhere it can only change basin assignment for ICs near separatrices. Verified by the pilot tolerance check (§13). |
| D8 | Trajectory length | tspan = 1000, 50 % burn-in, block diagnostics (default) |
| D9 | Expression outputs (2026-10-22) | Per-replicate mean and SD of every node, first 100 replicates per cell (9.7) |
| D10 | Continuous dynamics measure (2026-10-22) | Windowed log fold-change at lags 1–100, 2- and 4-fold (9.8), C1 and C3 only |
| D11 | LogAdditive noise mode (2026-10-22) | Piloted on TT (pilot 2) before deciding on production |
| D12 | Multiplicative noise (2026-11-06) | Geometric mode replaces the target-based rule in the paper; inhibitory floor 0.001 (1000-fold maximal repression) for this mode only: induced-fit binding and epigenetic locking are stronger than typical transcriptional repression. Run `2026-11-csbj-geometric` (config/production_geometric.toml); Additive and Stationary stay in `2026-10-csbj-v2` with floor 0.01 |

## 16. Changes from the old pipeline (for Methods / rebuttal)

1. Deterministic resampling: 10 → 1000 draws; λ drawn from the run's own recorded λ
   (analysis window) instead of global histograms sampled over t ∈ [50, 100]; 1 IC per
   draw instead of 100; measured as time occupancy like C1.
2. Multiplicative mode: one definition (target-based) for all networks. The old
   `MultiplicativeInvLambda` data used for the self-activation networks in Fig. 5 ran a
   constant-step additive walk on activation edges and is superseded.
3. New controls C3 (hold-resample ladder) and C4 (mean-λ).
4. Fixed derived seeds; paired ICs across all conditions.
5. All states written including zero occupancy; occupancy terminology; dwell times.
6. λ bounds [0.01, 1] / [1, 100] (old code: [0.001, 0.999] / [1.001, 100]); uniform
   selection within class (old: balance-weighted); tolerance unchanged (1e-4 / 1e-5).
7. Single code location (git), single config, validation gate.
8. Solver failures are errors. The old pipeline never checked the ODE return code, so a diverged
   trajectory (seen in production for TSSA, Multiplicative, Δt = 10) would have been discretized
   silently (NaN counts as "low"). The new code rejects negative-expression steps and errors on any
   failed solve. Production run `2026-10-csbj` was cancelled after this was found and restarted as
   `2026-10-csbj-v2` with the fix, so the whole run comes from one commit.
