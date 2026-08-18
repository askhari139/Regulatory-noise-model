# Methods vs. code audit

This document records discrepancies found between `main.tex`'s Methods section
(and `boolean_static.tex`'s second Methods block) and the actual simulation
code in `codes/julia/` and `codes/R/` (mirroring `RACIPEdata/`). Compiled from a
line-by-line audit run against the codebase as of 2026-08-18. Organized by
whether the discrepancy affects any reported figure/result.

## Affects reported results — manuscript text already updated for these

- **Burn-in fraction.** Methods said "the first 30% of each trajectory was
  discarded as transient." The code (`src/StochasticSimulations.jl`,
  `run_multiple_stochastic_simulations`, default `cut_fraction = 0.5`) actually
  discards the first **50%**, and `analyze_noise_effects` — the function
  actually called by `scripts/3_analyze_single_parameter.jl` for every
  production run — has no `cut_fraction` parameter at all, so this can't be
  overridden anywhere in the real call chain. **Fixed in `main.tex`** (now
  reads 50%). You're separately updating Fig 2B's crosshatch-region rendering
  to match.
- **Parameter-set count per stability class.** Methods said "a representative
  subset of 100 parameter sets was selected from each stability class." The
  current `final/` dataset actually used for every reported figure has up to
  1000 parameter sets per class (e.g. TS: 1000 monostable + 1000 bistable);
  "100" only describes the older `final_small/` dataset. **Fixed in `main.tex`**
  (now reads "up to 1000... fewer where a stability class had fewer available
  parameter sets"), in three places (Results line ~159, Methods line ~372).
- **Reachability claim (Fig 3C).** Methods/Results said "Reachability was
  satisfied across all sampled parameters." This was disproved by this
  session's `fill_mrt()`/`add_reachability()` fix: before the fix, ParamIDs
  that never visited a state were silently absent from the results instead of
  correctly showing MRT = 0, making it look like every parameter was
  reachable. Recomputed fresh against the current dataset: for TS, 89.0% of
  parameter sets are reachable for both nodes, 10.8% for only one, 0.3% for
  neither; every non-fully-reachable monostable parameter set has MRT = 0 for
  the all-high state (mean and max both exactly 0, confirmed directly).
  **Fixed in `main.tex`** (Results paragraph rewritten).

## Naming/documentation inconsistencies — code is correct, wording should match it (not yet fixed in `main.tex`; flagging for the authors)

- **"Stationary" (manuscript) vs. `"Fluctuating"` (code).** The manuscript's
  third noise mode, "Stationary" (`λ(t+dt) = λ(0) + N(0,σ_eff²)`), is
  implemented in code as noise mode `"Fluctuating"`
  (`scripts/lambda_sampler.jl:49-50`, `src/StochasticSimulations.jl:138-139`).
  The formula matches exactly — only the name differs. Either rename the
  manuscript's mode to "Fluctuating" or note the correspondence explicitly, to
  avoid confusion for anyone reading the code alongside the paper.
- **dt sweep undersold.** Methods states "dt = 0.01" as if it were the only
  value used. Production runs actually sweep dt ∈ {0.01, 1.0, 10.0} for every
  network/noise-mode combination (`final/sim_script.sh`, confirmed present in
  the results CSVs' `DT` column) — this is exactly what the existing
  Fig4S3/Fig5S3 supplementary figures show. The Methods prose should mention
  that the dt sweep exists (even briefly), rather than reading as dt=0.01-only.
- **Undocumented extra noise modes.** The code implements three noise modes
  beyond Additive/Multiplicative/Stationary(Fluctuating) that the manuscript
  never mentions: `"Lognormal"`, `"Jumping"`, `"Extreme"`
  (`scripts/lambda_sampler.jl:47,51-54`, `src/StochasticSimulations.jl:136,140-143`).
  Please confirm none of these were used for any reported figure — if so, no
  manuscript change is needed, but it's worth knowing they exist in case
  someone reruns the pipeline and picks one by accident.
- **Two different `get_lambda_dist`-equivalent implementations exist**:
  `src/StochasticSimulations.jl:318`'s `get_lambda_dist` (100 trajectories ×
  10,000 steps, burn-in at step 5000) and `src/simulate_racipe_resampled.jl`'s
  `build_lambda_dist` (50,000 steps, 10,000-step burn-in). Please confirm which
  one actually produced the resampled-λ distributions used in Figure 4/Figure
  5's deterministic-resampling comparisons, so the Methods' description of
  this step points at the right implementation.

## Internal-only inconsistencies — dead code paths, no effect on any reported figure

- `simulate_with_noise_tracked` (used only for λ-trajectory tracking/diagnostics,
  not on the path that produces any reported MRT figure) omits the explicit
  `reltol=1e-4, abstol=1e-5` that the rest of the pipeline sets, falling back
  to DifferentialEquations.jl's own Tsit5 defaults instead.
- The internal `u0` fallback inside `src/StochasticSimulations.jl` (never
  triggered in production, since `u0List` is always pre-populated by
  `scripts/3_analyze_single_parameter.jl` before this fallback would run) uses
  an extra `1.5×` factor not matching the manuscript's
  `u_i(0) = U(0,1) × (g_i/d_i)` formula — which the actual production
  initial-condition code does implement correctly.

## Not fully verified — flagged for follow-up, not asserted as fact

- **Boolean/large-network noise parameters.** Methods states σ² = 0.01 noise
  variance, and 1000 trajectories × 10,000 steps for biological networks vs.
  100 × 1,000 for artificial networks. A quick check of
  `Boolean/script_randCont.jl` found a *swept* range of noise values
  (0.001–0.05, not a single fixed σ²=0.01) and `nInit=1000, nIter=10000` in a
  script whose name suggests it might apply to *artificial* (random-topology)
  networks — which would contradict the manuscript's 100×1,000 claim for that
  category. This wasn't resolved within the audit's time budget (the
  biological-network-specific script wasn't located, and it's unclear whether
  "rand" in the filename refers to the manuscript's "artificial networks" or
  something else). Recommend a dedicated follow-up review of the `Boolean/`
  directory before trusting either number.

## Checked and confirmed matching (no discrepancy)

- RACIPE ensemble size (10,000 parameter sets): `src/RACIPErunner.jl:29`,
  `scripts/1_select_parameters.jl:87`, `data/TT4.cfg`.
- 100 random initial conditions for RACIPE's own steady-state search (external
  `racipemt` binary default, `NumberOfRIVs 100` in `data/TT4.cfg`).
- Integration window t ∈ [0, 1000]: `src/StochasticSimulations.jl:530`,
  `scripts/3_analyze_single_parameter.jl:148`.
- Tsit5 solver, reltol=1e-4/abstol=1e-5 on the production path:
  `src/StochasticSimulations.jl:157-162`.
- 100 stochastic trajectories per parameter set/noise level (code default is
  25, but every production run explicitly sets `NUM_SIMS=100`, matching the
  manuscript).
- Initial condition formula `u_i(0) = U(0,1) × (g_i/d_i)` on the production
  path: `scripts/3_analyze_single_parameter.jl:137-139`.
- Full noise-level set σ ∈ {0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1} (code
  default omits 0.5/1.0, but production runs explicitly set the full set).
- λ clamp bounds ([0.001, 0.999] inhibitory, [1.001, λmax] activatory):
  `src/StochasticSimulations.jl:97-149`, `src/simulate_racipe_resampled.jl:33-41`.
- Deterministic λ-resampling: 10 resampled draws per parameter set
  (`DET_ITERS=10` default, `scripts/3_analyze_single_parameter_det.jl:179`).
- Additive/Multiplicative noise-update formulas match the code exactly
  (`scripts/lambda_sampler.jl:36,38`, `src/StochasticSimulations.jl:117,119`),
  as does the σ_eff = σ·λ_max scaling (`src/StochasticSimulations.jl:114`).
- Figure 1C's specific percentages/counts (100%, 99.2%, 77.5%, 46.2%;
  9/12/97/119 distinct attractor combinations for TS/TSSA/TT/TTSA) were
  recomputed directly against `get_racipe_attractor_combo()`'s exact logic and
  match the manuscript exactly — these come from the full 10,000-parameter
  RACIPE ensemble (a pipeline stage independent of the 100→1000 stochastic
  subset-selection change), so they didn't need updating.
