# Production run summary — `2026-10-csbj-v2` (2026-10-24)

**Code:** commit `d17f3a0` (branch `regnoise-pilot`). **Config hash:** see any raw file's metadata.
**Scope:** TS, TSSA, TT, TTSA, TS4, TT4, TT4SA × Additive, Multiplicative, Stationary; up to 1,000
parameter sets per stability class (19,296 per mode); all controls C0–C4. One 128-core ctbp node per
network × mode; longest job 4 h 20 min (TT4SA Stationary).
**Data:** raw `/scratch/a.hari/regnoise/results/2026-10-csbj-v2/raw/` (74 GB);
collected `/scratch/a.hari/regnoise/results/2026-10-csbj-v2/collected/` (21 GB; the summary tables are
also copied locally, the per-replicate and fold-change tables are on scratch only).
**Figures:** `CSBJ_review/production_figures/`. **Full per-combination report:** `production_report_full_2026-10-24.md`.

A first attempt (`2026-10-csbj`, commit `fca1635`) was cancelled at ~92% after solver divergences in
self-activating networks (expression stepped below 0 with long λ holds). The fix (`isoutofdomain`,
rejecting negative-expression steps) was applied and the whole run repeated from one commit.

## Checks

| Check | Result |
|---|---|
| Validation (every condition cell at the expected replicate count, occupancies sum to 1, λ within bounds) | 21 / 21 PASS |
| Units | 57,888 / 57,888 |
| Failed ODE solves | 0 |
| C2 frozen draws: 95th-pct per-set error with 10 draws (old setting) | 0.06–0.22 |
| C2 with 100 draws | 0.005–0.07 |
| C2 bootstrap SE at 1,000 draws (worst cell) | ≤ 0.020 |
| Block drift in the analysis window, σ ≥ 0.005 | ≤ 0.011 (stationary) |
| Block drift at σ = 0.001, Additive | 0.02–0.08 (not stationary within t = 1000; report as transient) |

## All-high occupancy, Additive noise (mean over parameter sets, Δt = 0.01)

| Network | σ | C0 original | C4 mean λ | C2 frozen | C3 fast (T = 0.01) | C1 stochastic | fast − mean λ |
|---|---|---|---|---|---|---|---|
| TS | 0.01 / 0.1 / 1 | 0.066 | 0.78 / 0.78 / 0.78 | 0.62 / 0.59 / 0.40 | 0.78 / 0.78 / 0.78 | 0.65 / 0.76 / 0.78 | 0.000 / −0.001 / −0.001 |
| TSSA | 0.01 / 0.1 / 1 | 0.100 | 0.61 / 0.61 / 0.61 | 0.53 / 0.51 / 0.42 | 0.64 / 0.71 / 0.89 | 0.46 / 0.69 / 0.89 | +0.03 / +0.10 / +0.28 |
| TT | 0.01 / 0.1 / 1 | 0.007 | 0.81 / 0.81 / 0.81 | 0.51 / 0.42 / 0.14 | 0.81 / 0.81 / 0.80 | 0.58 / 0.78 / 0.80 | −0.001 / −0.001 / −0.001 |
| TTSA | 0.01 / 0.1 / 1 | 0.008 | 0.73 / 0.73 / 0.73 | 0.49 / 0.42 / 0.22 | 0.75 / 0.78 / 0.92 | 0.46 / 0.77 / 0.92 | +0.01 / +0.05 / +0.19 |
| TS4 | 0.01 / 0.1 / 1 | 0.004 | 0.82 / 0.82 / 0.81 | 0.54 / 0.42 / 0.11 | 0.82 / 0.82 / 0.81 | 0.62 / 0.80 / 0.81 | 0.000 / 0.000 / +0.001 |
| TT4 | 0.01 / 0.1 / 1 | 0.001 | 0.88 / 0.88 / 0.88 | 0.48 / 0.34 / 0.07 | 0.88 / 0.88 / 0.88 | 0.60 / 0.85 / 0.88 | −0.001 / +0.001 / 0.000 |
| TT4SA | 0.01 / 0.1 / 1 | 0.000 | 0.80 / 0.80 / 0.80 | 0.43 / 0.32 / 0.12 | 0.80 / 0.85 / 0.96 | 0.45 / 0.82 / 0.96 | +0.002 / +0.05 / +0.16 |

## Main findings

1. **Fast fluctuations of inhibitory interactions act exactly like their mean.** In the four
   inhibition-only networks, fast hold-resampling equals the fixed-mean-λ control to within 0.001 at every
   σ. The inhibitory production term H⁻ + λH⁺ is linear in λ, so time-averaging λ is the same as
   averaging the production. The large all-high increase under additive noise in these networks is the
   shift of mean λ (≈ 0.05 → 0.5, i.e. ~20-fold → ≤ 2-fold repression).
2. **Fluctuations of activating interactions have an effect beyond their mean.** In the three
   self-activating networks, fast fluctuations give more all-high than mean λ, by up to +0.28 (TSSA,
   σ = 1), growing with noise amplitude. The activating production term H⁻/λ + H⁺ is convex in λ (the
   1/λ term), so its time average exceeds its value at mean λ (Jensen). Excursions toward λ ≈ 1 (no
   regulation) give full production regardless of the activator, raising the floor of every node. This
   is the paper's "unconditional derepression vs conditional regulation" asymmetry, made precise: beyond
   the mean, it acts only through the nonlinearity of activating edges.
3. **Frozen heterogeneity and temporal fluctuation differ strongly.** With the same λ distribution,
   frozen draws (C2) give far less all-high than time-varying λ (C1, C3), and the hold-time ladder
   interpolates between them with a crossover at hold times ≈ 1–10, the expression relaxation time.
4. **The stochastic walk (C1) sits on that ladder according to its own speed.** At σ ≥ 0.1 it behaves
   like fast fluctuation; at σ = 0.01 it is slow (C1 between the T = 10 and T = 100 rungs), and in the
   self-activating networks it falls below the frozen level (e.g. TSSA 0.46 vs 0.53).
