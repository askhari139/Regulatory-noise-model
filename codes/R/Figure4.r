# figure4.r
#
# PURPOSE: Paper Figure 4 — lambda distribution shift and its consequences
#   Panel A: Lambda distributions (Inh/Act), faceted into one plot (shared
#            x axis, one strip per type) instead of two stacked plots
#   Panel B: Resampling effect alone (Original det vs Det-resampled), all-high
#            state, violin of the per-ParamID difference across every noise
#            level, networks arranged 2x2 (the single-sigma scatter version
#            now lives in Fig4S1)
#   Panel C (i): Added stochastic effect (Det-resampled vs Stochastic),
#            all-high state -- same violin, networks arranged 2x2
#   Panel C (ii): TS/TT-only scatter at sigma = 0.1, stacked 2x1, showing the
#            underlying sigmoidal MRT_det -> MRT relationship the violin's
#            signed difference can't show on its own (the full 4-network
#            scatter also lives in Fig4S1)
#   Panel D: All networks compared at sigma = 0.001/0.01/0.1, all three MRT
#            series, faceted by stability class (rows) x noise level (cols),
#            with paired Wilcoxon test significance brackets (resampled vs.
#            stochastic, resampled vs. original/point-lambda)
#
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpubr)
library(ggsignif)
.scriptDir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(fa) > 0) return(dirname(normalizePath(fa[1])))
    for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
    getwd()
})()
source(file.path(.scriptDir, "figure_common.r"))

# ── CONFIG ────────────────────────────────────────────────────────────────────
figDir        <- "figures/individual"
finalDir      <- "figures/final"
noiseType     <- "Additive"
dt            <- 0.01
nets          <- c("TS", "TSSA", "TT", "TTSA")
sigmaFixed    <- 0.1

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

# get_det_stoch_comparison_data(), plot_paramwise_comparison(), and
# plot_det_stoch_network_comparison() all come from figure_common.r (shared
# with Figure4_supp.r, which reuses them at additional noise levels).

# ── PANEL A: Lambda distributions, faceted by type (Inh/Act) instead of two
# separate stacked plots -- one shared x axis (text/title drawn once, on the
# bottom facet only) instead of duplicating it per type ─────────────────────
pA <- plot_lambda_density(
	noiseType = noiseType, dataFolder = dataFolder,
	noise_levels = c(0.001, 0.01, 0.1, 1.0), dt = dt,
	lambda_type = c("Inh", "Act")
) +
	theme(strip.text = element_text(size = rel(1.2)))
ggsave(file.path(figDir, "fig4_panelA_lambda.jpg"), pA, width = 5, height = 8)

# ── PANEL B: Resampling effect, violin across all noise levels, networks
# arranged 2x2 (all-high state only -- the per-state-class breakdown across
# all-high / single-high / double-high already lives in Fig4S2). The facet
# override replaces plot_class_diff_violin()'s default facet_grid(rows =
# StateClass, cols = Network) -- a no-op row split here since stateClasses
# is a single value -- with a plain 2-column network wrap. ─────────────────
p2 <- plot_class_diff_violin(
	nets, noiseType, dt, dataFolder, resultsFolder,
	compare = "resample", stateClasses = "all-high"
) + facet_wrap(vars(Network), ncol = 2, scales = "fixed") +
	theme(strip.text = element_text(size = rel(1.2)))
ggsave(file.path(figDir, "fig4_panelB_resample_violin.jpg"), p2, width = 8, height = 8)

# ── PANEL C (i)/(ii): Stochastic effect -- same violin (2x2), plus a TS/TT
# scatter at sigma = 0.1, stacked 2x1, to show the sigmoidal MRT_det -> MRT
# relationship directly (the violin's signed difference collapses that
# shape away) ────────────────────────────────────────────────────────────────
p3_violin <- plot_class_diff_violin(
	nets, noiseType, dt, dataFolder, resultsFolder,
	compare = "stochastic", stateClasses = "all-high"
) + facet_wrap(vars(Network), ncol = 2, scales = "fixed")
ggsave(file.path(figDir, "fig4_panelC_i_stochastic_violin.jpg"), p3_violin, width = 8, height = 8) +
	theme(strip.text = element_text(size = rel(1.2)))

p3_scatter <- plot_paramwise_comparison(
	nets = c("TS", "TT"), noiseType, dt, dataFolder, resultsFolder,
	compare = "stochastic", sigma = sigmaFixed
) + facet_wrap(vars(Network), ncol = 1) +
	theme(strip.text = element_text(size = rel(1.2)))
ggsave(file.path(figDir, "fig4_panelC_ii_stochastic_scatter_TS_TT.jpg"), p3_scatter, width = 4, height = 8)

p3 <- plot_grid(p3_violin, p3_scatter, nrow = 1, rel_widths = c(1.3, 1),
                 labels = c("C (i)", "(ii)"), label_size = LS)

# ── PANEL D: All networks compared at sigma = 0.001/0.01/0.1, faceted by
# stability class (rows) x noise level (cols), with paired Wilcoxon test
# significance brackets ─────────────────────────────────────────────────────
# "Resampled" = MRT_det (deterministic, resampled lambda). Two comparisons,
# paired by ParamID within each (Network, StabilityClass, NoiseLevel):
# resampled vs. stochastic (MRT_det vs. MRT), and resampled vs. original
# point-lambda (MRT_det vs. MRT_orig) -- the same two effects isolated
# separately in Panels B/C, now tested for significance per network, per
# noise level.
p4 <- plot_det_stoch_network_comparison(
	nets, noiseType, dt, dataFolder, resultsFolder,
	sigma = c(0.001, 0.01, 0.1),
	outFile = file.path(figDir, "fig4_panelD_allNets_byStabilityAndNoise.jpg"),
	width = 16, height = 8
) +
	theme(strip.text = element_text(size = rel(1.2)))

# ── ASSEMBLE COMBINED FIGURE ──────────────────────────────────────────────────
# p3 already carries its own internal "(i)/(ii)" labels, so its outer
# plot_grid slot only needs a plain "" placeholder (not a second C) to avoid
# a duplicate tag.
row1 <- plot_grid(pA, p2, nrow = 1, rel_widths = c(1, 1.3),
                   labels = c("A", "B"), label_size = LS)
row2 <- plot_grid(p3, nrow = 1)
row3 <- plot_grid(p4, nrow = 1, labels = c("D"), label_size = LS)

fig4 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1.3, 1.3, 1.1))

ggsave(file.path(finalDir, "Fig4.jpg"), fig4,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 16,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 20)

message("Figure 4 complete. Individual panels and combined figure saved to: ", figDir)