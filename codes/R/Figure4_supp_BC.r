# figure4_supp_BC.r
#
# PURPOSE: Supplement to paper Figure 4 -- now Fig4S1 (previously Fig4S2;
# the old Fig4S1 -- full-network paramwise scatter at every noise level for
# BOTH resample and stochastic effects, plus network-comparison-by-
# stability at every noise level -- was retired: Panel C below (the
# stochastic-effect scatter) is the one piece of it that wasn't redundant
# with the main figure or with Panels A/B here, so it moved in; the rest
# was dropped rather than carried forward).
#   Panel A: resampling effect, MRT_det - MRT_orig, violin across every
#            noise level at once, faceted by state class x network
#   Panel B: stochastic effect, MRT - MRT_det, same violin design
#   Panel C: stochastic effect paramwise scatter (MRT_det vs. MRT), one
#            column per noise level (0.001, 0.01, 0.1, 1.0), faceted by
#            network -- the full 4-network version of Fig4's Panel C(ii)
#            TS/TT-only scatter, across every noise level rather than the
#            main figure's single fixed sigma = 0.1
# 2-node networks (TS, TSSA) have no double-high states of their own (their
# 2-high state is already "all-high"), so that row comes back empty for
# those two networks in Panels A/B -- expected, not an error.
#
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpubr)
.scriptDir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(fa) > 0) return(dirname(normalizePath(fa[1])))
    for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
    getwd()
})()
source(file.path(.scriptDir, "figure_common.r"))

# ── CONFIG ────────────────────────────────────────────────────────────────────
figDir      <- "figures/individual"
finalDir    <- "figures/final"
noiseType   <- "Additive"
dt          <- 0.01
nets        <- c("TS", "TSSA", "TT", "TTSA")
suppSigmas  <- c(0.001, 0.01, 0.1, 1.0)

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

# plot_class_diff_violin() and plot_paramwise_comparison() come from
# figure_common.r.

# ── PANEL A: resampling effect (MRT_det - MRT_orig) ──────────────────────────
pA <- plot_class_diff_violin(
	nets, noiseType, dt, dataFolder, resultsFolder,
	compare = "resample",
	outFile = file.path(figDir, "fig4S1_panelA_resample_violin.jpg")
)

# ── PANEL B: stochastic effect (MRT - MRT_det) ───────────────────────────────
pB <- plot_class_diff_violin(
	nets, noiseType, dt, dataFolder, resultsFolder,
	compare = "stochastic",
	outFile = file.path(figDir, "fig4S1_panelB_stochastic_violin.jpg")
)

# ── PANEL C: stochastic effect paramwise scatter, one column per noise level ─
rowC <- map(suppSigmas, function(sigma) {
	plot_paramwise_comparison(
		nets, noiseType, dt, dataFolder, resultsFolder,
		compare = "stochastic", sigma = sigma,
		outFile = file.path(figDir, paste0("fig4S1_panelC_stochastic_sigma", sigma, ".jpg"))
	)
})
pC <- plot_grid(plotlist = rowC, nrow = 1, labels = c("C", "", "", ""), label_size = LS)

# ── ASSEMBLE COMBINED FIGURE ──────────────────────────────────────────────────
# Supplementary panels are labeled A, B, C... independently of the main
# figure's own B/C/D letters for these same comparisons.
fig4supp_BC <- plot_grid(pA, pB, pC, ncol = 1, labels = c("A", "B", ""), label_size = LS,
                          rel_heights = c(1, 1, 1))

ggsave(file.path(finalDir, "Fig4S1.jpg"), fig4supp_BC,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 20,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 16)

message("Figure 4 supplement (now Fig4S1: Panels A/B violin + C stochastic scatter) complete. Individual panels and combined figure saved to: ", figDir)
