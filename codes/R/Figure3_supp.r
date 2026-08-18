# figure3_supp.r
#
# PURPOSE: Supplement to paper Figure 3.
#   Row B: mrt_boxplot_stability() (Figure3.r Panel B), for TSSA, TT, TTSA --
#          each network facets automatically into however many of
#          monostable/bistable/tristable it actually has.
#   Row C: plot_reach_scatter() (Figure3.r Panel C), for:
#            - TS at sigma = 0.001, 0.01, 1.0 (0.1 is already the main figure)
#            - TSSA at sigma = 0.001, 0.01, 0.1, 1.0
#          TT/TTSA are skipped for this panel: it's an X-vs-Y scatter of two
#          nodes' margins, which doesn't generalize to their 3 nodes.
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
figDir      <- "figures/individual"
finalDir    <- "figures/final"
noiseType   <- "Additive"
dt          <- 0.01

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

# mrt_boxplot_stability() and plot_reach_scatter() come from figure_common.r
# (shared with Figure3.r, which uses them for TS at the main figure's fixed
# sigma = 0.1).

# ── ROW B: mrt_boxplot_stability(), one panel per network ───────────────────
rowB_nets <- c("TSSA", "TT", "TTSA")
rowB <- map(rowB_nets, function(net) {
	mrt_boxplot_stability(
		net, noiseType, resultsFolder, dt = dt,
		outFile = file.path(figDir, paste0("fig3supp_panelB_", net, ".jpg"))
	) + labs(title = net)
})

# ── ROW C: plot_reach_scatter(), TS at 3 extra sigmas + TSSA at 4 sigmas ────
rowC_TS <- map(c(0.001, 0.01, 1.0), function(sigma) {
	plot_reach_scatter(
		"TS", noiseType, dt, dataFolder, resultsFolder,
		state = "(1, 1)", sigma = sigma, param_type = "monostable",
		outFile = file.path(figDir, paste0("fig3supp_panelC_TS_sigma", sigma, ".jpg"))
	) + labs(title = paste0("TS, σ = ", sigma))
})

rowC_TSSA <- map(c(0.001, 0.01, 0.1, 1.0), function(sigma) {
	plot_reach_scatter(
		"TSSA", noiseType, dt, dataFolder, resultsFolder,
		state = "(1, 1)", sigma = sigma, param_type = "monostable",
		outFile = file.path(figDir, paste0("fig3supp_panelC_TSSA_sigma", sigma, ".jpg"))
	) + labs(title = paste0("TSSA, σ = ", sigma))
})

# ── ASSEMBLE ──────────────────────────────────────────────────────────────────
# Supplementary panels are labeled A, B, C... independently of whatever
# letter the analogous panel has in the main figure (main Figure3.r's B and
# C, here). Row 1: A (mrt_boxplot_stability, one net per row) and B (TS
# reach-scatter, one sigma per row) side by side, each stacked vertically.
# Row 2: C (TSSA reach-scatter) unchanged, side by side.
rowB_grid     <- plot_grid(plotlist = rowB, ncol = 1, labels = c("B (i)", "(ii)", "(iii)"), label_size = LS)
rowC_TS_grid  <- plot_grid(plotlist = rowC_TS, ncol = 1, labels = c("A", "", ""), label_size = LS)
rowC_TSSA_grid <- plot_grid(plotlist = rowC_TSSA, nrow = 1, labels = c("C", "", "", ""), label_size = LS)

row1 <- plot_grid(rowC_TS_grid, rowB_grid, nrow = 1, rel_widths = c(1, 2.7))

fig3supp <- row1

ggsave(file.path(finalDir, "Fig3S1.jpg"), fig3supp,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 13,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 15)

message("Figure 3 supplement complete. Individual panels and combined figure saved to: ", figDir)
