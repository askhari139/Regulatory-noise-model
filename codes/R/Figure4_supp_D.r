# figure4_supp_D.r
#
# PURPOSE: Supplement to paper Figure 4 Panel D -- now Fig4S2 (previously
# Fig4S3; renumbered when the old Fig4S1 was retired and Figure4_supp_BC.r
# became Fig4S1). Same paired Wilcoxon test
# networks-compared comparison as the main figure's all-high-state panel
# (fixed there at sigma = 0.1), generalized across state class (all-high /
# single-high / double-high) and repeated at all four noise levels (0.001,
# 0.01, 0.1, 1.0). One row per state class, one column per noise level.
# Each cell itself facets into Monostable/Multistable, same as the main
# figure panel. 2-node networks (TS, TSSA) have no double-high states of
# their own (their 2-high state is already "all-high"), so those networks
# come back missing from that row's x axis -- expected, not an error.
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
figDir       <- "figures/individual"
finalDir     <- "figures/final"
noiseType    <- "Additive"
dt           <- 0.01
nets         <- c("TS", "TSSA", "TT", "TTSA")
sigmas       <- c(0.001, 0.01, 0.1, 1.0)
stateClasses <- c("all-high", "single-high", "double-high")

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

# plot_det_stoch_network_comparison() comes from figure_common.r (shared with
# Figure4.r and Figure4_supp.r, which fix it at state_class = NULL, i.e.
# all-high).

# ── ONE ROW PER STATE CLASS, ONE COLUMN PER NOISE LEVEL ─────────────────────
rows <- map(stateClasses, function(cls) {
	map(sigmas, function(sigma) {
		plot_det_stoch_network_comparison(
			nets, noiseType, dt, dataFolder, resultsFolder,
			state_class = cls, sigma = sigma,
			outFile = file.path(figDir, paste0("fig4supp_D_", cls, "_sigma", sigma, ".jpg"))
		) + labs(title = paste0(cls, ", σ = ", sigma))
	})
})

# ── ASSEMBLE COMBINED FIGURE ──────────────────────────────────────────────────
# Rows are labeled A, B, C... (not the main figure's "D"); the state class
# itself is already on each subplot's own title, so the row letter doesn't
# need to repeat it.
row_letters <- LETTERS[seq_along(stateClasses)]
grid_rows <- map2(rows, row_letters, function(row_plots, lab) {
	labs_row <- c(lab, rep("", length(row_plots) - 1))
	plot_grid(plotlist = row_plots, nrow = 1, labels = labs_row, label_size = LS)
})

fig4supp_D <- plot_grid(plotlist = grid_rows, ncol = 1)

ggsave(file.path(finalDir, "Fig4S2.jpg"), fig4supp_D,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 30,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 20,
       limitsize = FALSE)

message("Figure 4 supplement (Panel D, by state class) complete. Individual panels and combined figure saved to: ", figDir)
