# figure5_supp2.r
#
# PURPOSE: Figure 5 supplement 2 (Fig5S2) -- continues from Fig5S1 with the
# rest of old Fig5S1's panels (formerly labeled E-I, renumbered A-G so this
# supplement's own panel labels start at A), plus two new resampling-effect
# counterparts to the two stochastic-effect violins (D, E):
#   Row 1 -- Panel A/B: single/double/all-high MRT vs noise,
#            TT4/TT4SA/TS4, faceted by network x stability class --
#            Additive (A) then Multiplicative (B)
#   Row 2 -- Panel C: state-wise MRT for TS4's double-high states
#            (Multiplicative); Panel D: stochastic-effect violin,
#            single-high, TSSA/TTSA/TT4/TS4, 2x2
#   Row 3 -- Panel E: stochastic-effect violin, all-high + double-high,
#            TS/TT/TT4; Panel F: same networks/state-class as D but the
#            resampling effect instead of the stochastic effect, 2x2;
#            Panel G: same as E but the resampling effect
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
figDir     <- "figures/individual"
finalDir   <- "figures/final"
noiseType  <- "Multiplicative"
dt         <- 0.01

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

# ════════════════════════════════════════════════════════════════════════════
# PANEL A / B: single/double/all-high MRT vs noise, TT4/TT4SA/TS4,
# faceted by network x stability class -- Additive (A) then Multiplicative (B)
# ════════════════════════════════════════════════════════════════════════════
pA <- plot_stateclass_by_stability(
	nets = c("TT4", "TT4SA", "TS4"), noiseType = "Additive", resultsFolder = resultsFolder,
	dt = dt,
	outFile = file.path(figDir, "fig5S2_panelA_stateclass_additive.jpg")
)

pB <- plot_stateclass_by_stability(
	nets = c("TT4", "TT4SA", "TS4"), noiseType = "Multiplicative", resultsFolder = resultsFolder,
	dt = dt,
	outFile = file.path(figDir, "fig5S2_panelB_stateclass_multiplicative.jpg")
)

# ════════════════════════════════════════════════════════════════════════════
# PANEL C: state-wise MRT for TS4's double-high states (Multiplicative) --
# every double-high state gets its own facet (topStates = NULL, i.e. no
# top-N truncation) to show they stay the same set across noise levels
# rather than being dominated by only a couple of them.
# ════════════════════════════════════════════════════════════════════════════
pC <- plot_mrt_trajectories(
	net = "TS4", noiseType = noiseType, resultsFolder = resultsFolder, dt = dt,
	facet_var = "State", stateClasses = "double-high", topStates = NULL, ncol = 3,
	sampleN = 100,
	outFile = file.path(figDir, "fig5S2_panelC_TS4_doubleHigh_statewise.jpg"),
	panelWidth = 4, panelHeight = 3.5
)

# ════════════════════════════════════════════════════════════════════════════
# PANEL D / F: stochastic- and resampling-effect violins, single-high,
# TSSA/TTSA/TT4/TS4, networks arranged 2x2
# ════════════════════════════════════════════════════════════════════════════
pD <- plot_class_diff_violin(
	nets = c("TSSA", "TTSA", "TT4", "TS4"), noiseType = noiseType, dt = dt,
	dataFolder = dataFolder, resultsFolder = resultsFolder,
	compare = "stochastic", stateClasses = "single-high"
) + facet_wrap(vars(Network), ncol = 2, scales = "free_y")
ggsave(file.path(figDir, "fig5S2_panelD_stochDiff_singleHigh.jpg"), pD, width = 8, height = 8)

pDresample <- plot_class_diff_violin(
	nets = c("TSSA", "TTSA", "TT4", "TS4"), noiseType = noiseType, dt = dt,
	dataFolder = dataFolder, resultsFolder = resultsFolder,
	compare = "resample", stateClasses = "single-high"
) + facet_wrap(vars(Network), ncol = 2, scales = "free_y")
ggsave(file.path(figDir, "fig5S2_panelDprime_resampleDiff_singleHigh.jpg"), pDresample, width = 8, height = 8)

# ════════════════════════════════════════════════════════════════════════════
# PANEL E / G: same as D / F, for all-high and double-high states,
# TS/TT/TT4 (2 state classes x 3 networks -- already a compact 2-row grid)
# ════════════════════════════════════════════════════════════════════════════
pE <- plot_class_diff_violin(
	nets = c("TS", "TT", "TT4"), noiseType = noiseType, dt = dt,
	dataFolder = dataFolder, resultsFolder = resultsFolder,
	compare = "stochastic", stateClasses = c("all-high", "double-high"),
	outFile = file.path(figDir, "fig5S2_panelE_stochDiff_allDoubleHigh.jpg")
)

pEresample <- plot_class_diff_violin(
	nets = c("TS", "TT", "TT4"), noiseType = noiseType, dt = dt,
	dataFolder = dataFolder, resultsFolder = resultsFolder,
	compare = "resample", stateClasses = c("all-high", "double-high"),
	outFile = file.path(figDir, "fig5S2_panelEprime_resampleDiff_allDoubleHigh.jpg")
)

# ════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ════════════════════════════════════════════════════════════════════════════
row1 <- plot_grid(pA, pB, nrow = 1, labels = c("A", "B"), label_size = LS)
row2 <- plot_grid(pC, pD + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
                   nrow = 1, labels = c("C", "D"), label_size = LS)
row3 <- plot_grid(pE + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
                   pDresample + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
                   pEresample + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
                   nrow = 1, labels = c("E", "F", "G"), label_size = LS)

fig5S2 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1, 1.4, 1))

ggsave(file.path(finalDir, "Fig5S2.jpg"), fig5S2,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 22,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 26,
       limitsize = FALSE)

message("Figure 5 supplement 2 complete. Panels + combined saved to: ", figDir)
