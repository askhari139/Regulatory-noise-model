# figure5_supp.r
#
# PURPOSE: Figure 5 supplement 1 (Fig5S1).
#   Row 1 -- Panel A: same trajectory plot as Fig5 Panel A, but for
#            TS / TSSA / TT
#   Row 2 -- Panel B (lambda distributions) + Panel C (TT4/TT4SA/TS4
#            network topologies)
#   Row 3 -- Panel D: corresponding deterministic state frequency for
#            Panel C's networks
#
# Old Fig5S1 Panels E-I (stateclass-by-stability, TS4 double-high
# statewise, and the two stochastic-effect violins) moved to
# Figure5_supp2.r, alongside two new resampling-effect counterparts.
#
# ============================================================

library(funcsKishore)
library(tidyverse)
library(ggforce)
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
# PANEL A: same as Fig5 Panel A, repeated for TS / TSSA / TT
# ════════════════════════════════════════════════════════════════════════════
# fill_mrt() labels a 2-node network's 2-high state "all-high" (Expression
# == n_nodes is checked before Expression == 2), never "double-high" -- so
# for TS/TSSA specifically, filtering on stateClasses = c("single-high",
# "double-high") as-is silently drops that facet entirely (no error, just
# a missing panel). Pull "all-high" instead for those two networks, then
# relabel its facet strip to "double-high" so it reads consistently with
# TT's genuine double-high facet -- same state (all nodes high, which for
# a 2-node network the two things are one and the same), just named per
# the network's own node count downstream vs. this panel's shared vocabulary.
pA_nets <- c("TS", "TSSA", "TT")
pA_list <- map(pA_nets, function(net) {
	is2node <- net %in% c("TS", "TSSA")
	classes <- if (is2node) c("single-high", "all-high") else c("single-high", "double-high")
	p <- plot_mrt_trajectories(
		net = net, noiseType = noise_type_for(net, noiseType), resultsFolder = resultsFolder, dt = dt,
		facet_var = "StateClass",
		stateClasses = classes,
		topStates = NULL, ncol = 1,
		sampleN = 100,
		panelWidth = 6, panelHeight = 4
	)
	if (is2node) {
		p <- p + facet_wrap(vars(StateClass), ncol = 1, scales = "free_y",
							 labeller = as_labeller(c("single-high" = "single-high",
													   "all-high"   = "double-high")))
	}
	p
})
pA <- plot_grid(plotlist = pA_list, nrow = 1,
                 labels = c("A (i)", "(ii)", "(iii)"), label_size = LS)
ggsave(file.path(figDir, "fig5S1_panelA_traj.jpg"), pA, width = 18, height = 4.5)

# ════════════════════════════════════════════════════════════════════════════
# PANEL B: lambda distributions (Multiplicative), like Fig4 Panel A -- one
# faceted plot (shared x axis) instead of two plots combined side by side,
# same as the Fig4 Panel A fix.
# ════════════════════════════════════════════════════════════════════════════
pB <- plot_lambda_density(
	noiseType = noiseType, dataFolder = dataFolder,
	noise_levels = c(0.001, 0.01, 0.1, 1.0), dt = dt,
	lambda_type = c("Inh", "Act")
)
ggsave(file.path(figDir, "fig5S1_panelB_lambda.jpg"), pB, width = 5, height = 8)

# ════════════════════════════════════════════════════════════════════════════
# PANEL C: TT4 / TT4SA / TS4 network topologies
# ════════════════════════════════════════════════════════════════════════════
pC_i   <- plot_network_diagram(file.path(dataFolder, "TT4.topo"))
pC_ii  <- plot_network_diagram(file.path(dataFolder, "TT4SA.topo"))
pC_iii <- plot_network_diagram(file.path(dataFolder, "TS4.topo"))

pC <- plot_grid(pC_i, pC_ii, pC_iii, nrow = 1,
                 labels = c("C (i)", "(ii)", "(iii)"), label_size = LS)
ggsave(file.path(figDir, "fig5S1_panelC_networks.jpg"), pC, width = 12, height = 4.2, bg = "white")

# ════════════════════════════════════════════════════════════════════════════
# PANEL D: corresponding top RACIPE attractor-combination frequency for
# Panel C's networks (per-ParamID "Original" combination, not per-state --
# see plot_original_frequency() in figure_common.r)
# ════════════════════════════════════════════════════════════════════════════
pD <- plot_original_frequency(
	nets = c("TT4", "TT4SA", "TS4"), noiseType = noiseType, resultsFolder = resultsFolder,
	dataFolder = dataFolder, dt = dt, topN = 10, ncol = 3,
	outFile = file.path(figDir, "fig5S1_panelD_stateFreq.jpg")
)

# ════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ════════════════════════════════════════════════════════════════════════════
row1 <- pA
row2 <- plot_grid(pB, pC, nrow = 1, rel_widths = c(1, 1.3), labels = c("B", ""), label_size = LS)
row3 <- plot_grid(pD, nrow = 1, labels = c("D"), label_size = LS)

fig5S1 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1.7, 1.6, 1))

ggsave(file.path(finalDir, "Fig5S1.jpg"), fig5S1,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 18,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 15)

message("Figure 5 supplement 1 complete. Panels + combined saved to: ", figDir)
