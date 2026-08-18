# figure1.r
#
# PURPOSE: Paper Figure 1 -- network schematic + RACIPE background.
#   Panel A: network topology diagrams (i)-(iv) -- TS, TSSA, TT, TTSA --
#            drawn directly from the .topo files, tikz-style, arranged 2x2.
#   Panel B (i): RACIPE_schematic.png (existing illustrative figure).
#   Panel B (ii): TS RACIPE steady-state solutions, A vs B on log scale,
#            threshold lines + discrete-state quadrant labels.
#   Panel C: deterministic (NoiseLevel == 0) attractor-combination frequency
#            for TS/TSSA/TT/TTSA, one subpanel per network.
#   Panel D: simulated lambda(t) trajectories under the three noise rules,
#            for an inhibition-range edge (lambda in [0, 1]).
#
# The general ODE/shifted-Hill and noise-update-rule equation panels that
# used to be this figure's Panels B/C are dropped here -- the new B/C
# defined above replace them; they aren't referenced anywhere in this
# script anymore.
# ============================================================

library(funcsKishore)
library(tidyverse)
library(ggforce)
library(cowplot)
library(patchwork)
.scriptDir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(fa) > 0) return(dirname(normalizePath(fa[1])))
    for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
    getwd()
})()
source(file.path(.scriptDir, "figure_common.r"))

figDir    <- "figures/individual"
finalDir  <- "figures/final"
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)
dataDir   <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/data"
noiseType <- "Additive"
dt        <- 0.01

# ════════════════════════════════════════════════════════════════════════════
# PANEL A: network diagrams, arranged 2x2
# ════════════════════════════════════════════════════════════════════════════
pA_i   <- plot_network_diagram(file.path(dataDir, "TS.topo"), node_text_size = 10)
pA_ii  <- plot_network_diagram(file.path(dataDir, "TSSA.topo"),node_text_size = 10)
pA_iii <- plot_network_diagram(file.path(dataDir, "TT.topo"),node_text_size = 10)
pA_iv  <- plot_network_diagram(file.path(dataDir, "TTSA.topo"),node_text_size = 10)

pA <- plot_grid(pA_i, pA_ii, pA_iii, pA_iv, nrow = 2,
				 labels = c("A (i)", "(ii)", "(iii)", "(iv)"), label_size = LS)
ggsave(file.path(figDir, "fig1_panelA_networks.jpg"), pA, width = 8, height = 8, bg = "white")

# ════════════════════════════════════════════════════════════════════════════
# PANEL B (i): existing illustrative schematic (figures/RACIPE_schematic.png)
# ════════════════════════════════════════════════════════════════════════════
schematic_img <- png::readPNG("/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/figures/RACIPE_schematic.png")
pBi <- ggplot() +
	annotation_custom(grid::rasterGrob(schematic_img, interpolate = TRUE),
					   xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
	theme_void()

# ════════════════════════════════════════════════════════════════════════════
# PANEL B (ii): TS RACIPE steady-state solutions -- A vs B, log scale,
# threshold lines (geometric mean of each node's solutions, i.e. mean in the
# log2 space read_solutions() returns, exponentiated back), quadrants
# labeled with their corresponding discrete state (see plots_script.r for
# the original, unlabeled version of this scatter).
# ════════════════════════════════════════════════════════════════════════════
ts_sol <- read_solutions(file.path(dataDir, "TS_solution.dat")) %>% arrange(A, B)
a_thresh <- 2^mean(ts_sol$A)
b_thresh <- 2^mean(ts_sol$B)

rangeA <- range(2^ts_sol$A); rangeB <- range(2^ts_sol$B)
# Quadrant label positions: geometric midpoint of {axis edge, threshold} in
# log space, i.e. halfway between the edge and the threshold line on the
# log-scaled axis actually drawn.
geo_mid <- function(lo, hi) 10^(mean(log10(c(lo, hi))))
xLow  <- geo_mid(rangeA[1], a_thresh); xHigh <- geo_mid(a_thresh, rangeA[2])
yLow  <- geo_mid(rangeB[1], b_thresh); yHigh <- geo_mid(b_thresh, rangeB[2])

quadLabels <- tibble(
	x     = c(xLow, xHigh, xLow, xHigh),
	y     = c(yLow, yLow, yHigh, yHigh),
	label = c("(0, 0)", "(1, 0)", "(0, 1)", "(1, 1)")
)

pBii <- ggplot(ts_sol, aes(x = 2^A, y = 2^B)) +
	geom_point(alpha = 0.4, size = 0.7, color = "steelblue") +
	geom_hline(yintercept = b_thresh, linetype = "dashed", color = "red", linewidth = 0.6) +
	geom_vline(xintercept = a_thresh, linetype = "dashed", color = "red", linewidth = 0.6) +
	geom_label(data = quadLabels, aes(x = x, y = y, label = label), inherit.aes = FALSE,
			   fontface = "bold", size = 5, fill = "white", alpha = 0.8, label.size = 0) +
	scale_x_log10() + scale_y_log10() +
	theme_Publication() +
	labs(x = "Node A logscale", y = "Node B logscale")

pB <- plot_grid(pBi, pBii, ncol = 1, rel_heights = c(1.4, 1),
				 labels = c("B (i)", "(ii)"), label_size = LS)
ggsave(file.path(figDir, "fig1_panelB_schematic_TSsolutions.jpg"), pB, width = 6, height = 10, bg = "white")

# ════════════════════════════════════════════════════════════════════════════
# PANEL C: deterministic attractor-combination frequency, TS/TSSA/TT/TTSA,
# 4 subpanels side by side -- sourced directly from RACIPE's own
# steady-state search (solution.dat, via plot_racipe_original_frequency())
# rather than the noise-simulation pipeline's NoiseLevel == 0 rows
# (plot_original_frequency(), still used by Fig5S1 Panel D) -- this panel
# is specifically about what RACIPE itself finds.
# ════════════════════════════════════════════════════════════════════════════
pC <- plot_racipe_original_frequency(
	nets = c("TS", "TSSA", "TT", "TTSA"), dataFolder = dataDir, topN = 10, ncol = 4,
	outFile = file.path(figDir, "fig1_panelC_stateFreq.jpg"), panelWidth = 4, panelHeight = 4
) +
theme(strip.text = element_text(size = rel(1.2)))

# ════════════════════════════════════════════════════════════════════════════
# PANEL D: simulated lambda(t) under the three noise rules, for an
# inhibition-range edge (lambda in [0, 1]) -- 10 individual trajectories
# drawn as black dashed lines, mean +/- SD ribbon computed over 100.
# plot_lambda_trajectories()/simulate_lambda() live in figure_common.r
# (shared with Figure1_supp.r's activation-range version).
# ════════════════════════════════════════════════════════════════════════════
pD <- plot_lambda_trajectories(
	lo = 0.001, hi = 0.999, n_stats = 100, n_display = 5, display_every = 1, display_alpha = 0.7,
) +
theme(strip.text = element_text(size = rel(1.2)))

# ════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ════════════════════════════════════════════════════════════════════════════
row1 <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1.3, 1))
row2 <- plot_grid(pC, nrow = 1, labels = c("C"), label_size = LS)
row3 <- plot_grid(pD, nrow = 1, labels = c("D"), label_size = LS)

fig1 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1.6, 1, 0.9))
ggsave(file.path(finalDir, "Fig1.jpg"), fig1,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 14,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 18,
       bg = "white")

message("Figure 1 complete. Individual panels and combined figure saved to: ", figDir)
