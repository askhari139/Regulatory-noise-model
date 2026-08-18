# figure6.r
#
# PURPOSE: Figure 6 (placeholder -- not yet assigned a final position in
# the manuscript). Built for TT4 and TT4SA (run twice, once per network),
# each producing its own Fig6_<net>.jpg. Three sections:
#   Section 1 (Panel A): how stable are the observed states under noise --
#     distribution of switching-event counts by noise level. Recorded
#     window already excludes a burn-in (see plot_switching_stability()'s
#     comment in figure_common.r for exactly how).
#   Section 2 (Panels B/C): stability patterns among double-high ParamIDs
#     -- frequency of the 3 deterministic stability classes (Panel B), then
#     each class's own multistability breakdown (how many fates it
#     accesses under noise, collapsed to Monostable/Bistable/Tristable(/
#     Tetrastable) and shown as a fraction-of-parameter-sets stack per
#     noise level), all 3 classes side by side (Panel C).
#   Section 3 (Panels D/E): what the transitions actually look like -- one
#     example monostable double-high attractor's own per-parameter MRT
#     distribution at noise = 0.01 (Panel D), then the pooled
#     monostable-double-high population as a single noise-level -> fate
#     alluvial, colored by the same multistability classification as
#     Panel C's "Monostable (double-high)" facet (Panel E).
#
# See Figure6_supp.r for the circular flow diagrams and the bistable
# mirror/non-mirror counterparts of Panel E.
#
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpubr)
library(ggalluvial)
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
sigmasABC  <- c(0.001, 0.01, 0.05, 0.1)
sigmasBC   <- c(0.01, 0.05, 0.1)  # drop 0.001 -- too few accessible fates to be informative

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

build_fig6 <- function(net) {

	tag <- function(name) paste0("fig6_", net, "_", name)

	accessData <- compute_access_breakdown(
		net = net, noiseType = noiseType, resultsFolder = resultsFolder,
		dataFolder = dataFolder, dt = dt, sigmas = sigmasABC
	)

	# ════════════════════════════════════════════════════════════════════
	# SECTION 1 -- PANEL A: switching-event stability, Additive vs.
	# Multiplicative, restricted to noise levels <= 0.1 (the far tail --
	# 0.5, 1 -- dominates the shared y range without adding much beyond
	# "switching keeps increasing")
	# ════════════════════════════════════════════════════════════════════
	pA <- plot_switching_stability_combined(net, c("Additive", "Multiplicative"), resultsFolder, dt = dt,
		max_sigma = 0.1, outFile = file.path(figDir, tag("panelA_switching.jpg")))

	# ════════════════════════════════════════════════════════════════════
	# SECTION 2 -- PANEL B/C: stability-class frequency + multistability
	# class summary, all 3 stability classes side by side in Panel C
	# ════════════════════════════════════════════════════════════════════
	freqData <- compute_stability_class_freq(net = net, noiseType = noiseType,
		resultsFolder = resultsFolder, dataFolder = dataFolder, dt = dt)
	pB <- plot_stability_class_freq(freqData,
		outFile = file.path(figDir, tag("panelB_stabilityFreq.jpg")))

	classDataAll <- compute_nfates_class_all(accessData)
	pC <- plot_nfates_class_summary_combined(classDataAll, sigmas = sigmasABC,
		outFile = file.path(figDir, tag("panelC_nfatesClassSummary.jpg")))

	# ════════════════════════════════════════════════════════════════════
	# SECTION 3 -- PANEL D/E: one example attractor's own MRT distribution,
	# then the pooled monostable-double-high population as a single
	# noise-level -> fate alluvial
	# ════════════════════════════════════════════════════════════════════
	topLabel <- accessData %>%
		filter(StabilityClass == "Monostable (double-high)") %>%
		group_by(OriginalLabel) %>% summarise(Tot = sum(MRT), .groups = "drop") %>%
		slice_max(Tot, n = 1, with_ties = FALSE) %>% pull(OriginalLabel)

	pD <- plot_example_mrt_violin(accessData, originalLabel = topLabel, sigma = 0.01,
		outFile = file.path(figDir, tag("panelD_exampleViolin.jpg")))

	roleMono  <- classify_double_high_role(accessData)
	classMono <- compute_nfates_class(roleMono, fateRoles = c("SH1_1", "SH1_2", "SH3"))

	pE <- plot_double_high_alluvial(roleMono, classMono, sigmas = sigmasABC,
		outFile = file.path(figDir, tag("panelE_noiseAlluvial.jpg")))

	# ════════════════════════════════════════════════════════════════════
	# ASSEMBLE
	# ════════════════════════════════════════════════════════════════════
	row1 <- plot_grid(pA, nrow = 1, labels = c("A"), label_size = LS)
	row2 <- plot_grid(pB + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
	                   pC + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
	                   nrow = 1, rel_widths = c(0.5, 1), labels = c("B", "C"), label_size = LS)
	row3 <- plot_grid(pD + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
	                   pE + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
	                   nrow = 1, rel_widths = c(0.7, 1), labels = c("D", "E"), label_size = LS)

	fig6 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(0.9, 0.9, 1.1))

	ggsave(file.path(finalDir, paste0("Fig6_", net, ".jpg")), fig6,
	       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 18,
	       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 22,
	       limitsize = FALSE)

	message("Figure 6 (", net, ") complete. Panels + combined saved to: ", figDir)
}

for (net in c("TT4", "TT4SA")) build_fig6(net)
