# figure6_supp.r
#
# PURPOSE: Figure 6 supplement 1 (Fig6S1). Built for TT4 and TT4SA (run
# twice, once per network), each producing its own Fig6S1_<net>.jpg, single
# column:
#   Panel A: circular double-high -> single-high flow diagram (monostable
#     double-high ParamIDs only, A/B/C/D fates only -- Self/Other dropped),
#     at noise = 0.01, via a circlize chord diagram. (A ggraph circular
#     network layout was also tried as an alternative rendering -- dropped
#     after comparison; plot_double_high_circular_alluvial() in
#     figure_common.r still implements it if wanted again later.)
#   Panel C/D: Figure6.r Panel E's noise-level -> fate alluvial, repeated
#     for the two bistable stability classes (mirror / non-mirror) instead
#     of monostable double-high. These use the 4 raw single-high letters as
#     fate buckets (classify_bistable_role()), not the monostable case's
#     Hamming-distance SH1/SH3 split, which doesn't generalize to a
#     compound bistable Original -- see classify_bistable_role()'s comment.
#
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpubr)
library(ggalluvial)
library(circlize)
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
circleSigma <- 0.01

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

build_fig6_supp <- function(net) {

	tag <- function(name) paste0("fig6S1_", net, "_", name)

	accessData <- compute_access_breakdown(
		net = net, noiseType = noiseType, resultsFolder = resultsFolder,
		dataFolder = dataFolder, dt = dt, sigmas = sigmasABC
	)

	# ════════════════════════════════════════════════════════════════════
	# PANEL A: circular flow diagram, monostable double-high -> A/B/C/D,
	# noise = 0.01. circlize draws with base graphics, not grid graphics --
	# grid::grid.grabExpr() (which only captures grid's own display list)
	# comes back blank for it, so the standalone file is instead read back
	# as a raster and wrapped in a grob to sit in the same cowplot grid as
	# the ggplot2 panels.
	# ════════════════════════════════════════════════════════════════════
	chordFile <- file.path(figDir, tag("panelA_chord.png"))
	plot_double_high_chord(accessData, sigma = circleSigma, outFile = chordFile)
	chordGrob <- grid::rasterGrob(png::readPNG(chordFile), interpolate = TRUE)
	pA <- cowplot::ggdraw() + cowplot::draw_grob(chordGrob)

	# ════════════════════════════════════════════════════════════════════
	# PANEL C/D: noise-level alluvial, Bistable (mirror) / (non-mirror)
	# ════════════════════════════════════════════════════════════════════
	roleMirror  <- classify_bistable_role(accessData, "Bistable (mirror)")
	classMirror <- compute_nfates_class(roleMirror, fateRoles = c("A", "B", "C", "D"))
	pC <- plot_double_high_alluvial(roleMirror, classMirror, sigmas = sigmasABC,
		roleLevels = c("A", "B", "C", "D", "Self", "Other"),
		outFile = file.path(figDir, tag("panelC_mirrorAlluvial.jpg")))

	roleNonMirror  <- classify_bistable_role(accessData, "Bistable (non-mirror)")
	classNonMirror <- compute_nfates_class(roleNonMirror, fateRoles = c("A", "B", "C", "D"))
	pD <- plot_double_high_alluvial(roleNonMirror, classNonMirror, sigmas = sigmasABC,
		roleLevels = c("A", "B", "C", "D", "Self", "Other"),
		outFile = file.path(figDir, tag("panelD_nonMirrorAlluvial.jpg")))

	# ════════════════════════════════════════════════════════════════════
	# ASSEMBLE -- single column
	# ════════════════════════════════════════════════════════════════════
	row1 <- plot_grid(pA, labels = c("A"), label_size = LS)
	row2 <- plot_grid(pC + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
	                   labels = c("C"), label_size = LS)
	row3 <- plot_grid(pD + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 30)),
	                   labels = c("D"), label_size = LS)

	fig6S1 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1, 0.8, 0.8))

	ggsave(file.path(finalDir, paste0("Fig6S1_", net, ".jpg")), fig6S1,
	       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 9,
	       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 16,
	       limitsize = FALSE)

	message("Figure 6 supplement 1 (", net, ") complete. Panels + combined saved to: ", figDir)
}

for (net in c("TT4", "TT4SA")) build_fig6_supp(net)
