# MRT_by_DT.r
#
# PURPOSE: Exploratory check of DT (integration/resampling time step)
# sensitivity -- MRT vs. noise level (0.001, 0.01, 0.1, 1.0), colored by DT
# (0.01, 1, 10), faceted by state class x stability class (rows) x network
# (cols), with paired Wilcoxon test significance brackets between adjacent
# DT values (paired by ParamID), plus a dashed horizontal line per facet
# marking the deterministic (NoiseLevel = 0) mean MRT as a reference --
# DT doesn't apply at zero noise (nothing to resample), so it's drawn as a
# single reference line rather than another dodged/bracketed x position.
# One standalone figure per noise mode
# (Additive, Multiplicative, Fluctuating -- "Stationary" was requested too,
# but that noise type only exists under a separate FoldChange/ experiment
# with an incompatible schema, no DT column at all; "Fluctuating" is the
# noise mode in the main results pipeline with a matching schema, used as
# the substitute).
#
# Depends on figure_common.r's plot_mrt_by_dt().
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
library(ggsignif)
.scriptDir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(fa) > 0) return(dirname(normalizePath(fa[1])))
    for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
    getwd()
})()
source(file.path(.scriptDir, "figure_common.r"))

finalDir <- "figures/final"
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

nets        <- c("TS", "TSSA", "TT", "TTSA")
noiseLevels <- c(0.001, 0.01, 0.1, 1.0)
dtValues    <- c(0.01, 1, 10)
stateClasses <- c("all-high", "single-high", "double-high")
noiseModes  <- c("Additive", "Multiplicative", "Fluctuating")

for (nm in noiseModes) {
	plot_mrt_by_dt(
		nets, noiseType = nm, resultsFolder = resultsFolder,
		noise_levels = noiseLevels, dt_values = dtValues,
		state_classes = stateClasses,
		outFile = file.path(finalDir, paste0("mrtByDT_", tolower(nm), ".jpg")),
		panelWidth = if (exists("FIG_WIDTH")) FIG_WIDTH else 3.2,
		panelHeight = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 3
	)
}

message("MRT-by-DT plots complete (one figure per noise mode). Saved to: ", finalDir)
