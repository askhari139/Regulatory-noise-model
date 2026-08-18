# MRT_vs_Switches.r
#
# PURPOSE: Reachability-style scatter of per-ParamID MRT (of a given state
# class) against the mean number of state switches per stochastic
# trajectory (normalized by the number of saved timesteps), faceted by
# network (rows) x noise level (cols, 0.001/0.01/0.1/1.0). One figure per
# (noise mode, state class) combination -- TS/TSSA have no genuine
# double-high state of their own (fill_mrt() folds a 2-node network's
# 2-high state into "all-high"), so those rows come back empty for the
# double-high figures -- expected, not an error. Also produces a 2D-density
# heatmap version of the same data (mrtVsSwitchesHeatmap_*.jpg) -- more
# legible than the scatter where a panel is dominated by thousands of
# overplotted points.
#
# Depends on figure_common.r's compute_mrt_switches_data() /
# plot_mrt_vs_switches() / plot_mrt_vs_switches_heatmap().
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
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

nets         <- c("TS", "TSSA", "TT", "TTSA")
noiseLevels  <- c(0.001, 0.01, 0.1, 1.0)
dt           <- 0.01
stateClasses <- c("all-high", "single-high", "double-high")
noiseModes   <- c("Additive", "Multiplicative", "Fluctuating")

for (nm in noiseModes) {
	for (cls in stateClasses) {
		clsTag <- gsub("-", "", cls)
		plot_mrt_vs_switches(
			nets, noiseType = nm, resultsFolder = resultsFolder, dt = dt,
			stateClass = cls, noise_levels = noiseLevels,
			outFile = file.path(finalDir, paste0("mrtVsSwitches_", tolower(nm), "_", clsTag, ".jpg")),
			panelWidth = 3.2, panelHeight = 3
		)
		plot_mrt_vs_switches_heatmap(
			nets, noiseType = nm, resultsFolder = resultsFolder, dt = dt,
			stateClass = cls, noise_levels = noiseLevels,
			outFile = file.path(finalDir, paste0("mrtVsSwitchesHeatmap_", tolower(nm), "_", clsTag, ".jpg")),
			panelWidth = 3.2, panelHeight = 3
		)
	}
}

message("MRT-vs-switches plots complete (scatter + heatmap, one of each per noise mode x state class). Saved to: ", finalDir)
