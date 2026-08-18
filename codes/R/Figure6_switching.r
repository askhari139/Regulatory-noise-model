# figure6_switching.r
#
# PURPOSE: Switching-event stability, generalized across noise types and
# networks in one panel -- Figure6.r Panel A (plot_switching_stability())
# repeated per network only for "Multiplicative" noise; this instead
# compares all 3 canonical noise types (Additive/Multiplicative/
# Fluctuating) at once, fill = NoiseType, faceted by Network (cols) x
# StabilityClass (rows). MeanSwitches is normalized into a per-timestep
# rate using the pipeline's actual post-burn-in timestep count (500, not
# 700 -- see POST_BURNIN_TIMESTEPS's comment in figure_common.r for how
# this was verified against the real simulation code and a submitted
# cluster job script).
#
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

# ── CONFIG ────────────────────────────────────────────────────────────────────
finalDir   <- "figures/final"
nets       <- c("TT4", "TT4SA")
noiseTypes <- c("Additive", "Multiplicative", "Fluctuating")
dt         <- 0.01

dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

p <- plot_switching_stability_combined(nets, noiseTypes, resultsFolder, dt = dt,
	outFile = file.path(finalDir, "Fig6_switchingByNoiseType.jpg"), width = 10, height = 8)

message("Switching-events-by-noise-type figure complete: figures/final/Fig6_switchingByNoiseType.jpg")
