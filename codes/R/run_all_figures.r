# run_all_figures.r
#
# PURPOSE: Master driver script. For each figure-generating script, sets LS
# (plot_grid() label size) and the final figure's width/height, then sources
# that script. figure_common.r (sourced once below, and again -- harmlessly
# -- by every Figure*.r script) only supplies LS/FIG_WIDTH/FIG_HEIGHT
# defaults when they aren't already set, so setting them here overrides the
# per-script defaults without needing to touch each script every time.
#
# Each Figure*.r script is still independently runnable on its own (Rscript
# scripts/FigureN.r) -- when run that way LS/FIG_WIDTH/FIG_HEIGHT don't
# exist yet, so every script falls back to its own hardcoded defaults.
#
# Panels saved to figures/individual/ keep their own fixed sizes regardless
# of this file -- only the combined, lettered output (Fig1.jpg, Fig2.jpg,
# ...) responds to FIG_WIDTH/FIG_HEIGHT. Figure2.r emits two combined
# figures (Fig2.jpg, Fig2S1.jpg); one LS/width/height pair applies to both.
# MRT_by_DT.r has no lettered panels (LS is ignored) and emits three files
# (one per noise mode) built from a single panelWidth/panelHeight -- FIG_
# WIDTH/FIG_HEIGHT double as that panel size for this one script.
#
# ============================================================

scriptDir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(fa) > 0) return(dirname(normalizePath(fa[1])))
    for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
    getwd()
})()
source(file.path(scriptDir, "figure_common.r"))  # defines default LS, resultsFolder, dataFolder

figureConfigs <- list(
    list(script = "Figure1.r",         LS = 20, width = 14, height = 18),
    list(script = "Figure1_supp.r",    LS = 20, width = 12, height = 8),
    list(script = "Figure2.r",         LS = 20, width = 20, height = 12),
    list(script = "Figure3.r",         LS = 20, width = 15, height = 13),
    list(script = "Figure3_supp.r",    LS = 20, width = 13, height = 15),
    list(script = "Figure4.r",         LS = 20, width = 16, height = 20),
    list(script = "Figure4_supp_BC.r", LS = 20, width = 20, height = 16),
    list(script = "Figure4_supp_D.r",  LS = 20, width = 30, height = 20),
    list(script = "Figure5.r",         LS = 20, width = 16, height = 14),
    list(script = "Figure5_supp.r",    LS = 20, width = 18, height = 15),
    list(script = "Figure5_supp2.r",   LS = 20, width = 22, height = 16),
    list(script = "Figure6.r",         LS = 20, width = 18, height = 20),
    list(script = "MRT_by_DT.r",       LS = 20, width = 3.2, height = 3)
)

for (cfg in figureConfigs) {
    LS         <- cfg$LS
    FIG_WIDTH  <- cfg$width
    FIG_HEIGHT <- cfg$height
    message("== Running ", cfg$script,
            " (LS = ", LS, ", ", FIG_WIDTH, " x ", FIG_HEIGHT, ") ==")
    source(file.path(scriptDir, cfg$script))
}

message("All figures regenerated (panels in figures/individual/, final figures in figures/final/).")
