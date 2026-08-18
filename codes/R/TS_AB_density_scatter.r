# TS_AB_density_scatter.r
#
# PURPOSE: Same scope as TS_noise_violin_reconstruction.r (TS network, 5
# noise levels, one figure per noise mode) but showing the JOINT (A, B)
# relationship instead of each node's marginal distribution: one point per
# ParamID at (Mean_A, Mean_B), with a 2D density contour overlay, faceted
# by noise level. Unlike the violin reconstruction, this uses the raw
# per-parameter Mean_A/Mean_B pairs directly (real data, no skew-normal
# fitting needed -- the "distribution" here is the actual empirical joint
# distribution of per-parameter means).
#
# Threshold lines (same per-node RACIPE deterministic-mean values as the
# violin figures) are drawn as a vertical (node A) + horizontal (node B)
# line, in each node's own color, marking the same quadrant-defining
# reference the violins showed on their single axis.
# ============================================================

library(funcsKishore)
suppressMessages(library(tidyverse))

resultsFolder <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final"
dataFolder    <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/data"
finalDir      <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/figures/final"

net         <- "TS"
dt          <- 0.01
noiseLevels <- c(0, 0.001, 0.01, 0.1, 1.0)
nodes       <- c("A", "B")

solutions  <- read_solutions(file.path(dataFolder, paste0(net, "_solution.dat")))
# TS_solution.dat's node columns are in log2 space (values go negative,
# e.g. -6.4, which raw expression can't) -- convert back to linear units
# via a basin-frequency-weighted mean in log2 space, then 2^(.), matching
# the Julia pipeline's own get_mean_expression() (src/RACIPEdata.jl) and
# figure_common.r's add_reachability() (same fix applied there).
solWeights <- solutions$basin / 100
thresholds <- setNames(sapply(nodes, function(nd) 2^(sum(solutions[[nd]] * solWeights) / sum(solWeights))), nodes)
nodeColors <- setNames(scales::hue_pal()(length(nodes)), nodes)

run_for_mode <- function(noiseType) {
    message("== ", noiseType, " ==")
    f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_stats.csv")
    # Pairing key beyond (ParamID, ParamType, NoiseLevel): the stats file
    # has multiple rows per (ParamID, Node) whenever there's more than one
    # thing to summarise separately -- multistable parameters at sigma=0
    # get one row per co-existing attractor (Count = how many of the
    # num_sims initial conditions landed there), while sigma>0 gets one
    # row per individual stochastic replicate (Count == 1 always). Both
    # cases are consistently interleaved A,B,A,B,... in file order (verified
    # directly against the raw CSV), so a plain within-group row_number()
    # pairs each node's rows correctly either way without needing to know
    # which case applies.
    d <- read_csv(f, show_col_types = FALSE) %>%
        filter(DT == dt, NoiseLevel %in% noiseLevels, Node %in% nodes) %>%
        group_by(ParamID, ParamType, NoiseLevel, Node) %>%
        mutate(Rep = row_number()) %>%
        ungroup() %>%
        select(ParamID, ParamType, NoiseLevel, Rep, Node, Mean) %>%
        pivot_wider(names_from = Node, values_from = Mean, names_prefix = "Mean_") %>%
        filter(!is.na(Mean_A), !is.na(Mean_B), Mean_A > 0, Mean_B > 0) %>%
        mutate(NoiseLevel = factor(NoiseLevel, levels = noiseLevels))

    p <- ggplot(d, aes(x = Mean_A, y = Mean_B)) +
        geom_point(alpha = 0.08, size = 0.4, color = "grey30") +
        geom_density_2d(color = "black", linewidth = 0.3) +
        geom_vline(xintercept = thresholds[["A"]], color = nodeColors[["A"]],
                   linetype = "dashed", linewidth = 0.7) +
        geom_hline(yintercept = thresholds[["B"]], color = nodeColors[["B"]],
                   linetype = "dashed", linewidth = 0.7) +
        scale_x_log10() + scale_y_log10() +
        facet_wrap(vars(NoiseLevel), nrow = 1) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
              plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 16, unit = "pt")) +
        labs(x = "Mean expression, Node A (log scale)", y = "Mean expression, Node B (log scale)",
             title = paste0("TS, ", noiseType, " noise: joint (A, B) mean-expression distribution"),
             subtitle = "One point per ParamID; contours = 2D density. Dashed lines = each node's threshold.")

    outFile <- file.path(finalDir, paste0("TS_AB_density_scatter_", noiseType, ".jpg"))
    ggsave(outFile, p, width = 18, height = 4.6, dpi = 300)
    message("saved: ", outFile)
}

for (nt in c("Additive", "Fluctuating", "Multiplicative")) run_for_mode(nt)
