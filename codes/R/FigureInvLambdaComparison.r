# FigureInvLambdaComparison.r
#
# PURPOSE: Compare Multiplicative vs. MultiplicativeInvLambda noise, for the
# self-activation networks (TSSA, TTSA, TT4SA, NFSA) plus DA -- the exact
# set covered by the MultiplicativeInvLambda run. Two separate figures:
#
#   FigInvLambda_MRT.jpg       -- raw MRT vs. noise level, single-high /
#                                  double-high / all-high state classes,
#                                  colored by Method, faceted by Network.
#   FigInvLambda_StochDiff.jpg -- stochastic effect (Stochastic - Resampled
#                                  MRT), same state classes/networks, violin
#                                  per noise level like Figure5.r's Panel D
#                                  (plot_class_diff_violin), but with two
#                                  dodged violins per x-position (filled by
#                                  Method) instead of one solid-color violin,
#                                  since that function only ever compares
#                                  stoch-vs-det within a single noise mode.
#
# Depends on figure_common.r (fill_mrt, mean_sd, get_det_stoch_comparison_data)
# and funcsKishore (theme_Publication).
# ============================================================

library(funcsKishore)
library(tidyverse)
.scriptDir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(fa) > 0) return(dirname(normalizePath(fa[1])))
    for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
    getwd()
})()
source(file.path(.scriptDir, "figure_common.r"))

figDir   <- "figures/individual"
finalDir <- "figures/final"
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

nets         <- c("DA", "TSSA", "TTSA", "TT4SA", "NFSA")
dt           <- 0.01
stateClasses <- c("single-high", "double-high", "all-high")
methods      <- c(Multiplicative = "Multiplicative", MultiplicativeInv = "MultiplicativeInvLambda")

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 1: raw MRT by state class, colored by Method, faceted by Network
# ════════════════════════════════════════════════════════════════════════════
mrtData <- map_dfr(names(methods), function(methodLabel) {
    noiseType <- methods[[methodLabel]]
    map_dfr(nets, function(net) {
        f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_results.csv")
        if (!file.exists(f)) { warning("Missing: ", noiseType, "/", net); return(NULL) }
        d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
        fill_mrt(d, unique(d$State)) %>%
            filter(StateClass %in% stateClasses) %>%
            group_by(ParamID, ParamType, NoiseLevel, StateClass) %>%
            summarise(MRT = sum(MRT), .groups = "drop") %>%
            mutate(Network = net, Method = methodLabel)
    })
})

mrtData <- mrtData %>%
    mutate(Network    = factor(Network, levels = nets),
           StateClass = factor(StateClass, levels = stateClasses),
           NoiseLevel = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
           Method     = factor(Method, levels = names(methods)))

pMRT <- ggplot(mrtData, aes(x = NoiseLevel, y = MRT, color = Method, fill = Method)) +
    stat_summary(fun.data = mean_sd, geom = "errorbar",
                 position = position_dodge(width = 0.6), width = 0.3) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 0.6), size = 2) +
    facet_grid(rows = vars(StateClass), cols = vars(Network), scales = "free_y") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
    labs(x = "Noise Level", y = "Average MRT",
         title = "Multiplicative vs. MultiplicativeInvLambda: MRT by state class")

ggsave(file.path(finalDir, "FigInvLambda_MRT.jpg"), pMRT,
       width = 4 * length(nets), height = 3.5 * length(stateClasses), limitsize = FALSE)

# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2: stochastic effect (Stochastic - Resampled MRT), violin per noise
# level, dodged/filled by Method -- same networks and state classes as above.
# ════════════════════════════════════════════════════════════════════════════
diffData <- map_dfr(names(methods), function(methodLabel) {
    noiseType <- methods[[methodLabel]]
    map_dfr(stateClasses, function(cls) {
        get_det_stoch_comparison_data(nets, noiseType, dt, dataFolder, resultsFolder,
                                       focal_class = cls) %>%
            mutate(StateClass = cls)
    }) %>% mutate(Method = methodLabel)
})

diffData <- diffData %>%
    mutate(Diff = MRT - MRT_det) %>%
    filter(NoiseLevel != 0) %>%
    mutate(Network    = factor(Network, levels = nets),
           StateClass = factor(StateClass, levels = stateClasses),
           NoiseLevel = factor(NoiseLevel, levels = sort(unique(NoiseLevel))),
           Method     = factor(Method, levels = names(methods)))

pDiff <- ggplot(diffData, aes(x = NoiseLevel, y = Diff, fill = Method)) +
    geom_violin(alpha = 0.6, scale = "width", color = NA,
                position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.4) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6, color = "black",
                 position = position_dodge(width = 0.8)) +
    facet_grid(rows = vars(StateClass), cols = vars(Network), scales = "free_y") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
    labs(x = "Noise Level", y = "Stochastic - Resampled MRT",
         title = "Stochastic effect: Multiplicative vs. MultiplicativeInvLambda")

ggsave(file.path(finalDir, "FigInvLambda_StochDiff.jpg"), pDiff,
       width = 4 * length(nets), height = 3.5 * length(stateClasses), limitsize = FALSE)

message("Done: FigInvLambda_MRT.jpg, FigInvLambda_StochDiff.jpg -> ", finalDir)
