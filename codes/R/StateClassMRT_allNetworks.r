# StateClassMRT_allNetworks.r
#
# PURPOSE: Diagnostic figure -- for every network, Multiplicative noise
# (MultiplicativeInvLambda substituted in for self-activation networks +
# DA, via figure_common.r's noise_type_for(); labels stay "Multiplicative"),
# plot mean MRT (summed within each StateClass, averaged across ParamIDs,
# +-SD errorbars) vs NoiseLevel for exactly 4 classes: single-high,
# double-high, all-high, all-low. Checks whether the all-low class
# visibly grows at high noise -- prompted by TS's density-scatter plot
# showing a real (if minority, ~16% of points) low-low cluster at sigma=1
# that TS's MRT-based Fig5 Panel B analysis doesn't surface on its own
# (that panel folds all-high into double-high and never plots all-low at
# all). For 2-node networks (TS/TSSA/NF/NFSA/DA), "double-high" and
# "all-high" are the same physical state -- fill_mrt() always classifies
# by n_nodes first, so for these networks the state is labeled
# "all-high" only, and "double-high" is simply empty/absent (not folded
# in, unlike Figure5.r's Panel B/C which folds it in as "double-high" for
# a different comparison's sake -- here every network's classes are kept
# separate and literal, since the whole point is to see all 4 classes at
# once).
# ============================================================

library(funcsKishore)
suppressMessages(library(tidyverse))
library(cowplot)
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

nets       <- c("TS", "TSSA", "TT", "TTSA", "TT4", "TT4SA", "TS4", "NF", "NFSA", "DA")
noiseType  <- "Multiplicative"
dt         <- 0.01
wantClasses <- c("single-high", "double-high", "all-high", "all-low")

dAll <- map_dfr(nets, function(net) {
    effectiveNoiseType <- noise_type_for(net, noiseType)
    f <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
    if (!file.exists(f)) { warning("Missing: ", net); return(NULL) }
    d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
    d <- fill_mrt(d, unique(d$State)) %>%
        filter(StateClass %in% wantClasses) %>%
        group_by(ParamID, ParamType, NoiseLevel, StateClass) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        mutate(Network = net)
})

dAll <- dAll %>%
    mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
           Network    = factor(Network, levels = nets),
           StateClass = factor(StateClass, levels = wantClasses))

p <- ggplot(dAll, aes(x = NoiseLevel, y = MRT, color = StateClass)) +
    stat_summary(fun.data = mean_sd, geom = "errorbar",
                 position = position_dodge(width = 0.6), width = 0.3) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 0.6), size = 1.6) +
    stat_summary(fun = mean, geom = "line", aes(group = StateClass),
                 position = position_dodge(width = 0.6), linewidth = 0.4, alpha = 0.5) +
    facet_wrap(vars(Network), ncol = 5, scales = "fixed") +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
    labs(x = "Noise Level", y = "Mean MRT (+/- SD across parameters)",
         color = "State class",
         title = "MRT by state class across noise levels, Multiplicative noise (InvLambda for SA networks + DA)")

outFile <- file.path(finalDir, "StateClassMRT_allNetworks.jpg")
ggsave(outFile, p, width = 20, height = 9)
message("saved: ", outFile)
