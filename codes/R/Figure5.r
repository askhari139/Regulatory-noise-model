# figure5.r
#
# PURPOSE: Paper Figure 5 — MULTIPLICATIVE noise (contrast to additive Figs 3-4).
#   Panel A: TTSA MRT trajectories — single+double high -> single high alone
#   Panel B: single-high vs double-high frequency, TS/TSSA/TT/TTSA, 2x2
#   Panel C: single-high vs double-high frequency, TT4/TS4, stacked 2x1
#   Panel D: stochastic - resampled MRT difference (violin), single-high
#            states, across all noise levels, TS/TSSA/TT/TTSA/TT4/TT4SA, 2x3
#
# The TT4 double-high-accessibility panels (old Panels E/F/G: nFates
# distribution, heatmap, ratio violin) moved out to Figure6.r; the
# TS/TSSA/TT trajectory + lambda/topology/state-frequency panels (old
# Fig5S1 Panels A-D) stay in Figure5_supp.r; old Fig5S1's remaining panels
# E-I moved to the new Figure5_supp2.r along with two new resample-effect
# counterparts.
#
# Depends on figure_common.r (fill_mrt, mean_sd, plot_mrt_trajectories,
# get_det_stoch_comparison_data, plot_class_diff_violin) and funcsKishore
# (theme_Publication, read_parameters, read_solutions).
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpubr)
.scriptDir <- (function() {
    a <- commandArgs(trailingOnly = FALSE)
    fa <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(fa) > 0) return(dirname(normalizePath(fa[1])))
    for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
    getwd()
})()
source(file.path(.scriptDir, "figure_common.r"))

# ── CONFIG ────────────────────────────────────────────────────────────────────
figDir        <- "figures/individual"
finalDir      <- "figures/final"
noiseType     <- "Multiplicative"
dt            <- 0.01
sigmaFixed    <- 0.1

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

# ════════════════════════════════════════════════════════════════════════════
# PANEL A: TTSA trajectories — single+double -> single alone
# ════════════════════════════════════════════════════════════════════════════
pA <- plot_mrt_trajectories(
	net = "TTSA", noiseType = noise_type_for("TTSA", noiseType), resultsFolder = resultsFolder, dt = dt,
	facet_var = "StateClass",
	stateClasses = c("single-high", "double-high"),
	topStates = NULL, ncol = 1,
	sampleN = 100,
	outFile = file.path(figDir, "fig5_panelA_TTSA_traj.jpg"),
	panelWidth = 6, panelHeight = 4
) +
	theme(strip.text = element_text(size = rel(1.2)))

# ════════════════════════════════════════════════════════════════════════════
# PANEL B / C: single-high vs double-high frequency across networks
# (SD errorbars) -- TS/TSSA all-high folded locally into "double-high" for
# this panel only; fill_mrt untouched.
# ════════════════════════════════════════════════════════════════════════════
plot_single_vs_double <- function(nets, noiseType, resultsFolder, dt = 0.01,
                                   param_type = NULL,
                                   outFile = NULL, ncol = 3, panelWidth = 3.5, panelHeight = 4) {

	dAll <- map_dfr(nets, function(net) {
		effectiveNoiseType <- noise_type_for(net, noiseType)
		f <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
		if (!file.exists(f)) { warning("Missing: ", net); return(NULL) }
		d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
		if (!is.null(param_type)) d <- d %>% filter(ParamType %in% param_type)
		# TT4/TT4SA were only ever simulated at DT = 1/10 (never 0.01), so a
		# DT == 0.01 filter leaves 0 rows for them -- d$State[1] would index
		# past the end of an empty vector (NA), and str_count(NA, ",") + 1
		# blows up the n_nodes == 2 check below with "missing value where
		# TRUE/FALSE needed". Skip the network instead of crashing, same as
		# the missing-file case above.
		if (nrow(d) == 0) { warning("No rows at DT=", dt, ": ", net); return(NULL) }
		d <- fill_mrt(d, unique(d$State))
		n_nodes <- str_count(d$State[1], ",") + 1
		if (n_nodes == 2) {
			d <- d %>% mutate(StateClass = ifelse(StateClass == "all-high",
												  "double-high", StateClass))
		}
		d %>% mutate(Network = net)
	})
	if (nrow(dAll) == 0) stop("No data found")

	dAll <- dAll %>%
		filter(StateClass %in% c("single-high", "double-high")) %>%
		group_by(ParamID, ParamType, NoiseLevel, Network, StateClass) %>%
		summarise(MRT = sum(MRT), .groups = "drop") %>%
		mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
			   Network    = factor(Network, levels = nets),
			   StateClass = factor(StateClass, levels = c("single-high", "double-high")))

	p <- ggplot(dAll, aes(x = NoiseLevel, y = MRT, color = StateClass)) +
		stat_summary(fun.data = mean_sd, geom = "errorbar",
					 position = position_dodge(width = 0.6), width = 0.3) +
		stat_summary(fun = mean, geom = "point",
					 position = position_dodge(width = 0.6), size = 2) +
		facet_wrap(vars(Network), ncol = ncol, scales = "fixed") +
		theme_Publication() +
		theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
		labs(x = "Noise Level", y = "Average MRT")

	if (!is.null(outFile)) {
		n_panels <- length(unique(dAll$Network))
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = panelWidth * min(n_panels, ncol),
			   height = panelHeight * ceiling(n_panels / ncol))
	}
	invisible(p)
}

# Panel B: 2x2
pB <- plot_single_vs_double(
	nets = c("TS", "TSSA", "TT", "TTSA"),
	noiseType = noiseType, resultsFolder = resultsFolder, dt = dt, ncol = 2,
	outFile = file.path(figDir, "fig5_panelB_single_vs_double.jpg")
) +
	theme(strip.text = element_text(size = rel(1.2)))

# Panel C: stacked 2x1
pC <- plot_single_vs_double(
	nets = c("TT4", "TS4"),
	noiseType = noiseType, resultsFolder = resultsFolder, dt = dt, ncol = 1,
	outFile = file.path(figDir, "fig5_panelC_single_vs_double_4node.jpg")
) +
	theme(strip.text = element_text(size = rel(1.2)))

# ════════════════════════════════════════════════════════════════════════════
# PANEL D: stochastic effect (Stochastic - Resampled MRT), single-high
# states, across all noise levels -- TS/TSSA/TT/TTSA/TT4/TT4SA (2x3),
# including the SA networks alongside the base ones in the same plot.
# ════════════════════════════════════════════════════════════════════════════
pD <- plot_class_diff_violin(
	nets = c("TS", "TSSA", "TT", "TTSA", "TT4", "TT4SA"), noiseType = noiseType, dt = dt,
	dataFolder = dataFolder, resultsFolder = resultsFolder,
	compare = "stochastic", stateClasses = "single-high"
) + facet_wrap(vars(Network), ncol = 3, scales = "fixed")
ggsave(file.path(figDir, "fig5_panelD_stochDiff_singleHigh.jpg"), pD, width = 12, height = 8)

# ════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ════════════════════════════════════════════════════════════════════════════
row1 <- plot_grid(pA, pB + theme(legend.position = "top"), nrow = 1,
                   rel_widths = c(1, 1.6), labels = c("A", "B"), label_size = LS)
row2 <- plot_grid(pC + theme(legend.position = "top"), pD, nrow = 1,
                   rel_widths = c(1, 1.8), labels = c("C", "D"), label_size = LS)

fig5 <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1.2))

ggsave(file.path(finalDir, "Fig5.jpg"), fig5,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 16,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 14)

message("Figure 5 complete. Panels + combined saved to: ", figDir)
