# figure3.r
#
# PURPOSE: Paper Figure 3 — additive noise reshapes the phenotypic landscape,
#          and reachability / threshold-margin explains the variation.
#   Panel A: TS all-high / single-high MRT trajectories (facet by State)
#   Panel B: TT / TTSA state-class transition (single -> double -> all-high)
#   Panel C: Monostable vs Bistable MRT boxplots by state class (TS)
#   Panel D: Reachability scatter (Max - Threshold), colored by MRT (TS monostable)
#   Panel E: Margin vs MRT correlation across noise levels (TS monostable)
#   Panel F: Margin-colored MRT trajectories (TS monostable, all-high)
#
# Depends on figure_common.r for: fill_mrt, get_top_original, mean_sd,
# add_reachability, plot_mrt_trajectories; and funcsKishore for:
# theme_Publication, read_parameters, read_solutions.
# ============================================================

library(funcsKishore)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggpubr)
library(ggsignif)
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
noiseType     <- "Additive"
dt            <- 0.01

dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(finalDir, recursive = TRUE, showWarnings = FALSE)

# Panel A uses the shared plot_mrt_trajectories() from figure_common.r
# (Panel F's margin-colored variant is commented out below in favor of the
# correlation heatmap, but would also use the shared version).

# ════════════════════════════════════════════════════════════════════════════
# PANEL B: state-class transition across networks (errorbar of mean per class)
# ════════════════════════════════════════════════════════════════════════════
plot_stateclass_transition <- function(nets, noiseType, resultsFolder, dt = 0.01,
                                        param_type = NULL,
                                        stateClasses = c("all-high", "double-high", "single-high"),
                                        outFile = NULL, ncol = 1, panelWidth = 5, panelHeight = 4,
                                        dodge_width = 0.6) {

	dAll <- map_dfr(nets, function(net) {
		f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_results.csv")
		if (!file.exists(f)) { warning("Missing: ", net); return(NULL) }
		d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
		if (!is.null(param_type)) d <- d %>% filter(ParamType %in% param_type)
		d <- fill_mrt(d, unique(d$State))
		d %>% mutate(Network = net)
	})
	if (nrow(dAll) == 0) stop("No data found")

	dAll <- dAll %>%
		filter(StateClass %in% stateClasses) %>%
		group_by(ParamID, ParamType, NoiseLevel, Network, StateClass) %>%
		summarise(MRT = sum(MRT), .groups = "drop") %>%
		mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
			   Network    = factor(Network, levels = nets),
			   StateClass = factor(StateClass, levels = c("single-high", "double-high", "multi-high", "all-high")))

	# Adjacent-pair paired Wilcoxon tests along the single -> double -> all-high
	# narrative (ggplot's discrete scales drop unused factor levels by
	# default, so "multi-high" -- never present here -- doesn't consume a
	# dodge slot; the 3 classes that ARE present dodge in declaration order).
	presentClasses <- levels(droplevels(dAll$StateClass))
	adjacent_pairs <- list(presentClasses[1:2], presentClasses[2:3])

	dWide <- dAll %>%
		select(ParamID, Network, NoiseLevel, StateClass, MRT) %>%
		pivot_wider(names_from = StateClass, values_from = MRT)

	dSig <- dWide %>%
		group_by(Network, NoiseLevel) %>%
		summarise(p_1 = paired_wilcox_p(.data[[adjacent_pairs[[1]][1]]], .data[[adjacent_pairs[[1]][2]]]),
				  p_2 = paired_wilcox_p(.data[[adjacent_pairs[[2]][1]]], .data[[adjacent_pairs[[2]][2]]]),
				  .groups = "drop") %>%
		mutate(stars_1 = stars_from_p(p_1), stars_2 = stars_from_p(p_2))

	dSumm <- dAll %>%
		group_by(Network, NoiseLevel, StateClass) %>%
		summarise(Mean = mean(MRT), SD = sd(MRT), .groups = "drop")

	off <- dodge_width * (seq_along(presentClasses) - (length(presentClasses) + 1) / 2) / length(presentClasses)
	names(off) <- presentClasses

	dBrackets <- dSumm %>%
		group_by(Network, NoiseLevel) %>%
		summarise(yTop = max(Mean + SD, na.rm = TRUE), .groups = "drop") %>%
		left_join(dSig, by = c("Network", "NoiseLevel")) %>%
		mutate(noise_x = as.numeric(NoiseLevel)) %>%
		group_by(Network) %>%
		mutate(yStep = diff(range(yTop, na.rm = TRUE)) %>%
				 { ifelse(. == 0, max(yTop, 0.05) * 0.1, . * 0.12) }) %>%
		ungroup() %>%
		mutate(y_position = yTop + yStep)

	# brack2 sits a step higher than brack1 -- with 8 noise levels and a 3-way
	# dodge, the two brackets' star labels are horizontally close enough that
	# sharing one height would run them together (e.g. "*** ***" reading as
	# one smear of stars); stacking them keeps each legible on its own row.
	brack1 <- dBrackets %>%
		transmute(Network, xmin = noise_x + off[[adjacent_pairs[[1]][1]]],
				  xmax = noise_x + off[[adjacent_pairs[[1]][2]]],
				  y_position, annotation = stars_1) %>%
		mutate(id = row_number())
	brack2 <- dBrackets %>%
		transmute(Network, xmin = noise_x + off[[adjacent_pairs[[2]][1]]],
				  xmax = noise_x + off[[adjacent_pairs[[2]][2]]],
				  y_position = y_position + yStep, annotation = stars_2) %>%
		mutate(id = row_number())

	p <- ggplot(dAll, aes(x = NoiseLevel, y = MRT, color = StateClass)) +
		stat_summary(fun.data = mean_sd, geom = "errorbar",
					 position = position_dodge(width = dodge_width), width = 0.3) +
		stat_summary(fun = mean, geom = "point",
					 position = position_dodge(width = dodge_width), size = 2) +
		geom_signif(data = brack1, aes(xmin = xmin, xmax = xmax, y_position = y_position,
									   annotations = annotation, group = id),
					manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, size = 0.3, textsize = 2.6) +
		geom_signif(data = brack2, aes(xmin = xmin, xmax = xmax, y_position = y_position,
									   annotations = annotation, group = id),
					manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, size = 0.3, textsize = 2.6) +
		facet_wrap(vars(Network), ncol = ncol) +
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

# mrt_boxplot_stability() (Panel B) and plot_reach_scatter() (Panel C) now
# live in figure_common.r, shared with Figure3_supp.r.

# ════════════════════════════════════════════════════════════════════════════
# PANEL E: margin vs MRT correlation across noise levels
# ════════════════════════════════════════════════════════════════════════════
plot_reach_correlation <- function(net, noiseType, dt, dataFolder, resultsFolder,
                                    state, sigma = c(0.001, 0.01, 0.1, 1),
                                    param_type = NULL, reachableOnly = FALSE,
                                    outFile = NULL, panelWidth = 5, panelHeight = 5) {

	df <- add_reachability(net, noiseType, dt, dataFolder, resultsFolder) %>%
		filter(NoiseLevel %in% sigma, State == state)
	diffCols <- names(df)[str_ends(names(df), "_diff")]
	df <- df %>% mutate(MinMargin = do.call(pmin, as.list(df[diffCols])))

	if (!is.null(param_type)) df <- df %>% filter(ParamType %in% param_type)
	if (reachableOnly) {
		reachCols <- names(df)[str_ends(names(df), "_reach")]
		df <- df %>% filter(if_all(all_of(reachCols), ~ . == 1))
	}
	if (nrow(df) == 0) stop("No rows remain after filtering")

	df <- df %>%
		mutate(ParamType  = str_to_sentence(ParamType),
			   NoiseLevel = factor(NoiseLevel, levels = sort(unique(NoiseLevel))))

	# "All-high" if every discretized dimension in the state tuple is 1
	# (e.g. "(1, 1)"), otherwise fall back to the literal state label.
	y_label <- if (!str_detect(state, "0")) "MRT of all-high state" else paste0("MRT of ", state)

	# Bigger font, no exact p-value text -- just a single "*" appended when
	# p < 0.05 (output.type = "text" so r.label is plain "R = 0.42", not a
	# plotmath expression, and can be pasted with the star safely).
	p <- ggplot(df, aes(x = MinMargin, y = MRT, color = NoiseLevel)) +
		geom_point(alpha = 0.6, size = 1) +
		ggpubr::stat_cor(aes(color = NoiseLevel,
							 label = paste0(after_stat(r.label),
							                ifelse(after_stat(p) < 0.05, " *", ""))),
						 method = "spearman", output.type = "text",
						 label.x.npc = 0.7, label.y.npc = 0.7, size = 4.2) +
		scale_color_viridis_d(name = "σ") +
		theme_Publication() +
		labs(x = "Min(Max - Threshold) across nodes", y = y_label)

	if (!is.null(param_type) && length(param_type) > 1) p <- p + facet_wrap(vars(ParamType))

	if (!is.null(outFile)) {
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = panelWidth, height = panelHeight)
	}
	invisible(p)
}

plot_correlation_heatmap <- function(nets, noiseType, dt, dataFolder, resultsFolder,
                                      state = "auto", sigma = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1),
                                      reachableOnly = FALSE, method = "spearman",
                                      show_values = TRUE,
                                      outFile = NULL, panelWidth = 6, panelHeight = 6) {

	dCorr <- map_dfr(nets, function(net) {
		st <- if (identical(state, "auto")) {
			parameters <- read_parameters(file.path(dataFolder, paste0(net, "_parameters.dat")))
			n_nodes <- sum(str_detect(colnames(parameters), "^Prod_of_"))
			paste0("(", paste(rep(1, n_nodes), collapse = ", "), ")")
		} else state

		d <- add_reachability(net, noiseType, dt, dataFolder, resultsFolder) %>%
			filter(NoiseLevel %in% sigma, State == st)
		if (nrow(d) == 0) { warning("No data for ", net); return(NULL) }

		diffCols <- names(d)[str_ends(names(d), "_diff")]
		d <- d %>% mutate(MinMargin = do.call(pmin, as.list(d[diffCols])))

		if (reachableOnly) {
			reachCols <- names(d)[str_ends(names(d), "_reach")]
			d <- d %>% filter(if_all(all_of(reachCols), ~ . == 1))
		}
		if (nrow(d) == 0) return(NULL)

		# collapse stability class: monostable vs multistable (bi + tri + ...)
		d <- d %>%
			mutate(StabilityClass = ifelse(ParamType == "monostable",
										   "Monostable", "Multistable"))

		# one correlation + p-value per (network, noise level, stability class)
		d %>%
			group_by(StabilityClass, NoiseLevel) %>%
			summarise(
				res = list(
					if (n() >= 3 && sd(MinMargin) > 0 && sd(MRT) > 0)
						suppressWarnings(cor.test(MinMargin, MRT, method = method))
					else NULL
				),
				nSets = n(),
				.groups = "drop"
			) %>%
			mutate(
				R = map_dbl(res, ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)),
				P = map_dbl(res, ~ if (is.null(.x)) NA_real_ else .x$p.value),
				Network = net
			) %>%
			select(-res)
	})

	if (nrow(dCorr) == 0) stop("No correlation data computed for any network")

	dCorr <- dCorr %>%
		mutate(Network        = factor(Network, levels = rev(nets)),
			   NoiseLevel     = factor(NoiseLevel, levels = sort(unique(NoiseLevel))),
			   StabilityClass = factor(StabilityClass, levels = c("Monostable", "Multistable")),
			   Stars = case_when(
				   is.na(P)   ~ "ns",
				   P < 0.05   ~ "*",
				   TRUE       ~ "ns"
			   ),
			   TileLabel = ifelse(is.na(R), "",
							 ifelse(show_values,
									paste0(sprintf("%.2f", R), "\n", Stars),
									Stars)))

	p <- ggplot(dCorr, aes(x = NoiseLevel, y = Network, fill = R)) +
		geom_tile(color = "white", linewidth = 0.5) +
		geom_text(aes(label = TileLabel), size = 3, color = "black", lineheight = 0.8) +
		scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
							  midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1),
							  name = paste0(str_to_title(method), "\ncorrelation")) +
		facet_wrap(vars(StabilityClass)) +
		theme_Publication() +
		theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
		labs(x = "Noise Level", y = NULL)

	if (!is.null(outFile)) {
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = panelWidth, height = panelHeight)
	}
	invisible(list(plot = p, data = dCorr))
}


# ════════════════════════════════════════════════════════════════════════════
# BUILD MARGIN DATA (for Panel F), then generate all panels
# ════════════════════════════════════════════════════════════════════════════
marginData_TS <- add_reachability("TS", noiseType, dt, dataFolder, resultsFolder) %>%
	distinct(ParamID, across(ends_with("_diff"))) %>%
	mutate(MinMargin = do.call(pmin, as.list(select(., ends_with("_diff"))))) %>%
	select(ParamID, MinMargin)

# ── Panel A ──
pA <- plot_mrt_trajectories("TS", noiseType, resultsFolder, dt = dt,
	facet_var = "State", topStates = 2, ncol = 1,
	sampleN = 100,
	outFile = file.path(figDir, "fig3_panelA_TS_traj.jpg"),
	panelWidth = 6, panelHeight = 4) +
	theme(strip.text = element_text(size = rel(1.2)))

# ── Panel B ──
pB <- mrt_boxplot_stability("TS", noiseType, resultsFolder, dt = dt,
	outFile = file.path(figDir, "fig3_panelC_boxplot.jpg")) +
	theme(strip.text = element_text(size = rel(1.2)))

# ── Panel C ──
# The real data only extends down to about -15/-18 on the two axes (99% of
# points are well above 0), so pushing the visible range out to -300 opens
# up an empty margin along the bottom-left with no data in it -- the
# inset legend moves there (was sitting at (0.22, 0.5), which overlapped a
# real diagonal band of points and was unreadable against them).
pC <- plot_reach_scatter("TS", noiseType, dt, dataFolder, resultsFolder,
	state = "(1, 1)", sigma = 0.1, param_type = "monostable",
	outFile = file.path(figDir, "fig3_panelD_reachScatter.jpg")) +
    coord_cartesian(xlim = c(-300, 1000), ylim = c(-300, 1000)) +
    theme(legend.position = c(0.15, 0.5), legend.direction = "vertical", legend.key.height = unit(0.8, "cm"), legend.background = element_rect(fill = NA), legend.key.width = unit(0.1, "cm")) +
	theme(strip.text = element_text(size = rel(1.2))) + labs(color = "MRT of\nall-high\nstate")

# ── Panel D ──
pD <- plot_reach_correlation("TS", noiseType, dt, dataFolder, resultsFolder,
	state = "(1, 1)", sigma = c(0.001, 0.01, 0.1, 1), param_type = "monostable",
	outFile = file.path(figDir, "fig3_panelE_correlation.jpg")) +
    scale_color_discrete() + labs(y = "All-high state MRT") +
	theme(axis.title.y = element_text(hjust = 1)) +
	theme(strip.text = element_text(size = rel(1.2)))

# ── Panel E ──
pE <- plot_stateclass_transition(c("TT", "TTSA"), noiseType, resultsFolder, dt = dt,
	param_type = "monostable", ncol = 1,
	outFile = file.path(figDir, "fig3_panelB_transition.jpg")) + 
    theme(legend.position = "top")

# ── Panel E ──
pF <- plot_correlation_heatmap(
	nets = c("TS", "TSSA", "TT", "TTSA"),
	noiseType = noiseType, dt = dt,
	dataFolder = dataFolder, resultsFolder = resultsFolder,
	show_values = F, 
	outFile = file.path(figDir, "fig3_panelF_corrHeatmap.jpg")
)$plot + theme(legend.key.width = unit(0.4, "cm"))

# # ── Panel F ──
# pF <- plot_mrt_trajectories("TS", noiseType, resultsFolder, dt = dt,
# 	param_type = "monostable", marginData = marginData_TS,
# 	facet_var = "StateClass", stateClasses = "all-high", topStates = NULL, ncol = 1,
# 	outFile = file.path(figDir, "fig3_panelF_marginTraj.jpg"),
# 	panelWidth = 6, panelHeight = 4)

# ── ASSEMBLE ──
middle_top <- plot_grid(pC, pD, ncol = 1, labels = c("C", "D"), rel_heights = c(1,1), label_size = LS)
left_top <- plot_grid(pA,middle_top, ncol = 2, labels = c("A", ""), rel_widths = c(1, 1), label_size = LS)
left_col  <- plot_grid(left_top, pB, ncol = 1, labels = c("", "B"), rel_heights = c(1.7, 1), label_size = LS)
# mid_col   <- plot_grid(pB, ncol = 1, labels = "B")
right_col <- plot_grid(pE, pF, ncol = 1, labels = c("E", "F"), rel_heights = c(2, 1), label_size = LS)

fig3 <- plot_grid(left_col, right_col, ncol = 2, rel_widths = c(1.8, 1))
ggsave(file.path(finalDir, "Fig3.jpg"), fig3,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 15,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 13)

message("Figure 3 complete. Panels + combined saved to: ", figDir)

### Figure 3 supplementary


pA <- plot_mrt_trajectories("TSSA", noiseType, resultsFolder, dt = dt,
	facet_var = "State", topStates = 2, ncol = 1,
	outFile = file.path(figDir, "fig3_panelA_TS_traj.jpg"),
	panelWidth = 6, panelHeight = 4)

# ── Panel B ──
pB <- mrt_boxplot_stability("TSSA", noiseType, resultsFolder, dt = dt,
	outFile = file.path(figDir, "fig3_panelC_boxplot.jpg"))
