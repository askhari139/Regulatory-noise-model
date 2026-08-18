library(funcsKishore)
library(ggh4x)
library(cowplot)
library(ggpattern)
library(gridExtra)
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

# Computes "nice" axis breaks over [0, max] only, ignoring whatever range
# the deterministic segment's negative PlotTime values span — so no numeric
# labels ever appear over the unlabeled deterministic portion of the axis.
stoch_breaks <- function(limits) {
    br <- scales::extended_breaks()(c(0, limits[2]))
    br[br >= 0]
}

plot_burnin_trajectories <- function(net, paramID, noiseLevel, noiseMode,
	save_path = NULL,
	width = 10,
	height = 8,
    det_time = 300, stoch_time = 300,
	base_path = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata",
	scripts_path = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/scripts",
	data_path = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/data",
	dpi = 300, nTraj = 0, cut_fraction = 0.5,
    top_labs = c("A(i)", "(ii)"), bot_labs = "B") {

	# Locate the trajectory file, generating it via Julia if it doesn't exist yet
	file_path <- paste0(base_path, "/figures/burnin_trajectories/", net, "/param", paramID, "_iter1_",
		noiseMode, "_sigma", noiseLevel %>% as.character() %>% str_replace("\\.", "n"), ".csv")

	if (!file.exists(file_path)) {
		cmd <- paste0("julia ", scripts_path, "/plot_traj.jl network=", net, " noise-mode=", noiseMode,
			" param-id=", paramID, " sigma=", noiseLevel, " iter=1")
		system(cmd)
	}

	# Load and prep trajectory data
	d <- read_csv(file_path, show_col_types = FALSE)

	nodes <- colnames(d)
	nodes <- nodes[4:(ncol(d) - 1)]

	d <- d %>%
		mutate(across(all_of(nodes), .fns = function(x) log2(x)))

	# Compute discretization thresholds from RACIPE solutions
	solution_path <- paste0(data_path, "/", net, "_solution.dat")
	racip_sol <- read_solutions(solution_path)
	thresholds <- sapply(nodes, function(nd) {
		mean(racip_sol[[nd]])
	})

	# Discretize each node against its threshold
	disc_cols <- paste0(nodes, "d")

	d <- d %>%
		mutate(across(
			all_of(nodes),
			.fns = function(x) ifelse(x > thresholds[[cur_column()]], 1, 0),
			.names = "{.col}d"
		))

	# Build the combined state label, e.g. "(1, 0, 1)"
	d <- d %>%
		mutate(State = paste0("(", do.call(paste, c(select(d, all_of(disc_cols)), sep = ", ")), ")"))

	if (nTraj != 0) {
		ics <- unique(d$IC) %>% sample(nTraj)
		d <- d %>% filter(IC %in% ics)
	}

	states <- d$State %>% sort() %>% unique()

	# Fixed color per state, shared between the State trace panel and the MRT bar
	state_colors <- setNames(scales::hue_pal()(length(states)), states)

	# ------------------------------------------------------------------
	# Identify the burn-in window from Phase, for shading + MRT windowing
	# ------------------------------------------------------------------
	phase_bounds <- d %>%
		group_by(Phase) %>%
		summarise(start = min(Time), end = max(Time), .groups = "drop")

	# ------------------------------------------------------------------
	# The two Phase values are not "burn-in vs. analysis" — they are
	# "Deterministic" (shown to illustrate the states that exist without
	# noise; kept fully visible, never excluded) vs. "Stochastic" (where
	# noise is on). The actual burn-in to discard is the transient at the
	# *start* of the stochastic segment, matching cut_fraction = 0.5 used
	# in the Julia pipeline (keep the last 50% of the stochastic run --
	# src/StochasticSimulations.jl's cut_fraction default, never overridden
	# by the systematic analyze_noise_effects()/3_analyze_single_parameter.jl
	# call chain that produces all_parameters_results.csv).
	# ------------------------------------------------------------------
	det_phase   <- phase_bounds$Phase[which.min(phase_bounds$start)]

	stoch_phase <- setdiff(phase_bounds$Phase, det_phase)[1]
	stoch_start <- phase_bounds$start[phase_bounds$Phase == stoch_phase]
	stoch_end   <- phase_bounds$end[phase_bounds$Phase == stoch_phase]
    det_start <- phase_bounds$start[phase_bounds$Phase == det_phase]
	det_end   <- phase_bounds$end[phase_bounds$Phase == det_phase]
    if (det_time > det_end) det_time <- det_end
    if (stoch_time > (stoch_end - stoch_start)) stoch_time <- stoch_end - stoch_start
    

	stoch_burnin_end <- stoch_start + (1 - cut_fraction) * (stoch_end - stoch_start)

	# Only the stochastic segment gets a labeled time origin. The
	# deterministic run is shown purely for context (which states exist
	# without noise), not meant to be read off in absolute time — so its
	# clock stays unlabeled rather than getting its own axis or origin.
	d <- d %>% mutate(PlotTime = Time - stoch_start)

	# Rows actually used for MRT / switch-counting: stochastic phase only,
	# past its own transient
	analysis_rows <- d %>% filter(Phase == stoch_phase, Time >= stoch_burnin_end)

	# ------------------------------------------------------------------
	# Pick one representative trajectory to highlight: the IC whose
	# switch count *within the analysis window* is closest to the median
	# across all ICs, rather than an arbitrary or cherry-picked one.
	# ------------------------------------------------------------------
	switch_counts <- analysis_rows %>%
		arrange(IC, Time) %>%
		group_by(IC) %>%
		summarise(n_switch = sum(State != lag(State), na.rm = TRUE), .groups = "drop")

	median_switch <- median(switch_counts$n_switch)
	highlight_ic  <- switch_counts$IC[which.min(abs(switch_counts$n_switch - median_switch))]

	d_bg   <- d %>% filter(IC != highlight_ic)
	d_high <- d %>% filter(IC == highlight_ic)

	highlight_color <- "black"

	# ------------------------------------------------------------------
	# Panel 1: per-node expression traces, one trajectory highlighted
	# ------------------------------------------------------------------
	node_plots <- lapply(nodes, function(nd) {
		thr <- thresholds[[nd]]

		y_range <- range(d[[nd]], na.rm = TRUE)
		y_pad   <- diff(y_range) * 0.05
		burnin_rect <- data.frame(
			xmin = 0, xmax = stoch_burnin_end - stoch_start,
			ymin = y_range[1] - y_pad, ymax = y_range[2] + y_pad
		)
		lbl <- data.frame(x = 10, y = thr + 0.5,
			label = paste0("theta[", nd, "]"))

		ggplot() +
			# geom_rect_pattern(data = burnin_rect,
			# 	aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
			# 	inherit.aes = FALSE, fill = NA, color = "grey60", linewidth = 0.2,
			# 	pattern = "crosshatch", pattern_fill = NA, pattern_color = "grey60",
			# 	pattern_angle = 45, pattern_density = 0.02, pattern_spacing = 0.03,
			# 	pattern_linewidth = 0.15) +
            geom_vline(xintercept = 0, color = "red", linewidth = 0.4) +
			geom_line(data = d_bg, aes(x = PlotTime, y = .data[[nd]], group = IC),
				alpha = 0.15, linewidth = 0.2, color = "#999999") +
			geom_line(data = d_high, aes(x = PlotTime, y = .data[[nd]], group = IC),
				alpha = 1, linewidth = 0.8, color = highlight_color) +
			geom_hline(yintercept = thr, linetype = "dashed", color = "black",
				linewidth = 0.5, alpha = 0.7) +
			geom_text(data = lbl, aes(x = x, y = y, label = label),
				inherit.aes = FALSE, parse = TRUE, hjust = 1, size = 5) +
			scale_x_continuous(breaks = stoch_breaks) +
			labs(x = "Time", y = nd) +
			theme_Publication() +
			theme(legend.position = "none", axis.text = element_text(size = rel(1.4)))
	})

	# ------------------------------------------------------------------
	# Panel 2: State trace. Background ICs shown flat grey (de-emphasized);
	# the highlighted IC alone is colored by state, so its transitions
	# read clearly against a quiet background instead of competing with
	# every other trajectory's switches.
	# ------------------------------------------------------------------
	state_hlines <- data.frame(yintercept = seq_along(states))
	state_burnin_rect <- data.frame(
		xmin = 0, xmax = stoch_burnin_end - stoch_start,
		ymin = 0.5, ymax = length(states) + 0.5
	)

	p_state <- ggplot() +
		geom_rect_pattern(data = state_burnin_rect,
			aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
			inherit.aes = FALSE, fill = NA, color = "grey60", linewidth = 0.2,
			pattern = "crosshatch", pattern_fill = NA, pattern_color = "grey60",
			pattern_angle = 45, pattern_density = 0.02, pattern_spacing = 0.03,
			pattern_linewidth = 0.15) +
        geom_vline(xintercept = 0, color = "red", linewidth = 0.4) +
		geom_hline(data = state_hlines, aes(yintercept = yintercept),
			inherit.aes = FALSE, color = "grey85", linewidth = 0.3) +
		geom_line(data = d_bg, aes(x = PlotTime, y = match(State, states), group = IC),
			color = "#CCCCCC", linewidth = 0.3, alpha = 0.5) +
		geom_line(data = d_high, aes(x = PlotTime, y = match(State, states), color = State, group = IC),
			linewidth = 1.1, alpha = 1) +
		scale_color_manual(values = state_colors) +
		scale_y_continuous(breaks = seq_along(states), labels = states) +
		scale_x_continuous(breaks = stoch_breaks) +
		labs(x = "Time", y = "State") +
		theme_Publication() +
		theme(legend.position = "none", axis.text = element_text(size = rel(1.4)))

	# ------------------------------------------------------------------
	# Per-state residence-time summary -- population MRT (all ICs) vs. the
	# single highlighted trajectory's own residence time -- written
	# directly onto the state panel's bottom-right corner instead of a
	# separate colored bar chart below it (reported as a confusing extra
	# panel: two differently-colored bars conveying the same states again).
	# ------------------------------------------------------------------
	frac_by_state <- function(rows) {
		full <- setNames(rep(0, length(states)), states)
		tab  <- rows %>% count(State, name = "n") %>% mutate(frac = n / sum(n))
		full[tab$State] <- tab$frac
		full
	}
	mrt_vec  <- frac_by_state(analysis_rows)
	traj_vec <- frac_by_state(analysis_rows %>% filter(IC == highlight_ic))

	fmt_line <- function(term, vals) {
		paste0(term, ":  ", paste0(names(vals), "=", sprintf("%.2f", vals), collapse = "  "))
	}
	corner_text <- paste0(fmt_line("Traj. residence time", traj_vec), "\n",
						   fmt_line("MRT", mrt_vec))

	# p_state <- p_state +
	# 	annotate("label", x = Inf, y = -Inf, label = corner_text,
	# 			 hjust = 1.02, vjust = -0.15, size = 5.6, fontface = "bold",
	# 			 fill = "white", alpha = 0.85, label.size = 0)

	# ------------------------------------------------------------------
	# Assemble: node panels across the top, State panel below
	# ------------------------------------------------------------------
	top_row <- plot_grid(plotlist = node_plots, ncol = length(node_plots), align = "h", labels = top_labs, label_size = LS)

	p <- plot_grid(top_row, p_state, ncol = 1, rel_heights = c(2, 2), labels = c("", bot_labs), label_size = LS)

	if (!is.null(save_path)) {
		ggsave(save_path, plot = p, width = width, height = height, dpi = dpi)
	}

	p
}

# Shared data prep for mrt_plot() and the Friedman-results table (built
# separately in the script body below): load a network's results, filter
# to dt/param_type, zero-fill via fill_mrt(), and optionally keep only the
# topOriginal/topStates subset. Pulled out of mrt_plot() so the table can
# reuse the exact same filtering/zero-filling logic instead of duplicating it.
load_mrt_data <- function(net, noiseType, resultsFolder, dt = 0.01,
                           param_type = NULL, topOriginal = NULL, topStates = NULL) {
    results_file <- file.path(resultsFolder, noiseType, net,
                               "results", "all_parameters_results.csv")
    if (!file.exists(results_file)) stop("Results file not found: ", results_file)

    dAll <- read_csv(results_file, show_col_types = FALSE) %>% filter(DT == dt)

    if (!is.null(param_type)) {
        dAll <- dAll %>% filter(ParamType %in% param_type)
        if (nrow(dAll) == 0) stop("No rows remain after filtering param_type")
    }

    if (!is.null(topOriginal)) {
        keepIDs <- get_top_original(dAll, topOriginal)
        dAll <- dAll %>%
            semi_join(keepIDs, by = c("ParamType", "ParamID"))
    }

    # Zero-fill across every state currently present, BEFORE ranking topStates,
    # so the ranking isn't itself biased by the same missing-row problem.
    candidateStates <- unique(dAll$State)
    dAll <- fill_mrt(dAll, candidateStates)

    if (!is.null(topStates)) {
        keepStates <- dAll %>%
            group_by(ParamType, State) %>%
            summarise(MeanMRT = mean(MRT), .groups = "drop") %>%
            group_by(ParamType) %>%
            slice_max(MeanMRT, n = topStates, with_ties = FALSE) %>%
            ungroup()

        dAll <- dAll %>%
            semi_join(keepStates, by = c("ParamType", "State")) %>%
            mutate(State = factor(State,
                                   levels = keepStates %>%
                                       arrange(State) %>%
                                       pull(State) %>% unique()))
    }

    dAll %>%
        mutate(NoiseLevel = factor(NoiseLevel,
                                    levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
               ParamType  = str_to_sentence(ParamType))
}

# Single-network MRT summary (mean +/- SE by default), faceted by ParamType --
# distinct from the shared multi-network plot_mrt_trajectories()/mrt_plot() in
# figure_common.r / Figure3.r, which use SD errorbars and facet by StateClass.
#
# Significance (Friedman test, MRT vs. noise level) is NOT annotated on this
# plot -- collapsing it to one caption/bracket per state was either too small
# to read or hid the low-noise-vs-high-noise structure (see friedman_by_range()
# in figure_common.r and the combined results table built in the script body
# below, which replaces what used to be an in-plot caption here).
mrt_plot <- function(net, noiseType, resultsFolder, dt = 0.01,
                      group_var = "State", param_type = NULL,
                      facet_var = "ParamType",
                      style = c("errorbar", "boxplot", "violin"),
                      topOriginal = NULL, topStates = NULL,
                      y_label = "Average MRT", log_y = FALSE, angle_x = 60,
                      outFile = NULL, width = 6, height = 5) {
    style <- match.arg(style)

    dAll <- load_mrt_data(net, noiseType, resultsFolder, dt = dt,
                           param_type = param_type, topOriginal = topOriginal,
                           topStates = topStates)

    p <- ggplot(dAll, aes(x = NoiseLevel, y = MRT,
                          color = .data[[group_var]], fill = .data[[group_var]]))

    p <- p + switch(style,
        errorbar = list(
            stat_summary(fun.data = mean_sd, geom = "errorbar",
                         position = position_dodge(width = 0.6), width = 0.3),
            stat_summary(fun = mean, geom = "point",
                         position = position_dodge(width = 0.6), size = 2)
        ),
        boxplot = geom_boxplot(outlier.size = 0.5, alpha = 0.4,
                               position = position_dodge(width = 0.7)),
        violin  = geom_violin(alpha = 0.4, scale = "width",
                              position = position_dodge(width = 0.7))
    )

    if (!is.null(facet_var)) p <- p + facet_wrap(vars(.data[[facet_var]]))

    p <- p +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = angle_x, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = y_label)

    if (log_y) p <- p + scale_y_log10()

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }

    invisible(p)
}


# Per-network lambda violin, straight from all_parameters_lambdas.csv --
# distinct from the shared plot_lambda_density() in figure_common.r, which
# reads the flat histogram cache instead (not split by network).
plot_lambda_density_by_net <- function(net, noiseType, resultsFolder, noise_levels, dt,
                                 lambda_type = c("Inh", "Act"),
                                 y_label = expression(lambda), log_y = FALSE,
                                 angle_x = 60, outFile = NULL,
                                 width = 7, height = 5) {
    lambda_type <- match.arg(lambda_type)
    lambda_file <- file.path(resultsFolder, noiseType, net,
                               "results_det", "all_parameters_lambdas.csv")
    if (!file.exists(lambda_file)) stop("Results file not found: ", lambda_file)
    d <- read_csv(lambda_file, show_col_types = FALSE) %>%
        filter(DT %in% dt, NoiseLevel %in% noise_levels)

    lambdaCols <- names(d) %>% str_subset(paste0("^lambda_", lambda_type, "_of_"))
    if (length(lambdaCols) == 0) {
        stop("No lambda columns found for type = ", lambda_type)
    }

    dLong <- d %>%
        select(ParamID, NoiseLevel, all_of(lambdaCols)) %>%
        pivot_longer(cols = all_of(lambdaCols), names_to = "Edge", values_to = "Lambda") %>%
        mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(NoiseLevel))))

    p <- ggplot(dLong, aes(x = NoiseLevel, y = Lambda)) +
        geom_violin(fill = "steelblue", alpha = 0.5, scale = "width", color = NA) +
        # geom_boxplot(width = 0.08, outlier.shape = NA, alpha = 0,
        #              color = "black", linewidth = 0.4) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = angle_x, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = y_label)#,
            #  title = paste0(lambda_type, " \u03bb densities"))

    if (log_y) p <- p + scale_y_log10()

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }

    invisible(p)
}

p1 <- plot_burnin_trajectories("TS", 136, 0.01, "Additive", save_path = NULL)
p2 <- plot_mrt_trajectories(
    net = "TS", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    sampleN = 100,
    outFile = NULL, ncol = 1
) + labs(title = "TS - All parameters")
# Fixed State -> color mapping shared across every panel that uses it, so a
# shared legend is actually correct (not just visually tidy). Both TT and
# TTSA are 3-node networks with the same 8 possible discrete states, so
# this mapping is valid across both regardless of which subset topStates
# happens to keep for each network individually.
allStates   <- c("(0, 0, 0)", "(0, 0, 1)", "(0, 1, 0)", "(0, 1, 1)",
                 "(1, 0, 0)", "(1, 0, 1)", "(1, 1, 0)", "(1, 1, 1)")
stateColors <- setNames(scales::hue_pal()(length(allStates)), allStates)

allStatesTS   <- c("(0, 0)", "(0, 1)", "(1, 0)", "(1, 1)")
stateColorsTS <- setNames(scales::hue_pal()(length(allStatesTS)), allStatesTS)
p3_full <- mrt_plot(net = "TT", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    style = "errorbar", topStates = 7, facet_var = NULL, param_type = "bistable", outFile = NULL) +
    scale_color_manual(values = stateColors, breaks = names(stateColors)) +
    scale_fill_manual(values = stateColors, breaks = names(stateColors))

p4_full <- mrt_plot(net = "TTSA", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    style = "errorbar", topStates = 7, facet_var = NULL, param_type = "bistable", outFile = NULL) +
    scale_color_manual(values = stateColors, breaks = names(stateColors)) +
    scale_fill_manual(values = stateColors, breaks = names(stateColors)) +
    theme(legend.direction = "vertical")

p5_full <- mrt_plot(net = "TSSA", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    style = "errorbar", topStates = 4, facet_var = NULL, param_type = "bistable", outFile = NULL) +
    scale_color_manual(values = stateColorsTS, breaks = names(stateColorsTS)) +
    scale_fill_manual(values = stateColorsTS, breaks = names(stateColorsTS)) +
    theme(legend.direction = "vertical")

# legend <- cowplot::get_legend(p4_full)
# legendTS <- cowplot::get_legend(p5_full)

# Strip x-axis text/title from every panel except the bottom-most one (E),
# same convention a facet_wrap column would give you for free.
p3 <- p3_full + theme(legend.position = "right", legend.direction = "vertical",
                       axis.text.x = element_blank(), axis.title.x = element_blank()) + labs(title = "TT - Bistable parameters")
p4 <- p4_full + theme(legend.position = "right", legend.direction = "vertical",) + 
                        labs(title = "TTSA - Bistable parameters")
p5 <- p5_full + theme(legend.position = "right", legend.direction = "vertical",
                       axis.text.x = element_blank(), axis.title.x = element_blank()) + labs(title = "TSSA - Bistable parameters")
# p5 <- plot_lambda_density(net = "TS", noiseType = "Fluctuating", resultsFolder = resultsFolder,
#     noise_levels = c(0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0),
#     dt = 0.01, lambda_type = "Inh", outFile = NULL)

d_col  <- plot_grid(p5, p3, p4, ncol = 1, align = "v",
                     labels = c("D", "E", "F"), label_size = LS,
                     rel_heights = c(1, 1, 1.2))

# pright <- plot_grid(d_col, legend, ncol = 2, rel_widths = c(4, 1))
pright <- d_col

pAll <- plot_grid(p1, p2, pright, ncol = 3, align = "h",
                   labels = c("", "C", ""), label_size = LS,
                   rel_widths = c(2, 1, 1.3))
jpeg(file.path(finalDir, "Fig2.jpg"),
     width = if (exists("FIG_WIDTH")) FIG_WIDTH else 20,
     height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 12,
     units = "in", res = 300)
print(pAll)
dev.off()



p1 <- mrt_plot(net = "TS", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    style="errorbar",
    param_type = NULL, outFile = NULL) + labs(title = "TS")
p2 <- mrt_plot(net = "TSSA", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    style="errorbar",
    param_type = c("monostable", "tristable"), outFile = NULL) + labs(title = "TSSA")
p3 <- mrt_plot(net = "TT", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    style="errorbar",
    param_type = c("monostable", "tristable"), outFile = NULL) + labs(title = "TT")
p4 <- mrt_plot(net = "TTSA", noiseType = "Fluctuating", resultsFolder = resultsFolder, dt = 0.01,
    style="errorbar",
    param_type = c("monostable", "tristable"), outFile = NULL) + labs(title = "TTSA")

# ------------------------------------------------------------------
# Panels E(i)-(iv): one Friedman-test results table per network (TS,
# TSSA, TT, TTSA), replacing what used to be a per-panel caption on
# mrt_plot() (too small to read, and misleading once collapsed to one
# p-value across the whole 0-1 noise range -- see friedman_by_range() in
# figure_common.r). Each table has one row per (Stability, SigmaRange),
# with columns Stability | SigmaRange | FigurePanel | one per state.
#
# The Bistable class isn't covered by Fig2S1's own param_type filters
# (TSSA/TT/TTSA restrict to monostable+tristable there) -- it's what
# Fig2's D/E/F panels show instead (now filtered to param_type =
# "bistable" to match their titles). So each network's table combines
# both sources: Bistable rows come from the Fig2 D/E/F data, Monostable/
# Tristable rows from the Fig2S1 A-D data -- two disjoint stability
# coverages, not competing versions of the same thing, so both go in.
# ------------------------------------------------------------------

# One row per (Stability, SigmaRange) for a single source dataset,
# States pivoted out to columns (rather than kept as their own long-
# format column) since the table is now organized around Stability/
# SigmaRange as the row identity.
stability_sigma_rows <- function(dSub, panelLabel, low_thresh) {
    stabilities <- unique(as.character(dSub$ParamType))
    lapply(stabilities, function(st) {
        sub <- dSub %>% filter(ParamType == st) %>% arrange(State)
        res <- friedman_by_range(sub, group_var = "State", low_thresh = low_thresh)
        lowWide  <- res %>% select(State, stars_low)  %>% pivot_wider(names_from = State, values_from = stars_low)
        highWide <- res %>% select(State, stars_high) %>% pivot_wider(names_from = State, values_from = stars_high)
        bind_rows(
            bind_cols(data.frame(Stability = st, SigmaRange = paste0("<=", low_thresh),
                                  FigurePanel = panelLabel, check.names = FALSE), lowWide),
            bind_cols(data.frame(Stability = st, SigmaRange = paste0(">", low_thresh),
                                  FigurePanel = panelLabel, check.names = FALSE), highWide)
        )
    }) %>% bind_rows()
}

# Builds one network's full table: Bistable (from the Fig2 D/E/F-style
# call, no topStates truncation -- the table wants the complete state
# set, not the plot's decluttered top-N subset) unioned with Monostable/
# Tristable (from the Fig2S1-style call).
build_network_sig_table <- function(net, low_thresh, bistable_panel, other_panel, other_types) {
    dBi    <- load_mrt_data(net, "Fluctuating", resultsFolder, dt = 0.01, param_type = "bistable")
    dOther <- load_mrt_data(net, "Fluctuating", resultsFolder, dt = 0.01, param_type = other_types)
    bind_rows(
        stability_sigma_rows(dBi, bistable_panel, low_thresh),
        stability_sigma_rows(dOther, other_panel, low_thresh)
    )
}

# TS never appears in Fig2 D/E/F -- its Fig2S1 A call (param_type = NULL)
# already covers whatever stability classes exist for it (Bistable,
# Monostable) in one dataset, so it doesn't need the bistable/other split.
dTS       <- load_mrt_data("TS", "Fluctuating", resultsFolder, dt = 0.01, param_type = NULL)
sigTableTS   <- stability_sigma_rows(dTS, "Fig2S1 A", low_thresh = 0.1)
sigTableTSSA <- build_network_sig_table("TSSA", low_thresh = 0.1,
                    bistable_panel = "Fig2 D", other_panel = "Fig2S1 B", other_types = c("monostable", "tristable"))
sigTableTT   <- build_network_sig_table("TT", low_thresh = 0.01,
                    bistable_panel = "Fig2 E", other_panel = "Fig2S1 C", other_types = c("monostable", "tristable"))
sigTableTTSA <- build_network_sig_table("TTSA", low_thresh = 0.01,
                    bistable_panel = "Fig2 F", other_panel = "Fig2S1 D", other_types = c("monostable", "tristable"))

sigTheme <- ttheme_minimal(base_size = 11, core = list(padding = unit(c(4, 3), "mm")))
grobTS   <- tableGrob(sigTableTS,   rows = NULL, theme = sigTheme)
grobTSSA <- tableGrob(sigTableTSSA, rows = NULL, theme = sigTheme)
grobTT   <- tableGrob(sigTableTT,   rows = NULL, theme = sigTheme)
grobTTSA <- tableGrob(sigTableTTSA, rows = NULL, theme = sigTheme)

pleft <- plot_grid(p1, p2, ncol = 1, labels = c("A", "B"), label_size = LS)
pright <- plot_grid(p3, p4, ncol = 1, labels = c("C", "D"), label_size = LS)
pTop <- plot_grid(pleft, pright, nrow = 1, labels = c("", ""), label_size = LS)

# Each table needs real vertical space of its own (header + one row per
# Stability x SigmaRange combo) rather than an arbitrary fixed fraction.
# Mirror A-D's two-column layout: TS/TSSA tables stacked under the left
# column, TT/TTSA stacked under the right column.
rowHeight  <- 0.16
tsHeight   <- 0.6 + nrow(sigTableTS)   * rowHeight
tssaHeight <- 0.6 + nrow(sigTableTSSA) * rowHeight
ttHeight   <- 0.6 + nrow(sigTableTT)   * rowHeight
ttsaHeight <- 0.6 + nrow(sigTableTTSA) * rowHeight

eLeft  <- plot_grid(grobTS, grobTSSA, ncol = 1, labels = c("E(i)", "E(ii)"),
                     label_size = LS, rel_heights = c(tsHeight, tssaHeight))
eRight <- plot_grid(grobTT, grobTTSA, ncol = 1, labels = c("E(iii)", "E(iv)"),
                     label_size = LS, rel_heights = c(ttHeight, ttsaHeight))
eRow <- plot_grid(eLeft, eRight, nrow = 1)

gridHeight  <- if (exists("FIG_HEIGHT")) FIG_HEIGHT else 12
tableHeight <- max(tsHeight + tssaHeight, ttHeight + ttsaHeight)
pAll <- plot_grid(pTop, eRow, ncol = 1, rel_heights = c(gridHeight, tableHeight))
jpeg(file.path(finalDir, "Fig2S1.jpg"),
     width = if (exists("FIG_WIDTH")) FIG_WIDTH else 20,
     height = gridHeight + tableHeight,
     units = "in", res = 300)
print(pAll)
dev.off()
# parD <- data.frame(
#     Network = c(rep("TS", 6)), 
#     ParamIDs = c(77, 77, 136, 4367, 1914, 582), 
#     Sigmas = c(0.5, 0.1, 0.01, 0.1, 0.005, 0.5),
#     NoiseTypes = c("Fluctuating", "Additive", "Additive", "Additive", "Additive", "Multiplicative")
# )

# apply(parD, 1, function(x) {
#     net <- x["Network"]
#     parID <- x["ParamIDs"] %>% as.integer
#     noise <- x["Sigmas"]
#     noiseLab <- noise %>% str_replace("\\.", "n") 
#     noise <- noise %>% as.numeric()
#     nT <- x["NoiseTypes"]
#     nL <- substr(nT, 1, 1)
#     savePath <- paste0("figures/", net, "/", parID, "_", noiseLab, "_", nL, ".jpg")
#     plot_burnin_trajectories(net, parID, noise, nT, save_path = savePath, width = 9, height = 8.5)
# })

# plot_mrt_trajectories(
#     net = "TS", noiseType = "Fluctuating", resultsFolder = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final", dt = 0.01,
#     outFile = file.path("figures", "individual", "Fluctuating_TS_topStates_traj.jpg"), ncol = 1
# )
# mrt_plot(net = "TS", noiseType = "Fluctuating", resultsFolder = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final", dt = 0.01,
#     style="errorbar",  
#     facet_var = NULL, param_type = "bistable", outFile = file.path("figures", "individual", "Fluctuating_TS_topStates_MRT_bistable.jpg"), width = 7, height = 7.5)
# mrt_plot(net = "TS", noiseType = "Fluctuating", resultsFolder = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final", dt = 0.01,
#     style="errorbar",  
#     param_type = NULL, outFile = file.path("figures", "individual", "Fluctuating_TS_topStates_MRT_all.jpg"), width = 13, height = 7.5)
# mrt_plot(net = "TSSA", noiseType = "Fluctuating", resultsFolder = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final", dt = 0.01,
#     style="errorbar",  
#     param_type = NULL, outFile = file.path("figures", "individual", "Fluctuating_TSSA_topStates_MRT_all.jpg"), width = 13, height = 7.5)
# mrt_plot(net = "TT", noiseType = "Fluctuating", resultsFolder = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final", dt = 0.01,
#     style="errorbar", topStates = 7,
#     facet_var = NULL, param_type = "bistable", outFile = file.path("figures", "individual", "Fluctuating_TT_topStates_MRT_bistable.jpg"), width = 6, height = 5.5)
# mrt_plot(net = "TTSA", noiseType = "Fluctuating", resultsFolder = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final", dt = 0.01,
#     style="errorbar", topStates = 7, 
#     facet_var = NULL, param_type = "bistable", outFile = file.path("figures", "individual", "Fluctuating_TTSA_topStates_MRT_bistable.jpg"), width = 6, height = 5.5)
# plot_lambda_density(net = "TS", noiseType = "Fluctuating", resultsFolder = "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final",
#                      noise_levels = c(0, 0.01, 0.1, 0.5, 1.0),
#                      dt = 0.01, lambda_type = "Inh", 
#                      outFile = file.path("figures", "individual", "Fluctuating_Inh_lambdas.jpg"), width = 6, height = 5.5)
