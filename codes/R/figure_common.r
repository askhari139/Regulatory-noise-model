# figure_common.r
#
# Shared helpers for the Figure2-5 scripts: MRT zero-filling / state
# classification, reachability/margin computation, the deterministic- vs.
# stochastic-MRT comparison loader, and the two plotting functions
# (MRT trajectories, lambda density) reused across multiple figures.
#
# Each figure script sources this after its own library() calls:
#   source("/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/scripts/figure_common.r")
# ============================================================

resultsFolder <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final"
dataFolder    <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/data"
# LS (plot_grid label size) is left alone if a caller (run_all_figures.r) has
# already set it -- lets the master script pick a per-figure LS while each
# Figure*.r script still gets a sane default when run standalone.
if (!exists("LS")) LS <- 20

# Mean +/- SD, for stat_summary(fun.data = mean_sd) errorbars (as opposed to
# stat_summary(fun.data = mean_se), which some panels use deliberately instead).
mean_sd <- function(x) {
    data.frame(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
}

# Zero-fills every (ParamID, ParamType, NoiseLevel, State) combo for the
# states actually present in df. Any ParamID that never visits a given state
# at a given noise level gets MRT = 0 for that combo instead of being
# silently absent -- otherwise every downstream mean (ranking, plotting,
# whatever) is computed over a shifting, state-dependent subset of ParamIDs,
# which biases everything. Also attaches:
#   - StateClass: single-high / double-high / all-high / multi-high, by
#     count of high ("1") nodes in the state label
#   - Original: the compound deterministic attractor label per ParamID
#     (every state visited at NoiseLevel == 0, sorted and concatenated)
fill_mrt <- function(df, states) {
    dOrig <- df %>% filter(NoiseLevel == 0) %>% group_by(ParamID) %>%
        summarise(Original = State %>% unique %>% sort %>% paste0(collapse = "_"))

    full_grid <- df %>%
        distinct(ParamID, ParamType, NoiseLevel) %>%
        tidyr::crossing(State = states)
    n_nodes <- str_count(df$State[1], ",") + 1

    df <- df %>%
        select(ParamID, ParamType, NoiseLevel, State, MRT) %>%
        right_join(full_grid, by = c("ParamID", "ParamType", "NoiseLevel", "State")) %>%
        mutate(MRT = tidyr::replace_na(MRT, 0)) %>%
        mutate(Expression = str_count(State, "1")) %>%
        mutate(StateClass = case_when(
            Expression == 0        ~ "all-low",
            Expression == n_nodes  ~ "all-high",
            Expression == 1        ~ "single-high",
            Expression == 2        ~ "double-high",
            .default               = "multi-high"
        ))

    full_join(df, dOrig, by = "ParamID")
}

# Compound deterministic attractor labels per ParamID (combine every state
# visited at NoiseLevel == 0), ranked by frequency within each ParamType.
# Used to keep only the n most common deterministic outcomes.
get_top_original <- function(df, n) {
    origMap <- df %>%
        filter(NoiseLevel == 0) %>%
        group_by(ParamType, ParamID) %>%
        summarise(Original = State %>% unique() %>% sort() %>% paste0(collapse = "_"),
                  .groups = "drop")

    keepOriginal <- origMap %>%
        count(ParamType, Original, name = "Freq") %>%
        group_by(ParamType) %>%
        slice_max(Freq, n = n, with_ties = FALSE) %>%
        ungroup()

    origMap %>%
        inner_join(keepOriginal %>% select(ParamType, Original),
                   by = c("ParamType", "Original"))
}

# Networks with an activatory lambda (self-activation loops, or -- for DA --
# mutual activation edges) that have a MultiplicativeInvLambda counterpart
# dataset. Same list as scripts/FigureInvLambdaComparison.r.
sa_networks <- c("DA", "TSSA", "TTSA", "TT4SA", "NFSA")

# Redirects a network to its MultiplicativeInvLambda data when the caller
# asked for "Multiplicative" and the network actually has an invertible
# activatory lambda; every other network/noiseType passes through unchanged.
# Labels/titles built from the original `noiseType` argument at the call
# site are unaffected -- this only changes which data folder gets read.
noise_type_for <- function(net, noiseType) {
    if (noiseType == "Multiplicative" && net %in% sa_networks) "MultiplicativeInvLambda" else noiseType
}

# Adds per-node reachability (ceiling > threshold) and margin (ceiling -
# threshold) to a results file. Threshold = mean deterministic expression
# per node, from the RACIPE solution file; ceiling = Prod/Deg per node,
# from the parameter file.
#
# The solution file's node columns are LOG2 expression (RACIPE's native
# convention -- values go negative, confirmed against the raw .dat and
# against src/RACIPEdata.jl's read_solutions() docstring: "4+: Log2
# expression values for each node"), while ceiling (Prod/Deg) is a raw
# linear-scale steady-state concentration. Comparing a linear ceiling
# against a raw log2 mean is comparing different units -- the threshold
# needs converting back to linear first. This mirrors exactly what the
# Julia pipeline's own get_mean_expression() does (src/RACIPEdata.jl):
# a basin-frequency-weighted mean in log2 space, then 2^(.) to return to
# linear units -- basin = % of initial conditions converging to that
# solution row, the same weighting Julia's racipe_thresholds already use
# for state discretization/attractor ID. This function previously used a
# plain, unweighted mean(solutions[[node]]) with no 2^ conversion at all,
# which both used the wrong scale (log2 instead of linear) and ignored
# basin weighting -- caught because a plotted threshold line landed at
# ~4 instead of the expected ~16-19 for TS.
add_reachability <- function(net, noiseType, dt, dataFolder, resultsFolder) {
    results_file <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_results.csv")
    results_df <- read_csv(results_file, show_col_types = FALSE) %>% filter(DT == dt)
    # The raw CSV only ever has a row for a (ParamID, State) combo that was
    # actually visited -- a ParamID that never reaches a given state (e.g.
    # because it fails the reachability check below for one node) has ZERO
    # rows there, not an MRT=0 row. Any caller filtering this by State ==
    # would then silently drop that ParamID instead of correctly showing it
    # at MRT=0. Zero-fill before any State-based use downstream.
    #
    # fill_mrt() drops every column except ParamID/ParamType/NoiseLevel/
    # State/MRT (and the Original/StateClass columns it adds) -- including
    # DT -- so DT must be filtered to a single value BEFORE calling it,
    # not after (previously `dt` was accepted but never actually used here;
    # every caller filtered DT downstream instead, which no longer works
    # once DT no longer exists as a column post-fill_mrt).
    results_df <- fill_mrt(results_df, unique(results_df$State))

    parameters <- read_parameters(file.path(dataFolder, paste0(net, "_parameters.dat")))
    solutions  <- read_solutions(file.path(dataFolder, paste0(net, "_solution.dat")))
    nodes <- colnames(parameters)
    nodes <- nodes[str_detect(nodes, "Prod_of")] %>% str_remove("Prod_of_")

    w <- solutions$basin / 100
    means <- sapply(nodes, function(node) 2^(sum(solutions[[node]] * w) / sum(w)))
    names(means) <- nodes

    reaches <- parameters %>% select(ParamID)
    for (node in nodes) {
        node_max <- parameters[[paste0("Prod_of_", node)]] / parameters[[paste0("Deg_of_", node)]]
        reaches[[paste0(node, "_reach")]] <- as.integer(node_max > means[node])
        reaches[[paste0(node, "_diff")]]  <- node_max - means[node]
        reaches[[paste0(node, "_max")]]   <- node_max
    }

    results_df %>% left_join(reaches, by = "ParamID")
}

# Deterministic (resampled-lambda) vs. stochastic MRT, paramwise, for one
# focal state (focal_state = "auto" resolves to the all-high state for that
# network's node count) or one focal StateClass (MRT summed across every
# state in the class). Provide exactly one of focal_state / focal_class.
get_det_stoch_comparison_data <- function(nets, noiseType, dt, dataFolder, resultsFolder,
                                            focal_state = NULL, focal_class = NULL,
                                            sigma = NULL) {

    if (is.null(focal_state) && is.null(focal_class))
        stop("Provide exactly one of focal_state or focal_class")

    map_dfr(nets, function(net) {
        effectiveNoiseType <- noise_type_for(net, noiseType)
        stoch_file <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
        det_file   <- file.path(resultsFolder, effectiveNoiseType, net, "results_det", "all_parameters_results.csv")
        if (!file.exists(stoch_file) || !file.exists(det_file)) {
            warning("Missing stoch or det-resampled results, skipping: ", net)
            return(NULL)
        }

        st <- if (!is.null(focal_state) && identical(focal_state, "auto")) {
            parameters <- read_parameters(paste0(dataFolder, "/", net, "_parameters.dat"))
            n_nodes <- sum(str_detect(colnames(parameters), "^Prod_of_"))
            paste0("(", paste(rep(1, n_nodes), collapse = ", "), ")")
        } else focal_state

        dStoch <- read_csv(stoch_file, show_col_types = FALSE) %>% filter(DT == dt)
        dDet   <- read_csv(det_file,   show_col_types = FALSE) %>% filter(DT == dt)

        candidateStates <- unique(c(dStoch$State, dDet$State))
        dStoch <- fill_mrt(dStoch, candidateStates)
        dDet   <- fill_mrt(dDet, candidateStates)

        collapse <- function(df, valcol) {
            if (!is.null(focal_class)) {
                df %>%
                    filter(StateClass == focal_class) %>%
                    group_by(ParamID, ParamType, NoiseLevel) %>%
                    summarise(!!valcol := sum(MRT), .groups = "drop")
            } else {
                df %>%
                    filter(State == st) %>%
                    group_by(ParamID, ParamType, NoiseLevel) %>%
                    summarise(!!valcol := sum(MRT), .groups = "drop")
            }
        }

        sStoch <- collapse(dStoch, "MRT")
        sDet   <- collapse(dDet,   "MRT_det")

        # MRT_orig: sigma = 0 row of det-resampled data = original point-lambda
        # deterministic outcome (no resampling occurs at sigma = 0).
        dOrig <- sDet %>% filter(NoiseLevel == 0) %>% select(ParamID, MRT_orig = MRT_det)

        full_join(sStoch, sDet, by = c("ParamID", "ParamType", "NoiseLevel")) %>%
            mutate(MRT     = tidyr::replace_na(MRT, 0),
                   MRT_det = tidyr::replace_na(MRT_det, 0)) %>%
            left_join(dOrig, by = "ParamID") %>%
            mutate(MRT_orig = tidyr::replace_na(MRT_orig, 0),
                   Network  = net) %>%
            { if (!is.null(sigma)) filter(., NoiseLevel %in% sigma) else . }
    })
}

# MRT-vs-noise trajectories, one line per ParamID, faceted by State (or by
# StateClass, with topStates = NULL, to sum MRT within each class first).
# Optionally colored by reachability-margin quartile (marginData) instead of
# plotted as a flat grey background.
plot_mrt_trajectories <- function(net, noiseType, resultsFolder, dt = 0.01,
                                   param_type = NULL, facet_var = "State",
                                   topOriginal = NULL, topStates = 3,
                                   stateClasses = NULL, marginData = NULL,
                                   meta_cols = NULL, outFile = NULL,
                                   ncol = 3, panelWidth = 7.5, panelHeight = 6,
                                   sampleN = NULL) {

    results_file <- file.path(resultsFolder, noiseType, net,
                               "results", "all_parameters_results.csv")
    if (!file.exists(results_file)) stop("Results file not found: ", results_file)

    df_raw <- read_csv(results_file, show_col_types = FALSE) %>%
        filter(DT == dt) %>%
        mutate(NoiseLevel = as.numeric(as.character(NoiseLevel)))

    if (!is.null(param_type)) {
        df_raw <- df_raw %>% filter(ParamType %in% param_type)
        if (nrow(df_raw) == 0) stop("No rows remain after filtering param_type")
    }

    if (!is.null(topOriginal)) {
        keepIDs <- get_top_original(df_raw, topOriginal)
        df_raw <- df_raw %>% semi_join(keepIDs, by = c("ParamType", "ParamID"))
    }

    candidateStates <- unique(df_raw$State)
    df <- fill_mrt(df_raw, candidateStates)

    if (!is.null(stateClasses)) {
        df <- df %>% filter(StateClass %in% stateClasses)
        if (nrow(df) == 0) stop("No rows remain after filtering stateClasses")

        if (facet_var == "StateClass") {
            groupCols <- intersect(c("ParamID", "ParamType", "NoiseLevel", "StateClass",
                                      "Original", meta_cols), names(df))
            df <- df %>%
                group_by(across(all_of(groupCols))) %>%
                summarise(MRT = sum(MRT), .groups = "drop") %>%
                mutate(State = StateClass)
        }
    }

    if (!is.null(topStates)) {
        keepStates <- df %>%
            group_by(State) %>%
            summarise(MeanMRT = mean(MRT), .groups = "drop") %>%
            slice_max(MeanMRT, n = topStates, with_ties = FALSE)

        df <- df %>%
            semi_join(keepStates, by = "State") %>%
            mutate(State = factor(State,
                                   levels = keepStates %>%
                                       arrange(desc(MeanMRT)) %>%
                                       pull(State)))
    }

    df <- df %>%
        mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(NoiseLevel))))

    if (!is.null(meta_cols)) {
        meta <- df_raw %>% distinct(ParamID, across(all_of(meta_cols)))
        df <- df %>% left_join(meta, by = "ParamID")
    }

    if (!is.null(marginData)) {
        df <- df %>% left_join(marginData %>% select(ParamID, MinMargin), by = "ParamID")
        breaks <- quantile(marginData$MinMargin, probs = seq(0, 1, 0.25), na.rm = TRUE)
        df <- df %>%
            mutate(MarginBin = cut(MinMargin, breaks = breaks, include.lowest = TRUE, dig.lab = 3))
    }

    # Thin per-ParamID lines are drawn from a (optionally sampled) subset,
    # but every stat_summary() mean/aggregate below still inherits the full
    # df set at the ggplot() call -- so the sample only thins the visual
    # spaghetti, it never changes what the mean line represents.
    dfLines <- df
    if (!is.null(sampleN)) {
        set.seed(1)
        idPool <- df %>% distinct(ParamID)
        keepIDs <- idPool %>% slice_sample(n = min(sampleN, nrow(idPool)))
        dfLines <- df %>% semi_join(keepIDs, by = "ParamID")
    }

    p <- ggplot(df, aes(x = NoiseLevel, y = MRT))

    if (!is.null(marginData)) {
        p <- p +
            geom_line(data = dfLines, aes(group = ParamID, color = MarginBin), alpha = 0.25, linewidth = 0.3) +
            stat_summary(aes(group = MarginBin, color = MarginBin), fun = mean, geom = "line",
                        linewidth = 1.3, linetype = "dashed") +
            scale_color_viridis_d(name = "Margin (Max - Threshold)")
    } else {
        p <- p +
            geom_line(data = dfLines, aes(group = ParamID), alpha = 0.15, linewidth = 0.3) +
            stat_summary(aes(group = 1), fun = mean, geom = "line",
                        linewidth = 1.3, linetype = "dashed", color = "black")
    }

    p <- p +
        facet_wrap(vars(.data[[facet_var]]), ncol = ncol, scales = "free_y") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = "MRT")

    if (!is.null(outFile)) {
        n_panels <- length(unique(df[[facet_var]]))
        this_ncol <- min(n_panels, ncol)
        this_nrow <- ceiling(n_panels / ncol)
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = panelWidth * this_ncol, height = panelHeight * this_nrow)
    }

    invisible(p)
}

# Violin of lambda distribution (Inh or Act edge multiplier) across noise
# levels, from the pre-computed flat histogram cache (one CSV per dt).
plot_lambda_density <- function(noiseType, dataFolder, noise_levels, dt,
                                 lambda_type = c("Inh", "Act"),
                                 y_label = expression(lambda), log_y = FALSE,
                                 angle_x = 60, outFile = NULL,
                                 width = 7, height = 5) {
    # several.ok lets the caller pass both types at once to get one faceted
    # plot back (see below) instead of assembling two separate plots with
    # plot_grid() -- avoids the duplicated x-axis text/title that stacking
    # two full standalone plots produces.
    lambda_type <- match.arg(lambda_type, several.ok = TRUE)

    d <- map_dfr(lambda_type, function(this_type) {
        ltype_key <- tolower(this_type)   # "inh" or "act" -- matches LambdaType column
        if (ltype_key == "inh") {
            # Inh_of_AToB is naturally in [0, 1] in steps of 0.01, matching the
            # histogram cache's "inh" LambdaInit grid at that resolution.
            p_d <- read_parameters(file.path(dataFolder, "TS_parameters.dat"))
            lambdas <- p_d %>% pull(Inh_of_AToB) %>% round(2)
        } else {
            # Act_of_AToA is continuous over ~[1, 100], but the cache's "act"
            # LambdaInit grid is only the 100 whole numbers 1..100 -- round to
            # the nearest integer (not 2 decimals) so the join below actually
            # matches instead of missing almost everything.
            p_d <- read_parameters(file.path(dataFolder, "TSSA_parameters.dat"))
            lambdas <- p_d %>% pull(Act_of_AToA) %>% round(0)
        }
        lambdas_df <- data.frame(LambdaInit = lambdas %>% sample(100))

        # "Act" (activatory) lambda has an InvLambda counterpart cache;
        # "Inh" (inhibitory) doesn't -- lambda_sampler.jl's InvLambda rule
        # only kicks in for max_lambda > 1 (activation-range edges), falling
        # back to plain multiplicative noise for inhibition-range edges.
        histNoiseType <- if (ltype_key == "act") noise_type_for("TSSA", noiseType) else noiseType

        # One CSV per dt (histogram cache isn't split by network -- see note below)
        map_dfr(dt, function(this_dt) {
            hist_file <- file.path(dataFolder,
                                    paste0("lambda_hist_flat_", histNoiseType, "_dt", this_dt, ".csv"))
            if (!file.exists(hist_file)) {
                warning("Histogram file not found, skipping: ", hist_file)
                return(NULL)
            }
            read_csv(hist_file, show_col_types = FALSE) %>%
                mutate(DT = this_dt) %>% right_join(lambdas_df, by = "LambdaInit")
        }) %>%
            filter(LambdaType == ltype_key)
    })
    if (nrow(d) == 0) stop("No histogram files found for dt = ", paste(dt, collapse = ", "))

    d <- d %>%
        filter(NoiseLevel %in% noise_levels, Count > 0) %>%
        mutate(NoiseLevel  = factor(NoiseLevel, levels = sort(unique(NoiseLevel))),
               LambdaType  = factor(LambdaType, levels = c("inh", "act"),
                                     labels = c("lambda[Inh]", "lambda[Act]")))

    p <- ggplot(d, aes(x = NoiseLevel, y = BinMid, weight = Count)) +
        geom_violin(fill = "steelblue", alpha = 0.5, scale = "width", color = NA) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = angle_x, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = y_label)

    if (length(lambda_type) > 1) {
        # Facet instead of the y-axis title carrying a single label -- the
        # strip text (parsed as plotmath) distinguishes Inh/Act instead, and
        # the shared, non-free x axis means its text/title is drawn once
        # (bottom row only), not once per type.
        p <- p +
            facet_wrap(vars(LambdaType), ncol = 1, scales = "free_y", labeller = label_parsed) +
            labs(y = NULL)
    }

    if (log_y) p <- p + scale_y_log10()

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }

    invisible(p)
}

# ── Det/stoch paramwise + network-comparison panels (Figure 4 + its
#    supplement) -- these two need library(ggsignif) loaded by the caller ──

stars_from_p <- function(p) {
    case_when(
        is.na(p)  ~ "ns",
        p < 0.05  ~ "*",
        TRUE      ~ "ns"
    )
}

# Friedman test (MRT vs. NoiseLevel, blocked by ParamID -- fill_mrt()
# already zero-fills every ParamID at every NoiseLevel, so the design is
# complete/balanced), run separately on the low (<= low_thresh) and high
# (>= low_thresh) NoiseLevel ranges instead of once across the whole range.
# One row per group_var level, with both ranges' p-value/stars side by
# side -- used to build a results table (Figure2.r) rather than an
# in-plot annotation, so there's no need to collapse per-state results
# into a single label.
friedman_by_range <- function(dAll, group_var = "State", low_thresh = 0.1) {
    dAll <- dAll %>% mutate(NoiseLevelNum = as.numeric(as.character(NoiseLevel)))

    one_range <- function(sub) {
        sub %>%
            group_by(.data[[group_var]]) %>%
            summarise(p = if (length(unique(NoiseLevelNum)) >= 2) tryCatch(
                          stats::friedman.test(y = MRT, groups = factor(NoiseLevel),
                                                blocks = factor(ParamID))$p.value,
                          error = function(e) NA_real_
                      ) else NA_real_,
                      .groups = "drop") %>%
            mutate(stars = stars_from_p(p))
    }

    low  <- one_range(dAll %>% filter(NoiseLevelNum <= low_thresh)) %>%
        rename(stars_low = stars, p_low = p)
    high <- one_range(dAll %>% filter(NoiseLevelNum >= low_thresh)) %>%
        rename(stars_high = stars, p_high = p)

    full_join(low, high, by = group_var)
}

# Paired Wilcoxon signed-rank test. A mean-based paired z/t-test assumes the
# sampling distribution of the mean difference is ~normal; checked directly
# against this data (get_det_stoch_comparison_data() output, sigma = 0.1):
# group sizes here are only ~100-200 ParamIDs, and several (Network,
# StabilityClass) groups have visibly skewed differences (skewness up to
# ~1.8, e.g. TT/Multistable) -- not surprising since MRT is bounded in
# [0, 1] and piles up near the boundaries. At that n/skew combination a
# z-test (which additionally assumes the *population* SD is known, not
# estimated) measurably overstates significance relative to a t-test in the
# tails -- and the nonparametric, rank-based Wilcoxon test sidesteps the
# normality assumption entirely rather than just approximating around it.
paired_wilcox_p <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < 2 || all(x == y)) return(NA_real_)
    suppressWarnings(wilcox.test(x, y, paired = TRUE)$p.value)
}

# Paramwise scatter of one MRT series against another, at a fixed sigma,
# faceted by network: compare = "resample" is MRT_orig vs. MRT_det (the
# effect of resampling lambda alone), "stochastic" is MRT_det vs. MRT (the
# effect of adding noise on top of the resampled-lambda deterministic run).
# Previously colored by a reachability-margin statistic (Max - Threshold
# across nodes); dropped -- across every noise level and margin_stat tried,
# the color was too uniform within a panel to convey anything.
plot_paramwise_comparison <- function(nets, noiseType, dt, dataFolder, resultsFolder,
                                       compare = c("resample", "stochastic"),
                                       state = "auto", sigma = 0.1, param_type = NULL,
                                       facet_by = "Network",
                                       outFile = NULL, panelWidth = 3.5, panelHeight = 3.5) {

	compare <- match.arg(compare)
	d <- get_det_stoch_comparison_data(nets, noiseType, dt, dataFolder, resultsFolder,
										state, sigma = sigma)
	if (nrow(d) == 0) stop("No data found")
	if (!is.null(param_type)) d <- d %>% filter(ParamType %in% param_type)

	if (compare == "resample") {
		xcol <- "MRT_orig"; ycol <- "MRT_det"
		xlab <- "Original (point λ) MRT"; ylab <- "Deterministic (resampled λ) MRT"
	} else {
		xcol <- "MRT_det"; ycol <- "MRT"
		xlab <- "Deterministic (resampled λ) MRT"; ylab <- "Stochastic MRT"
	}

	d <- d %>%
		mutate(Network = factor(Network, levels = nets),
			   ParamType = str_to_sentence(ParamType))

	p <- ggplot(d, aes(x = .data[[xcol]], y = .data[[ycol]])) +
		geom_point(alpha = 0.5, size = 1.3, color = "#e41a1c") +
		geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.6) +
		theme_Publication() +
		labs(x = xlab, y = ylab,
			 title = paste0(if (compare == "resample") "Resampling effect" else "Stochastic effect",
							"  (σ = ", sigma, ")"))

	if (!is.null(facet_by)) p <- p + facet_wrap(vars(.data[[facet_by]]))

	if (!is.null(outFile)) {
		n_panels <- if (!is.null(facet_by)) length(unique(d[[facet_by]])) else 1
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = panelWidth * min(n_panels, 3),
			   height = panelHeight * ceiling(n_panels / 3))
	}
	invisible(p)
}

# All networks compared at a fixed sigma, all three MRT series (Original
# point-lambda / Deterministic resampled-lambda / Stochastic), split by
# stability class (mono- vs multistable) -- with paired Wilcoxon test significance
# brackets, paired by ParamID within each (Network, StabilityClass):
# "Resampled" (MRT_det) vs. stochastic (MRT), and resampled vs. original
# point-lambda (MRT_orig) -- the same two effects isolated separately by
# plot_paramwise_comparison(), now tested for significance per network.
# stability = NULL shows both classes as facets (the main-figure panel);
# stability = "Monostable" or "Multistable" shows just that one (for the
# supplementary figure, where each class gets its own labeled panel).
# state_class, if given (one of fill_mrt()'s StateClass values -- "all-high",
# "single-high", "double-high"), sums MRT across every state in that class
# instead of using the single focal state (default: the all-high state, as
# in the main figure) -- for 2-node networks "double-high" is empty (their
# 2-high states are already "all-high"), so those facets come back blank,
# which is expected rather than an error.
plot_det_stoch_network_comparison <- function(nets, noiseType, dt, dataFolder, resultsFolder,
                                                state = "auto", state_class = NULL, sigma,
                                                param_type = NULL,
                                                stability = NULL, dodge_width = 0.6,
                                                outFile = NULL, width = 10, height = 5) {

	typeLevels <- c("Original (point λ)", "Deterministic (resampled λ)", "Stochastic")

	d <- if (!is.null(state_class)) {
		get_det_stoch_comparison_data(nets, noiseType, dt, dataFolder, resultsFolder,
									   focal_class = state_class, sigma = sigma)
	} else {
		get_det_stoch_comparison_data(nets, noiseType, dt, dataFolder, resultsFolder,
									   focal_state = state, sigma = sigma)
	}
	if (!is.null(param_type)) d <- d %>% filter(ParamType %in% param_type)
	if (nrow(d) == 0) stop("No rows remain -- check state_class exists for these networks")

	d <- d %>%
		mutate(Network       = factor(Network, levels = nets),
			   StabilityClass = ifelse(ParamType == "monostable", "Monostable", "Multistable") %>%
			   					 factor(levels = c("Monostable", "Multistable")))

	if (!is.null(stability)) d <- d %>% filter(StabilityClass %in% stability)
	if (nrow(d) == 0) stop("No rows remain after filtering stability class")

	# sigma may now be a vector -- when it resolves to more than one distinct
	# NoiseLevel actually present in the data, every grouping/bracket step
	# below also groups by NoiseLevel (harmless extra grouping level when
	# there's only one, so the sigma = <scalar> callers -- Figure4.r's main
	# panel, Figure4_supp.r's per-sigma loop -- are unaffected) so multiple
	# noise levels are never pooled together into one paired test or one
	# Mean/SD bar, and get their own facet column instead.
	d <- d %>% mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(NoiseLevel))))
	multiSigma <- length(levels(d$NoiseLevel)) > 1

	# Paired significance tests, per (Network, StabilityClass, NoiseLevel), on
	# the raw (non-summarized) MRT triples -- NOT on the Mean/SD summary below.
	dSig <- d %>%
		group_by(Network, StabilityClass, NoiseLevel) %>%
		summarise(p_resample_vs_stoch = paired_wilcox_p(MRT_det, MRT),
				  p_resample_vs_orig  = paired_wilcox_p(MRT_det, MRT_orig),
				  .groups = "drop") %>%
		mutate(stars_stoch = stars_from_p(p_resample_vs_stoch),
			   stars_orig  = stars_from_p(p_resample_vs_orig))

	dLong <- d %>%
		select(ParamID, Network, StabilityClass, NoiseLevel, MRT_orig, MRT_det, MRT) %>%
		pivot_longer(c(MRT_orig, MRT_det, MRT), names_to = "Type", values_to = "Value") %>%
		mutate(Type = recode(Type,
							  MRT_orig = "Original (point λ)",
							  MRT_det  = "Deterministic (resampled λ)",
							  MRT      = "Stochastic") %>%
					   factor(levels = typeLevels))

	dSumm <- dLong %>%
		group_by(Network, StabilityClass, NoiseLevel, Type) %>%
		summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop")

	# Dodge offsets for 3 evenly-spaced groups at width `dodge_width`, matching
	# what position_dodge(width = dodge_width) places them at: level i (of n)
	# sits at (i - (n+1)/2) * dodge_width / n relative to the integer x position.
	off <- dodge_width * (seq_along(typeLevels) - (length(typeLevels) + 1) / 2) / length(typeLevels)
	names(off) <- typeLevels

	# Bracket geometry + label, one row per (Network, StabilityClass,
	# NoiseLevel): the Original-vs-Deterministic bracket sits left of center,
	# the Deterministic-vs-Stochastic bracket sits right of center. y_position
	# is placed just above the tallest errorbar top in that (Network,
	# StabilityClass, NoiseLevel) facet.
	dBrackets <- dSumm %>%
		group_by(Network, StabilityClass, NoiseLevel) %>%
		summarise(yTop = max(Mean + SD, na.rm = TRUE), .groups = "drop") %>%
		left_join(dSig, by = c("Network", "StabilityClass", "NoiseLevel")) %>%
		mutate(net_x = as.numeric(Network)) %>%
		group_by(StabilityClass, NoiseLevel) %>%
		mutate(yStep = diff(range(yTop, na.rm = TRUE)) %>% { ifelse(. == 0, yTop[1] * 0.1, . * 0.12) }) %>%
		ungroup() %>%
		mutate(y_pos = yTop + yStep)

	# geom_signif's stat collapses same-group rows to a single annotation, so
	# each bracket needs its own group -- otherwise every network's label
	# gets merged into one and mis-placed (all n rows share the default
	# group = -1 unless told otherwise).
	brack_orig  <- dBrackets %>%
		transmute(StabilityClass, NoiseLevel, xmin = net_x + off[["Original (point λ)"]],
				  xmax = net_x + off[["Deterministic (resampled λ)"]],
				  y_position = y_pos, annotation = stars_orig) %>%
		mutate(id = row_number())
	brack_stoch <- dBrackets %>%
		transmute(StabilityClass, NoiseLevel, xmin = net_x + off[["Deterministic (resampled λ)"]],
				  xmax = net_x + off[["Stochastic"]],
				  y_position = y_pos, annotation = stars_stoch) %>%
		mutate(id = row_number())

	p <- ggplot(dSumm, aes(x = Network, y = Mean, color = Type)) +
		geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
					  position = position_dodge(width = dodge_width), width = 0.3) +
		geom_point(position = position_dodge(width = dodge_width), size = 2.5) +
		geom_signif(data = brack_orig, aes(xmin = xmin, xmax = xmax, y_position = y_position,
										   annotations = annotation, group = id),
					manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, size = 0.4, textsize = 3.2) +
		geom_signif(data = brack_stoch, aes(xmin = xmin, xmax = xmax, y_position = y_position,
											annotations = annotation, group = id),
					manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, size = 0.4, textsize = 3.2) +
		theme_Publication() +
		labs(x = NULL,
			 y = if (multiSigma) paste0("MRT of ", state_class %||% "all-high", " state")
			     else paste0("MRT of ", state_class %||% "all-high", " state (σ = ", sigma, ")"),
			 color = "")

	# Single sigma: unchanged facet_wrap(StabilityClass) behavior. Multiple:
	# NoiseLevel becomes its own facet column (rows = StabilityClass) instead
	# of being folded into the y-axis label, which can't hold more than one
	# sigma value.
	p <- p + if (multiSigma) {
		facet_grid(rows = vars(StabilityClass), cols = vars(NoiseLevel), scales = "free_y")
	} else {
		facet_wrap(vars(StabilityClass), scales = "free_y")
	}

	if (!is.null(outFile)) {
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = width, height = height)
	}
	invisible(p)
}

# Violin of the paramwise MRT difference -- the resampling effect
# (MRT_det - MRT_orig, i.e. panel B collapsed to one number per ParamID) or
# the stochastic effect (MRT - MRT_det, panel C collapsed the same way) --
# across every noise level at once (x axis), faceted by StateClass (rows) x
# Network (cols). This is the state-class generalization of
# plot_paramwise_comparison()'s single-sigma, single-state (all-high)
# scatter: same two underlying comparisons, but summarized as a signed
# difference so single-high/double-high/all-high can share one figure
# instead of needing a scatter panel per class. NoiseLevel == 0 (the
# deterministic point-lambda baseline itself) is dropped since the
# difference there is 0 by construction. 2-node networks have no
# "double-high" states (see plot_det_stoch_network_comparison's note), so
# those facets come back empty -- expected, not an error.
plot_class_diff_violin <- function(nets, noiseType, dt, dataFolder, resultsFolder,
                                    compare = c("resample", "stochastic"),
                                    stateClasses = c("all-high", "single-high", "double-high"),
                                    param_type = NULL,
                                    outFile = NULL, panelWidth = 3, panelHeight = 3) {

	compare <- match.arg(compare)

	d <- map_dfr(stateClasses, function(cls) {
		get_det_stoch_comparison_data(nets, noiseType, dt, dataFolder, resultsFolder,
									   focal_class = cls) %>%
			mutate(StateClass = cls)
	})
	if (!is.null(param_type)) d <- d %>% filter(ParamType %in% param_type)

	if (compare == "resample") {
		d <- d %>% mutate(Diff = MRT_det - MRT_orig)
		ylab <- "Resampled - Original MRT"
		ttl  <- "Resampling effect"
	} else {
		d <- d %>% mutate(Diff = MRT - MRT_det)
		ylab <- "Stochastic - Resampled MRT"
		ttl  <- "Stochastic effect"
	}

	d <- d %>%
		filter(NoiseLevel != 0) %>%
		mutate(Network    = factor(Network, levels = nets),
			   StateClass = factor(StateClass, levels = stateClasses),
			   NoiseLevel = factor(NoiseLevel, levels = sort(unique(NoiseLevel))))
	if (nrow(d) == 0) stop("No rows remain -- check requested stateClasses exist for these networks")

	p <- ggplot(d, aes(x = NoiseLevel, y = Diff)) +
		geom_violin(fill = "steelblue", alpha = 0.5, scale = "width", color = NA) +
		geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.4) +
		stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6,
					 fill = "white", color = "black") +
		facet_grid(rows = vars(StateClass), cols = vars(Network), scales = "free_y") +
		theme_Publication() +
		theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
		labs(x = "Noise Level", y = ylab, title = ttl)

	if (!is.null(outFile)) {
		n_col <- length(unique(d$Network)); n_row <- length(unique(d$StateClass))
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = panelWidth * n_col, height = panelHeight * n_row)
	}
	invisible(p)
}

# ── Figure 3 (+ its supplement) panels -- these two need library(ggsignif)
#    loaded by the caller ──

# Mono/bi/tri-stable MRT boxplots by state class (all-high vs. single-high
# by default), one facet per ParamType actually present in the data, with a
# paired Wilcoxon test significance bracket between the two StateClasses at every
# (ParamType, NoiseLevel) x-position.
mrt_boxplot_stability <- function(net, noiseType, resultsFolder, dt = 0.01,
                                   stateClasses = c("all-high", "single-high"),
                                   outFile = NULL, panelWidth = 5, panelHeight = 4.5,
                                   dodge_width = 0.7) {

	f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_results.csv")
	d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
	d <- fill_mrt(d, unique(d$State)) %>%
		filter(StateClass %in% stateClasses) %>%
		group_by(ParamID, ParamType, NoiseLevel, StateClass) %>%
		summarise(MRT = sum(MRT), .groups = "drop") %>%
		mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
			   ParamType  = str_to_sentence(ParamType),
			   StateClass = factor(StateClass, levels = stateClasses))

	# Paired Wilcoxon test between the two StateClasses (paired by ParamID -- fill_mrt
	# already zero-filled both classes for every ParamID), per (ParamType,
	# NoiseLevel) facet/x-position.
	dWide <- d %>%
		select(ParamID, ParamType, NoiseLevel, StateClass, MRT) %>%
		pivot_wider(names_from = StateClass, values_from = MRT)

	dSig <- dWide %>%
		group_by(ParamType, NoiseLevel) %>%
		summarise(p = paired_wilcox_p(.data[[stateClasses[1]]], .data[[stateClasses[2]]]),
				  yTop = max(c(.data[[stateClasses[1]]], .data[[stateClasses[2]]]), na.rm = TRUE),
				  .groups = "drop") %>%
		mutate(stars = stars_from_p(p),
			   noise_x = as.numeric(NoiseLevel))

	off <- dodge_width * (seq_along(stateClasses) - (length(stateClasses) + 1) / 2) / length(stateClasses)

	dBrackets <- dSig %>%
		group_by(ParamType) %>%
		mutate(yStep = diff(range(yTop, na.rm = TRUE)) %>%
				 { ifelse(. == 0, max(yTop, 0.05) * 0.1, . * 0.12) }) %>%
		ungroup() %>%
		mutate(y_position = yTop + yStep,
			   xmin = noise_x + off[1], xmax = noise_x + off[2]) %>%
		mutate(id = row_number())

	p <- ggplot(d, aes(x = NoiseLevel, y = MRT, color = StateClass, fill = StateClass)) +
		geom_boxplot(outlier.size = 0.5, alpha = 0.4, position = position_dodge(width = dodge_width)) +
		stat_summary(aes(fill = StateClass), fun = mean, geom = "point",
					 shape = 23, size = 2, color = "black",
					 position = position_dodge(width = dodge_width)) +
		geom_signif(data = dBrackets, aes(xmin = xmin, xmax = xmax, y_position = y_position,
										  annotations = stars, group = id),
					manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, size = 0.3, textsize = 2.6) +
		facet_wrap(vars(ParamType)) +
		theme_Publication() +
		theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
		labs(x = "Noise Level", y = "MRT")

	if (!is.null(outFile)) {
		n_panels <- length(unique(d$ParamType))
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = panelWidth * n_panels, height = panelHeight)
	}
	invisible(p)
}

# Node letters (A, B, C, ...) in the exact positional order used inside a
# State tuple like "(1, 0, 1, 0)" -- taken from the RACIPE solution file's
# own node columns rather than assumed from the .topo file's node order,
# since that's the file that actually defines the tuple's column order.
get_node_letters <- function(net, dataFolder) {
    sol <- read_solutions(file.path(dataFolder, paste0(net, "_solution.dat")))
    setdiff(colnames(sol), c("ParamID", "nStates", "basin"))
}

# "(1, 0, 1, 0)" -> "AC" (the letters of every node that's "on"); an
# all-off tuple (no "1"s) becomes the literal string "none" rather than ""
# so it still reads as a real axis label instead of a blank bar.
state_to_letters <- function(state_str, node_letters) {
    vec  <- as.integer(str_split(str_remove_all(state_str, "[()]"), ",\\s*")[[1]])
    code <- paste0(node_letters[vec == 1], collapse = "")
    if (code == "") "none" else code
}

# Same, for a compound "Original" label (multiple states joined by "_" for
# a multistable ParamID) -- "(1, 0, 1, 0)_(0, 1, 0, 1)" -> "AC_BD".
original_to_letters <- function(original_str, node_letters) {
    str_split(original_str, "_")[[1]] %>%
        map_chr(state_to_letters, node_letters = node_letters) %>%
        paste0(collapse = "_")
}

# Deterministic (NoiseLevel == 0) attractor-COMBINATION frequency, using
# fill_mrt()'s compound "Original" label (one per ParamID: the sorted,
# underscore-joined set of every state that ParamID visits at zero noise)
# instead of per-state frequency -- so a bistable ParamID's two states are
# counted together as the single combination it actually is, not as two
# independent state tallies. MRT is summed within each ParamID's own
# Original group first (== 1 by construction, since MRT partitions that
# ParamID's ICs/time across exactly the states making up its Original),
# then summed across ParamIDs sharing the same Original and normalized per
# network -- i.e. the fraction of parameter sets landing on each
# attractor combination, computed via the MRT sum rather than a raw
# ParamID count. Labels are relabeled to on-node letters (e.g. "AC",
# "AC_BD") via state_to_letters()/original_to_letters() so multistable
# combinations stay readable instead of cluttered with full tuples.
plot_original_frequency <- function(nets, noiseType, resultsFolder, dataFolder, dt = 0.01,
                                     param_type = NULL, topN = 10,
                                     outFile = NULL, panelWidth = 4, panelHeight = 4, ncol = 3) {

    dAll <- map_dfr(nets, function(net) {
        effectiveNoiseType <- noise_type_for(net, noiseType)
        f <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
        if (!file.exists(f)) { warning("Missing: ", net); return(NULL) }
        d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
        if (!is.null(param_type)) d <- d %>% filter(ParamType %in% param_type)
        d <- fill_mrt(d, unique(d$State)) %>% filter(NoiseLevel == 0, !is.na(Original))

        node_letters <- get_node_letters(net, dataFolder)
        d %>%
            group_by(ParamID, Original) %>%
            summarise(MRT_total = sum(MRT), .groups = "drop") %>%
            mutate(Label   = map_chr(Original, original_to_letters, node_letters = node_letters),
                   Network = net)
    })
    if (nrow(dAll) == 0) stop("No data found")

    dFreqFull <- dAll %>%
        group_by(Network, Label) %>%
        summarise(TotalMRT = sum(MRT_total), .groups = "drop") %>%
        group_by(Network) %>%
        mutate(Frac = TotalMRT / sum(TotalMRT)) %>%
        ungroup()

    render_label_frequency_plot(dFreqFull, nets, dataFolder, topN, outFile, panelWidth, panelHeight, ncol)
}

# Same "attractor-combination frequency" plot as plot_original_frequency(),
# but sourced directly from RACIPE's OWN steady-state search (solution.dat
# via read_solutions()) instead of the noise-simulation pipeline's
# NoiseLevel == 0 rows -- these are two different deterministic-attractor
# derivations that happen to usually agree, not guaranteed-identical, so
# "get it from RACIPE directly" is a real, distinct data source. Each row
# of solution.dat is one steady state RACIPE actually found for that
# ParamID (one basin); discretized against that node's own mean (the same
# convention used everywhere else in this project, e.g. Figure1.r's own
# Panel B(ii) thresholds) and aggregated per ParamID into the same sorted,
# underscore-joined compound label fill_mrt()'s "Original" builds (e.g. a
# bistable ParamID's two solution.dat rows combine into "(0, 1)_(1, 0)").
get_racipe_attractor_combo <- function(net, dataFolder) {
    solutions <- read_solutions(file.path(dataFolder, paste0(net, "_solution.dat")))
    nodes <- setdiff(colnames(solutions), c("ParamID", "nStates", "basin"))

    solutions %>%
        mutate(across(all_of(nodes), ~ if_else(.x > mean(.x), 1L, 0L))) %>%
        unite(col = "State", all_of(nodes), sep = ", ") %>%
        mutate(State = paste0("(", State, ")")) %>%
        group_by(ParamID) %>%
        summarise(Original = paste0(sort(unique(State)), collapse = "_"), .groups = "drop")
}

plot_racipe_original_frequency <- function(nets, dataFolder, topN = 10,
                                            outFile = NULL, panelWidth = 4, panelHeight = 4, ncol = 3) {

    dFreqFull <- map_dfr(nets, function(net) {
        node_letters <- get_node_letters(net, dataFolder)
        get_racipe_attractor_combo(net, dataFolder) %>%
            mutate(Label = map_chr(Original, original_to_letters, node_letters = node_letters)) %>%
            count(Label, name = "n") %>%
            mutate(Frac = n / sum(n), Network = net)
    })
    if (nrow(dFreqFull) == 0) stop("No data found")

    render_label_frequency_plot(dFreqFull, nets, dataFolder, topN, outFile, panelWidth, panelHeight, ncol)
}

# Shared renderer behind plot_original_frequency() and
# plot_racipe_original_frequency() -- both just differ in how dFreqFull
# (Network, Label, Frac; not yet topN-filtered) gets built.
render_label_frequency_plot <- function(dFreqFull, nets, dataFolder, topN = 10,
                                         outFile = NULL, panelWidth = 4, panelHeight = 4, ncol = 3) {

    dFreq <- dFreqFull %>%
        group_by(Network) %>%
        slice_max(Frac, n = topN, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(Network = factor(Network, levels = nets))

    # In-panel key mapping each single-node letter to its binary tuple (e.g.
    # "A: (1, 0)") for that network's own node count -- letter codes alone
    # (state_to_letters()/original_to_letters()) don't self-explain which
    # node is which position without this. Positioned in DATA coordinates,
    # not the usual x = -Inf/Inf corner trick -- with coord_flip() on a
    # discrete axis, Inf leaves almost no padding beyond the last category
    # and gets clipped to an invisible sliver. x = 1.5 sits just above the
    # bottom (smallest-Frac, i.e. shortest-bar) row instead, with y at ~55%
    # of that facet's own max Frac -- comfortably inside the empty space to
    # the right of those short bars, left-justified so it doesn't matter
    # exactly how long each network's key text is.
    keyText <- map_dfr(levels(dFreq$Network), function(net) {
        node_letters <- get_node_letters(net, dataFolder)
        # Every unique atomic letter-code actually shown for this network
        # (splitting each displayed compound label, e.g. "A_AB", back into
        # its individual pieces "A" and "AB") -- not every possible 2^n
        # combination, just the ones that appear -- one line each, no
        # additional generic text.
        atomicCodes <- dFreq$Label[dFreq$Network == net] %>%
            str_split("_") %>% unlist() %>% unique() %>% sort()
        lines <- map_chr(atomicCodes, function(code) {
            v <- as.integer(node_letters %in% strsplit(code, "")[[1]])
            paste0(code, ": (", paste(v, collapse = ", "), ")")
        })
        tibble(Network = net, label = paste(lines, collapse = "\n"))
    }) %>% mutate(Network = factor(Network, levels = levels(dFreq$Network)))

    # Anchored off the GLOBAL max Frac (shared across every facet, not each
    # facet's own) -- facet_wrap's scales = "free_y" frees the post-coord_flip
    # axis, which is the already-per-panel-free discrete state-label axis, not
    # the continuous Frac axis, so every panel actually shares the same Frac
    # range regardless of that network's own bars being much shorter (e.g.
    # TT/TTSA vs. TS/TSSA). A per-network anchor would then sit far short of
    # where the visible bars actually end in the shared range, and (right-
    # justified) grow back left across them.
    keyText <- keyText %>% mutate(y = 1.15 * max(dFreq$Frac))

    # reorder(Label, Frac) computes ONE global order from the mean Frac of
    # every row sharing that exact Label string across ALL networks at
    # once -- not each facet's own descending order -- so a given
    # network's bars can come out visibly out of order (e.g. a bigger bar
    # sitting below a smaller one) whenever other networks' rows for
    # differently-valued occurrences of a shared label string shift the
    # global ranking. Compositing Network into the sort key (and dropping
    # it again for display via scale_x_discrete's labeller) keeps each
    # facet's ordering fully independent, matching what "sorted per panel"
    # actually requires.
    dFreq <- dFreq %>%
        mutate(LabelKey = paste(Network, Label, sep = "\r")) %>%
        group_by(Network) %>%
        mutate(LabelKey = factor(LabelKey, levels = LabelKey[order(Frac)])) %>%
        ungroup()

    p <- ggplot(dFreq, aes(x = LabelKey, y = Frac)) +
        geom_col(fill = "steelblue") +
        geom_text(data = keyText, aes(x = 1.2, y = y, label = label), inherit.aes = FALSE,
                  hjust = 1, vjust = 0, size = 3.6, fontface = "italic", lineheight = 1) +
        coord_flip() +
        facet_wrap(vars(Network), scales = "free_y", ncol = ncol) +
        scale_x_discrete(labels = function(x) sub("^.*\r", "", x)) +
        theme_Publication() +
        labs(x = NULL, y = "Fraction of parameter sets (deterministic)")

    if (!is.null(outFile)) {
        n_panels <- length(unique(dFreq$Network))
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = panelWidth * min(n_panels, ncol),
               height = panelHeight * ceiling(n_panels / ncol))
    }
    invisible(p)
}

# MRT vs. noise level, colored by StateClass, faceted by Network (cols) x
# StabilityClass (rows) -- the two-way generalization of Figure3's
# plot_stateclass_transition() (which only facets by Network, one
# param_type at a time) for panels that need mono- vs multistable shown
# side by side rather than picked with a single param_type filter.
plot_stateclass_by_stability <- function(nets, noiseType, resultsFolder, dt = 0.01,
                                          stateClasses = c("single-high", "double-high", "all-high"),
                                          outFile = NULL, panelWidth = 3.2, panelHeight = 3.2,
                                          dodge_width = 0.6) {

    dAll <- map_dfr(nets, function(net) {
        f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_results.csv")
        if (!file.exists(f)) { warning("Missing: ", net); return(NULL) }
        d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
        d <- fill_mrt(d, unique(d$State))
        d %>% mutate(Network = net)
    })
    if (nrow(dAll) == 0) stop("No data found")

    dAll <- dAll %>%
        filter(StateClass %in% stateClasses) %>%
        group_by(ParamID, ParamType, NoiseLevel, Network, StateClass) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        mutate(NoiseLevel     = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
               Network        = factor(Network, levels = nets),
               StateClass     = factor(StateClass, levels = stateClasses),
               StabilityClass = ifelse(ParamType == "monostable", "Monostable", "Multistable") %>%
                                 factor(levels = c("Monostable", "Multistable")))

    p <- ggplot(dAll, aes(x = NoiseLevel, y = MRT, color = StateClass)) +
        stat_summary(fun.data = mean_sd, geom = "errorbar",
                     position = position_dodge(width = dodge_width), width = 0.3) +
        stat_summary(fun = mean, geom = "point",
                     position = position_dodge(width = dodge_width), size = 2) +
        facet_grid(rows = vars(StabilityClass), cols = vars(Network), scales = "free_y") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = "Average MRT")

    if (!is.null(outFile)) {
        n_col <- length(unique(dAll$Network))
        n_row <- length(unique(droplevels(dAll$StabilityClass)))
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = panelWidth * n_col, height = panelHeight * n_row)
    }
    invisible(p)
}

# MRT of one focal StateClass (default all-high) vs. noise level, colored
# by DT (the integration/resampling time step) instead of by state --
# facet_grid(rows = StabilityClass, cols = Network), one call per noise
# type (Additive/Multiplicative are separate calls, not folded into one
# facet dimension, since they're different noise-update mechanisms
# entirely, not just another factor level of the same thing).
# fill_mrt() is called separately per DT value, not once on the combined
# multi-DT data -- it drops the DT column internally (select() only keeps
# ParamID/ParamType/NoiseLevel/State/MRT before its zero-filling join), so
# feeding it multiple DTs' rows at once for the same (ParamID, NoiseLevel,
# State) key would silently scramble which MRT belongs to which DT.
# Needs library(ggsignif) loaded by the caller (same convention as the
# other bracket-drawing functions in this file).
plot_mrt_by_dt <- function(nets, noiseType, resultsFolder,
                            noise_levels = c(0.001, 0.01, 0.1, 1.0),
                            dt_values = c(0.01, 1, 10),
                            state_classes = c("all-high", "single-high", "double-high"),
                            outFile = NULL, panelWidth = 3.2, panelHeight = 3.2,
                            dodge_width = 0.7) {

    dAll <- map_dfr(nets, function(net) {
        f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_results.csv")
        if (!file.exists(f)) { warning("Missing: ", net); return(NULL) }
        d_raw <- read_csv(f, show_col_types = FALSE)
        map_dfr(dt_values, function(this_dt) {
            d <- d_raw %>% filter(DT == this_dt)
            if (nrow(d) == 0) return(NULL)
            fill_mrt(d, unique(d$State)) %>% mutate(DT = this_dt)
        }) %>% mutate(Network = net)
    })
    if (nrow(dAll) == 0) stop("No data found")

    # Deterministic (NoiseLevel = 0) baseline -- DT is irrelevant at zero
    # noise (nothing to resample, so all DT values collapse to the same
    # trajectory), so instead of dodging it in as another x-position it's
    # averaged across whatever DT rows exist and drawn as one dashed
    # reference line per facet below.
    #
    # State classes with more than one member state (e.g. "single-high" =
    # (1,0) or (0,1) for TS) need their per-state MRT summed within each
    # (ParamID, DT) BEFORE averaging across ParamIDs -- averaging the raw
    # per-state rows directly (one row per member state) mixes in the
    # unvisited member's MRT = 0 row for every reaching ParamID, roughly
    # halving the true value for 2-state classes (and worse for classes
    # with more members). "all-high" has only one member state so this
    # bug was invisible there.
    dBaseline <- dAll %>%
        filter(StateClass %in% state_classes, NoiseLevel == 0) %>%
        mutate(Network        = factor(Network, levels = nets),
               StateClass     = factor(StateClass, levels = state_classes),
               StabilityClass = ifelse(ParamType == "monostable", "Monostable", "Multistable") %>%
                                 factor(levels = c("Monostable", "Multistable"))) %>%
        group_by(ParamID, Network, StateClass, StabilityClass, DT) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        group_by(Network, StateClass, StabilityClass) %>%
        summarise(MRT = mean(MRT), .groups = "drop")

    dAll <- dAll %>%
        filter(StateClass %in% state_classes, NoiseLevel %in% noise_levels) %>%
        group_by(ParamID, ParamType, NoiseLevel, Network, DT, StateClass) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        mutate(NoiseLevel     = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
               Network        = factor(Network, levels = nets),
               DT             = factor(DT, levels = sort(dt_values)),
               StateClass     = factor(StateClass, levels = state_classes),
               StabilityClass = ifelse(ParamType == "monostable", "Monostable", "Multistable") %>%
                                 factor(levels = c("Monostable", "Multistable")))
    if (nrow(dAll) == 0) stop("No rows remain -- check state_classes exist for these networks")

    # Paired Wilcoxon test between ADJACENT DT values (paired by ParamID --
    # fill_mrt() already zero-fills every ParamID at every DT), per (Network,
    # StateClass, StabilityClass, NoiseLevel) x-position -- same
    # adjacent-pair-of-dodged-groups bracket pattern as
    # plot_det_stoch_network_comparison()'s Original/Deterministic/Stochastic
    # brackets, just for however many DT values are dodged here (2 brackets
    # for the default 3 DT values: 0.01-vs-1, 1-vs-10).
    dtLevels <- levels(dAll$DT)
    dWide <- dAll %>%
        select(ParamID, Network, StateClass, StabilityClass, NoiseLevel, DT, MRT) %>%
        pivot_wider(names_from = DT, values_from = MRT)

    dBrackets <- NULL
    if (length(dtLevels) >= 2) {
        off <- dodge_width * (seq_along(dtLevels) - (length(dtLevels) + 1) / 2) / length(dtLevels)
        names(off) <- dtLevels

        dSig <- map2_dfr(dtLevels[-length(dtLevels)], dtLevels[-1], function(dt1, dt2) {
            dWide %>%
                group_by(Network, StateClass, StabilityClass, NoiseLevel) %>%
                summarise(p = paired_wilcox_p(.data[[dt1]], .data[[dt2]]), .groups = "drop") %>%
                mutate(stars = stars_from_p(p), dt1 = dt1, dt2 = dt2)
        }) %>%
            mutate(pairIdx = match(dt1, dtLevels))

        dTop <- dAll %>%
            group_by(Network, StateClass, StabilityClass, NoiseLevel, DT) %>%
            summarise(Mean = mean(MRT), SD = sd(MRT), .groups = "drop") %>%
            group_by(Network, StateClass, StabilityClass, NoiseLevel) %>%
            summarise(yTop = max(Mean + SD, na.rm = TRUE), .groups = "drop")

        dBrackets <- dSig %>%
            left_join(dTop, by = c("Network", "StateClass", "StabilityClass", "NoiseLevel")) %>%
            mutate(noise_x = as.numeric(NoiseLevel)) %>%
            group_by(Network, StateClass, StabilityClass) %>%
            mutate(yStep = diff(range(yTop, na.rm = TRUE)) %>%
                     { ifelse(. == 0, max(yTop, 0.05, na.rm = TRUE) * 0.1, . * 0.12) }) %>%
            ungroup() %>%
            mutate(y_position = yTop + yStep * pairIdx,
                   xmin = noise_x + off[dt1], xmax = noise_x + off[dt2]) %>%
            mutate(id = row_number())
    }

    p <- ggplot(dAll, aes(x = NoiseLevel, y = MRT, color = DT)) +
        stat_summary(fun.data = mean_sd, geom = "errorbar",
                     position = position_dodge(width = dodge_width), width = 0.3) +
        stat_summary(fun = mean, geom = "point",
                     position = position_dodge(width = dodge_width), size = 2)

    if (!is.null(dBrackets)) {
        p <- p + geom_signif(data = dBrackets, aes(xmin = xmin, xmax = xmax, y_position = y_position,
                                                    annotations = stars, group = id),
                              manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, size = 0.3, textsize = 2.2)
    }

    if (nrow(dBaseline) > 0) {
        p <- p + geom_hline(data = dBaseline, aes(yintercept = MRT),
                             inherit.aes = FALSE, linetype = "dashed",
                             color = "grey30", linewidth = 0.5)
    }

    p <- p +
        facet_grid(rows = vars(StateClass, StabilityClass), cols = vars(Network), scales = "free_y") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = "MRT", title = noiseType)

    if (!is.null(outFile)) {
        n_col <- length(unique(dAll$Network))
        n_row <- length(unique(droplevels(dAll$StateClass))) * length(unique(droplevels(dAll$StabilityClass)))
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = panelWidth * n_col, height = panelHeight * n_row)
    }
    invisible(p)
}

# Iteratively drops the smallest-MRT entry from a named MRT vector until
# the remaining min/max ratio clears ratio_cutoff (or only one entry is
# left) -- the "meaningfully accessible" survivors at that threshold. A
# ratio_cutoff of 0.1 means a state only counts as accessible once its MRT
# is at least 10% of the top state's -- a fairly generous bar, not a
# ~50/50 "balanced" one. All-zero input returns an empty vector (nFates = 0).
prune_accessible <- function(mrt_named, ratio_cutoff = 0.1) {
    repeat {
        if (length(mrt_named) <= 1) break
        rmin <- min(mrt_named); rmax <- max(mrt_named)
        if (rmax == 0) { mrt_named <- mrt_named[0]; break }
        if (rmin / rmax >= ratio_cutoff) break
        mrt_named <- mrt_named[-which.min(mrt_named)[1]]
    }
    mrt_named
}

# Classifies a double-high-only "Original" label (fill_mrt()'s compound
# per-ParamID attractor combination) into exactly 3 stability classes, or
# NA if it's something else (tristable, or mixed with a non-double-high
# state) -- "mirror" bistable pairs are two double-high states sharing NO
# active node (true complements, e.g. AB/CD on a 4-node network);
# "non-mirror" pairs share exactly one (e.g. AB/AD).
classify_double_high_origin <- function(original_str) {
    parts  <- str_split(original_str, "_")[[1]]
    nHighs <- str_count(parts, "1")
    if (length(parts) == 1 && nHighs[1] == 2) return("Monostable (double-high)")
    if (length(parts) == 2 && all(nHighs == 2)) {
        actives <- map(parts, ~ which(as.integer(str_split(str_remove_all(.x, "[()]"), ",\\s*")[[1]]) == 1))
        shared  <- length(intersect(actives[[1]], actives[[2]]))
        if (shared == 0) return("Bistable (mirror)")
        if (shared == 1) return("Bistable (non-mirror)")
    }
    NA_character_
}

# For every ParamID whose deterministic Original is purely double-high (one
# of the 3 classify_double_high_origin() classes), breaks its occupancy at
# each of several noise levels into: each single-high state (letter-coded
# via state_to_letters()), "Original" (occupancy remaining in whichever
# double-high state(s) made up its own Original), and "Other" (everything
# else -- all-high, a different double-high state never part of Original,
# etc.). Long format, one row per (ParamID, NoiseLevel, Role) -- the shared
# data source behind the nFates / heatmap / ratio panels below.
compute_access_breakdown <- function(net, noiseType, resultsFolder, dataFolder, dt = 0.01,
                                      sigmas = c(0.001, 0.01, 0.05, 0.1)) {

    effectiveNoiseType <- noise_type_for(net, noiseType)
    f <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
    d <- read_csv(f, show_col_types = FALSE) %>%
        filter(DT == dt) %>%
        mutate(NoiseLevel = as.numeric(as.character(NoiseLevel)))
    d <- fill_mrt(d, unique(d$State)) %>% filter(!is.na(Original))

    node_letters <- get_node_letters(net, dataFolder)

    origMap <- d %>% distinct(ParamID, Original) %>%
        mutate(StabilityClass = map_chr(Original, classify_double_high_origin),
               OriginalLabel  = map_chr(Original, original_to_letters, node_letters = node_letters)) %>%
        filter(!is.na(StabilityClass))
    if (nrow(origMap) == 0) stop("No monostable/bistable double-high ParamIDs found for ", net)

    singleStates <- d %>% filter(StateClass == "single-high") %>% distinct(State) %>% pull(State)
    singleLabels <- setNames(map_chr(singleStates, state_to_letters, node_letters = node_letters), singleStates)

    d %>%
        filter(ParamID %in% origMap$ParamID, NoiseLevel %in% sigmas) %>%
        left_join(origMap %>% select(ParamID, StabilityClass, OriginalLabel), by = "ParamID") %>%
        mutate(Role = case_when(
            State %in% names(singleLabels)     ~ unname(singleLabels[State]),
            str_detect(Original, fixed(State)) ~ "Original",
            TRUE                                ~ "Other"
        )) %>%
        group_by(ParamID, NoiseLevel, StabilityClass, OriginalLabel, Role) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        mutate(Network = net)
}

# Panel E-style: distribution of the number of single-high states
# meaningfully accessed (prune_accessible() survivors) under noise, split
# by the 3 stability classes from compute_access_breakdown(), at several
# fixed noise levels shown side by side.
plot_access_nFates <- function(accessData, sigmas, ratio_cutoff = 0.1,
                                stabilityLevels = c("Monostable (double-high)", "Bistable (mirror)", "Bistable (non-mirror)"),
                                outFile = NULL, width = 10, height = 4) {

    node_letters <- setdiff(unique(accessData$Role), c("Original", "Other"))

    dN <- accessData %>%
        filter(Role %in% node_letters, NoiseLevel %in% sigmas) %>%
        group_by(ParamID, NoiseLevel, StabilityClass) %>%
        summarise(nFates = {
            v <- setNames(MRT, Role)
            length(prune_accessible(v, ratio_cutoff))
        }, .groups = "drop") %>%
        mutate(NoiseLevel     = factor(NoiseLevel, levels = sigmas),
               StabilityClass = factor(StabilityClass, levels = stabilityLevels),
               nFates         = factor(nFates, levels = as.character(0:length(node_letters))))

    dSumm <- dN %>%
        count(StabilityClass, NoiseLevel, nFates, .drop = FALSE) %>%
        group_by(StabilityClass, NoiseLevel) %>%
        mutate(Frac = n / sum(n)) %>%
        ungroup()

    p <- ggplot(dSumm, aes(x = nFates, y = Frac, fill = NoiseLevel)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        facet_wrap(vars(StabilityClass)) +
        theme_Publication() +
        labs(x = "Number of single-high states accessed", y = "Fraction of parameter sets", fill = "Noise Level")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Panel F-style: heatmap of mean occupancy (across ParamIDs sharing the
# same OriginalLabel) in each possible "role" -- every single-high state,
# plus "Original" (still sitting in the double-high state itself) and
# "Other" (anything else) -- at each noise level. Rows ordered by overall
# ParamID frequency (most common attractor combination at the top).
plot_access_heatmap <- function(accessData, sigmas, outFile = NULL, width = 18, height = 7) {

    node_letters <- setdiff(unique(accessData$Role), c("Original", "Other"))
    roleLevels   <- c(node_letters, "Original", "Other")

    labelOrder <- accessData %>% distinct(ParamID, OriginalLabel) %>%
        count(OriginalLabel) %>% arrange(desc(n)) %>% pull(OriginalLabel)

    dF <- accessData %>%
        filter(NoiseLevel %in% sigmas) %>%
        mutate(Role          = factor(Role, levels = roleLevels),
               NoiseLevel    = factor(NoiseLevel, levels = sigmas),
               OriginalLabel = factor(OriginalLabel, levels = rev(labelOrder))) %>%
        group_by(NoiseLevel, OriginalLabel, Role) %>%
        summarise(MeanMRT = mean(MRT), .groups = "drop")

    # nrow = 1 (rather than facet_wrap's default near-square grid) so every
    # one of the ~20 OriginalLabel rows only needs to fit once vertically,
    # not be legible twice over in a 2-row grid.
    p <- ggplot(dF, aes(x = Role, y = OriginalLabel, fill = MeanMRT)) +
        geom_tile(color = "white") +
        facet_wrap(vars(NoiseLevel), nrow = 1) +
        scale_fill_viridis_c(name = "Mean MRT") +
        theme_Publication() +
        theme(axis.text.y = element_text(size = 9)) +
        labs(x = "Accessed state", y = "Original (deterministic) attractor")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Panel G-style: for the same OriginalLabel rows as plot_access_heatmap(),
# the max/min MRT ratio among just the ACCESSED single-high states
# (prune_accessible() survivors) per ParamID -- how skewed access is even
# among the fates that clear the accessibility bar, not just whether they
# clear it (that's Panel E). Ratio == 1 when only one fate survives.
plot_access_ratio <- function(accessData, sigmas, ratio_cutoff = 0.1,
                               outFile = NULL, width = 18, height = 7) {

    node_letters <- setdiff(unique(accessData$Role), c("Original", "Other"))

    labelOrder <- accessData %>% distinct(ParamID, OriginalLabel) %>%
        count(OriginalLabel) %>% arrange(desc(n)) %>% pull(OriginalLabel)

    dG <- accessData %>%
        filter(Role %in% node_letters, NoiseLevel %in% sigmas) %>%
        group_by(ParamID, NoiseLevel, OriginalLabel) %>%
        summarise(
            nFates = { v <- setNames(MRT, Role); length(prune_accessible(v, ratio_cutoff)) },
            Ratio  = { v <- setNames(MRT, Role); s <- prune_accessible(v, ratio_cutoff)
                       if (length(s) == 0) NA_real_ else max(s) / min(s) },
            .groups = "drop"
        ) %>%
        filter(!is.na(Ratio)) %>%
        mutate(NoiseLevel    = factor(NoiseLevel, levels = sigmas),
               OriginalLabel = factor(OriginalLabel, levels = labelOrder),
               nFates        = factor(nFates))

    p <- ggplot(dG, aes(x = OriginalLabel, y = Ratio, fill = nFates)) +
        geom_violin(scale = "width", alpha = 0.7, color = NA,
                    position = position_dodge(width = 0.8)) +
        facet_wrap(vars(NoiseLevel)) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
        labs(x = "Original (deterministic) attractor", y = "Max / Min MRT ratio (accessed single-high states)",
             fill = "# single-high\nstates accessed")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Re-labels compute_access_breakdown()'s Role for monostable double-high
# ParamIDs by Hamming distance from the double-high Original to each
# single-high fate, instead of by raw node letter -- collapses the
# network's specific node identities down to a distance-based scheme so
# different OriginalLabels (AB, AC, ...) become directly comparable:
#   - "Original":  Role == "Original" (still occupying the double-high state)
#   - "SH1_1"/"SH1_2": the two single-high fates at Hamming distance 1 (the
#     two nodes already active in Original) -- split by the alphabetical
#     order of OriginalLabel's own two letters (state_to_letters() always
#     emits them in node order/alphabetical, e.g. "BD" not "DB"), so the
#     split is consistent across every ParamID sharing that OriginalLabel
#   - "SH3": the single-high fate(s) at Hamming distance 3 (the nodes that
#     were off in Original) -- pooled into one bucket rather than split;
#     for a 4-node network there are exactly 2 and they aren't kept distinct
#   - "Other": Role == "Other", unchanged
classify_double_high_role <- function(accessData) {
    accessData %>%
        filter(StabilityClass == "Monostable (double-high)") %>%
        mutate(
            Role2 = case_when(
                Role == "Original"              ~ "Self",
                Role == "Other"                  ~ "Other",
                Role == str_sub(OriginalLabel, 1, 1) ~ "SH1_1",
                Role == str_sub(OriginalLabel, 2, 2) ~ "SH1_2",
                TRUE                              ~ "SH3"
            )
        )
}

# For ParamIDs of a given stability class (classify_double_high_role() or
# classify_bistable_role() output), classifies each (ParamID, NoiseLevel)
# by how many of the given fateRoles buckets survive prune_accessible() --
# i.e. how many distinct single-high fates that ParamID meaningfully
# spreads its post-noise occupancy across. Feeds both
# plot_double_high_alluvial()'s ribbon color and its own summary panel,
# plot_nfates_class_summary(). nFates == 0 (all-zero MRT among the fate
# buckets) is floored to 1 ("Monostable") since there's no spread to speak
# of either way. fateRoles defaults to the monostable-double-high case's 3
# Hamming-distance buckets (SH1_1/SH1_2/SH3); pass c("A","B","C","D") for
# the bistable classes (classify_bistable_role() doesn't collapse any of
# the 4 raw single-high letters together, unlike the monostable case's SH3
# pooling), which can therefore reach nFates == 4 ("Tetrastable").
compute_nfates_class <- function(roleData, fateRoles = c("SH1_1", "SH1_2", "SH3"), ratio_cutoff = 0.1) {
    classLabels <- setNames(NFATES_TIER_NAMES[seq_along(fateRoles)], as.character(seq_along(fateRoles)))
    roleData %>%
        filter(Role2 %in% fateRoles) %>%
        # Collapses any Role2 buckets that pool more than one physical
        # single-high state (e.g. monostable's SH3) *before* pruning, so
        # the named vector below always has exactly length(fateRoles)
        # entries, never more from a repeated name.
        group_by(ParamID, NoiseLevel, Role2) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        group_by(ParamID, NoiseLevel) %>%
        summarise(nFates = {
            v <- setNames(MRT, Role2)
            length(prune_accessible(v, ratio_cutoff))
        }, .groups = "drop") %>%
        mutate(nFates = pmax(nFates, 1),
               ClassLabel = factor(classLabels[as.character(nFates)],
                                    levels = NFATES_TIER_NAMES[seq_along(fateRoles)]))
}

# Role-prep for the two bistable stability classes (mirror / non-mirror),
# analogous to classify_double_high_role() but without the Hamming-distance
# SH1/SH3 split: a compound bistable Original is 2 double-high states, and
# each can independently put a *different* single-high letter at distance
# 1 (e.g. Original "AB_CD": A and B are distance 1 from AB but distance 3
# from CD, while C and D are the reverse) -- there's no single fixed
# "close pair" to split on the way there is for one monostable double-high
# source. Keeps all 4 raw single-high letters as their own Role2 (instead
# of the monostable case's 2-populated-of-4 split), and relabels
# Role == "Original" to "Self" so it shares a name with
# classify_double_high_role()'s output.
classify_bistable_role <- function(accessData, stabilityClass) {
    accessData %>%
        filter(StabilityClass == stabilityClass) %>%
        mutate(Role2 = if_else(Role == "Original", "Self", Role))
}

# Shared multistability color scale, keyed by compute_nfates_class()'s
# ClassLabel -- 3 tiers for the monostable case (fateRoles has 3 entries),
# 4 for the bistable cases (fateRoles has 4, so Tetrastable can occur).
NFATES_TIER_NAMES  <- c("Monostable", "Bistable", "Tristable", "Tetrastable")
NFATES_CLASS_COLORS <- setNames(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"), NFATES_TIER_NAMES)

# Plot D-style: ONE combined alluvial (not faceted) with NoiseLevel itself
# as the first axis and the 5 destination buckets (classify_double_high_role()'s
# Role2) as the second -- a real two-box-per-stage flow (like a classic
# Origin -> Destination alluvial) instead of 4 small facets each fanning out
# from a single trivial "Double-high" source box, which had nothing to flow
# from/to and so read as stacked bars rather than a flow. Ribbon width is
# MRT (summed across ParamIDs within each (NoiseLevel, Role2, ClassLabel)
# cell); color is compute_nfates_class()'s per-ParamID multistability class
# -- a single ParamID's occupancy typically spans more than one destination
# bucket at once, which is the whole point of the classification.
plot_double_high_alluvial <- function(roleData, classData, sigmas,
                                       roleLevels = c("SH1_1", "SH1_2", "SH3", "Self", "Other"),
                                       classColors = NFATES_CLASS_COLORS,
                                       outFile = NULL, width = 10, height = 6) {

    # Each NoiseLevel stratum is normalized to its OWN total (Frac sums to
    # 1 within every noise level) rather than plotted as raw summed MRT --
    # raw MRT sums are only comparable *within* one noise level's own
    # ParamID population; stacking 4 noise levels' raw totals on one shared
    # axis (as this used to) produces a combined running total (e.g. up to
    # ~4x a single level's total) that mixes 4 separate experiments'
    # magnitudes into one number with no real meaning. The absolute Frac
    # values along the y axis are still just a stacking position, not a
    # meaningful quantity on their own -- only relative ribbon thickness
    # is -- so the axis is left unlabeled with numbers hidden entirely.
    d <- roleData %>%
        filter(NoiseLevel %in% sigmas) %>%
        left_join(classData, by = c("ParamID", "NoiseLevel")) %>%
        group_by(NoiseLevel, Role2, ClassLabel) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        filter(MRT > 0) %>%
        group_by(NoiseLevel) %>%
        mutate(Frac = MRT / sum(MRT)) %>%
        ungroup() %>%
        mutate(Role2      = factor(Role2, levels = roleLevels),
               NoiseLevel = factor(NoiseLevel, levels = sigmas))

    p <- ggplot(d, aes(axis1 = NoiseLevel, axis2 = Role2, y = Frac)) +
        ggalluvial::geom_alluvium(aes(fill = ClassLabel), alpha = 0.8) +
        ggalluvial::geom_stratum(fill = "grey90", color = "grey30") +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5) +
        scale_x_discrete(limits = c("Noise Level", "Accessed state"), expand = c(0.12, 0.12)) +
        scale_fill_manual(values = classColors, name = "Multistability\n(# fates accessed)") +
        theme_Publication(base_size = 18) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.y = element_blank(), axis.ticks.y = element_blank())

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Plot E-style: fraction of monostable-double-high ParamIDs in each
# compute_nfates_class() multistability class, by noise level -- the same
# red/blue/green classification driving plot_double_high_alluvial()'s
# ribbon colors, shown here as its own population-level summary.
plot_nfates_class_summary <- function(classData, sigmas,
                                       classColors = NFATES_CLASS_COLORS,
                                       outFile = NULL, width = 7, height = 5) {

    d <- classData %>%
        filter(NoiseLevel %in% sigmas) %>%
        mutate(NoiseLevel = factor(NoiseLevel, levels = sigmas)) %>%
        count(NoiseLevel, ClassLabel, .drop = FALSE) %>%
        group_by(NoiseLevel) %>%
        mutate(Frac = n / sum(n)) %>%
        ungroup()

    p <- ggplot(d, aes(x = NoiseLevel, y = Frac, fill = ClassLabel)) +
        geom_col(position = "stack", width = 0.7) +
        scale_fill_manual(values = classColors, name = "Multistability\n(# fates accessed)") +
        theme_Publication() +
        labs(x = "Noise Level", y = "Fraction of parameter sets")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Distribution of MeanSwitches (mean, across the noise simulation's
# replicate ICs, of the number of discrete-state transitions during the
# recorded window) by noise level, faceted by ParamType -- "how stable are
# the observed states" companion to the MRT-based panels above. The
# recorded window already excludes a burn-in: the systematic pipeline
# (scripts/3_analyze_single_parameter.jl) runs each simulation over
# tspan = (0, 1000) and discards the first cut_fraction = 0.5 (never
# overridden) before counting anything, so every switch count here already
# reflects only the back half of the run, not the full trajectory.
# MeanSwitches/StdSwitches are written once per (ParamID, NoiseLevel, DT)
# -- identical across every State row for that group (they're not
# state-specific) -- so distinct() first, or every violin would silently
# overweight ParamIDs by however many State rows they happen to have.
plot_switching_stability <- function(net, noiseType, resultsFolder, dt = 0.01,
                                      outFile = NULL, width = 10, height = 4) {

    effectiveNoiseType <- noise_type_for(net, noiseType)
    f <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
    d <- read_csv(f, show_col_types = FALSE) %>%
        filter(DT == dt) %>%
        distinct(ParamID, ParamType, NoiseLevel, MeanSwitches, StdSwitches) %>%
        mutate(NoiseLevel = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
               ParamType  = str_to_sentence(ParamType))

    p <- ggplot(d, aes(x = NoiseLevel, y = MeanSwitches)) +
        geom_violin(fill = "steelblue", alpha = 0.5, scale = "width", color = NA) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6, fill = "white", color = "black") +
        facet_wrap(vars(ParamType)) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = "Mean # switching events")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Number of saved timesteps surviving the pipeline's burn-in cut, for
# normalizing MeanSwitches into a rate. tspan = (0, 1000), saveat = 1.0
# (scripts/3_analyze_single_parameter.jl:148) -> 1000 saved timesteps;
# cut_fraction = 0.5 (src/StochasticSimulations.jl:491,542, never
# overridden by analyze_noise_effects() or any CLI arg in
# scripts/2_submit_jobs.jl or a real submitted cluster array script,
# final/Additive/NF/logs/submit_array_1.sh) discards the first half ->
# 500 remain. (Not 700/30% -- no 0.3 cut_fraction exists anywhere in this
# codebase; checked directly against the above cluster job script.)
POST_BURNIN_TIMESTEPS <- 500

# plot_switching_stability(), generalized across networks AND noise types
# in one panel: fill = NoiseType (e.g. Additive/Multiplicative/Fluctuating,
# this project's 3 canonical types -- MultiplicativeInvLambda is not a 4th
# type, just the corrected data source noise_type_for() substitutes in for
# self-activation networks under "Multiplicative", applied automatically
# here the same way), facet_grid rows = StabilityClass (ParamType),
# cols = Network (dropped when only one network is requested -- a single
# redundant facet column showing just that one network's name adds nothing
# over a plain facet_wrap by StabilityClass alone). MeanSwitches is
# normalized by n_timesteps (the number of saved post-burn-in timesteps,
# POST_BURNIN_TIMESTEPS by default) into a per-timestep switching rate, so
# networks/noise types with different switch-count scales are directly
# comparable on one shared y axis. max_sigma optionally restricts to noise
# levels <= that value (e.g. dropping the far tail -- 0.5, 1 -- to focus
# on the more informative low/mid range).
plot_switching_stability_combined <- function(nets, noiseTypes, resultsFolder, dt = 0.01,
                                               n_timesteps = POST_BURNIN_TIMESTEPS, max_sigma = NULL,
                                               outFile = NULL, width = 12, height = 8) {

    d <- map_dfr(nets, function(net) {
        map_dfr(noiseTypes, function(nt) {
            effectiveNoiseType <- noise_type_for(net, nt)
            f <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
            if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
            read_csv(f, show_col_types = FALSE) %>%
                filter(DT == dt) %>%
                distinct(ParamID, ParamType, NoiseLevel, MeanSwitches) %>%
                mutate(Network = net, NoiseType = nt)
        })
    })
    if (nrow(d) == 0) stop("No data found for requested nets/noiseTypes")

    if (!is.null(max_sigma)) d <- d %>% filter(as.numeric(as.character(NoiseLevel)) <= max_sigma)

    d <- d %>%
        mutate(SwitchRate      = MeanSwitches / n_timesteps,
               NoiseLevel      = factor(NoiseLevel, levels = sort(unique(as.numeric(as.character(NoiseLevel))))),
               StabilityClass  = str_to_sentence(ParamType),
               Network         = factor(Network, levels = nets),
               NoiseType       = factor(NoiseType, levels = noiseTypes))

    p <- ggplot(d, aes(x = NoiseLevel, y = SwitchRate, fill = NoiseType)) +
        geom_violin(alpha = 0.5, scale = "width", color = NA, position = position_dodge(width = 0.8)) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 1.2, fill = "white", color = "black",
                     position = position_dodge(width = 0.8)) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
        labs(x = "Noise Level", y = "Switching events per timestep", fill = "Noise Type")

    p <- p + if (length(nets) == 1) {
        facet_wrap(vars(StabilityClass), scales = "free_y")
    } else {
        facet_grid(rows = vars(StabilityClass), cols = vars(Network), scales = "free_y")
    }

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Frequency of the 3 classify_double_high_origin() classes among
# double-high-only deterministic (NoiseLevel == 0) ParamIDs -- the
# population-level companion to compute_access_breakdown(), answering "how
# common is each stability pattern in the first place" ahead of
# plot_access_nFates_combined()'s "how many fates does each pattern access
# under noise".
compute_stability_class_freq <- function(net, noiseType, resultsFolder, dataFolder, dt = 0.01) {
    effectiveNoiseType <- noise_type_for(net, noiseType)
    f <- file.path(resultsFolder, effectiveNoiseType, net, "results", "all_parameters_results.csv")
    d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
    d <- fill_mrt(d, unique(d$State)) %>% filter(!is.na(Original))

    d %>%
        distinct(ParamID, Original) %>%
        mutate(StabilityClass = map_chr(Original, classify_double_high_origin)) %>%
        filter(!is.na(StabilityClass)) %>%
        count(StabilityClass, name = "n") %>%
        mutate(Frac = n / sum(n), Network = net)
}

plot_stability_class_freq <- function(freqData,
                                       stabilityLevels = c("Monostable (double-high)", "Bistable (mirror)", "Bistable (non-mirror)"),
                                       outFile = NULL, width = 5, height = 4) {

    d <- freqData %>% mutate(StabilityClass = factor(StabilityClass, levels = stabilityLevels))

    p <- ggplot(d, aes(x = StabilityClass, y = Frac)) +
        geom_col(fill = "steelblue") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
        labs(x = NULL, y = "Fraction of double-high ParamIDs\n(deterministic)")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# compute_nfates_class(), run for all 3 stability classes at once and
# stacked into one long table -- the monostable class uses
# classify_double_high_role()'s Hamming-distance SH1/SH3 buckets (3 fate
# roles, so up to "Tristable"), the two bistable classes use
# classify_bistable_role()'s 4 raw single-high letters (so up to
# "Tetrastable"). Feeds plot_nfates_class_summary_combined(); each
# stability class keeps its own natural fate-bucket count rather than
# forcing all 3 onto a shared scheme.
compute_nfates_class_all <- function(accessData, ratio_cutoff = 0.1) {
    roleMono <- classify_double_high_role(accessData)
    classMono <- compute_nfates_class(roleMono, fateRoles = c("SH1_1", "SH1_2", "SH3"), ratio_cutoff = ratio_cutoff) %>%
        mutate(StabilityClass = "Monostable (double-high)")

    bind_rows(lapply(c("Bistable (mirror)", "Bistable (non-mirror)"), function(sc) {
        roleData <- classify_bistable_role(accessData, sc)
        compute_nfates_class(roleData, fateRoles = c("A", "B", "C", "D"), ratio_cutoff = ratio_cutoff) %>%
            mutate(StabilityClass = sc)
    })) %>%
        bind_rows(classMono)
}

# plot_nfates_class_summary(), combined across all 3 stability classes
# (one facet each) instead of just the pooled monostable-double-high
# population -- same red/blue/green(/purple) multistability classification,
# now directly comparable across every stability pattern, not just the one
# feeding the noise-level alluvial.
plot_nfates_class_summary_combined <- function(classDataAll, sigmas,
                                                stabilityLevels = c("Monostable (double-high)", "Bistable (mirror)", "Bistable (non-mirror)"),
                                                classColors = NFATES_CLASS_COLORS,
                                                outFile = NULL, width = 12, height = 4) {

    d <- classDataAll %>%
        filter(NoiseLevel %in% sigmas) %>%
        mutate(NoiseLevel = factor(NoiseLevel, levels = sigmas),
               StabilityClass = factor(StabilityClass, levels = stabilityLevels)) %>%
        count(StabilityClass, NoiseLevel, ClassLabel, .drop = FALSE) %>%
        group_by(StabilityClass, NoiseLevel) %>%
        mutate(Frac = n / sum(n)) %>%
        ungroup()

    p <- ggplot(d, aes(x = NoiseLevel, y = Frac, fill = ClassLabel)) +
        geom_col(position = "stack", width = 0.7) +
        facet_wrap(vars(StabilityClass), nrow = 1) +
        scale_fill_manual(values = classColors, name = "Multistability\n(# fates accessed)", drop = FALSE) +
        theme_Publication() +
        labs(x = "Noise Level", y = "Fraction of parameter sets")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Distribution of per-ParamID MRT at each destination Role, for ONE example
# specific double-high attractor (OriginalLabel) at one noise level --
# "what do the transitions for a single example monostable double-high
# attractor actually look like", shown as the actual per-parameter-set MRT
# distribution (a violin per destination) rather than a flow-summed
# aggregate, so the numbers stay directly interpretable (MRT in [0, 1] per
# ParamID) regardless of how many ParamIDs happen to share that attractor.
plot_example_mrt_violin <- function(accessData, originalLabel, sigma = 0.01,
                                     outFile = NULL, width = 7, height = 4.5) {

    d <- accessData %>%
        filter(StabilityClass == "Monostable (double-high)",
               OriginalLabel == originalLabel, NoiseLevel == sigma) %>%
        mutate(Role = if_else(Role == "Original", "Self", Role))

    roleOrder <- c(setdiff(unique(d$Role), c("Self", "Other")), "Self", "Other")
    d <- d %>% mutate(Role = factor(Role, levels = roleOrder))

    p <- ggplot(d, aes(x = Role, y = MRT)) +
        geom_violin(fill = "steelblue", alpha = 0.5, scale = "width", color = NA) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6, fill = "white", color = "black") +
        theme_Publication() +
        labs(x = "Accessed state", y = "MRT (per parameter set)",
             title = paste0(originalLabel, " (noise level = ", sigma, ")"))

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# Circular chord diagram (circlize) of MRT flow from each specific
# monostable double-high attractor to each single-high fate (A/B/C/D
# only -- "Self"/"Other" dropped, since this diagram is specifically about
# the double-high -> single-high transition, not the "stayed put"/"went
# elsewhere" residual), at one fixed noise level. circlize places every
# sector on a single ring with a gap between the double-high group and the
# single-high group (the standard bipartite chord-diagram convention) --
# not on two literal concentric radii. Requires library(circlize) loaded
# by the caller; draws via base graphics (not ggplot2), so outFile is
# written with png()/dev.off() rather than ggsave() -- PNG rather than
# JPEG so a caller can read it back losslessly (e.g. via png::readPNG()) to
# embed it as a grob in a larger cowplot/ggplot composite, since circlize's
# base-graphics output can't be captured directly via grid::grid.grabExpr().
plot_double_high_chord <- function(accessData, sigma = 0.01, outFile = NULL, width = 8, height = 8) {

    d <- accessData %>%
        filter(StabilityClass == "Monostable (double-high)", NoiseLevel == sigma,
               Role %in% c("A", "B", "C", "D")) %>%
        group_by(OriginalLabel, Role) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        filter(MRT > 0)

    origOrder <- d %>% group_by(OriginalLabel) %>% summarise(Tot = sum(MRT), .groups = "drop") %>%
        arrange(desc(Tot)) %>% pull(OriginalLabel)
    roleOrder <- c("A", "B", "C", "D")
    sectorOrder <- c(origOrder, roleOrder)

    grid.col <- setNames(
        c(scales::hue_pal()(length(origOrder)), rep("grey70", length(roleOrder))),
        sectorOrder
    )

    render <- function() {
        circos.clear()
        chordDiagram(as.data.frame(d), order = sectorOrder, grid.col = grid.col,
                     annotationTrack = "grid", preAllocateTracks = 1, big.gap = 15)
        circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
            xlim <- get.cell.meta.data("xlim")
            ylim <- get.cell.meta.data("ylim")
            sector.name <- get.cell.meta.data("sector.index")
            circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
                        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
        })
    }

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        png(outFile, width = width, height = height, units = "in", res = 300, bg = "white")
        render()
        dev.off()
    } else {
        render()
    }
    invisible(NULL)
}

# Same double-high -> single-high MRT flow as plot_double_high_chord(),
# rendered as a circular network layout (ggraph/tidygraph) instead of a
# chord diagram: nodes (6 double-high + 4 single-high) placed evenly
# around one circle, edges drawn as arcs with width proportional to MRT.
# Like the chord diagram, this places every node on one ring (grouped),
# not on two literal concentric radii -- ggraph's circular layouts don't
# support that directly. Requires library(ggraph) + library(tidygraph)
# loaded by the caller.
plot_double_high_circular_alluvial <- function(accessData, sigma = 0.01, outFile = NULL, width = 8, height = 8) {

    d <- accessData %>%
        filter(StabilityClass == "Monostable (double-high)", NoiseLevel == sigma,
               Role %in% c("A", "B", "C", "D")) %>%
        group_by(OriginalLabel, Role) %>%
        summarise(MRT = sum(MRT), .groups = "drop") %>%
        filter(MRT > 0)

    origOrder <- d %>% group_by(OriginalLabel) %>% summarise(Tot = sum(MRT), .groups = "drop") %>%
        arrange(desc(Tot)) %>% pull(OriginalLabel)
    roleOrder <- c("A", "B", "C", "D")

    nodes <- tibble(name = c(origOrder, roleOrder),
                     Type = c(rep("Double-high", length(origOrder)), rep("Single-high", length(roleOrder))))
    edges <- d %>% transmute(from = match(OriginalLabel, nodes$name), to = match(Role, nodes$name), MRT)

    g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

    p <- ggraph(g, layout = "linear", circular = TRUE) +
        geom_edge_arc(aes(width = MRT, color = MRT), alpha = 0.7, strength = 0.6, lineend = "round") +
        scale_edge_width(range = c(0.2, 4), guide = "none") +
        scale_edge_color_viridis(name = "MRT", guide = guide_edge_colorbar()) +
        geom_node_point(aes(color = Type), size = 4) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
        scale_color_manual(values = c("Double-high" = "#e41a1c", "Single-high" = "grey40"), name = NULL) +
        coord_fixed() +
        theme_void()

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# ── Illustrative lambda(t) trajectories (Figure1 Panel D + Figure1_supp
#    Panel B) -- self-contained simulation (no RACIPE run needed): these
#    update rules are exactly the three defined in Panel C / used throughout
#    StochasticSimulations.jl, just applied here to a bare scalar lambda
#    instead of an edge weight inside a full ODE integration. ──

# lambda0 defaults to 20% of the way from 0 to hi -- "start from 0.2" for
# an inhibition-range call (hi ~ 1) is lambda0 ~ 0.2 unchanged; the same
# rule automatically lands at a valid, comparable starting point (~20) for
# an activation-range call (hi ~ 100), rather than 0.2 itself, which would
# be clamped to lo on step 1 and never move.
# effective_sigma mirrors this project's real noise scaling (documented in
# Panel C: sigma_eff = sigma * max_lambda for Act-type edges, whose lambda
# ranges over ~[1, 100], vs. sigma itself for Inh-type edges, whose lambda
# is already in [0, 1]) -- so the same nominal sigma (e.g. 0.06) produces
# comparable relative drift whether lo/hi describe an inhibition or an
# activation edge, with no extra argument needed at the call site.
simulate_lambda <- function(mode, lambda0 = 0.2 * hi, sigma = 0.06, effective_sigma = TRUE,
                             t_max = 100, steps_per_unit = 100, n_traj = 100,
                             lo = 0.001, hi = 0.999, seed = 1, useInvLambda = FALSE) {
    set.seed(seed)
    sigma_use <- if (effective_sigma && hi > 1) sigma * hi else sigma
    dt <- 1 / steps_per_unit
    times <- seq(0, t_max, by = dt)
    map_dfr(seq_len(n_traj), function(traj) {
        lam <- numeric(length(times)); lam[1] <- lambda0
        for (k in 2:length(times)) {
            lam[k] <- switch(mode,
                Stationary     = lambda0 + rnorm(1, 0, sigma_use),
                Additive       = lam[k - 1] + rnorm(1, 0, sigma_use),
                # Mirrors lambda_sampler.jl's "MultiplicativeInvLambda" branch:
                # for activation-range edges (hi > 1) the noise kick is
                # divided by the CURRENT lambda instead of applied flat, so
                # it shrinks once lambda has already strengthened.
                Multiplicative = if (useInvLambda && hi > 1) {
                    lam[k - 1] * (1 + rnorm(1, 0, sigma_use) / lam[k - 1])
                } else {
                    lam[k - 1] * (1 + rnorm(1, 0, sigma_use))
                })
            lam[k] <- min(max(lam[k], lo), hi)
        }
        tibble(Trajectory = traj, Time = times, Lambda = lam)
    }) %>% mutate(Mode = mode)
}

# Runs simulate_lambda() for each noise mode with n_stats trajectories (the
# sample the mean/SD ribbon is computed from), but only draws n_display of
# them as individual black dashed lines -- showing every one of 100
# trajectories individually would just be visual noise, but the ribbon
# still needs the full sample to be a meaningful summary.
plot_lambda_trajectories <- function(lambda0 = 0.2 * hi, sigma = 0.06, lo = 0.001, hi = 0.999,
                                      n_stats = 100, n_display = 5, display_every = 1,
                                      display_alpha = 0.7,
                                      modes = c("Stationary", "Multiplicative", "Additive"),
                                      seeds = seq_along(modes),
                                      outFile = NULL, width = 11, height = 4,
                                      useInvLambda = FALSE) {

    dD <- map2_dfr(modes, seeds, function(m, s) {
        simulate_lambda(m, lambda0 = lambda0, sigma = sigma, lo = lo, hi = hi,
                         n_traj = n_stats, seed = s, useInvLambda = useInvLambda)
    }) %>% mutate(Mode = factor(Mode, levels = modes))

    dD_summ <- dD %>% group_by(Mode, Time) %>%
        summarise(Mean = mean(Lambda), SD = sd(Lambda), .groups = "drop")

    # display_every thins the individual trajectory lines down to one point
    # per `display_every` time units (default: every time unit, not every
    # simulation step) -- otherwise, at steps_per_unit = 100, each displayed
    # line is 10000 points of fine-grained noise that reads as a solid inked
    # block rather than a legible trajectory.
    dD_display <- dD %>% filter(Trajectory <= n_display)
    if (!is.null(display_every)) {
        dD_display <- dD_display %>%
            filter(abs(Time / display_every - round(Time / display_every)) < 1e-6)
    }

    p <- ggplot(dD, aes(x = Time, y = Lambda)) +
        geom_line(data = dD_display, aes(group = Trajectory), color = "black",
                  linewidth = 0.3, linetype = "21", alpha = display_alpha) +
        geom_ribbon(data = dD_summ, aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Mode), alpha = 0.35) +
        geom_line(data = dD_summ, aes(y = Mean, color = Mode), linewidth = 1.1, linetype = "21") +
        facet_wrap(vars(Mode), nrow = 1) +
        scale_fill_manual(values = c(Stationary = "#D62728", Multiplicative = "#2CA02C", Additive = "#1F77B4")) +
        scale_color_manual(values = c(Stationary = "#D62728", Multiplicative = "#2CA02C", Additive = "#1F77B4")) +
        labs(x = "Time (A.U.)", y = expression(lambda)) +
        theme_Publication() +
        theme(legend.position = "none")

    if (!is.null(outFile)) {
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = width, height = height)
    }
    invisible(p)
}

# ── Network topology diagrams (Figure1 Panel A + Figure5_supp Panel C) --
#    tikz-style circles + bar-headed inhibition edges + curved
#    self-activation loops, drawn straight from a .topo file. Layout only
#    defined for 2/3/4 nodes -- the only sizes used in this project's
#    topologies (toggle switch / triad / TT4-TT4SA-TS4's 4-node ring or
#    complete graph). Inhibition (Type 2) edges are drawn as a segment per
#    directed pair with a perpendicular bar near the target end, stopping
#    AT the bar rather than continuing on to touch the target node itself.
#    Self-activation (Type 1, Source == Target) is a small curved loop
#    bowing OUTWARD from the node's center -- bowing inward gets silently
#    hidden behind the node circle's opaque fill, drawn on top (verified by
#    checking rendered pixel colors directly, since the geoms otherwise
#    give no error or warning when this happens).
plot_network_diagram <- function(topoFile, node_radius = 0.28, bar_len = 0.16,
                                  bar_offset = 0.9, edge_color = "red", edge_width = 1.8,
                                  edge_offset = 0.16,
                                  loop_curvature = 1.3, loop_half_span = 40,
                                  node_text_size = 10) {
    topo <- read.table(topoFile, header = TRUE, stringsAsFactors = FALSE)
    nodes <- unique(c(topo$Source, topo$Target))
    n <- length(nodes)

    layout <- if (n == 2) {
        tibble(node = nodes, x = c(-1, 1), y = c(0, 0))
    } else if (n == 3) {
        tibble(node = nodes, x = c(0, -0.87, 0.87), y = c(1, -0.5, -0.5))
    } else if (n == 4) {
        tibble(node = nodes, x = c(-1, 1, 1, -1), y = c(1, 1, -1, -1))
    } else stop("plot_network_diagram: layout only defined for 2, 3, or 4 nodes")

    pos <- setNames(map(seq_len(n), ~ c(layout$x[.x], layout$y[.x])), layout$node)

    edge_point <- function(from, to, r) {
        d <- to - from; len <- sqrt(sum(d^2))
        from + d / len * r
    }

    # ---- inhibition edges: one directed edge per (source, target), running
    # from the source node's boundary and stopping AT the bar -- not
    # continuing on to actually touch the target node -- so the bar reads as
    # a real terminal arrowhead with a visible gap in front of the target,
    # the same way a normal arrow doesn't touch the box it points to. ----
    inhib <- topo %>% filter(Type == 2, Source != Target)
    pairKey <- function(a, b) paste(pmin(a, b), pmax(a, b))
    pairs <- inhib %>% mutate(key = pairKey(Source, Target)) %>% distinct(key) %>% pull(key)

    # Mutually-inhibiting pairs get two separate parallel edges (offset
    # perpendicular to the node-node axis) instead of one shared line with a
    # bar at both ends, which reads as one ambiguous double-headed edge
    # instead of two distinct inhibitions.
    seg_df <- list(); bar_df <- list()
    make_edge <- function(key, p_src, p_tgt, frac) {
        dir_local  <- p_tgt - p_src; dir_local <- dir_local / sqrt(sum(dir_local^2))
        perp_local <- c(-dir_local[2], dir_local[1])
        bar_pos <- p_src + (p_tgt - p_src) * frac
        seg_df[[key]] <<- tibble(x = p_src[1], y = p_src[2], xend = bar_pos[1], yend = bar_pos[2])
        bar_df[[key]] <<- tibble(
            x = bar_pos[1] - perp_local[1] * bar_len / 2, y = bar_pos[2] - perp_local[2] * bar_len / 2,
            xend = bar_pos[1] + perp_local[1] * bar_len / 2, yend = bar_pos[2] + perp_local[2] * bar_len / 2)
    }

    for (k in pairs) {
        ab <- strsplit(k, " ")[[1]]; a <- ab[1]; b <- ab[2]
        base_a <- edge_point(pos[[a]], pos[[b]], node_radius)
        base_b <- edge_point(pos[[b]], pos[[a]], node_radius)

        dir  <- (base_b - base_a); dir <- dir / sqrt(sum(dir^2))
        perp <- c(-dir[2], dir[1])

        has_ab <- nrow(inhib %>% filter(Source == a, Target == b)) > 0
        has_ba <- nrow(inhib %>% filter(Source == b, Target == a)) > 0

        if (has_ab && has_ba) {
            off <- perp * edge_offset / 2
            make_edge(paste0(k, "_ab"), base_a + off, base_b + off, bar_offset)  # a -> b
            make_edge(paste0(k, "_ba"), base_b - off, base_a - off, bar_offset)  # b -> a
        } else if (has_ab) {
            make_edge(paste0(k, "_ab"), base_a, base_b, bar_offset)
        } else if (has_ba) {
            make_edge(paste0(k, "_ba"), base_b, base_a, bar_offset)
        }
    }
    seg_df <- bind_rows(seg_df); bar_df <- bind_rows(bar_df)

    # ---- self-activation loops: two points +/- loop_half_span degrees
    # around the node's outward-facing direction, placed slightly OUTSIDE
    # the node's own circle (node_radius + loop_gap, not right on the
    # boundary) so the whole loop -- including its arrowhead -- sits fully
    # clear of the node's opaque fill (drawn on top of the loop) instead of
    # being partly covered by it, with a closed arrowhead pointing back
    # toward the node. Needs a wide angular span + high curvature or the
    # loop reads as merged into the node outline instead of a visibly
    # separate loop. "Outward-facing" is which side of the diagram the node
    # sits on: right-side nodes loop to their right, left-side nodes loop to
    # their left, and top nodes (x ~ 0, e.g. the triad's apex) keep looping
    # upward -- so loops on either side point away from the other nodes
    # instead of all defaulting to the top regardless of layout.
    self_acts <- topo %>% filter(Type == 1, Source == Target) %>% distinct(Source) %>% pull(Source)
    loop_gap <- 0.05
    self_df <- map_dfr(self_acts, function(nd) {
        p <- pos[[nd]]
        r <- node_radius + loop_gap
        center_ang <- if (p[1] < -0.1) 180 else if (p[1] > 0.1) 0 else 90
        lo <- center_ang - loop_half_span; hi <- center_ang + loop_half_span
        tibble(x    = p[1] + r * cos(lo * pi / 180),
               y    = p[2] + r * sin(lo * pi / 180),
               xend = p[1] + r * cos(hi * pi / 180),
               yend = p[2] + r * sin(hi * pi / 180))
    })

    node_df <- tibble(x = layout$x, y = layout$y, label = layout$node)

    p <- ggplot()
    if (nrow(seg_df) > 0) p <- p + geom_segment(data = seg_df, aes(x = x, y = y, xend = xend, yend = yend),
                                                 color = edge_color, linewidth = edge_width)
    if (nrow(bar_df) > 0) p <- p + geom_segment(data = bar_df, aes(x = x, y = y, xend = xend, yend = yend),
                                                 color = edge_color, linewidth = edge_width)
    if (nrow(self_df) > 0) p <- p + geom_curve(data = self_df, aes(x = x, y = y, xend = xend, yend = yend),
                                                curvature = loop_curvature, color = edge_color, linewidth = edge_width,
                                                arrow = arrow(length = unit(0.15, "cm"), type = "closed"))
    p <- p +
        ggforce::geom_circle(data = node_df, aes(x0 = x, y0 = y, r = node_radius),
                              fill = "white", color = edge_color, linewidth = edge_width) +
        geom_text(data = node_df, aes(x = x, y = y, label = label), size = node_text_size, fontface = "bold") +
        coord_fixed(clip = "off") +
        theme_void() +
        theme(plot.margin = margin(10, 10, 10, 10))
    p
}

# Reachability scatter (Max - Threshold per node), colored by MRT of the
# given state. 2-node networks only (needs exactly two "_diff" columns).
plot_reach_scatter <- function(net, noiseType, dt, dataFolder, resultsFolder,
                                state, sigma = 0.1, param_type = NULL,
                                reachableOnly = FALSE,
                                outFile = NULL, panelWidth = 5, panelHeight = 5, ncol = 3) {

	df <- add_reachability(net, noiseType, dt, dataFolder, resultsFolder) %>%
		filter(NoiseLevel == sigma, State == state)
	if (!is.null(param_type)) df <- df %>% filter(ParamType %in% param_type)
	if (reachableOnly) {
		reachCols <- names(df)[str_ends(names(df), "_reach")]
		df <- df %>% filter(if_all(all_of(reachCols), ~ . == 1))
	}
	if (nrow(df) == 0) stop("No rows remain after filtering")

	nodes <- names(df)[str_ends(names(df), "_diff")] %>% str_remove("_diff")
	if (length(nodes) != 2) stop("Scatter supports 2-node networks only; got: ", paste(nodes, collapse = ", "))
	xcol <- paste0(nodes[1], "_diff"); ycol <- paste0(nodes[2], "_diff")

	df <- df %>% mutate(ParamType = str_to_sentence(ParamType))

	p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]], color = MRT)) +
		geom_point() +
		geom_hline(yintercept = 0, color = "red") +
		geom_vline(xintercept = 0, color = "red") +
		scale_color_viridis_c(breaks = scales::breaks_pretty(n = 3)) +
		scale_x_continuous(trans = scales::pseudo_log_trans(base = 2), breaks = c(-50, 0, 10, 100, 1000)) +
		scale_y_continuous(trans = scales::pseudo_log_trans(base = 2), breaks = c(-50, 0, 10, 100, 1000)) +
		facet_wrap(vars(ParamType), ncol = ncol) +
		theme_Publication() +
		theme(legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.4, "cm")) +
		labs(x = paste0("Max ", nodes[1], " - Threshold ", nodes[1]),
			 y = paste0("Max ", nodes[2], " - Threshold ", nodes[2]),
			 color = paste0("MRT of ", state))

	if (!is.null(outFile)) {
		n_panels <- length(unique(df$ParamType))
		dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
		ggsave(outFile, p, width = panelWidth * min(n_panels, ncol),
			   height = panelHeight * ceiling(n_panels / ncol))
	}
	invisible(p)
}

# Per-ParamID MRT (summed over the member states of one state class) vs.
# mean number of state switches per stochastic trajectory, normalized by the
# number of saved timesteps. all_parameters_transitions.csv's Count column is
# already the mean count of one specific (from, to) transition pair across
# the num_sims replicate trajectories (aggregate_transitions() in
# StochasticSimulations.jl averages, doesn't sum, across replicates) -- so
# summing Count over every transition pair for a ParamID gives the mean
# TOTAL number of switches per trajectory (by linearity of expectation, this
# holds regardless of how the transition types are correlated within a
# single run). n_timesteps defaults to 1000, matching this project's actual
# simulation settings: analyze_noise_effects()/3_analyze_single_parameter.jl
# hardcode tspan = (0.0, 1000.0) and saveat = 1.0 regardless of the (unused)
# "saveat" CLI arg threaded through 0_run_all.jl -- i.e. 1000 saved
# timesteps per trajectory.
compute_mrt_switches_data <- function(nets, noiseType, resultsFolder, dt = 0.01,
                                       stateClass = c("all-high", "single-high", "double-high"),
                                       n_timesteps = 1000) {
    stateClass <- match.arg(stateClass)

    dMRT <- map_dfr(nets, function(net) {
        f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_results.csv")
        if (!file.exists(f)) { warning("Missing results: ", net); return(NULL) }
        d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt)
        fill_mrt(d, unique(d$State)) %>%
            filter(StateClass == stateClass) %>%
            group_by(ParamID, ParamType, NoiseLevel) %>%
            summarise(MRT = sum(MRT), .groups = "drop") %>%
            mutate(Network = net)
    })
    if (nrow(dMRT) == 0) stop("No MRT rows found for state class: ", stateClass)

    dSwitch <- map_dfr(nets, function(net) {
        f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_transitions.csv")
        if (!file.exists(f)) { warning("Missing transitions: ", net); return(NULL) }
        read_csv(f, show_col_types = FALSE) %>%
            filter(DT == dt) %>%
            group_by(ParamID, NoiseLevel) %>%
            summarise(Switches = sum(Count), .groups = "drop") %>%
            mutate(Network = net)
    })
    if (nrow(dSwitch) == 0) stop("No transition rows found")

    # A handful of ParamIDs (varies by network/noise mode, worst under
    # Fluctuating noise) report more switches than a single 1000-step
    # trajectory can physically contain -- confirmed against the project's
    # actual simulation settings (tspan = (0, 1000), saveat = 1.0 uniformly
    # across noise modes) that this is a data artifact, not a differing
    # timestep count, so NormSwitches is clipped at 1.0 rather than plotted
    # (or trusted) at face value.
    dMRT %>%
        left_join(dSwitch, by = c("ParamID", "NoiseLevel", "Network")) %>%
        mutate(Switches      = replace_na(Switches, 0),
               NormSwitches  = pmin(Switches / n_timesteps, 1))
}

# Reachability-style scatter: MRT (x) vs. normalized mean switches (y), one
# point per ParamID, faceted by network (rows) x noise level (cols). One
# state class and noise mode per call/output file -- MRT_vs_Switches.r loops
# over both to produce the full set.
plot_mrt_vs_switches <- function(nets, noiseType, resultsFolder, dt = 0.01,
                                  stateClass = c("all-high", "single-high", "double-high"),
                                  noise_levels = c(0.001, 0.01, 0.1, 1.0),
                                  n_timesteps = 1000,
                                  outFile = NULL, panelWidth = 3.2, panelHeight = 3.2) {
    stateClass <- match.arg(stateClass)

    d <- compute_mrt_switches_data(nets, noiseType, resultsFolder, dt = dt,
                                    stateClass = stateClass, n_timesteps = n_timesteps) %>%
        filter(NoiseLevel %in% noise_levels) %>%
        mutate(Network    = factor(Network, levels = nets),
               NoiseLevel = factor(NoiseLevel, levels = sort(noise_levels)))
    if (nrow(d) == 0) stop("No rows remain -- check state class exists for these networks/noise levels")

    p <- ggplot(d, aes(x = MRT, y = NormSwitches)) +
        geom_point(alpha = 0.4, size = 1, color = "#e41a1c") +
        facet_grid(rows = vars(Network), cols = vars(NoiseLevel)) +
        theme_Publication() +
        labs(x = paste0("MRT (", stateClass, ")"), y = "Mean switches per parameter set (normalized)",
             title = paste0(noiseType, " noise — ", stateClass))

    if (!is.null(outFile)) {
        n_col <- length(unique(droplevels(d$NoiseLevel)))
        n_row <- length(unique(droplevels(d$Network)))
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = panelWidth * n_col, height = panelHeight * n_row)
    }
    invisible(p)
}

# Same data as plot_mrt_vs_switches(), rendered as a 2D density heatmap
# (geom_bin2d, fill = per-panel bin density) instead of a point scatter --
# more legible when a facet is dominated by thousands of overplotted
# points. bins is the number of bins along each axis (both MRT and
# NormSwitches run over a fixed [0, 1] range, so a shared bin count keeps
# bin size comparable between axes).
plot_mrt_vs_switches_heatmap <- function(nets, noiseType, resultsFolder, dt = 0.01,
                                          stateClass = c("all-high", "single-high", "double-high"),
                                          noise_levels = c(0.001, 0.01, 0.1, 1.0),
                                          n_timesteps = 1000, bins = 30,
                                          outFile = NULL, panelWidth = 3.2, panelHeight = 3.2) {
    stateClass <- match.arg(stateClass)

    d <- compute_mrt_switches_data(nets, noiseType, resultsFolder, dt = dt,
                                    stateClass = stateClass, n_timesteps = n_timesteps) %>%
        filter(NoiseLevel %in% noise_levels) %>%
        mutate(Network    = factor(Network, levels = nets),
               NoiseLevel = factor(NoiseLevel, levels = sort(noise_levels)))
    if (nrow(d) == 0) stop("No rows remain -- check state class exists for these networks/noise levels")

    # log10-scaled fill: raw bin density is dominated by a handful of huge
    # spikes right at MRT/NormSwitches ~ 0 or 1 (most parameter sets barely
    # switch at all), which on a linear scale saturates the palette and
    # flattens every other bin to the same dark color. log1p barely helped
    # since per-bin density values are already small fractions (log1p(x) ~
    # x for small x); log10 (every plotted bin has count >= 1, so density >
    # 0 -- no zero-handling issue) spreads the whole multi-order-of-
    # magnitude range out instead.
    p <- ggplot(d, aes(x = MRT, y = NormSwitches)) +
        geom_bin2d(aes(fill = after_stat(density)), bins = bins) +
        scale_fill_viridis_c(name = "Density", trans = "log10",
                              breaks = scales::breaks_log(n = 4),
                              labels = scales::label_number(accuracy = 0.01)) +
        facet_grid(rows = vars(Network), cols = vars(NoiseLevel)) +
        theme_Publication() +
        theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(1.2, "cm"),
              legend.text = element_text(size = rel(0.7))) +
        labs(x = paste0("MRT (", stateClass, ")"), y = "Mean switches per parameter set (normalized)",
             title = paste0(noiseType, " noise — ", stateClass))

    if (!is.null(outFile)) {
        n_col <- length(unique(droplevels(d$NoiseLevel)))
        n_row <- length(unique(droplevels(d$Network)))
        dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)
        ggsave(outFile, p, width = panelWidth * n_col, height = panelHeight * n_row)
    }
    invisible(p)
}
