# Distribution across parameter sets of one node's mean expression under every condition (C0-C4),
# from the collected expression summaries (per parameter set and condition cell: mean of the
# trajectories' time-averaged expression, and mean within-trajectory SD).
#   Rscript expression_summaries.r <collected_dir> <racipe_data_dir> <output_dir> [network] [node]

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
collected <- args[1]; dataDir <- args[2]; outDir <- args[3]
net <- if (length(args) >= 4) args[4] else "TT"
node <- if (length(args) >= 5) args[5] else "A"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# --- threshold and maximum level of the node, from the RACIPE files ------------------
prs <- fread(file.path(dataDir, paste0(net, ".prs")))$Parameter
nodes <- sub("Prod_of_", "", grep("^Prod_of_", prs, value = TRUE))
sol <- fread(file.path(dataDir, paste0(net, "_solution.dat")), header = FALSE)
setnames(sol, c("ParamID", "NumStates", "Frequency", nodes))
theta <- 2^(sum(sol$Frequency * sol[[node]]) / sum(sol$Frequency))   # weighted geometric mean
par <- fread(file.path(dataDir, paste0(net, "_parameters.dat")), header = FALSE)
setnames(par, c("ParamID", "NumStates", prs))
gk <- par[, .(ParamID, GK = get(paste0("Prod_of_", node)) / get(paste0("Deg_of_", node)))]
cat(sprintf("%s node %s: threshold θ = %.3f\n", net, node, theta))

# --- expression summaries -------------------------------------------------------------
files <- list.files(collected, pattern = paste0("^", net, "_(Additive|Multiplicative|Stationary)_expression\\.csv\\.gz$"), full.names = TRUE)
e <- rbindlist(lapply(files, fread))[Node == node]
e <- merge(e, gk, by = "ParamID")
ctrlSigmas <- sort(unique(e[Condition == "C3", Sigma]))
poolDT <- unique(e[Condition == "C3", DT])
e <- e[Condition == "C0" | (Sigma %in% ctrlSigmas & DT == poolDT)]
# as.character, not format(): format() pads a vector to a common width ("  0.01")
e[, Cond := fifelse(Condition == "C3", paste0("C3 T=", as.character(THold)), Condition)]
c0 <- e[Condition == "C0"]
e <- rbind(e[Condition != "C0"], rbindlist(lapply(ctrlSigmas, function(s) copy(c0)[, Sigma := s])))

condLevels <- c("C0", "C4", "C2", "C3 T=100", "C3 T=10", "C3 T=1", "C3 T=0.1", "C3 T=0.01", "C1")
condLabels <- c("C0 original", "C4 mean λ", "C2 frozen", "C3 T=100", "C3 T=10", "C3 T=1",
				"C3 T=0.1", "C3 T=0.01", "C1 stochastic")
famColors <- c(Reference = "#8a8984", Frozen = "#eb6834", Hold = "#1baf7a", Stochastic = "#2a78d6")
e[, Cond := factor(Cond, levels = condLevels, labels = condLabels)]
e[, Family := fcase(Cond == "C0 original", "Reference", Cond %in% c("C4 mean λ", "C2 frozen"), "Frozen",
					Cond == "C1 stochastic", "Stochastic", default = "Hold")]
e[, SigmaLab := factor(paste0("σ = ", Sigma), levels = paste0("σ = ", ctrlSigmas))]
e[, RelThr := log2(pmax(Mean, 1e-9) / theta)]
e[, FracMax := Mean / GK]
e[, SDFrac := MeanSD / GK]

theme_dist <- function() theme_bw(base_size = 10) +
	theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
		  panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
		  strip.background = element_rect(fill = "grey96", colour = "grey80"),
		  strip.text.y.left = element_text(angle = 0, hjust = 1), strip.placement = "outside",
		  axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",
		  plot.title = element_text(face = "bold"), plot.caption = element_text(colour = "grey35", hjust = 0))

ridge <- function(x, var, xlab, ref, xlim, title, cap, bnd = c(-Inf, Inf)) {
	out <- mean(x[[var]] < xlim[1] | x[[var]] > xlim[2])
	if (out > 0) cap <- paste0(cap, sprintf("\n%.1f%% of values fall outside the plotted range (density computed on all values).", 100 * out))
	ggplot(x, aes(.data[[var]], fill = Family)) +
		ref +
		geom_density(aes(y = after_stat(scaled)), alpha = 0.75, colour = NA, bounds = bnd) +
		facet_grid(Cond ~ SigmaLab, switch = "y") +
		scale_fill_manual(values = famColors) +
		coord_cartesian(xlim = xlim) +
		labs(x = xlab, y = NULL, title = title, caption = cap) +
		theme_dist()
}

for (mode in unique(e$NoiseMode)) {
	x <- e[NoiseMode == mode]
	n <- uniqueN(x$ParamID)
	tag <- paste0(net, "_node", node, "_", tolower(substr(mode, 1, 4)))
	base_cap <- paste0("Distribution over ", n, " parameter sets of node ", node, "'s mean expression (time average over t = 500–1000, ",
					   "averaged over replicates). Pilot run; densities scaled to the same height.")
	plots <- list(
		rel = ridge(x, "RelThr", paste0("log2(mean expression of ", node, " / threshold)   (0 = high/low threshold)"),
					geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4), c(-10, 7),
					paste0(net, " node ", node, " — ", mode, " noise: mean expression relative to threshold"), base_cap),
		frac = ridge(x, "FracMax", paste0("Mean expression of ", node, " as a fraction of its maximum (g/k)"),
					 geom_vline(xintercept = 0.505^2, colour = "grey30", linewidth = 0.4, linetype = "dashed"), c(0, 1),
					 paste0(net, " node ", node, " — ", mode, " noise: mean expression as a fraction of maximum"),
					 paste0(base_cap, "\nDashed line: 0.505² = 0.255, the lowest level ", node, " can reach with both inhibitory λ at the additive pool mean."),
					 bnd = c(0, 1)),
		sd = ridge(x, "SDFrac", paste0("Within-trajectory SD of ", node, " as a fraction of its maximum (g/k)"),
				   geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3), c(0, 0.4),
				   paste0(net, " node ", node, " — ", mode, " noise: size of temporal fluctuations"),
				   paste0("Distribution over ", n, " parameter sets of the mean within-trajectory SD (0 for constant expression)."),
				   bnd = c(0, Inf))
	)
	for (k in names(plots)) {
		ggsave(file.path(outDir, paste0(tag, "_", k, ".pdf")), plots[[k]], width = 10, height = 9, device = cairo_pdf)
		ggsave(file.path(outDir, paste0(tag, "_", k, ".png")), plots[[k]], width = 10, height = 9, dpi = 200)
		cat("wrote", paste0(tag, "_", k), "\n")
	}
	s <- x[, .(MedianRelThr = median(RelThr), FracAboveThr = mean(RelThr > 0), MedianFracMax = median(FracMax),
				MedianSDFrac = median(SDFrac)), by = .(Sigma, Cond)][order(Sigma, Cond)]
	cat("\n", mode, " summary (node ", node, "):\n", sep = ""); print(s, digits = 3)
}
