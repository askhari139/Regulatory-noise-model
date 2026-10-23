# Per-replicate expression distributions and windowed log fold-change for one node, all conditions.
# Inputs: collected replicates / lfc tables (codes/regnoise/scripts/collect.jl) and the RACIPE files.
#   Rscript replicate_expression_figures.r <collected_dir> <racipe_data_dir> <output_dir> [network] [node]
# Each parameter set is weighted equally (its replicates share weight 1).

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
collected <- args[1]; dataDir <- args[2]; outDir <- args[3]
net <- if (length(args) >= 4) args[4] else "TT"
node <- if (length(args) >= 5) args[5] else "A"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

prs <- fread(file.path(dataDir, paste0(net, ".prs")))$Parameter
nodes <- sub("Prod_of_", "", grep("^Prod_of_", prs, value = TRUE))
sol <- fread(file.path(dataDir, paste0(net, "_solution.dat")), header = FALSE)
setnames(sol, c("ParamID", "NumStates", "Frequency", nodes))
theta <- 2^(sum(sol$Frequency * sol[[node]]) / sum(sol$Frequency))
par <- fread(file.path(dataDir, paste0(net, "_parameters.dat")), header = FALSE)
setnames(par, c("ParamID", "NumStates", prs))
pk <- par[, .(ParamID, GK = get(paste0("Prod_of_", node)) / get(paste0("Deg_of_", node)),
			  Relax = 1 / get(paste0("Deg_of_", node)))]

famColors <- c(Reference = "#8a8984", Frozen = "#eb6834", Hold = "#1baf7a", Stochastic = "#2a78d6")
condLevels <- c("C0", "C4", "C2", "C3 T=100", "C3 T=10", "C3 T=1", "C3 T=0.1", "C3 T=0.01", "C1")
condLabels <- c("C0 original", "C4 mean λ", "C2 frozen", "C3 T=100", "C3 T=10", "C3 T=1",
				"C3 T=0.1", "C3 T=0.01", "C1 stochastic")
label_conditions <- function(x, ctrl) {
	x <- x[Condition == "C0" | (Sigma %in% ctrl & DT == min(DT[Condition == "C3"], na.rm = TRUE))]
	x[, Cond := fifelse(Condition == "C3", paste0("C3 T=", as.character(THold)), Condition)]
	c0 <- x[Condition == "C0"]
	if (nrow(c0)) x <- rbind(x[Condition != "C0"], rbindlist(lapply(ctrl, function(s) copy(c0)[, Sigma := s])))
	x[, Cond := factor(Cond, levels = condLevels, labels = condLabels)]
	x[, Family := fcase(Cond == "C0 original", "Reference", Cond %in% c("C4 mean λ", "C2 frozen"), "Frozen",
						Cond == "C1 stochastic", "Stochastic", default = "Hold")]
	x[, SigmaLab := factor(paste0("σ = ", Sigma), levels = paste0("σ = ", ctrl))]
	x
}
theme_r <- function() theme_bw(base_size = 10) +
	theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
		  panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
		  strip.background = element_rect(fill = "grey96", colour = "grey80"),
		  strip.text.y.left = element_text(angle = 0, hjust = 1), strip.placement = "outside",
		  axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none",
		  plot.title = element_text(face = "bold"), plot.caption = element_text(colour = "grey35", hjust = 0))
save2 <- function(p, name, w, h) {
	ggsave(file.path(outDir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
	ggsave(file.path(outDir, paste0(name, ".png")), p, width = w, height = h, dpi = 200)
	cat("wrote", name, "\n")
}

rf <- list.files(collected, pattern = paste0("^", net, "_\\w+_replicates\\.csv\\.gz$"), full.names = TRUE)
lf <- list.files(collected, pattern = paste0("^", net, "_\\w+_lfc\\.csv\\.gz$"), full.names = TRUE)
for (f in rf) {
	r <- merge(fread(f), pk, by = "ParamID")
	mode <- r$NoiseMode[1]; tag <- paste0(net, "_node", node, "_", mode)
	ctrl <- sort(unique(r[Condition == "C3", Sigma]))
	r <- label_conditions(r, ctrl)
	r[, Mean := get(paste0("Mean_", node))]
	r[, W := 1 / .N, by = .(ParamID, Cond, Sigma)]                      # equal weight per parameter set
	r[, RelThr := log2(pmax(Mean, 1e-9) / theta)]
	# expression cannot exceed g/k; values above 1 are Float32 rounding of nodes at their maximum (≤ 1e-5)
	stopifnot(max(r$Mean / r$GK) < 1 + 1e-5)
	r[, FracMax := pmin(Mean / GK, 1)]
	n <- uniqueN(r$ParamID)
	for (v in c("RelThr", "FracMax")) {
		xlim <- if (v == "RelThr") c(-10, 7) else c(0, 1)
		ref <- if (v == "RelThr") geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) else
			geom_vline(xintercept = 0.505^2, colour = "grey30", linewidth = 0.4, linetype = "dashed")
		xlab <- if (v == "RelThr") paste0("log2(expression of ", node, " / threshold)   (0 = high/low threshold)") else
			paste0("Expression of ", node, " as a fraction of its maximum (g/k)")
		p <- ggplot(r, aes(.data[[v]], weight = W, fill = Family)) + ref +
			geom_density(aes(y = after_stat(scaled)), alpha = 0.75, colour = NA,
						 bounds = if (v == "FracMax") c(0, 1) else c(-Inf, Inf)) +
			facet_grid(Cond ~ SigmaLab, switch = "y") + scale_fill_manual(values = famColors) +
			coord_cartesian(xlim = xlim) +
			labs(x = xlab, y = NULL, title = paste0(net, " node ", node, " — ", mode, " noise: per-replicate expression"),
				 caption = paste0("Each replicate's time-averaged expression (t = 500–1000); ", n, " parameter sets, equal weight per set; ",
								  "first 100 replicates per condition (C2: first 100 frozen draws).")) +
			theme_r()
		save2(p, paste0(tag, "_replicates_", tolower(v)), 10, 9)
	}
}
for (f in lf) {
	l <- merge(fread(f), pk, by = "ParamID")
	mode <- l$NoiseMode[1]; tag <- paste0(net, "_node", node, "_", mode)
	ctrl <- sort(unique(l[Condition == "C3", Sigma]))
	l <- label_conditions(l, ctrl)
	cols <- grep(paste0("^", node, "_Frac_lag"), names(l), value = TRUE)
	m <- melt(l, id.vars = c("ParamID", "Cond", "Family", "SigmaLab", "Relax"), measure.vars = cols,
			  variable.name = "Col", value.name = "U")
	m[, Frac := U / 65535]
	m[, Lag := as.numeric(sub(".*_lag([0-9.]+)_.*", "\\1", Col))]
	m[, Fold := sub(".*_([0-9]+fold)$", "\\1", Col)]
	s <- m[, .(Frac = mean(Frac)), by = .(ParamID, Cond, Family, SigmaLab, Lag, Fold)]       # per set
	s <- s[, .(Mean = mean(Frac), Lo = quantile(Frac, 0.25), Hi = quantile(Frac, 0.75)), by = .(Cond, Family, SigmaLab, Lag, Fold)]
	p <- ggplot(s[Fold == "2fold"], aes(Lag, Mean, colour = Family, group = Cond)) +
		geom_ribbon(aes(ymin = Lo, ymax = Hi, fill = Family), alpha = 0.15, colour = NA) +
		geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
		facet_grid(Cond ~ SigmaLab, switch = "y") +
		scale_x_log10(breaks = sort(unique(s$Lag))) + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
		scale_colour_manual(values = famColors) + scale_fill_manual(values = famColors) +
		labs(x = "Window τ (time units, log scale)", y = NULL,
			 title = paste0(net, " node ", node, " — ", mode, " noise: how often expression shifts ≥ 2-fold within τ"),
			 caption = paste0("Fraction of time points with |log2 x(t+τ) − log2 x(t)| ≥ 1; mean over parameter sets (band: interquartile range). ",
							  "Node ", node, " relaxation time 1/k: median ", round(median(pk[ParamID %in% l$ParamID, Relax]), 1), " time units.")) +
		theme_r() + theme(axis.text.y = element_text(size = 7), axis.ticks.y = element_line())
	save2(p, paste0(tag, "_lfc"), 10, 9)
}
