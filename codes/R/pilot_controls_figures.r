# Stochastic vs deterministic comparisons for all control conditions (C0-C4), pilot run.
# Input: state-class summary from codes/regnoise/scripts/state_class_summary.jl
#   Rscript pilot_controls_figures.r <pilot_state_classes.csv> <output_dir>
# Additive noise is shown for the all-high state class (as in Fig. 4), Multiplicative noise
# for the single-high class (as in Fig. 5). Means are unweighted across parameter sets.

suppressPackageStartupMessages({
	library(data.table); library(ggplot2); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
d <- fread(args[1])
outDir <- args[2]
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# --- encodings --------------------------------------------------------------------
# Condition families (identity is also on the x axis / direct labels):
famColors <- c(Reference = "#8a8984", Frozen = "#eb6834", Hold = "#1baf7a", Stochastic = "#2a78d6")
classColors <- c(monostable = "#2a78d6", bistable = "#eb6834", tristable = "#1baf7a")
classShapes <- c(monostable = 16, bistable = 17, tristable = 15)
targetClass <- c(Additive = "all-high", Multiplicative = "single-high", LogAdditive = "all-high", Stationary = "all-high")
ctrlSigmas <- c(0.01, 0.1, 1)
poolDT <- 0.01

theme_pilot <- function(base = 11) {
	theme_bw(base_size = base) +
		theme(panel.grid.minor = element_blank(),
			  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
			  strip.background = element_rect(fill = "grey96", colour = "grey80"),
			  legend.position = "bottom",
			  plot.title = element_text(face = "bold", size = base + 1),
			  plot.caption = element_text(colour = "grey35", hjust = 0))
}

save_fig <- function(p, name, w, h) {
	# cairo_pdf: the default pdf device cannot encode λ, σ, −, — and prints dots instead
	ggsave(file.path(outDir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
	ggsave(file.path(outDir, paste0(name, ".png")), p, width = w, height = h, dpi = 200)
	cat("wrote", name, "\n")
}

d[, StabilityClass := factor(StabilityClass, levels = c("monostable", "bistable", "tristable"))]

# One row per parameter set and condition cell for the target state class.
target <- function(mode) d[NoiseMode == mode & StateClass == targetClass[[mode]]]

# Wide table of the conditions at the control sigmas (pool dt), one row per ParamID x sigma.
wide_controls <- function(x) {
	c0 <- x[Condition == "C0", .(Network, ParamID, StabilityClass, C0 = Occupancy)]
	ctl <- x[Condition != "C0" & Sigma %in% ctrlSigmas & DT == poolDT]
	ctl[, key := fifelse(Condition == "C3", paste0("C3_", THold), Condition)]
	w <- dcast(ctl, Network + ParamID + StabilityClass + Sigma ~ key, value.var = "Occupancy")
	merge(w, c0, by = c("Network", "ParamID", "StabilityClass"))
}

condLevels <- c("C0", "C4", "C2", "C3_100", "C3_10", "C3_1", "C3_0.1", "C3_0.01", "C1")
condLabels <- c("C0\norig.", "C4\nmean λ", "C2\nfrozen", "C3\nT=100", "C3\nT=10",
				"C3\nT=1", "C3\nT=0.1", "C3\nT=0.01", "C1\nstoch.")
condFamily <- c(C0 = "Reference", C4 = "Frozen", C2 = "Frozen", C1 = "Stochastic")

# --- Figure 1: condition comparison (Fig. 4D analogue) ------------------------------
fig_conditions <- function(mode, net) {
	w <- wide_controls(target(mode))[Network == net]
	long <- melt(w, id.vars = c("Network", "ParamID", "StabilityClass", "Sigma"),
				 measure.vars = condLevels, variable.name = "Cond", value.name = "Occ")
	s <- long[, .(Mean = mean(Occ), SD = sd(Occ), N = .N), by = .(StabilityClass, Sigma, Cond)]
	s[, Cond := factor(Cond, levels = condLevels)]
	s[, Family := fifelse(startsWith(as.character(Cond), "C3"), "Hold", condFamily[as.character(Cond)])]
	s[, SigmaLab := factor(paste0("σ = ", Sigma), levels = paste0("σ = ", ctrlSigmas))]
	ggplot(s, aes(Cond, Mean, colour = Family)) +
		geom_line(data = s[Family == "Hold"], aes(group = 1), linewidth = 0.5) +
		geom_linerange(aes(ymin = pmax(0, Mean - SD), ymax = pmin(1, Mean + SD)), linewidth = 0.5) +
		geom_point(size = 2.4) +
		facet_grid(StabilityClass ~ SigmaLab) +
		# limits pins the order: the first layer only holds C3, and ggplot would otherwise
		# train the axis on those levels first and append the rest after them
		scale_x_discrete(limits = condLevels, labels = condLabels) +
		scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
		scale_colour_manual(values = famColors, name = NULL,
							labels = c(Reference = "Original λ", Frozen = "Frozen λ (deterministic)",
									   Hold = "Hold-resample (C3)", Stochastic = "Stochastic (C1)")) +
		labs(x = NULL, y = paste0("Occupancy of ", targetClass[[mode]], " states"),
			 title = paste0(net, " — ", mode, " noise: occupancy under each condition"),
			 caption = paste0("Mean ± SD over parameter sets (", s$N[1], " per class). C2–C4 use the λ values of the ",
							  "stochastic run at the same σ (Δt = 0.01).")) +
		theme_pilot() + theme(axis.text.x = element_text(size = 7.5))
}

# --- Figure 2: effect distributions (Fig. 4B / 4C(i) analogue) ----------------------
fig_effects <- function(mode) {
	x <- target(mode)[DT == poolDT | Condition == "C0"]
	c0 <- x[Condition == "C0", .(Network, ParamID, C0 = Occupancy)]
	w <- dcast(x[Condition %in% c("C1", "C2")], Network + ParamID + Sigma ~ Condition, value.var = "Occupancy")
	w <- merge(w, c0, by = c("Network", "ParamID"))
	e <- rbind(w[, .(Network, Sigma, Effect = "Resampling effect (C2 − C0)", Value = C2 - C0)],
			   w[, .(Network, Sigma, Effect = "Stochastic effect (C1 − C2)", Value = C1 - C2)])
	e[, Effect := factor(Effect, levels = c("Resampling effect (C2 − C0)", "Stochastic effect (C1 − C2)"))]
	ggplot(e, aes(factor(Sigma), Value)) +
		geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
		geom_violin(fill = "grey85", colour = "grey45", linewidth = 0.3, scale = "width") +
		stat_summary(fun = mean, geom = "point", shape = 3, size = 2.5, stroke = 0.9, colour = "black") +
		facet_grid(Network ~ Effect) +
		scale_y_continuous(limits = c(-1, 1)) +
		labs(x = "σ", y = paste0("Difference in ", targetClass[[mode]], " occupancy"),
			 title = paste0(mode, " noise: resampling and stochastic effects across σ (Δt = 0.01)"),
			 caption = "Violins: distribution over parameter sets; + marks the mean. C2 uses 1,000 frozen draws per parameter set.") +
		theme_pilot()
}

fig_decomposition <- function(mode) {
	w <- wide_controls(target(mode))
	e <- rbind(w[, .(Network, Sigma, Step = "Mean\nC4 − C0", Value = C4 - C0)],
			   w[, .(Network, Sigma, Step = "Spread\nC2 − C4", Value = C2 - C4)],
			   w[, .(Network, Sigma, Step = "Persistence\nC3(0.01) − C2", Value = `C3_0.01` - C2)],
			   w[, .(Network, Sigma, Step = "Memory\nC1 − C3(0.01)", Value = C1 - `C3_0.01`)])
	e[, Step := factor(Step, levels = unique(Step))]
	e[, SigmaLab := factor(paste0("σ = ", Sigma), levels = paste0("σ = ", ctrlSigmas))]
	ggplot(e, aes(Step, Value)) +
		geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
		geom_violin(fill = "grey85", colour = "grey45", linewidth = 0.3, scale = "width") +
		stat_summary(fun = mean, geom = "point", shape = 3, size = 2.5, stroke = 0.9, colour = "black") +
		facet_grid(Network ~ SigmaLab) +
		scale_y_continuous(limits = c(-1, 1)) +
		labs(x = NULL, y = paste0("Change in ", targetClass[[mode]], " occupancy"),
			 title = paste0(mode, " noise: step-by-step decomposition of the noise effect"),
			 caption = paste("Each step changes one aspect of λ: its mean, its spread, its persistence time, and its temporal correlation.",
							 "C3(0.01) = hold-resample with T_hold = 0.01.")) +
		theme_pilot() + theme(axis.text.x = element_text(size = 7.5))
}

# --- Figure 3: per-parameter scatter (Fig. 4C(ii) analogue) --------------------------
fig_scatter <- function(mode, sigma = 0.1) {
	w <- wide_controls(target(mode))[Sigma == sigma]
	s <- rbind(w[, .(Network, StabilityClass, X = C2, Y = C1, Panel = "Stochastic (C1)")],
			   w[, .(Network, StabilityClass, X = C2, Y = `C3_0.01`, Panel = "Fast hold-resample (C3, T = 0.01)")])
	s[, Panel := factor(Panel, levels = c("Stochastic (C1)", "Fast hold-resample (C3, T = 0.01)"))]
	ggplot(s, aes(X, Y, colour = StabilityClass, shape = StabilityClass)) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
		geom_point(size = 1.6, alpha = 0.7) +
		facet_grid(Network ~ Panel) +
		coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
		scale_colour_manual(values = classColors, name = NULL) +
		scale_shape_manual(values = classShapes, name = NULL) +
		labs(x = paste0("Frozen-λ occupancy, C2 (", targetClass[[mode]], ")"),
			 y = paste0("Occupancy (", targetClass[[mode]], ")"),
			 title = paste0(mode, " noise: frozen vs time-varying λ per parameter set (σ = ", sigma, ", Δt = 0.01)"),
			 caption = "Each point is one parameter set. Dashed line: y = x.") +
		theme_pilot()
}

# --- Figure 4: hold-time ladder (new) ------------------------------------------------
fig_ladder <- function(mode) {
	w <- wide_controls(target(mode))
	holds <- c(0.01, 0.1, 1, 10, 100)
	lad <- rbindlist(lapply(holds, function(h) w[, .(Network, Sigma, THold = h, Occ = get(paste0("C3_", h)))]))
	lad <- lad[, .(Mean = mean(Occ), SE = sd(Occ) / sqrt(.N)), by = .(Network, Sigma, THold)]
	ref <- melt(w[, .(C4 = mean(C4), C2 = mean(C2), C1 = mean(C1)), by = .(Network, Sigma)],
				id.vars = c("Network", "Sigma"), variable.name = "Cond", value.name = "Mean")
	ref[, Label := c(C4 = "mean λ (C4)", C2 = "frozen (C2)", C1 = "stochastic (C1)")[as.character(Cond)]]
	ref[, Family := c(C4 = "Frozen", C2 = "Frozen", C1 = "Stochastic")[as.character(Cond)]]
	ref[, LT := c(C4 = "dashed", C2 = "dotted", C1 = "solid")[as.character(Cond)]]
	lab <- function(x) factor(paste0("σ = ", x), levels = paste0("σ = ", ctrlSigmas))
	lad[, SigmaLab := lab(Sigma)]; ref[, SigmaLab := lab(Sigma)]
	# direct labels: keep at least `gap` between label positions within a panel
	gap <- 0.07
	setorder(ref, Network, Sigma, Mean)
	ref[, LabY := {y <- Mean; for (i in seq_along(y)[-1]) y[i] <- max(y[i], y[i - 1] + gap)
				   y - max(0, max(y) - 0.98)}, by = .(Network, Sigma)]   # keep the stack inside the panel
	ggplot(lad, aes(THold, Mean)) +
		geom_hline(data = ref, aes(yintercept = Mean, colour = Family, linetype = LT), linewidth = 0.6) +
		geom_label(data = ref, aes(x = 160, y = LabY, label = Label, colour = Family),
				   hjust = 0, size = 2.7, linewidth = 0, fill = "white", label.padding = unit(0.08, "lines"),
				   show.legend = FALSE) +
		geom_ribbon(aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE), fill = famColors[["Hold"]], alpha = 0.2) +
		geom_line(colour = famColors[["Hold"]], linewidth = 0.8) +
		geom_point(colour = famColors[["Hold"]], size = 2.2) +
		facet_grid(Network ~ SigmaLab) +
		scale_x_log10(breaks = holds, labels = c("0.01", "0.1", "1", "10", "100"), expand = expansion(mult = c(0.05, 0.6))) +
		scale_y_continuous(limits = c(0, 1)) +
		scale_colour_manual(values = famColors, guide = "none") +
		scale_linetype_identity() +
		labs(x = "Hold time T_hold (λ redrawn from the pool every T_hold)",
			 y = paste0("Mean ", targetClass[[mode]], " occupancy"),
			 title = paste0(mode, " noise: occupancy vs the persistence time of λ"),
			 caption = "Green: hold-resample (C3), mean ± 95% CI over parameter sets. Lines: references at the same σ (Δt = 0.01).") +
		theme_pilot()
}

warnCounts <- list()
for (mode in intersect(c("Additive", "Multiplicative", "LogAdditive", "Stationary"), unique(d$NoiseMode))) {
	withCallingHandlers({
	tag <- c(Additive = "addi", Multiplicative = "mult", LogAdditive = "logadd", Stationary = "stat")[[mode]]
	for (net in sort(unique(d[NoiseMode == mode, Network]))) save_fig(fig_conditions(mode, net), paste0("fig1_conditions_", tag, "_", net), 11, 8)
	save_fig(fig_effects(mode), paste0("fig2a_effects_", tag), 10, 6.5)
	save_fig(fig_decomposition(mode), paste0("fig2b_decomposition_", tag), 11, 6.5)
	save_fig(fig_scatter(mode), paste0("fig3_scatter_", tag), 9, 8.5)
	save_fig(fig_ladder(mode), paste0("fig4_ladder_", tag), 11, 6.5)
	}, warning = function(w) { msg <- conditionMessage(w); if (!grepl("built under", msg)) warnCounts[[msg]] <<- (if (is.null(warnCounts[[msg]])) 0 else warnCounts[[msg]]) + 1; invokeRestart("muffleWarning") })
}
if (length(warnCounts)) { cat("warnings:\n"); for (m in names(warnCounts)) cat(sprintf("  [%d×] %s\n", warnCounts[[m]], m)) }
