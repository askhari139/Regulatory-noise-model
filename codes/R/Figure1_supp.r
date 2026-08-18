# figure1_supp.r
#
# PURPOSE: Figure 1 supplement.
#   Panel A: TS RACIPE solution heatmap (ID x node, log2-scaled expression)
#            -- the heatmap companion to Fig1 Panel B(ii)'s scatter, from
#            plots_script.r's original p1.
#   Panel B: simulated lambda(t) under the three noise rules, for an
#            activation-range edge (lambda in [1, 100]) -- same design as
#            Fig1 Panel D, just a different lo/hi (simulate_lambda()'s
#            lambda0/effective-sigma scaling handle the rest automatically).
#
# ============================================================

library(funcsKishore)
library(tidyverse)
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
dataDir <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/data"

# ════════════════════════════════════════════════════════════════════════════
# PANEL A: TS solution heatmap (ID x node, log2-scaled expression)
# ════════════════════════════════════════════════════════════════════════════
ts_sol <- read_solutions(file.path(dataDir, "TS_solution.dat")) %>% arrange(A, B)
a_mean <- mean(ts_sol$A); b_mean <- mean(ts_sol$B)
log_breaks    <- pretty(c(ts_sol$A, ts_sol$B), n = 6)
linear_breaks <- 2^log_breaks

pA <- ts_sol %>% mutate(ID = row_number()) %>%
	select(ID, A, B) %>%
	pivot_longer(names_to = "Node", values_to = "Expression", -ID) %>%
	ggplot(aes(x = ID, y = Node, fill = Expression)) +
		geom_tile() +
		scale_fill_gradient2(low = "red", high = "blue", mid = "white",
							  midpoint = (a_mean + b_mean) / 2,
							  breaks = log_breaks, labels = round(linear_breaks, 2)) +
		labs(x = "", y = "") +
		theme_Publication() +
		theme(legend.position = "top", legend.key.width = unit(1.8, "cm")) +
		scale_x_continuous(expand = c(0, 0)) +
		scale_y_discrete(expand = c(0, 0))
ggsave(file.path(figDir, "fig1S1_panelA_TSheatmap.jpg"), pA, width = 10, height = 4, bg = "white")

# ════════════════════════════════════════════════════════════════════════════
# PANEL B: simulated lambda(t), activation-range edge (lambda in [1, 100])
# ════════════════════════════════════════════════════════════════════════════
pB <- plot_lambda_trajectories(
	lo = 1, hi = 100, n_stats = 100, n_display = 5, display_every = 1, display_alpha = 0.7,
	useInvLambda = TRUE,
	outFile = file.path(figDir, "fig1S1_panelB_lambda_act.jpg")
)

# ════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ════════════════════════════════════════════════════════════════════════════
fig1S1 <- plot_grid(pA, pB, ncol = 1, labels = c("A", "B"), label_size = LS,
					 rel_heights = c(1, 0.9))
ggsave(file.path(finalDir, "Fig1S1.jpg"), fig1S1,
       width = if (exists("FIG_WIDTH")) FIG_WIDTH else 12,
       height = if (exists("FIG_HEIGHT")) FIG_HEIGHT else 8,
       bg = "white")

message("Figure 1 supplement complete. Individual panels and combined figure saved to: ", figDir)
