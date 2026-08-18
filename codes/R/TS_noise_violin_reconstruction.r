# TS_noise_violin_reconstruction.r
#
# PURPOSE: Reconstruct TS's per-node expression distribution from
# all_parameters_stats.csv's per-(ParamID, Node) summary moments (Mean,
# StDev, MeanSkewness), at 5 fixed noise levels (0, 0.001, 0.01, 0.1, 1),
# and render as violins (one figure per noise mode: Additive, Fluctuating,
# Multiplicative). A horizontal threshold line per node (RACIPE's mean
# deterministic expression, same definition as figure_common.r's
# add_reachability()) is drawn in that node's own fill color.
#
# METHOD (same as the earlier point-density reconstruction, generalized to
# multiple noise levels + rendered as violins instead of overlaid curves):
#   - NoiseLevel > 0: per (ParamID, Node), fit a skew-normal matching
#     (Mean, StDev, Skewness) via a vectorized method-of-moments lookup
#     (skewness -> delta, precomputed on a fine grid once, then linear
#     interpolation -- equivalent to sn::cp2dp() but ~1000x faster at this
#     scale: millions of rows across 3 modes x 5 levels x 2 nodes). Mixture
#     density = equal-weight average of every parameter's fitted density,
#     evaluated on a shared per-(mode, level, node) log-spaced grid.
#   - NoiseLevel == 0: ~98-99% of rows have StDev == 0 (no noise term, so
#     each parameter sits at its own fixed point with no within-run
#     spread) -- a skew-normal fit needs SD > 0, so this level instead
#     uses a standard Gaussian KDE (small bandwidth) directly on the Mean
#     values, i.e. the population distribution of deterministic fixed
#     points rather than a stochastic-spread reconstruction.
#   - Skewness is clipped to +-0.99 (skew-normal's valid range is
#     +-0.9953) before fitting -- see TS_stationary_dist script/discussion
#     for why raw skewness estimates go far outside that range for a
#     handful of near-degenerate (small-StDev) parameters.
#
# Violins are drawn manually (geom_polygon), not geom_violin(), since the
# input here is a precomputed mixture density curve, not raw draws.
# ============================================================

library(funcsKishore)
suppressMessages(library(tidyverse))

resultsFolder <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final"
dataFolder    <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/data"
finalDir      <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/figures/final"

net          <- "TS"
dt           <- 0.01
noiseLevels  <- c(0, 0.001, 0.01, 0.1, 1.0)
skewClip     <- 0.99
gridN        <- 500
halfwidth    <- 0.22
dodgeWidth   <- 0.45

# ── Vectorized skew-normal moment-matching (skewness -> delta via a
#    precomputed lookup, replacing sn::cp2dp()'s per-call root-finding) ──
delta_grid  <- seq(-0.9999, 0.9999, length.out = 20001)
gamma1_of_delta <- function(delta) ((4 - pi) / 2) * (delta * sqrt(2/pi))^3 / (1 - 2*delta^2/pi)^1.5
gamma_grid  <- gamma1_of_delta(delta_grid)
skew_to_delta <- function(gamma1) approx(gamma_grid, delta_grid, xout = gamma1, rule = 2)$y

dsn_vec <- function(x, xi, omega, alpha) {
    z <- (x - xi) / omega
    (2 / omega) * dnorm(z) * pnorm(alpha * z)
}

# Node thresholds (RACIPE mean deterministic expression per node) -- same
# for every noise mode/level, drawn as one horizontal line per node.
solutions <- read_solutions(file.path(dataFolder, paste0(net, "_solution.dat")))
nodes     <- c("A", "B")
# TS_solution.dat's node columns are in log2 space (values go negative,
# e.g. -6.4, which raw expression can't) -- convert back to linear units
# via a basin-frequency-weighted mean in log2 space, then 2^(.), matching
# the Julia pipeline's own get_mean_expression() (src/RACIPEdata.jl) and
# figure_common.r's add_reachability() (same fix applied there).
solWeights <- solutions$basin / 100
thresholds <- setNames(sapply(nodes, function(nd) 2^(sum(solutions[[nd]] * solWeights) / sum(solWeights))), nodes)

# ── Per (Node, NoiseLevel) mixture-density curve ──
build_density <- function(d, nd, nl) {
    dsub <- d %>% filter(Node == nd, NoiseLevel == nl)

    if (nl == 0) {
        # Deterministic: reconstruct the population distribution of fixed
        # points directly from the Mean values (Silverman bandwidth KDE).
        # Fit in LOG space -- these values are plotted on a log axis, and a
        # bandwidth chosen on the raw linear scale is wildly too coarse
        # once log-transformed (it over-smooths the small-value end into a
        # near-uniform block instead of tapering).
        logMean <- log(pmax(dsub$Mean, 1e-3))
        dens <- density(logMean, n = gridN)
        return(tibble(y = exp(dens$x), density = dens$y))
    }

    dsub <- dsub %>% filter(StDev > 0) %>%
        mutate(SkewClipped = pmax(pmin(MeanSkewness, skewClip), -skewClip))
    delta <- skew_to_delta(dsub$SkewClipped)
    alpha <- delta / sqrt(1 - delta^2)
    omega <- dsub$StDev / sqrt(1 - 2*delta^2/pi)
    xi    <- dsub$Mean - omega * delta * sqrt(2/pi)

    lo <- max(min(dsub$Mean - 5*dsub$StDev), 1e-3)
    hi <- max(dsub$Mean + 5*dsub$StDev)
    ygrid <- exp(seq(log(lo), log(hi), length.out = gridN))
    dens  <- vapply(ygrid, function(y) mean(dsn_vec(y, xi, omega, alpha)), numeric(1))
    tibble(y = ygrid, density = dens)
}

# ── Manual violin polygon: mirror the density curve left/right around a
#    node-specific dodge position offset from the noise level's x-position
#    (scale = "width": max density maps to a fixed half-width, matching
#    this codebase's other violins' scale="width" convention). Each node
#    gets its own fully mirrored violin at (centerX + dodgeOffset), the
#    same as ggplot's position_dodge() -- NOT a split/half-violin sharing
#    one center line. ──
make_violin_polygon <- function(densDf, centerX, dodgeOffset, nd, nl) {
    densDf <- densDf %>% filter(density > 0)
    w <- densDf$density / max(densDf$density) * halfwidth
    x0 <- centerX + dodgeOffset
    tibble(
        x = c(x0 - w, rev(x0 + w)),
        y = c(densDf$y, rev(densDf$y)),
        Node = nd, NoiseLevel = nl
    )
}

run_for_mode <- function(noiseType) {
    message("== ", noiseType, " ==")
    f <- file.path(resultsFolder, noiseType, net, "results", "all_parameters_stats.csv")
    d <- read_csv(f, show_col_types = FALSE) %>% filter(DT == dt, NoiseLevel %in% noiseLevels)

    xpos <- setNames(seq_along(noiseLevels), as.character(noiseLevels))

    violinData <- map_dfr(noiseLevels, function(nl) {
        map_dfr(nodes, function(nd) {
            dens <- build_density(d, nd, nl)
            dodgeOffset <- if (nd == nodes[1]) -dodgeWidth/2 else dodgeWidth/2
            make_violin_polygon(dens, xpos[[as.character(nl)]], dodgeOffset, nd, nl)
        })
    }) %>%
        mutate(NoiseLevel = factor(NoiseLevel, levels = noiseLevels),
               Node = factor(Node, levels = nodes))

    nodeColors <- setNames(scales::hue_pal()(length(nodes)), nodes)

    threshDf <- tibble(Node = factor(nodes, levels = nodes), Threshold = thresholds[nodes])

    p <- ggplot() +
        geom_polygon(data = violinData, aes(x = x, y = y, fill = Node, group = interaction(Node, NoiseLevel)),
                     color = NA, alpha = 0.7) +
        geom_hline(data = threshDf, aes(yintercept = Threshold, color = Node),
                   linetype = "dashed", linewidth = 0.8) +
        scale_x_continuous(breaks = seq_along(noiseLevels), labels = as.character(noiseLevels)) +
        scale_y_log10() +
        scale_fill_manual(values = nodeColors) +
        scale_color_manual(values = nodeColors, guide = "none") +
        theme_Publication() +
        labs(x = "Noise Level", y = "Expression (log scale)",
             title = paste0("TS, ", noiseType, " noise: reconstructed expression distribution"),
             subtitle = "Skew-normal mixture per parameter (sigma>0); KDE on fixed points at sigma=0. Dashed lines = node threshold.")

    outFile <- file.path(finalDir, paste0("TS_violin_dist_", noiseType, ".jpg"))
    ggsave(outFile, p, width = 9, height = 6, dpi = 300)
    message("saved: ", outFile)
}

for (nt in c("Additive", "Fluctuating", "Multiplicative")) run_for_mode(nt)
