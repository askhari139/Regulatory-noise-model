# Hybrid ("all-medium") states: is the all-high state an enrichment of intermediate expression?
# Answers R3.7. Every expression value is taken relative to its own parameter set's maximum g/k,
# and every comparison is made on individual trajectories (never on set means, which would make
# a bistable set look intermediate).
#
#   Rscript hybrid_state_figures.r <run_dir> <racipe_data_dir> <output_dir> [nets=TS,TT] [mode=Additive]
#
# run_dir holds collected/ and selection/ of one regnoise run, e.g. RACIPEdata/regnoise/2026-10-csbj-v2.
# Writes, per network:
#   <NET>_<mode>_A_distributions   per-trajectory expression / (g/k), all nodes, every condition
#   <NET>_<mode>_B_joint           A vs B per trajectory
#   <NET>_<mode>_B_threshold       A / (g/k) vs A / ensemble threshold (log) per trajectory
#   <NET>_<mode>_B_composition     what threshold-called states are made of (near zero / intermediate / near max)
#   <NET>_<mode>_C_paired          per-set change vs C0, mean with 95% bootstrap CI, by stability class
#   <NET>_<mode>_D_steadiness      within-trajectory CV and windowed 2-fold change of intermediate nodes
#   <NET>_<mode>_E_floor           expression vs the weak-repression floor (product of mean λ on in-edges)
#   hybrid_summary_<mode>.csv, hybrid_paired_<mode>.csv, hybrid_bands_<mode>.csv

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
runDir <- args[1]; dataDir <- args[2]; outDir <- args[3]
nets <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else c("TS", "TT")
mode <- if (length(args) >= 5) args[5] else "Additive"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261113)

MID   <- c(0.15, 0.85)                                   # intermediate band, fraction of g/k
BANDS <- list(c(0.10, 0.90), c(0.15, 0.85), c(0.20, 0.80))
SIG   <- c(0.01, 0.1, 1)
DT0   <- 0.01
LAGS  <- c(1, 3, 10, 30, 100)
NBOOT <- 1000

condLevels <- c("C0", "C4", "C2", "C3 T=100", "C3 T=10", "C3 T=1", "C3 T=0.1", "C3 T=0.01", "C1")
condLabels <- c("C0 original", "C4 mean λ", "C2 frozen", "C3 T=100", "C3 T=10", "C3 T=1",
                "C3 T=0.1", "C3 T=0.01", "C1 stochastic")
famColors  <- c(Reference = "#8a8984", Frozen = "#eb6834", Hold = "#1baf7a", Stochastic = "#2a78d6")
classColors <- c(all = "#222222", monostable = "#2a78d6", bistable = "#eb6834", tristable = "#1baf7a")
levelColors <- c("all intermediate" = "#1baf7a", "≥1 near max" = "#2a78d6",
                 "≥1 near zero" = "#8a8984", "near max and near zero" = "#eb6834")

theme_h <- function() theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey96", colour = "grey80"),
          strip.text.y.left = element_text(angle = 0, hjust = 1), strip.placement = "outside",
          plot.title = element_text(face = "bold"), plot.caption = element_text(colour = "grey35", hjust = 0))
wrapcap <- function(x, w) paste(strwrap(x, width = w), collapse = "\n")
save2 <- function(p, name, w, h) {
    if (!is.null(p$labels$caption)) p$labels$caption <- wrapcap(p$labels$caption, round(w * 17))
    ggsave(file.path(outDir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(outDir, paste0(name, ".png")), p, width = w, height = h, dpi = 200)
    cat("  wrote", name, "\n")
}

# Read only Δt ∈ {0, DT0} and σ ∈ {0} ∪ SIG; checks that columns 8/9 are Sigma/DT first.
awkFilter <- sprintf("NR==1 || (($9==0 || $9==%s) && ($8==0 || %s))", DT0, paste0("$8==", SIG, collapse = " || "))
read_filtered <- function(f, select) {
    hdr <- names(fread(cmd = sprintf("gzcat %s | head -2", shQuote(f))))
    stopifnot(hdr[8] == "Sigma", hdr[9] == "DT", all(select %in% hdr))
    fread(cmd = sprintf("gzcat %s | awk -F, '%s'", shQuote(f), awkFilter), select = select)
}
# Condition label per row; C0 (σ = 0) is copied into every σ column so it sits beside each noise level.
label_conditions <- function(x) {
    x[, Cond := fifelse(Condition == "C3", paste0("C3 T=", THold), Condition)]
    c0 <- x[Condition == "C0"]
    if (nrow(c0)) x <- rbind(x[Condition != "C0"], rbindlist(lapply(SIG, function(s) copy(c0)[, Sigma := s])))
    x <- x[Sigma %in% SIG & Cond %in% condLevels]
    x[, Cond := factor(Cond, levels = condLevels, labels = condLabels)]
    x[, Family := fcase(Cond == "C0 original", "Reference", Cond %in% c("C4 mean λ", "C2 frozen"), "Frozen",
                        Cond == "C1 stochastic", "Stochastic", default = "Hold")]
    x[, SigmaLab := factor(paste0("σ = ", Sigma), levels = paste0("σ = ", SIG))]
    x[]
}
# Bimodality coefficient (Pfister et al. 2013): > 5/9 suggests bi- or multimodality.
bimod_coef <- function(x) {
    n <- length(x); m <- mean(x); s <- sqrt(mean((x - m)^2))
    if (n < 4 || s == 0) return(NA_real_)
    g <- mean((x - m)^3) / s^3; k <- mean((x - m)^4) / s^4 - 3
    (g^2 + 1) / (k + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}
boot_ci <- function(x) {
    x <- x[is.finite(x)]; n <- length(x)
    if (n < 2) return(list(Lo = NA_real_, Hi = NA_real_))
    b <- colMeans(matrix(x[sample.int(n, n * NBOOT, replace = TRUE)], nrow = n))
    q <- quantile(b, c(0.025, 0.975), names = FALSE)
    list(Lo = q[1], Hi = q[2])
}

summaryAll <- list(); pairedAll <- list(); bandsAll <- list()
for (net in nets) {
    cat(net, mode, "\n")
    tag <- paste0(net, "_", mode)
    collected <- file.path(runDir, "collected")

    # ---- inputs: nodes, g/k and 1/k per set, thresholds, topology -----------------------------
    prs <- fread(file.path(dataDir, paste0(net, ".prs")))$Parameter
    nodes <- sub("Prod_of_", "", grep("^Prod_of_", prs, value = TRUE)); N <- length(nodes)
    par <- fread(file.path(dataDir, paste0(net, "_parameters.dat")), header = FALSE)
    setnames(par, c("ParamID", "NumStates", prs))
    gkw <- par[, c(list(ParamID = ParamID),
                   setNames(lapply(nodes, function(n) get(paste0("Prod_of_", n)) / get(paste0("Deg_of_", n))), paste0("GK_", nodes)),
                   setNames(lapply(nodes, function(n) 1 / get(paste0("Deg_of_", n))), paste0("Relax_", nodes)))]
    thr <- fread(file.path(runDir, "selection", paste0(net, "_thresholds.csv")))
    theta <- setNames(thr$Threshold, thr$Node)[nodes]

    # ---- per-trajectory expression -----------------------------------------------------------
    r <- read_filtered(file.path(collected, paste0(tag, "_replicates.csv.gz")),
                       c("ParamID", "StabilityClass", "Condition", "Sigma", "THold", "Replicate",
                         paste0(rep(c("Mean_", "SD_"), each = N), nodes)))
    r <- label_conditions(r)
    r <- merge(r, gkw, by = "ParamID")
    nSets <- uniqueN(r$ParamID)
    cat("  ", nrow(r), "trajectory rows,", nSets, "parameter sets\n")
    for (n in nodes) {
        fm <- r[[paste0("Mean_", n)]] / r[[paste0("GK_", n)]]
        # expression cannot exceed g/k; small excesses (≤ 1e-3) are solver overshoot at λ jumps
        # (reltol 1e-4; seen only in C3 at large σ) and are capped at 1. Anything larger is an input mismatch.
        stopifnot(max(fm) < 1 + 1e-3)
        if (any(fm > 1 + 1e-5)) cat("   ", n, ":", sum(fm > 1 + 1e-5), "trajectories above g/k (max", signif(max(fm), 6), "), capped at 1\n")
        set(r, j = paste0("F_", n), value = pmin(fm, 1))
        set(r, j = paste0("CV_", n), value = r[[paste0("SD_", n)]] / pmax(r[[paste0("Mean_", n)]], 1e-12))
    }
    lvl <- function(lo, hi) lapply(nodes, function(n) { f <- r[[paste0("F_", n)]]; fcase(f < lo, -1L, f > hi, 1L, default = 0L) })
    L <- lvl(MID[1], MID[2])
    r[, NMid := Reduce(`+`, lapply(L, function(v) v == 0L))]
    r[, NHi  := Reduce(`+`, lapply(L, function(v) v == 1L))]
    r[, NLo  := Reduce(`+`, lapply(L, function(v) v == -1L))]
    r[, ThrAllHigh := Reduce(`&`, lapply(nodes, function(n) r[[paste0("Mean_", n)]] > theta[[n]]))]
    r[, Composition := factor(fcase(NMid == N, "all intermediate", NLo == 0, "≥1 near max",
                                    NHi == 0, "≥1 near zero", default = "near max and near zero"),
                              levels = names(levelColors))]

    # ---- A: per-trajectory expression distributions (all nodes pooled) ----------------------
    fl <- melt(r, id.vars = c("ParamID", "Cond", "Family", "SigmaLab"), measure.vars = paste0("F_", nodes),
               variable.name = "Node", value.name = "F")
    ann <- fl[, .(MidFrac = mean(F >= MID[1] & F <= MID[2]), LoFrac = mean(F < MID[1]), HiFrac = mean(F > MID[2]),
               BC = bimod_coef(F)), by = .(Cond, Family, SigmaLab)]
    hist <- fl[, .(n = .N), by = .(Cond, Family, SigmaLab, Bin = pmin(floor(F * 50), 49) / 50 + 0.01)]
    hist[, Dens := n / max(n), by = .(Cond, SigmaLab)]
    pA <- ggplot(hist, aes(Bin, Dens, fill = Family)) +
        annotate("rect", xmin = MID[1], xmax = MID[2], ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.6) +
        geom_col(width = 0.02, alpha = 0.85) +
        geom_text(data = ann, aes(x = 0.5, y = 0.93, label = sprintf("near 0: %.0f%% · intermediate: %.0f%% · near max: %.0f%%", 100 * LoFrac, 100 * MidFrac, 100 * HiFrac)),
                  inherit.aes = FALSE, size = 2.6, colour = "grey15") +
        facet_grid(Cond ~ SigmaLab, switch = "y") + scale_fill_manual(values = famColors, guide = "none") +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) + scale_y_continuous(breaks = NULL, limits = c(0, 1.05)) +
        labs(x = "Expression as a fraction of the node's maximum g/k (one value per trajectory and node)", y = NULL,
             title = paste0(net, " — ", mode, " noise: expression relative to each parameter set's maximum"),
             caption = paste0("Time-averaged expression over t = 500–1000, all ", N, " nodes pooled; ", nSets,
                              " parameter sets × 100 trajectories (C2: first 100 frozen draws); Δt = ", DT0, ". ",
                              "Grey band: intermediate (", MID[1], "–", MID[2], " of g/k); near 0 and near max are the bins outside it.")) +
        theme_h() + theme(axis.ticks.y = element_blank())
    save2(pA, paste0(tag, "_A_distributions"), 10, 10)

    # ---- B: joint levels (A vs B) and composition of threshold-called states -----------------
    showB <- condLabels[c(1, 2, 3, 4, 6, 9)]
    jb <- r[Cond %in% showB, .(n = .N), by = .(Cond, SigmaLab,
                                               X = pmin(floor(F_A * 40), 39) / 40 + 0.0125,
                                               Y = pmin(floor(F_B * 40), 39) / 40 + 0.0125)]
    jb[, Frac := n / sum(n), by = .(Cond, SigmaLab)]
    jAnn <- r[Cond %in% showB, .(Both = mean(F_A >= MID[1] & F_A <= MID[2] & F_B >= MID[1] & F_B <= MID[2])), by = .(Cond, SigmaLab)]
    pB <- ggplot(jb, aes(X, Y, fill = Frac)) + geom_tile(width = 0.025, height = 0.025) +
        annotate("rect", xmin = MID[1], xmax = MID[2], ymin = MID[1], ymax = MID[2], fill = NA, colour = "grey20", linetype = "dashed", linewidth = 0.3) +
        geom_text(data = jAnn, aes(x = 0.98, y = 0.98, label = sprintf("%.0f%% in box", 100 * Both)),
                  inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.6) +
        facet_grid(Cond ~ SigmaLab, switch = "y") + coord_equal() +
        scale_fill_viridis_c(trans = "log10", name = "Fraction of\ntrajectories", option = "magma", direction = -1) +
        scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        labs(x = "A / (g/k)", y = "B / (g/k)", title = paste0(net, " — ", mode, " noise: joint expression of A and B per trajectory"),
             caption = paste0("Dashed box: both nodes intermediate (", MID[1], "–", MID[2], " of g/k). ", nSets, " parameter sets × 100 trajectories.")) +
        theme_h()
    save2(pB, paste0(tag, "_B_joint"), 8, 11)

    # A relative to its own maximum (x) vs A relative to the ensemble threshold (y, log): where do
    # "high" calls sit on the g/k scale? Axis range from the central 99.8% of values, rounded out to a decade.
    jt <- r[Cond %in% showB, .(Cond, SigmaLab, F = F_A, T = log10(pmax(Mean_A, 1e-12) / theta[["A"]]))]
    yl <- c(floor(quantile(jt$T, 0.001)), ceiling(quantile(jt$T, 0.999)))
    nyb <- 40; dy <- diff(yl) / nyb
    jtb <- jt[, .(n = .N), by = .(Cond, SigmaLab,
                                  X = pmin(floor(F * 40), 39) / 40 + 0.0125,
                                  Y = yl[1] + (pmin(pmax(floor((T - yl[1]) / dy), 0), nyb - 1) + 0.5) * dy)]
    jtb[, Frac := n / sum(n), by = .(Cond, SigmaLab)]
    tAnn <- jt[, .(Above = mean(T > 0), MidGivenAbove = if (any(T > 0)) mean(F[T > 0] >= MID[1] & F[T > 0] <= MID[2]) else NA_real_),
               by = .(Cond, SigmaLab)]
    pBt <- ggplot(jtb, aes(X, Y, fill = Frac)) +
        annotate("rect", xmin = MID[1], xmax = MID[2], ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5) +
        geom_tile(width = 0.025, height = dy) +
        geom_hline(yintercept = 0, colour = "grey20", linetype = "dashed", linewidth = 0.3) +
        geom_text(data = tAnn, aes(x = 0.02, y = yl[2] - 0.2 * dy,
                                   label = sprintf("%.0f%% above threshold\n%.0f%% of those intermediate", 100 * Above, 100 * MidGivenAbove)),
                  inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.3, lineheight = 0.9) +
        facet_grid(Cond ~ SigmaLab, switch = "y") +
        scale_fill_viridis_c(trans = "log10", name = "Fraction of\ntrajectories", option = "magma", direction = -1) +
        scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(add = 0.02)) +
        scale_y_continuous(limits = yl, breaks = seq(yl[1], yl[2]), labels = function(b) parse(text = paste0("10^", b)),
                           expand = expansion(mult = 0.04)) +
        labs(x = "A / (g/k): fraction of the set's own maximum",
             y = "A / ensemble threshold (log scale)",
             title = paste0(net, " — ", mode, " noise: node A vs its own maximum and vs the threshold"),
             caption = paste0("One value per trajectory (time-averaged A, t = 500–1000). Dashed line: the ensemble threshold used to call high/low (",
                              signif(theta[["A"]], 3), "); above it A is called high. Grey band: intermediate (", MID[1], "–", MID[2], " of g/k). ",
                              "Values outside 10^", yl[1], "–10^", yl[2], " are piled at the edges. ", nSets, " parameter sets × 100 trajectories.")) +
        theme_h() + theme(panel.spacing.y = unit(0.6, "lines"))
    save2(pBt, paste0(tag, "_B_threshold"), 8.5, 12)
    rm(jt, jtb)

    comp <- rbind(r[, .(Set = "All trajectories", .N), by = .(Cond, SigmaLab, Composition)],
                  r[ThrAllHigh == TRUE, .(Set = "Trajectories called all-high (threshold)", .N), by = .(Cond, SigmaLab, Composition)])
    comp[, Frac := N / sum(N), by = .(Set, Cond, SigmaLab)]
    nAH <- r[, .(AH = mean(ThrAllHigh), n = sum(ThrAllHigh)), by = .(Cond, SigmaLab)]
    pBc <- ggplot(comp, aes(Cond, Frac, fill = Composition)) + geom_col(width = 0.8) +
        geom_text(data = nAH, aes(Cond, 1.04, label = sprintf("%.0f%%", 100 * AH)), inherit.aes = FALSE, size = 2.3,
                  colour = "grey30") +
        facet_grid(Set ~ SigmaLab) + scale_fill_manual(values = levelColors, name = NULL, drop = FALSE) +
        scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.08)) +
        labs(x = NULL, y = "Fraction of trajectories",
             title = paste0(net, " — ", mode, " noise: what the threshold-defined states are made of"),
             caption = paste0("Near zero < ", MID[1], ", intermediate ", MID[1], "–", MID[2], ", near max > ", MID[2], " of g/k, per node. ",
                              "Numbers above bars: fraction of trajectories whose time-averaged expression is above threshold for every node.")) +
        theme_h() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save2(pBc, paste0(tag, "_B_composition"), 10, 7)

    # ---- C: per-set paired change against C0 --------------------------------------------------
    ps <- fl[, .(MidFrac = mean(F >= MID[1] & F <= MID[2])), by = .(ParamID, Cond, SigmaLab)]
    ps <- merge(ps, r[, .(AllMidFrac = mean(NMid == N), ThrAllHighFrac = mean(ThrAllHigh)), by = .(ParamID, StabilityClass, Cond, SigmaLab)],
                by = c("ParamID", "Cond", "SigmaLab"))
    pl <- melt(ps, id.vars = c("ParamID", "StabilityClass", "Cond", "SigmaLab"),
               measure.vars = c("MidFrac", "AllMidFrac", "ThrAllHighFrac"), variable.name = "Metric")
    base <- pl[Cond == "C0 original", .(ParamID, SigmaLab, Metric, Base = value)]
    pl <- merge(pl[Cond != "C0 original"], base, by = c("ParamID", "SigmaLab", "Metric"))
    pl[, Diff := value - Base]
    pl <- rbind(pl, copy(pl)[, StabilityClass := "all"])
    paired <- pl[, c(list(Mean = mean(Diff), Median = median(Diff), FracUp = mean(Diff > 0), FracDown = mean(Diff < 0), NSets = .N),
                     boot_ci(Diff)), by = .(StabilityClass, Metric, Cond, SigmaLab)]
    paired[, Metric := factor(Metric, levels = c("MidFrac", "AllMidFrac", "ThrAllHighFrac"),
                              labels = c("Node-level: fraction\nintermediate", "Trajectory-level: every\nnode intermediate",
                                         "Trajectory-level: called\nall-high by threshold"))]
    paired[, StabilityClass := factor(StabilityClass, levels = names(classColors))]
    pC <- ggplot(paired, aes(Cond, Mean, colour = StabilityClass)) +
        geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
        geom_pointrange(aes(ymin = Lo, ymax = Hi), position = position_dodge(width = 0.6), size = 0.25) +
        facet_grid(Metric ~ SigmaLab, scales = "free_y", switch = "y") +
        scale_colour_manual(values = classColors, name = "Stability class") +
        labs(x = NULL, y = NULL, title = paste0(net, " — ", mode, " noise: paired change per parameter set, relative to C0"),
             caption = paste0("Per set: fraction of its 100 trajectories (× nodes for the top row) in the category, minus the same set's C0 value. ",
                              "Points: mean over sets; bars: 95% bootstrap CI (", NBOOT, " resamples).")) +
        theme_h() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save2(pC, paste0(tag, "_C_paired"), 10, 8)

    # ---- D: steadiness of intermediate nodes ---------------------------------------------------
    cvl <- melt(r, id.vars = c("ParamID", "Cond", "Family", "SigmaLab", "Replicate"),
                measure.vars = list(paste0("F_", nodes), paste0("CV_", nodes)), variable.name = "NodeIdx", value.name = c("F", "CV"))
    cvl[, Node := nodes[as.integer(NodeIdx)]][, NodeIdx := NULL]
    cvl <- cvl[F >= MID[1] & F <= MID[2]]
    qs <- function(x) as.list(setNames(quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9), names = FALSE), c("q10", "q25", "q50", "q75", "q90")))
    dCV <- cvl[, c(list(Measure = "Within-trajectory CV (SD / mean)", n = .N), qs(CV)), by = .(Cond, Family, SigmaLab)]

    lf <- read_filtered(file.path(collected, paste0(tag, "_lfc.csv.gz")),
                        c("ParamID", "Condition", "Sigma", "THold", "Replicate",
                          as.vector(outer(nodes, LAGS, function(n, l) paste0(n, "_Frac_lag", l, "_2fold")))))
    lf <- label_conditions(lf)
    lf <- merge(lf, gkw[, c("ParamID", paste0("Relax_", nodes)), with = FALSE], by = "ParamID")
    lfl <- rbindlist(lapply(nodes, function(n) {
        m <- as.matrix(lf[, paste0(n, "_Frac_lag", LAGS, "_2fold"), with = FALSE])
        k <- max.col(-abs(outer(log(lf[[paste0("Relax_", n)]]), log(LAGS), `-`)), ties.method = "first")
        lf[, .(ParamID, Cond, Family, SigmaLab, Replicate, Node = n, Lfc = m[cbind(seq_len(nrow(m)), k)] / 65535)]
    }))
    lfl <- merge(lfl, cvl[, .(ParamID, Cond, SigmaLab, Replicate, Node)], by = c("ParamID", "Cond", "SigmaLab", "Replicate", "Node"))
    dL <- lfl[, c(list(Measure = "Fraction of time with ≥ 2-fold change\nover one relaxation time (1/k)", n = .N), qs(Lfc)), by = .(Cond, Family, SigmaLab)]
    dD <- rbind(dCV, dL)
    pD <- ggplot(dD, aes(Cond, q50, colour = Family)) +
        geom_linerange(aes(ymin = q10, ymax = q90), linewidth = 0.4) +
        geom_crossbar(aes(ymin = q25, ymax = q75), width = 0.5, linewidth = 0.3, fill = "white") +
        facet_grid(Measure ~ SigmaLab, scales = "free_y", switch = "y") + scale_colour_manual(values = famColors, guide = "none") +
        labs(x = NULL, y = NULL, title = paste0(net, " — ", mode, " noise: are intermediate nodes steady or flickering?"),
             caption = paste0("Only node-trajectories at an intermediate level (", MID[1], "–", MID[2], " of g/k). Box: interquartile range and median; line: 10th–90th percentile. ",
                              "Windowed change is stored for C1 and C3 only (deterministic conditions are at a fixed point); lag nearest to each node's 1/k.")) +
        theme_h() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save2(pD, paste0(tag, "_D_steadiness"), 10, 6.5)
    rm(cvl, lf, lfl); gc(verbose = FALSE)

    # ---- E: the weak-repression floor ----------------------------------------------------------
    # min of H⁻ + λH⁺ over the regulator is λ (inhibition); min of (H⁻ + λH⁺)/λ is 1/λ (activation).
    lam <- read_filtered(file.path(collected, paste0(tag, "_lambda.csv.gz")), c("ParamID", "Condition", "Sigma", "Edge", "EdgeType", "Mean"))
    lam <- lam[Condition == "C1"]
    lam[, Target := sub("^.*To", "", Edge)]
    fl0 <- lam[, .(Floor = prod(fifelse(EdgeType == "inh", Mean, 1 / Mean))), by = .(ParamID, Sigma, Target)]
    fl0[, SigmaLab := factor(paste0("σ = ", Sigma), levels = paste0("σ = ", SIG))]
    showE <- condLabels[c(2, 8, 9)]
    ev <- melt(r[Cond %in% showE], id.vars = c("ParamID", "Cond", "SigmaLab", "Replicate"), measure.vars = paste0("F_", nodes),
               variable.name = "Target", value.name = "F")
    ev[, Target := sub("^F_", "", Target)]
    ev <- merge(ev, fl0[, .(ParamID, SigmaLab, Target, Floor)], by = c("ParamID", "SigmaLab", "Target"))
    ev[, Ratio := F / Floor]
    eAnn <- ev[, .(Near = mean(Ratio < 1.25), Below = mean(Ratio < 0.95)), by = .(Cond, SigmaLab)]
    eh <- ev[, .(n = .N), by = .(Cond, SigmaLab, Bin = pmin(pmax(floor(log2(pmax(Ratio, 1e-6)) * 20) / 20, -2), 2) + 0.025)]
    eh[, Dens := n / max(n), by = .(Cond, SigmaLab)]
    eFl <- fl0[, .(Med = median(Floor)), by = SigmaLab]
    pE <- ggplot(eh, aes(Bin, Dens)) + geom_col(width = 0.05, fill = "#2a78d6", alpha = 0.85) +
        geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.4) +
        geom_text(data = eAnn, aes(x = 1.95, y = 0.95, label = sprintf("%.0f%% within 1.25× of floor\n%.1f%% below floor", 100 * Near, 100 * Below)),
                  inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.6) +
        facet_grid(Cond ~ SigmaLab, switch = "y") + scale_y_continuous(breaks = NULL) +
        scale_x_continuous(limits = c(-2.05, 2.05), breaks = -2:2, labels = c("¼×", "½×", "1×", "2×", "4×")) +
        labs(x = "Expression ÷ floor (log scale; values beyond ¼× and 4× are piled at the edges)", y = NULL,
             title = paste0(net, " — ", mode, " noise: expression sits on the floor set by weak repression"),
             caption = paste0("One value per trajectory and node. Floor: the lowest expression (fraction of g/k) a node can reach when all its regulators are fully on, ",
                              "= product over in-edges of mean λ (inhibition) or 1/mean λ (activation), mean λ from the set's own C1 pool. ",
                              "Median floor: ", paste(sprintf("%s %.2f", eFl$SigmaLab, eFl$Med), collapse = ", "), ".")) +
        theme_h() + theme(axis.ticks.y = element_blank())
    save2(pE, paste0(tag, "_E_floor"), 9, 5.5)

    # ---- tables --------------------------------------------------------------------------------
    s <- r[, .(Trajectories = .N, AllIntermediate = mean(NMid == N), ThrAllHigh = mean(ThrAllHigh),
               AllIntermediateGivenThrAllHigh = if (any(ThrAllHigh)) mean(NMid[ThrAllHigh] == N) else NA_real_,
               AnyNearMaxGivenThrAllHigh = if (any(ThrAllHigh)) mean(NHi[ThrAllHigh] > 0) else NA_real_),
           by = .(Cond, SigmaLab)]
    s <- merge(s, ann[, .(Cond, SigmaLab, NodeIntermediate = MidFrac, BimodalityCoef = BC)], by = c("Cond", "SigmaLab"))
    s <- merge(s, eAnn[, .(Cond, SigmaLab, WithinFloor1.25 = Near)], by = c("Cond", "SigmaLab"), all.x = TRUE)
    s[, `:=`(Network = net, NoiseMode = mode, NSets = nSets)]
    summaryAll[[net]] <- s
    pairedAll[[net]] <- paired[, `:=`(Network = net, NoiseMode = mode)]
    bandsAll[[net]] <- rbindlist(lapply(BANDS, function(b) {
        Lb <- lvl(b[1], b[2]); nm <- Reduce(`+`, lapply(Lb, function(v) v == 0L))
        r[, .(Band = sprintf("%.2f-%.2f", b[1], b[2]), AllIntermediate = mean(nm[.I] == N)), by = .(Cond, SigmaLab)]
    }))[, `:=`(Network = net, NoiseMode = mode)]

    cat("  key numbers (Δt = ", DT0, "):\n", sep = "")
    print(s[Cond %in% condLabels[c(1, 2, 3, 4, 8, 9)],
            .(Cond, SigmaLab, NodeInt = round(NodeIntermediate, 3), AllInt = round(AllIntermediate, 3),
              ThrAllHigh = round(ThrAllHigh, 3), AllIntGivenAH = round(AllIntermediateGivenThrAllHigh, 3),
              BC = round(BimodalityCoef, 2))][order(SigmaLab, Cond)], row.names = FALSE)
    rm(r, fl, ps, pl, ev); gc(verbose = FALSE)
}
fwrite(rbindlist(summaryAll, use.names = TRUE), file.path(outDir, paste0("hybrid_summary_", mode, ".csv")))
fwrite(rbindlist(pairedAll, use.names = TRUE), file.path(outDir, paste0("hybrid_paired_", mode, ".csv")))
fwrite(rbindlist(bandsAll, use.names = TRUE), file.path(outDir, paste0("hybrid_bands_", mode, ".csv")))
cat("done ->", outDir, "\n")
