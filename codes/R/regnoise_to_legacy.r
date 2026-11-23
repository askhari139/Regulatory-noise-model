# regnoise_to_legacy.r
#
# PURPOSE: Converts the collected tables of a regnoise production run
# (codes/regnoise, SPEC §11) into the per-network CSV layout the existing
# figure scripts read, so Figure*.r can be rerun on the new data unchanged:
#
#   <outDir>/<Mode>/<NET>/results/all_parameters_{results,transitions,stats}.csv
#   <outDir>/<Mode>/<NET>/results_det/all_parameters_{results,lambdas}.csv
#
# Mapping (new -> old):
#   C1 (stochastic)          -> results/       (NoiseLevel = Sigma)
#   C2 (frozen pool draws)   -> results_det/
#   C0 (deterministic)       -> NoiseLevel 0 rows of both, repeated for every DT
#   Occupancy Mean           -> MRT            (full window: Block 0, ThresholdMult 1)
#   State "HL"               -> "(1, 0)"
#   Stationary               -> Fluctuating    (old folder name)
#   Multiplicative           -> Multiplicative, plus a MultiplicativeInvLambda
#                               symlink for self-activating networks, so
#                               noise_type_for() in figure_common.r still resolves
#
# Zero-occupancy states are dropped: the old files only had visited states,
# and fill_mrt() builds Original from every State row at NoiseLevel 0.
#
# The new run stores per-edge quantiles of the λ pool, not draws. The old
# results_det lambdas file had 10 draws per (ParamID, NoiseLevel, DT); here
# those 10 values are read off the pool's quantile function at
# p = 0.05, 0.15, ..., 0.95 (linear between the stored quantiles).
#
# C3, C4, dwell, convergence, replicates and lfc have no legacy counterpart
# and are read directly from collected/ by the new analyses.
#
# Usage: Rscript regnoise_to_legacy.r <runDir> <outDir> [NET ...]
# ============================================================

suppressPackageStartupMessages(library(data.table))

args   <- commandArgs(trailingOnly = TRUE)
runDir <- if (length(args) >= 1) args[1] else
    "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/regnoise/2026-10-csbj-v2"
outDir <- if (length(args) >= 2) args[2] else file.path(runDir, "legacy")
nets   <- if (length(args) >= 3) args[-(1:2)] else
    c("TS", "TSSA", "TT", "TTSA", "TS4", "TT4", "TT4SA")
modes  <- c(Additive = "Additive", Multiplicative = "Multiplicative", Stationary = "Fluctuating",
            Geometric = "Multiplicative")   # Geometric (multiplicative rule of the revision) fills the old Multiplicative slot
sa_networks <- c("TSSA", "TTSA", "TT4SA")
collected <- file.path(runDir, "collected")

letters_to_state <- function(x) {
    vapply(strsplit(x, ""), function(ch) paste0("(", paste(ifelse(ch == "H", 1, 0), collapse = ", "), ")"),
           character(1))
}

# Reads one collected table, optionally keeping only rows that pass an awk
# condition (cheaper than loading the full block/threshold breakdown).
read_tbl <- function(net, mode, table, awk_cond = NULL) {
    f <- file.path(collected, paste0(net, "_", mode, "_", table, ".csv.gz"))
    if (!file.exists(f)) stop("Missing collected table: ", f)
    if (is.null(awk_cond)) return(fread(cmd = paste("gzip -dc", shQuote(f))))
    fread(cmd = sprintf("gzip -dc %s | awk -F, 'NR==1 || (%s)'", shQuote(f), awk_cond))
}

# C0 has DT = 0; the old files repeat the noise-free rows under every DT.
expand_c0 <- function(d, dts) {
    rbindlist(lapply(dts, function(x) copy(d)[, DT := x]))
}

convert_one <- function(net, mode) {
    oldMode <- modes[[mode]]
    message("== ", net, " / ", mode, " -> ", oldMode)

    # Occupancy: column 7 = Condition, 13 = ThresholdMult, 14 = Block
    occ <- read_tbl(net, mode, "occupancy",
                    '$7 ~ /^C[012]$/ && $13 == 1.0 && $14 == 0 && $15 > 0')
    cells <- read_tbl(net, mode, "cells", '$7 ~ /^C[012]$/')
    dts <- sort(unique(occ[Condition != "C0", DT]))

    occ[, State := letters_to_state(State)]
    occ <- merge(occ, cells[, .(Condition, Sigma, DT, ParamID, MeanSwitches, SDSwitches)],
                 by = c("Condition", "Sigma", "DT", "ParamID"), all.x = TRUE)

    to_old <- function(d) {
        d[, .(ParamID, ParamType = StabilityClass, NoiseLevel = Sigma, DT, State,
              MRT = Mean, RelativeMRT = NA_real_,
              MeanSwitches = fcoalesce(MeanSwitches, 0), StdSwitches = fcoalesce(SDSwitches, 0))]
    }
    c0 <- expand_c0(occ[Condition == "C0"], dts)
    stoch <- to_old(rbind(c0, occ[Condition == "C1"]))
    det   <- to_old(rbind(c0, occ[Condition == "C2"]))[, RelativeMRT := NULL]
    setorder(stoch, ParamID, DT, NoiseLevel, State)
    setorder(det, ParamID, DT, NoiseLevel, State)

    # Transitions (C1): old Count was the per-trajectory mean
    tr <- read_tbl(net, mode, "transitions", '$7 == "C1"')
    tr <- tr[, .(ParamID, ParamType = StabilityClass, NoiseLevel = Sigma, DT,
                 Transition = paste0("(", letters_to_state(From), ", ", letters_to_state(To), ")"),
                 Count = Count / NReplicates)]

    # Node expression (C0 + C1), replicate-averaged
    ex <- read_tbl(net, mode, "expression", '$7 ~ /^C[01]$/')
    ex <- rbind(expand_c0(ex[Condition == "C0"], dts), ex[Condition == "C1"])
    ex <- ex[, .(ParamID, ParamType = StabilityClass, NoiseLevel = Sigma, DT, Node,
                 Mean, StDev = MeanSD, Count = NReplicates)]

    # λ pool -> 10 quantile-function draws per (ParamID, Sigma, DT, Edge)
    lam <- read_tbl(net, mode, "lambda")
    qp  <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
    qc  <- c("Q01", "Q05", "Q25", "Q50", "Q75", "Q95", "Q99")
    pd  <- seq(0.05, 0.95, by = 0.1)
    qmat <- as.matrix(lam[, ..qc])
    draws <- t(apply(qmat, 1, function(q) approx(qp, q, xout = pd, ties = "ordered")$y))
    lam_long <- data.table(lam[, .(ParamID, ParamType = StabilityClass, NoiseLevel = Sigma, DT,
                                   Col = paste0("lambda_", ifelse(EdgeType == "inh", "Inh", "Act"),
                                                "_of_", Edge))][rep(seq_len(nrow(lam)), each = length(pd))],
                           Draw = rep(seq_along(pd), nrow(lam)), Value = as.vector(t(draws)))
    lam_wide <- dcast(lam_long, ParamID + ParamType + NoiseLevel + DT + Draw ~ Col, value.var = "Value")
    lam_wide[, Draw := NULL]

    resDir <- file.path(outDir, oldMode, net, "results")
    detDir <- file.path(outDir, oldMode, net, "results_det")
    dir.create(resDir, recursive = TRUE, showWarnings = FALSE)
    dir.create(detDir, recursive = TRUE, showWarnings = FALSE)
    fwrite(stoch,    file.path(resDir, "all_parameters_results.csv"))
    fwrite(tr,       file.path(resDir, "all_parameters_transitions.csv"))
    fwrite(ex,       file.path(resDir, "all_parameters_stats.csv"))
    fwrite(det,      file.path(detDir, "all_parameters_results.csv"))
    fwrite(lam_wide, file.path(detDir, "all_parameters_lambdas.csv"))

    if (oldMode == "Multiplicative" && net %in% sa_networks) {
        link <- file.path(outDir, "MultiplicativeInvLambda", net)
        dir.create(dirname(link), recursive = TRUE, showWarnings = FALSE)
        if (!file.exists(link)) file.symlink(file.path("..", "Multiplicative", net), link)
    }
    message("   ", uniqueN(stoch$ParamID), " ParamIDs, ", nrow(stoch), " C0/C1 rows, ",
            nrow(det), " C0/C2 rows")
}

# Only the modes this run actually has (a Geometric run has no Additive/Stationary tables).
has_mode <- function(net, mode) file.exists(file.path(collected, paste0(net, "_", mode, "_occupancy.csv.gz")))
for (net in nets) for (mode in names(modes)) if (has_mode(net, mode)) convert_one(net, mode)
message("Legacy tables written to ", outDir)
