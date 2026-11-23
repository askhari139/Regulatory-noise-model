# compare_legacy_runs.r
#
# PURPOSE: Old pipeline (RACIPEdata/final) vs. regnoise run converted by
# regnoise_to_legacy.r, on the ParamIDs both runs simulated. For every
# (Network, NoiseType, run kind, sigma, DT, ParamID) the two state
# distributions are compared by total variation distance,
# TVD = 0.5 * sum_s |p_old(s) - p_new(s)|, and by the all-high and
# single-high occupancies. NoiseLevel 0 (no noise; the two runs differ only in
# their random ICs) is the floor that IC sampling alone produces.
#
# Old Multiplicative data for self-activating networks is the
# MultiplicativeInvLambda folder (what the submitted Fig. 5 used); both old
# folders are compared against the new target-based Multiplicative run.
#
# Usage: Rscript compare_legacy_runs.r <oldFolder> <newFolder> <outDir>
# ============================================================

suppressPackageStartupMessages(library(data.table))

args      <- commandArgs(trailingOnly = TRUE)
oldFolder <- if (length(args) >= 1) args[1] else "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final"
newFolder <- if (length(args) >= 2) args[2] else
    "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/regnoise/2026-10-csbj-v2/legacy"
outDir    <- if (length(args) >= 3) args[3] else
    "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/Regulatory-noise-model/CSBJ_review/production_match"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

nets <- c("TS", "TSSA", "TT", "TTSA", "TS4", "TT4", "TT4SA")
# (old folder, new folder) pairs
pairs <- list(c("Additive", "Additive"), c("Fluctuating", "Fluctuating"),
              c("Multiplicative", "Multiplicative"), c("MultiplicativeInvLambda", "Multiplicative"))
kinds <- c(stoch = "results", det = "results_det")

read_mrt <- function(folder, mode, net, kind) {
    f <- file.path(folder, mode, net, kinds[[kind]], "all_parameters_results.csv")
    if (!file.exists(f)) return(NULL)
    # old files have a blank line after every ParamID block
    d <- fread(f, select = c("ParamID", "NoiseLevel", "DT", "State", "MRT"), blank.lines.skip = TRUE)
    d <- d[!is.na(ParamID)]
    d[, .(MRT = sum(MRT)), by = .(ParamID, NoiseLevel, DT, State)]
}

n_high <- function(s) lengths(regmatches(s, gregexpr("1", s)))

compare_one <- function(net, oldMode, newMode, kind) {
    dOld <- read_mrt(oldFolder, oldMode, net, kind)
    dNew <- read_mrt(newFolder, newMode, net, kind)
    if (is.null(dOld) || is.null(dNew)) return(NULL)
    shared <- intersect(unique(dOld$ParamID), unique(dNew$ParamID))
    if (length(shared) == 0) return(NULL)
    dOld <- dOld[ParamID %in% shared]
    dNew <- dNew[ParamID %in% shared]
    # Normalise per cell: old det MRT can sum slightly off 1 after averaging iterations
    dOld[, MRT := MRT / sum(MRT), by = .(ParamID, NoiseLevel, DT)]
    dNew[, MRT := MRT / sum(MRT), by = .(ParamID, NoiseLevel, DT)]
    m <- merge(dOld, dNew, by = c("ParamID", "NoiseLevel", "DT", "State"),
               all = TRUE, suffixes = c("_old", "_new"))
    m[is.na(MRT_old), MRT_old := 0][is.na(MRT_new), MRT_new := 0]
    nn <- n_high(m$State[1]) + lengths(regmatches(m$State[1], gregexpr("0", m$State[1])))
    m[, k := n_high(State)]
    per <- m[, .(TVD        = 0.5 * sum(abs(MRT_old - MRT_new)),
                 AllHigh_old = sum(MRT_old[k == nn]),  AllHigh_new = sum(MRT_new[k == nn]),
                 Single_old  = sum(MRT_old[k == 1]),   Single_new  = sum(MRT_new[k == 1])),
             by = .(ParamID, NoiseLevel, DT)]
    per[, `:=`(Network = net, OldMode = oldMode, NewMode = newMode, Kind = kind)]
    per
}

per <- rbindlist(lapply(nets, function(net) rbindlist(lapply(pairs, function(p)
    rbindlist(lapply(names(kinds), function(k) compare_one(net, p[1], p[2], k))))))
)
fwrite(per, file.path(outDir, "paramwise_old_vs_new.csv.gz"))

summ <- per[, .(NShared        = uniqueN(ParamID),
                MedianTVD      = median(TVD),
                FracTVDlt0.1   = mean(TVD < 0.1),
                AllHigh_old    = mean(AllHigh_old), AllHigh_new = mean(AllHigh_new),
                AllHigh_cor    = suppressWarnings(cor(AllHigh_old, AllHigh_new)),
                Single_old     = mean(Single_old),  Single_new  = mean(Single_new)),
            by = .(Network, OldMode, NewMode, Kind, DT, NoiseLevel)]
setorder(summ, Network, OldMode, Kind, DT, NoiseLevel)
fwrite(summ, file.path(outDir, "summary_old_vs_new.csv"))
message("Wrote ", nrow(per), " paramwise rows and ", nrow(summ), " summary rows to ", outDir)
