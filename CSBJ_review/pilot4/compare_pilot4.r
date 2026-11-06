# Pilot 4: Geometric with (activation cap, inhibitory floor) = (100, 0.01) and (1000, 0.001), TT4 + TT4SA,
# 100 sets per class (a subset of production). Compared on the same ParamIDs with pilot 3's (100, 0.001),
# the old Multiplicative runs (submitted figures) and the production target-based run.
#   1. C1 per set: single-high and all-high occupancy, Fates (single-high states with occupancy > 1%)
#   2. C1 vs C2 vs C4 (stochastic effect)
#   3. C2 per frozen draw (10 ICs each): number of distinct discrete states reached, i.e. whether a
#      frozen λ vector leaves the circuit multistable
#
# Rscript compare_pilot4.r <tables dir> <outDir>
#   tables dir holds <run_id>/<NET>_Geometric_{occupancy,replicates}.csv.gz (collect.jl output) and
#   <run_id>/<NET>_thresholds.csv (the run's selection/<NET>_thresholds.csv)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
TAB  <- args[1]; OUT <- args[2]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
OLD  <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/final"
NEW  <- "/Users/kishorehari/Desktop/PostDoc/Abhay_Lakshmi/RACIPEdata/regnoise/2026-10-csbj-v2/legacy"
runs <- c(a100_f0.01   = "pilot4-2026-11-geometric-act100-floor01",
          a1000_f0.001 = "pilot4-2026-11-geometric-act1000-floor001",
          a100_f0.001  = "pilot3-2026-11-geometric-floor001")
pilot4 <- names(runs)[1:2]
nets <- c("TT4", "TT4SA")
N_DET_ICS <- 10

zcat <- function(f) fread(cmd = paste("gzip -dc", shQuote(f)))
nh_letters <- function(s) nchar(gsub("L", "", s))
nh_tuple   <- function(s) lengths(regmatches(s, gregexpr("1", s)))
per_set <- function(d, nn) d[, .(Single = sum(Occ[k == 1]), AllHigh = sum(Occ[k == nn]),
                                 Fates = sum(Occ[k == 1] > 0.01)),
                             by = .(Source, Condition, ParamID, Sigma, DT)]

read_geo <- function(net, label, ids = NULL) {
    f <- file.path(TAB, runs[[label]], paste0(net, "_Geometric_occupancy.csv.gz"))
    if (!file.exists(f)) return(NULL)
    d <- zcat(f)[Block == 0 & ThresholdMult == 1]
    if (!is.null(ids)) d <- d[ParamID %in% ids]
    nn <- nchar(d$State[1])
    d <- d[, .(Source = label, Condition, ParamID, Sigma, DT, State, Occ = Mean, k = nh_letters(State))]
    dts <- unique(d[Condition != "C0", DT])
    d <- rbind(d[Condition != "C0"], rbindlist(lapply(dts, function(x) d[Condition == "C0"][, DT := x])))
    per_set(d, nn)
}
read_legacy <- function(folder, mode, net, label, ids) {
    f <- file.path(folder, mode, net, "results", "all_parameters_results.csv")
    d <- fread(f, select = c("ParamID", "NoiseLevel", "DT", "State", "MRT"), blank.lines.skip = TRUE)[!is.na(ParamID) & ParamID %in% ids]
    nn <- nh_tuple(d$State[1]) + lengths(regmatches(d$State[1], gregexpr("0", d$State[1])))
    d <- d[, .(Occ = sum(MRT)), by = .(ParamID, Sigma = NoiseLevel, DT, State)]
    d[, Occ := Occ / sum(Occ), by = .(ParamID, Sigma, DT)]
    d[, `:=`(Source = label, Condition = "C1", k = nh_tuple(State))]
    per_set(d, nn)
}

all <- rbindlist(lapply(nets, function(net) {
    g <- rbindlist(lapply(pilot4, read_geo, net = net))
    if (nrow(g) == 0) return(NULL)
    ids <- unique(g$ParamID)
    oldMode <- if (net == "TT4SA") "MultiplicativeInvLambda" else "Multiplicative"
    rbind(g, read_geo(net, "a100_f0.001", ids), read_legacy(OLD, oldMode, net, "Old", ids),
          read_legacy(NEW, "Multiplicative", net, "Prod_target", ids))[, Network := net]
}))
fwrite(all, file.path(OUT, "pilot4_perset.csv.gz"))

nsrc <- all[Condition == "C1", .(n = uniqueN(Source)), by = .(Network, ParamID)]
common <- nsrc[, .SD[n == max(n)], by = Network][, .(Network, ParamID)]
a <- merge(all, common, by = c("Network", "ParamID"))
summ <- a[Condition == "C1", .(N = uniqueN(ParamID), Single = round(mean(Single), 3),
                               AllHigh = round(mean(AllHigh), 3), Fates = round(mean(Fates), 2)),
          by = .(Network, DT, Sigma, Source)]
wide <- dcast(summ, Network + DT + Sigma ~ Source, value.var = c("Single", "AllHigh", "Fates"))
setorder(wide, Network, DT, Sigma)
fwrite(wide, file.path(OUT, "pilot4_summary.csv"))
options(width = 250)
for (net in unique(wide$Network)) {
    cat("\n=====", net, "(", summ[Network == net, max(N)], "shared parameter sets ) C1\n")
    print(wide[Network == net, .SD, .SDcols = c("DT", "Sigma", grep("^(Single|AllHigh|Fates)", names(wide), value = TRUE))])
}

cmp <- dcast(a[Condition == "C1"], Network + ParamID + Sigma + DT ~ Source, value.var = "Single")
geo <- intersect(names(runs), names(cmp))
cat("\nMedian |single-high − old| per parameter set, σ > 0\n")
print(cmp[Sigma > 0, c(list(N = .N), lapply(.SD, function(x) round(median(abs(x - Old), na.rm = TRUE), 3))),
          by = .(Network, DT), .SDcols = c(geo, "Prod_target")][order(Network, DT)])

ctl <- a[Source %in% geo & Condition %in% c("C0", "C1", "C2", "C4")]
eff <- ctl[, .(Single = round(mean(Single), 3)), by = .(Network, Source, DT, Sigma, Condition)]
eff <- dcast(eff, Network + Source + DT + Sigma ~ Condition, value.var = "Single")
fwrite(eff, file.path(OUT, "pilot4_controls.csv"))
cat("\nMean single-high by condition (C4 only at Δt = 0.01, σ ∈ {0.01, 0.1, 1})\n")
print(eff[DT == 0.01][order(Network, Source, Sigma)])

# C2 per frozen draw: classify each deterministic replicate by its time-mean vs the thresholds
draws <- rbindlist(lapply(nets, function(net) rbindlist(lapply(pilot4, function(label) {
    f <- file.path(TAB, runs[[label]], paste0(net, "_Geometric_replicates.csv.gz"))
    th <- file.path(TAB, runs[[label]], paste0(net, "_thresholds.csv"))
    if (!file.exists(f) || !file.exists(th)) return(NULL)
    th <- fread(th)
    r <- zcat(f)[Condition == "C2"]
    M <- as.matrix(r[, paste0("Mean_", th$Node), with = FALSE])
    hi <- sweep(M, 2, th$Threshold, ">")
    r <- r[, .(ParamID, StabilityClass, Sigma, DT, Draw = (Replicate - 1) %/% N_DET_ICS + 1,
               State = apply(hi, 1, function(z) paste(as.integer(z), collapse = "")), K = rowSums(hi))]
    r[, .(NStates = uniqueN(State), NSingle = uniqueN(State[K == 1]), AnySingle = any(K == 1)),
      by = .(ParamID, StabilityClass, Sigma, DT, Draw)][, `:=`(Network = net, Source = label)]
}))))
if (nrow(draws)) {
    fwrite(draws, file.path(OUT, "pilot4_c2_draws.csv.gz"))
    cat("\nC2 frozen draws (10 ICs each): fraction multistable, mean # distinct states / single-high states\n")
    print(draws[, .(NDraws = .N, FracMultistable = round(mean(NStates > 1), 3),
                    MeanStates = round(mean(NStates), 2), MeanSingle = round(mean(NSingle), 2)),
                by = .(Network, Source, DT, Sigma)][order(Network, Source, DT, Sigma)])
    cat("\nsame, by stability class at Δt = 0.01\n")
    print(draws[DT == 0.01, .(FracMultistable = round(mean(NStates > 1), 3), MeanStates = round(mean(NStates), 2)),
                by = .(Network, Source, StabilityClass, Sigma)][order(Network, Source, StabilityClass, Sigma)])
}
