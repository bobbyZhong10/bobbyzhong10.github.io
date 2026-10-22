# ------------------------------------------------------------------
# 02_fuzzy.R -- Block A (fuzzy take-up): Wald-DID, Wald-TC and Wald-CIC with
# city-cluster bootstrap confidence intervals, on tessera_A (a stable
# control-group take-up rate, dCDH 2018 Assumption 2) and, for the
# partial-identification discussion, on variant A' (tessera_A_rising, where
# control-group take-up rises).
#
# Estimators: the hand-coded functions in fuzzy_did_fns.R (authoritative;
# de Chaisemartin & D'Haultfoeuille 2018, checked in Task 10).
#
# Bootstrap: B = 999 resamples of CITIES (the cluster). Each resample draws
# the 80 treatment-group cities and the 80 control-group cities separately,
# with replacement (the group is a feature of the design, so strata keep the
# four (g, t) cells non-empty), takes all couriers of the drawn cities, and
# recomputes the three estimators. Intervals are percentile intervals, the 2.5%
# and 97.5% quantiles of the 999 draws (quantile type 7); se is the standard
# deviation of the draws. The same seed is used for A and A', so both
# variants see the same 999 city resamples.
#
# Output res_fuzzy.rds:
#   $did, $tc, $cic   list(est, lo, hi, se)   block A (tessera_A)
#   $rising           list(did, tc, cic)       same fields, variant A'
#   $parts            first-stage and building blocks: DID_Y, DID_D, E(D_gt),
#                     delta_0, delta_1 (Wald-TC), for A and A'
#   $cells            cell means (n, take-up rate, mean y) by (g, t), A and A'
#   $boot             999 x 3 matrices of bootstrap draws (A, A')
#   $meta             B, seed, stratification, runtime, warnings
# ------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(code_dir, "fuzzy_did_fns.R"))

B_BOOT <- 999L
SEED   <- 20201L

A  <- readRDS(file.path(code_dir, "tessera_A.rds"))
AR <- readRDS(file.path(code_dir, "tessera_A_rising.rds"))
tr <- readRDS(file.path(code_dir, "truth20.rds"))
stopifnot(identical(A$city, AR$city), identical(A$unit, AR$unit))   # same couriers, so the same resamples apply

warn_log <- character()
with_log <- function(expr) withCallingHandlers(expr, warning = function(w) {
  warn_log <<- c(warn_log, conditionMessage(w)); invokeRestart("muffleWarning") })

est3 <- function(X) c(did = wald_did(X, "y", "d", "g", "t")$estimate,
                      tc  = wald_tc (X, "y", "d", "g", "t")$estimate,
                      cic = wald_cic(X, "y", "d", "g", "t")$estimate)

## ---- point estimates ----------------------------------------------------------------
point <- function(X) {
  wd <- wald_did(X, "y", "d", "g", "t"); wt <- wald_tc(X, "y", "d", "g", "t"); wc <- wald_cic(X, "y", "d", "g", "t")
  list(est = c(did = wd$estimate, tc = wt$estimate, cic = wc$estimate),
       parts = list(did_y = wd$did_y, did_d = wd$did_d, first_stage_treated = wt$first_stage,
                    delta_0 = unname(wt$delta["0"]), delta_1 = unname(wt$delta["1"])))
}
pA <- point(A); pAR <- point(AR)

## ---- city-cluster bootstrap -----------------------------------------------------------
## row indices of each city, and the cities of each group
rows_by_city <- split(seq_len(nrow(A)), A$city)
cities_by_g  <- lapply(split(unique(A[, .(city, g)]), by = "g", keep.by = FALSE), function(z) as.character(z$city))
cities_by_g  <- cities_by_g[c("0", "1")]

boot_both <- function(B, seed) {
  set.seed(seed)
  outA <- outAR <- matrix(NA_real_, B, 3, dimnames = list(NULL, c("did", "tc", "cic")))
  for (b in seq_len(B)) {
    drawn <- unlist(lapply(cities_by_g, function(cc) sample(cc, length(cc), replace = TRUE)), use.names = FALSE)
    rows  <- unlist(rows_by_city[drawn], use.names = FALSE)
    outA[b, ]  <- with_log(est3(A[rows]))
    outAR[b, ] <- with_log(est3(AR[rows]))
  }
  list(A = outA, AR = outAR)
}
t0 <- Sys.time()
bt <- boot_both(B_BOOT, SEED)
runtime_boot <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
stopifnot(all(is.finite(bt$A)), all(is.finite(bt$AR)))     # no degenerate resample

summ <- function(est, draws) {
  q <- quantile(draws, c(0.025, 0.975), names = FALSE)
  list(est = unname(est), lo = q[1], hi = q[2], se = sd(draws))
}
pack <- function(p, draws) setNames(lapply(c("did", "tc", "cic"), function(k) summ(p$est[[k]], draws[, k])),
                                    c("did", "tc", "cic"))
resA  <- pack(pA, bt$A)
resAR <- pack(pAR, bt$AR)

## ---- cell means -------------------------------------------------------------------------
cells <- function(X, label) {
  z <- X[, .(n = .N, take_up = mean(d), mean_y = mean(y)), keyby = .(g, t)]
  z[, variant := label]; setcolorder(z, "variant"); as.data.frame(z)
}

res_fuzzy <- list(
  did = resA$did, tc = resA$tc, cic = resA$cic,
  rising = resAR,
  parts = list(A = pA$parts, rising = pAR$parts),
  cells = rbind(cells(A, "A"), cells(AR, "A_rising")),
  boot = bt,
  meta = list(B = B_BOOT, seed = SEED, ci = "percentile 95%, city bootstrap stratified by group",
              n_cities = c(treated = length(cities_by_g[["1"]]), control = length(cities_by_g[["0"]])),
              n_couriers = uniqueN(A$unit), runtime_boot_sec = runtime_boot,
              warnings = warn_log)
)
saveRDS(res_fuzzy, file.path(code_dir, "res_fuzzy.rds"))

## ---- report -----------------------------------------------------------------------------
fmt <- function(r) sprintf("%8.3f  [%7.3f, %7.3f]  se %.3f", r$est, r$lo, r$hi, r$se)
cat(sprintf("== Block A: true LATE %.4f (population %.4f); B = %d, seed %d, %.0f s ==\n",
            tr$A$late, tr$A$late_pop, B_BOOT, SEED, runtime_boot))
for (k in c("did", "tc", "cic"))
  cat(sprintf("  Wald-%-4s %s   covers truth: %s\n", toupper(k), fmt(res_fuzzy[[k]]),
              res_fuzzy[[k]]$lo <= tr$A$late && tr$A$late <= res_fuzzy[[k]]$hi))
cat(sprintf("\n== Variant A': true LATE %.4f (same treatment group) ==\n", tr$A_rising$late))
for (k in c("did", "tc", "cic"))
  cat(sprintf("  Wald-%-4s %s   covers truth: %s\n", toupper(k), fmt(resAR[[k]]),
              resAR[[k]]$lo <= tr$A_rising$late && tr$A_rising$late <= resAR[[k]]$hi))
cat("\nCell means:\n"); print(res_fuzzy$cells, digits = 4, row.names = FALSE)
cat("\nBuilding blocks (A):\n"); str(pA$parts)
cat("\nBuilding blocks (A'):\n"); str(pAR$parts)
cat(sprintf("\nWarnings logged: %d\n", length(warn_log)))
cat("Saved res_fuzzy.rds\n")
