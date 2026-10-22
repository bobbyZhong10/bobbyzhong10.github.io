# tests/test_dgp.R
# Checks on the simulated Tessera world written by 01_dgp.R.
# Run from the code directory:  Rscript tests/test_dgp.R   (about 20 seconds:
# the last block runs did_multiplegt_dyn on the noise-free effect column)
suppressPackageStartupMessages(library(data.table))
source("fuzzy_did_fns.R")

A <- readRDS("tessera_A.rds"); B <- readRDS("tessera_B.rds"); H <- readRDS("tessera_B_had.rds"); C <- readRDS("tessera_C.rds")
AR <- readRDS("tessera_A_rising.rds")
tr <- readRDS("truth20.rds")

## ---- brief: design requirements -----------------------------------------------
rate <- function(g, t) mean(A$d[A$g == g & A$t == t])
stopifnot(rate(1, 0) <= 0.10, rate(1, 1) >= 0.45, abs(rate(0, 1) - rate(0, 0)) < 1e-12, rate(0, 0) > 0.05, rate(0, 0) < 0.95, abs(rate(0, 0) - rate(1, 0)) >= 0.15)  # 2026-10-20: stable control-group treatment rate (dCDH 2018 Assumption 2)
stopifnot(mean(B$dose[B$post == 1] == 0) >= 0.20, min(H$dose[H$post == 1]) > 0)
sw <- C[, .(n_switch = sum(diff(dose) != 0)), by = unit]
stopifnot(mean(sw$n_switch >= 2) >= 0.30, setequal(unique(C$dose), c(0, 1, 2)))

## ---- truth20 fields --------------------------------------------------------------
fin1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
stopifnot(fin1(tr$A$late), fin1(tr$A_rising$late), fin1(tr$B$att_overall), fin1(tr$C$avg_effect))
stopifnot(is.data.frame(tr$B$att_curve), all(c("dose", "att") %in% names(tr$B$att_curve)),
          nrow(tr$B$att_curve) >= 20, all(is.finite(tr$B$att_curve$att)))
stopifnot(is.numeric(tr$C$delta), length(tr$C$delta) == 8, all(is.finite(tr$C$delta)))

## ---- block A: fuzzy take-up, stable control group -----------------------------------
for (X in list(A, AR)) {
  stopifnot(is.data.table(X), all(c("unit", "city", "g", "t", "d", "y", "type_true", "tau_true") %in% names(X)))
  # balanced two-period panel of couriers; one group per courier and per city
  stopifnot(all(X[, .N, by = unit]$N == 2), all(X[, uniqueN(g), by = unit]$V1 == 1),
            all(X[, uniqueN(g), by = city]$V1 == 1), all(X[, uniqueN(city), by = g]$V1 >= 40))
  # take-up is monotone (Assumption 3: D = 1{V >= v_gt}, V time-invariant)
  wideX <- dcast(X, unit + g ~ t, value.var = "d")
  stopifnot(all(wideX[["1"]] >= wideX[["0"]]), all(is.finite(X$tau_true)))
}
# block A: no control-group switcher (D(0) = D(1) in g = 0), all three types in g = 1
stopifnot(!any(A[g == 0, type_true] == "switcher"), setequal(unique(A[g == 1, type_true]), c("always", "switcher", "never")))
# stored truth is the mean individual effect of treatment-group switchers at t = 1
stopifnot(abs(tr$A$late - A[g == 1 & t == 1 & type_true == "switcher", mean(tau_true)]) < 1e-10)
# population properties (closed form, and checked on a large auxiliary sample):
# Wald-TC and Wald-CIC recover the switchers' LATE (dCDH 2018 Theorems 2, 3);
# Wald-DID is biased because the effect of units treated at t = 0 changes
# (Assumption 5 fails): W_DID - Delta = (pA1 - pA0) lambda / pS1
stopifnot(abs(tr$A$pop_check$w_tc - tr$A$late_pop) < 0.15,
          abs(tr$A$pop_check$w_cic - tr$A$late_pop) < 0.15,
          abs(tr$A$pop_check$w_did - tr$A$late_pop) > 1.5,
          abs(tr$A$pop_check$w_did - tr$A$w_did_pop) < 0.25,
          abs(tr$A$w_did_pop - tr$A$late_pop -
              with(tr$params$A, (pA1 - pA0) * lambda / pS1)) < 1e-10)
# the hand-coded estimators run on the main sample
est <- c(wald_did(A, "y", "d", "g", "t")$estimate, wald_tc(A, "y", "d", "g", "t")$estimate,
         wald_cic(A, "y", "d", "g", "t")$estimate)
stopifnot(all(is.finite(est)))

## ---- variant A': rising control-group take-up, switchers self-selected on level -----
# same couriers, outcomes and treatment group as block A; only g = 0 differs
setkey(A, unit, t); setkey(AR, unit, t)
stopifnot(identical(A$unit, AR$unit), identical(A[g == 1], AR[g == 1]), identical(A$tau_true, AR$tau_true),
          isTRUE(all.equal(A$y - A$d * A$tau_true, AR$y - AR$d * AR$tau_true)))   # same Y(0) draws
rateR <- function(g, t) mean(AR$d[AR$g == g & AR$t == t])
stopifnot(rateR(0, 1) - rateR(0, 0) >= 0.05, rateR(0, 1) - rateR(0, 0) <= 0.15,
          abs(rateR(0, 0) - rate(0, 0)) < 1e-12, any(AR[g == 0, type_true] == "switcher"))
stopifnot(abs(tr$A_rising$late - AR[g == 1 & t == 1 & type_true == "switcher", mean(tau_true)]) < 1e-10)
# without Assumption 2, dCDH (2018) Theorem 4 only bounds Delta (the TC bounds
# are infinite with an unbounded outcome): the Wald-TC and Wald-CIC estimands
# miss the switchers' LATE by more than 2 orders a week (large-sample values)
stopifnot(abs(tr$A_rising$pop_check$w_tc - tr$A_rising$late_pop) > 2,
          abs(tr$A_rising$pop_check$w_cic - tr$A_rising$late_pop) > 2)
# the bias formula for W_TC (control switchers vs always-takers in Y(1) and vs
# never-takers in Y(0), evaluated with large-sample moments) matches the
# large-sample Wald-TC
stopifnot(length(tr$A_rising$tc_bias_parts) == 2,
          abs(tr$A_rising$pop_check$late + sum(tr$A_rising$tc_bias_parts) - tr$A_rising$pop_check$w_tc) < 0.25)

## ---- block B: continuous dose, with and without stayers -----------------------------
stopifnot(is.data.table(B), all(c("unit", "t", "post", "g", "dose", "d", "y", "tau_true") %in% names(B)))
for (X in list(B, H)) {
  stopifnot(all(X[, .N, by = unit]$N == 2), setequal(unique(X$t), 1:2), all(X$post == (X$t == 2)),
            all(X[, uniqueN(dose), by = unit]$V1 == 1),          # contdid: dose constant within unit
            all(X$d == X$dose * X$post),                         # d = treatment actually received
            all((X$g == 2) == (X$dose > 0)), all(X$g %in% c(0, 2)),
            all(X[post == 0, tau_true] == 0))
}
# same world: same units and pre-period outcomes; HAD doses equal B doses where B > 0
setkey(B, unit, t); setkey(H, unit, t)
stopifnot(identical(B$unit, H$unit), isTRUE(all.equal(B[t == 1, y], H[t == 1, y])),
          all(H[B$dose > 0, dose] == B[B$dose > 0, dose]))
# quasi-untreated units exist in the HAD variant (doses near zero)
stopifnot(mean(H$dose[H$post == 1] < 0.1) >= 0.05)
# stored overall ATT is the mean realized effect among treated cities
stopifnot(abs(tr$B$att_overall - B[post == 1 & dose > 0, mean(tau_true)]) < 1e-10)
# selection on gains: ATT(d|d) rises faster than the strong-PT curve ATT(d)
ac <- tr$B$att_curve
stopifnot(all(c("att_strong", "acrt", "slope_att") %in% names(ac)),
          ac$att[nrow(ac)] > ac$att_strong[nrow(ac)], ac$att[1] < ac$att_strong[1])

## ---- block C: switching, non-binary treatment ---------------------------------------
stopifnot(is.data.table(C), all(c("unit", "t", "g", "dose", "y", "tau_true") %in% names(C)))
setkey(C, unit, t)
nT <- uniqueN(C$t)
stopifnot(all(C[, .N, by = unit]$N == nT), all(C[t == 1, dose] == 0))     # everyone starts at 0
stopifnot(mean(C[t == 1, g] == 0) >= 0.15)                                # never-switchers
first_sw <- C[, .(F = if (any(dose != 0)) as.numeric(min(t[dose != 0])) else 0), by = unit]
stopifnot(all(first_sw$F == C[t == 1, g]), all(first_sw[F > 0, F] <= nT - 7))
stopifnot(all(C[g == 0 | t < g, tau_true] == 0))                          # no anticipation
stopifnot(abs(tr$C$avg_effect - C[dose > 0, sum(tau_true) / sum(dose)]) < 1e-10)
# did_multiplegt_dyn on the noise-free effect column reproduces the stored delta_l:
# the truth is defined exactly as the estimator aggregates switchers
suppressPackageStartupMessages({ library(polars); library(DIDmultiplegtDYN) })
nf <- as.data.frame(C[, .(unit, t, dose, tau_true)])
res <- suppressWarnings(did_multiplegt_dyn(nf, outcome = "tau_true", group = "unit", time = "t",
                                           treatment = "dose", effects = 8, graph_off = TRUE))
stopifnot(max(abs(res$results$Effects[, "Estimate"] - tr$C$delta)) < 1e-8)
resn <- suppressWarnings(did_multiplegt_dyn(nf, outcome = "tau_true", group = "unit", time = "t",
                                            treatment = "dose", effects = 8, normalized = TRUE, graph_off = TRUE))
stopifnot(max(abs(resn$results$Effects[, "Estimate"] - tr$C$delta_norm)) < 1e-8)

cat("test_dgp: OK\n")
