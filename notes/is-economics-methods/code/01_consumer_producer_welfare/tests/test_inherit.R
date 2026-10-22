# ------------------------------------------------------------------
# tests/test_inherit.R -- block W0 reproduces every number of the legacy chapter-1 script
# (snapshot: dev markdown/upgrade_2026-10-19/sdd_partI/legacy/code/01_consumer_producer_welfare/01_welfare.R)
# to the precision the legacy script and text print, and the card's list for W0: 27.73, 31.95, 24.21;
# the Slutsky decomposition -10 = -6 + (-4); the Slutsky matrix and its eigenvalues; the Gorman numbers
# 70, 95.5, 88.0; CE 100, risk premium 4 and its approximation 3.846; SOSD 10.196 against 10.000;
# cost 40 and z = (5, 20).
# Run from the chapter code directory after 01_dgp.R (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("welfare_fns.R")
W0 <- readRDS("truth01.rds")$W0
at <- function(x, v, d) length(x) == length(v) && all(abs(round(x, d) - v) < 1e-9)   # equal at d printed decimals

## ---- legacy 6.1: UMP, demand and indirect utility ----------------------------------
stopifnot(W0$a == 2 / 5, W0$w == 100, W0$p2 == 1, W0$p1_0 == 2, W0$p1_1 == 4)
stopifnot(at(W0$x1_0, 20, 4), at(W0$x2_0, 60, 4), at(W0$x1_1, 10, 4), at(W0$x2_1, 60, 4),
          at(W0$u0, 38.6636, 4), at(W0$u1, 29.3016, 4), at(W0$u0, 38.664, 3), at(W0$u1, 29.302, 3),
          at(W0$spend_0, 40, 6), at(W0$spend_1, 40, 6))

## ---- legacy 6.2: expenditure function and Shephard's lemma --------------------------
stopifnot(at(W0$e_check_0, 100, 6), at(W0$e_check_1, 100, 6), at(W0$h1_0, 20, 4))

## ---- legacy 6.3: Slutsky decomposition -10 = -6 + (-4) -------------------------------
stopifnot(at(W0$slutsky_total, -10, 9), at(W0$slutsky_subst, -6, 9), at(W0$slutsky_income, -4, 9),
          abs(W0$slutsky_total - (W0$slutsky_subst + W0$slutsky_income)) < 1e-12)

## ---- legacy 6.4: the three welfare numbers, EV < dCS < CV ------------------------------
stopifnot(at(W0$ev, 24.2142, 4), at(W0$dcs, 27.7259, 4), at(W0$cv, 31.9508, 4),
          at(W0$ev, 24.214, 3), at(W0$dcs, 27.726, 3), at(W0$cv, 31.951, 3),
          at(W0$ev, 24.21, 2), at(W0$dcs, 27.73, 2), at(W0$cv, 31.95, 2), W0$ev < W0$dcs, W0$dcs < W0$cv)
# the same numbers from the chapter's functions
e_cd <- function(p1, u) u * (p1 / 0.4)^0.4 * (1 / 0.6)^0.6
ce <- cv_ev_exact(e_cd, 2, 4, 100)
stopifnot(abs(ce[["cv"]] - W0$cv) < 1e-9, abs(ce[["ev"]] - W0$ev) < 1e-9,
          abs(dcs_numeric(function(p1, m) 0.4 * m / p1, 2, 4, 100) - W0$dcs) < 1e-10)

## ---- legacy 6.5: producer duality and the square-root lottery -----------------------
stopifnot(at(W0$cost, 40, 9), at(W0$z1, 5, 9), at(W0$z2, 20, 9), at(W0$cost_check, 40, 9))
stopifnot(at(W0$Ew, 104, 9), at(W0$Eu, 10, 9), at(W0$ce, 100, 9), at(W0$rp, 4, 9),
          at(W0$rp_approx, 3.8462, 4), at(W0$rp_approx, 3.846, 3), at(W0$A_Ew, 0.0048, 4), at(W0$varw, 1600, 9))

## ---- legacy 6.6: Slutsky matrix, Gorman aggregation, SOSD ---------------------------
stopifnot(identical(dim(W0$S), c(2L, 2L)), at(as.vector(W0$S), c(-6, 12, 12, -24), 9),
          at(W0$S_sym, 0, 9), at(W0$S_p1, 0, 9), at(W0$S_eig, c(0, -30), 9))
stopifnot(at(c(W0$x1_A_100, W0$x1_B_200, W0$x1_A_150, W0$x1_B_150), c(21, 49, 31, 39), 9),
          at(W0$X1, 70, 9), at(W0$X1_redis, 70, 9), at(W0$X1_b07, 95.5, 1), at(W0$X1_b07_redis, 88.0, 1))
stopifnot(at(W0$mean_F, 104, 4), at(W0$mean_G, 104, 4), at(W0$eu_F, 10.1962, 4), at(W0$eu_F, 10.196, 3),
          at(W0$eu_G, 10.0000, 4), at(W0$eu_G, 10.000, 3), W0$eu_F > W0$eu_G)

## ---- legacy 4.4 check: CRRA local risk premia rise with rho ----------------------------
stopifnot(at(W0$crra_rp_local, c(0.24, 0.481, 0.962), 3), all(diff(W0$crra_rp_local) > 0))

## ---- the card's W0 Afriat check: the two observed choices satisfy GARP ------------------
stopifnot(isTRUE(W0$garp_pass))
