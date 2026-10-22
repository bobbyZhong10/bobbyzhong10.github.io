# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 1.
# The first block is the chapter card, verbatim (comments translated, Ruling P6), with one field
# renamed by Ruling P7: spread_newgood -> spread_outside (the loss from raising p1 from 2 to 12 USD,
# far outside the experiment's 2 to 4 USD support; Stone-Geary demand has no finite choke price).
# The second block locks the Ruling P13 check of the Willig formula and the block W3 calibration.
# Run from the chapter code directory after 01 to 04 (tests/run_tests.R does this).
# ------------------------------------------------------------------

## ---- chapter card (verbatim, Ruling P7 field name) --------------------------------------
source("../../tools/R/chapter_kit.R"); tr <- readRDS("truth01.rds")
w1 <- readRDS("res_welfare.rds"); rp <- readRDS("res_rp.rds")
stopifnot(all(tr$W1$ev_i <= tr$W1$dcs_i + 1e-12), all(tr$W1$dcs_i <= tr$W1$cv_i + 1e-12))  # ordering, subscriber by subscriber
stopifnot(abs(w1$willig - tr$W1$cv_mean) < 0.2 * abs(w1$dcs - tr$W1$cv_mean))             # the Willig correction removes most of the CS error
stopifnot(abs(w1$cv_rep_g - tr$W1g$cv_mean) < 0.01 * tr$W1g$cv_mean)                       # Gorman holds: the representative consumer is unbiased
stopifnot(abs(w1$cv_rep - tr$W1$cv_mean) > 3 * abs(w1$cv_rep_g - tr$W1g$cv_mean))          # Gorman fails: clearly biased
stopifnot(w1$spread_in_range < 0.03, w1$spread_outside > 0.10)                              # stable within the support, unstable far outside it
stopifnot(rp$pass_rational == 1, rp$ccei_rational_min > 1 - 1e-8, rp$pass_noisy < 0.9, rp$pass_bronars < 0.5)

## ---- Ruling P13: the Willig formula CV ~ dCS + eta dCS^2 / (2 w) under the loss-positive convention ----
wc <- w1$willig_check
# W0 (Cobb-Douglas, eta = 1): the correction removes most of the CS error, and a coefficient of 1
# instead of 1/2 does worse. (The 01.1 assert "the opposite sign does worse than dCS" was dropped: with
# eta > 0 and EV < dCS < CV it held automatically. The opposite sign is tested where it can fail: as the
# first-order EV, below.)
stopifnot(abs(wc$W0$willig - wc$W0$cv) < 0.2 * abs(wc$W0$dcs - wc$W0$cv),
          abs(wc$W0$willig_coef1 - wc$W0$cv) > abs(wc$W0$willig - wc$W0$cv))
# with the opposite sign the same formula approximates EV (MWG's ordering EV < dCS < CV)
stopifnot(abs(wc$W0$willig_minus - wc$W0$ev) < 0.2 * abs(wc$W0$dcs - wc$W0$ev))
# W1 truth, subscriber by subscriber with each subscriber's own eta at (p0, w_i): every subscriber's
# Willig value is closer to the exact CV than dCS is, and on average it removes most of the CS error;
# the mirror formula (opposite sign) removes most of the CS error relative to the exact mean EV
stopifnot(all(abs(tr$W1$willig_i - tr$W1$cv_i) < abs(tr$W1$dcs_i - tr$W1$cv_i)),
          abs(wc$W1$willig_mean - tr$W1$cv_mean) < 0.2 * abs(tr$W1$dcs_mean - tr$W1$cv_mean),
          abs(wc$W1$willig_coef1_mean - tr$W1$cv_mean) > abs(wc$W1$willig_mean - tr$W1$cv_mean),
          abs(wc$W1$willig_minus_mean - tr$W1$ev_mean) < 0.2 * abs(tr$W1$dcs_mean - tr$W1$ev_mean))

## ---- W1: the representative consumer the analyst fits is the population's representative consumer ----
# (the fitted income-term form is a Gorman form: its CV at mean income is the population LES with
# b = mean beta, and it misses the mean CV by the heterogeneity term the cross-section cannot reveal)
stopifnot(abs(w1$cv_rep - tr$W1$cv_rep_pop) < 0.25 * abs(tr$W1$cv_rep_pop - tr$W1$cv_mean),
          sign(w1$cv_rep - tr$W1$cv_mean) == sign(tr$W1$cv_rep_pop - tr$W1$cv_mean))

## ---- W3: buying protection at 20 USD against a 12 USD expected loss implies rejecting every 50-50 bet
## that loses 100 USD, whatever it wins (CARA) ------------------------------------------------------
rk <- readRDS("res_risk.rds")
stopifnot(rk$A_interval[1] < rk$A_interval[2], rk$A_interval[1] > log(2) / 100, rk$A_rabin > log(2) / 100)

## ---- W2: the pass/fail verdicts of garp_check agree with revealedPrefs::checkGarp for every user and
## every Bronars chooser (cross-check only), and each failing CCEI sits where revealedPrefs flips -----------
W2 <- readRDS("plover_W2.rds"); W2b <- readRDS("plover_W2_bronars.rds")
cross <- function(d, id, res) {
  for (g in split(d, d[[id]])) {
    X <- cbind(g$x1, g$x2); Pm <- cbind(g$p1, g$p2); r <- res[res$id == g[[id]][1], ]
    stopifnot(nrow(r) == 1L, r$pass == !revealedPrefs::checkGarp(X, Pm)$violation)
    if (!r$pass) stopifnot(!revealedPrefs::checkGarp(X, Pm, afriat.par = r$ccei - 1e-6)$violation,
                           revealedPrefs::checkGarp(X, Pm, afriat.par = r$ccei + 1e-6)$violation)
  }
}
cross(W2, "user", rp$users); cross(W2b, "chooser", rp$bronars)
