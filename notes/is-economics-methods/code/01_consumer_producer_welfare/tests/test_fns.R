# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of welfare_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# The first block is the chapter card's assertions, verbatim (comments translated, Ruling P6);
# the rest checks the functions against closed forms (Cobb-Douglas, Stone-Geary), against
# numerical derivatives (numDeriv) and, for GARP and the Afriat efficiency index, against the
# revealedPrefs package (a cross-check only; the chapter's numbers come from welfare_fns.R).
# ------------------------------------------------------------------
source("welfare_fns.R")

## ---- chapter card (verbatim) -------------------------------------------------------
e_cd <- function(p1, u) u * (p1 / 0.4)^0.4 * (1 / 0.6)^0.6
x_cd <- function(p1, m) 0.4 * m / p1
stopifnot(abs(vartia_cv(x_cd, 2, 4, 100) - cv_ev_exact(e_cd, 2, 4, 100)[["cv"]]) < 1e-6)   # the ODE recovers the exact CV
P <- rbind(c(1, 1), c(1, 2)); X <- rbind(c(2, 0), c(0, 1.5))                                 # two periods, each strictly revealed preferred to the other
stopifnot(!garp_check(P, X)$pass, abs(ccei(P, X) - 0.75) < 1e-6)
stopifnot(garp_check(rbind(c(2, 1), c(4, 1)), rbind(c(20, 60), c(10, 60)))$pass)             # W0's two choices are consistent
stopifnot(cara_reject_threshold(10, 11) > log(2) / 100)                                       # from here on every bet that loses 100 is rejected

## ---- cv_ev_exact and vartia_cv: Cobb-Douglas closed forms (block W0) ------------------
cv0 <- 100 * (2^0.4 - 1); ev0 <- 100 * (1 - 2^-0.4)
ce <- cv_ev_exact(e_cd, 2, 4, 100)
stopifnot(identical(names(ce), c("cv", "ev")), abs(ce[["cv"]] - cv0) < 1e-9, abs(ce[["ev"]] - ev0) < 1e-9)
stopifnot(abs(vartia_cv(x_cd, 2, 4, 100) - cv0) < 1e-9,
          abs(-vartia_cv(x_cd, 4, 2, 100) - ev0) < 1e-9)            # reversed path: m(p0) - w = e(p0, u1) - w = -EV
stopifnot(max(abs(vartia_cv(x_cd, 2, 4, c(50, 100, 200)) - c(0.5, 1, 2) * cv0)) < 1e-9)   # vectorised over w
# loss-positive convention: a price fall has a negative CV and EV
ce_fall <- cv_ev_exact(e_cd, 4, 2, 100)
stopifnot(ce_fall[["cv"]] < 0, ce_fall[["ev"]] < 0, abs(ce_fall[["cv"]] + ev0) < 1e-9, abs(ce_fall[["ev"]] + cv0) < 1e-9)
# dcs_numeric: the area under Marshallian demand, 40 log 2 for W0
stopifnot(abs(dcs_numeric(x_cd, 2, 4, 100) - 40 * log(2)) < 1e-10,
          max(abs(dcs_numeric(x_cd, 2, 4, c(100, 200)) - c(40, 80) * log(2))) < 1e-10)

## ---- Stone-Geary helpers (block W1): closed forms against the generic functions -------
b <- 0.3; g1 <- 1.2; g2 <- 7; w <- 110
e_sg <- function(p1, u) sg_expenditure(p1, u, b, g1, g2)
x_sg <- function(p1, m) sg_demand(p1, m, b, g1, g2)
cf <- sg_welfare(2, 4, w, b, g1, g2)
ex <- cv_ev_exact(e_sg, 2, 4, w)
stopifnot(all(c("cv", "ev", "dcs", "eta") %in% names(cf)),
          abs(cf$cv - ex[["cv"]]) < 1e-9, abs(cf$ev - ex[["ev"]]) < 1e-9,
          abs(cf$cv - vartia_cv(x_sg, 2, 4, w)) < 1e-9, abs(cf$ev + vartia_cv(x_sg, 4, 2, w)) < 1e-9,
          abs(cf$dcs - dcs_numeric(x_sg, 2, 4, w)) < 1e-9,
          cf$ev < cf$dcs, cf$dcs < cf$cv)                              # normal good, price rise: EV < dCS < CV
# duality: e(p, v(p, w)) = w; Shephard: de/dp1 = Hicksian demand = Marshallian demand at income e(p, u)
u0 <- sg_indirect(2, w, b, g1, g2)
stopifnot(abs(e_sg(2, u0) - w) < 1e-10,
          abs(numDeriv::grad(function(p) e_sg(p, u0), 3) - x_sg(3, e_sg(3, u0))) < 1e-7)
# eta: the income elasticity of good 1 at (p0, w)
stopifnot(abs(cf$eta - numDeriv::grad(function(m) x_sg(2, m), w) * w / x_sg(2, w)) < 1e-8)
# Cobb-Douglas is the special case g1 = g2 = 0
cd <- sg_welfare(2, 4, 100, 0.4, 0, 0)
stopifnot(abs(cd$cv - cv0) < 1e-12, abs(cd$ev - ev0) < 1e-12, abs(cd$dcs - 40 * log(2)) < 1e-12, abs(cd$eta - 1) < 1e-12)

## ---- willig_cv --------------------------------------------------------------------------
stopifnot(abs(willig_cv(40 * log(2), 1, 100) - (40 * log(2) + (40 * log(2))^2 / 200)) < 1e-12,
          willig_cv(10, 0, 100) == 10,
          max(abs(willig_cv(c(1, 2), c(0.5, 1), c(100, 50)) - c(1 + 0.5 / 200, 2 + 4 / 100))) < 1e-12)
# first-order accuracy (Ruling P13): the Willig error shrinks faster than the CS error as the price rise
# shrinks (Cobb-Douglas, eta = 1, p1 from 2 to 2 + d)
rel_err <- sapply(c(2, 0.5, 0.1), function(d) {
  cv <- 100 * ((1 + d / 2)^0.4 - 1); dcs <- 40 * log(1 + d / 2)
  (cv - willig_cv(dcs, 1, 100)) / (cv - dcs)
})
stopifnot(all(rel_err > 0), all(diff(rel_err) < 0), rel_err[3] < 0.02)

## ---- garp_check and ccei -----------------------------------------------------------------
g <- garp_check(P, X)
stopifnot(identical(names(g), c("pass", "n_viol", "R")), g$n_viol == 2, is.logical(g$R), all(dim(g$R) == c(2, 2)), all(g$R))
stopifnot(ccei(rbind(c(2, 1), c(4, 1)), rbind(c(20, 60), c(10, 60))) == 1)    # consistent data: CCEI is exactly 1
# tie case: the direct relation R0 is weak. p1 . x1 = p1 . x2 = 2, so x1 R0 x2 only through the tie, and
# p2 . x2 = 6 > p2 . x1 = 4, so x2 P0 x1 strictly: GARP fails (it would pass if R0 were strict)
Pt <- rbind(c(1, 1), c(3, 1)); Xt <- rbind(c(1, 1), c(2, 0))
stopifnot(!garp_check(Pt, Xt)$pass, garp_check(Pt, Xt)$R[1, 2], revealedPrefs::checkGarp(Xt, Pt)$violation)
# Warshall closure: 1 R 2 and 2 R 3 directly, 1 R 3 only through the chain
P3 <- rbind(c(1, 10), c(1, 0.1), c(10, 1)); X3 <- rbind(c(3.95, 0.005), c(3.9, 0), c(0, 30))
g3 <- garp_check(P3, X3)
stopifnot(g3$pass, g3$R[1, 2], g3$R[2, 3], g3$R[1, 3], sum(P3[1, ] * X3[3, ]) > sum(P3[1, ] * X3[1, ]))
# cross-check against revealedPrefs::checkGarp on random data (2 or 3 goods, 3 to 8 budgets)
stopifnot(revealedPrefs::checkGarp(X, P, afriat.par = 0.75)$violation, !revealedPrefs::checkGarp(X, P, afriat.par = 0.7499)$violation)
set.seed(1011)
for (r in 1:300) {
  Tn <- sample(3:8, 1); K <- sample(2:3, 1)
  Pr <- matrix(runif(Tn * K, 1, 3), Tn); m <- runif(Tn, 50, 100)
  S <- matrix(rexp(Tn * K), Tn); Xr <- (S / rowSums(S)) * m / Pr
  gr <- garp_check(Pr, Xr)
  stopifnot(gr$pass == !revealedPrefs::checkGarp(Xr, Pr)$violation)
  ci <- ccei(Pr, Xr)
  stopifnot(ci > 0, ci <= 1, (ci == 1) == gr$pass)
  if (ci < 1) stopifnot(!revealedPrefs::checkGarp(Xr, Pr, afriat.par = ci - 1e-6)$violation,
                        revealedPrefs::checkGarp(Xr, Pr, afriat.par = ci + 1e-6)$violation)
}

## ---- cara_reject_threshold ---------------------------------------------------------------
A <- cara_reject_threshold(10, 11)
eu <- function(A, L, G, q = 0.5) -q * exp(A * L) - (1 - q) * exp(-A * G)       # CARA expected utility of the bet; 0 USD gives -1
stopifnot(A > 0, abs(0.5 * exp(10 * A) + 0.5 * exp(-11 * A) - 1) < 1e-12,
          eu(1.01 * A, 10, 11) < -1, eu(0.99 * A, 10, 11) > -1,              # rejected above the threshold, accepted below
          eu(A, 100, 1e6) < -1)                                                # lose 100 against any gain: rejected
# the printed bound ln2 / A = 76.30 USD at A = A_rabin: a 50-50 bet losing 76.31 is rejected even against a
# gain of 1e9, one losing 76.29 is accepted against it
stopifnot(abs(log(2) / A - 76.30) < 0.005, eu(A, 76.31, 1e9) < -1, eu(A, 76.29, 1e9) > -1)
stopifnot(cara_reject_threshold(11, 10) == 0, cara_reject_threshold(10, 10) == 0)   # unfavourable or fair: every risk averter rejects
# prob = probability of the loss: the protection decision of block W3 (decline = take the bet
# "lose L - price with probability q, keep the price otherwise")
Ai <- cara_reject_threshold(150 - 18, 18, prob = 0.08)
stopifnot(Ai > 0, abs(log(1 - 0.08 + 0.08 * exp(150 * Ai)) / Ai - 18) < 1e-9)

## ---- W1 demand forms: the income-term form is exact on noise-free Stone-Geary data with a common beta ----
dd <- expand.grid(p1 = c(2, 2.5, 3, 3.5, 4), w = c(60, 90, 120, 200))
dd$x1 <- sg_demand(dd$p1, dd$w, 0.35, 1, 5)
fm <- fit_demand_forms(dd)
stopifnot(identical(names(fm), c("linear", "loglinear", "income")),
          all(vapply(fm, function(f) is.function(f$x1_fun) && is.function(f$eta_fun), logical(1))))
stopifnot(max(abs(fm$income$x1_fun(c(2, 3, 7, 12), 150) - sg_demand(c(2, 3, 7, 12), 150, 0.35, 1, 5))) < 1e-10,
          abs(fm$income$eta_fun(2, 150) - sg_welfare(2, 4, 150, 0.35, 1, 5)$eta) < 1e-10)
# each form's eta_fun is the elasticity of its own x1_fun
for (f in fm) stopifnot(abs(f$eta_fun(2, 120) - numDeriv::grad(function(m) f$x1_fun(2, m), 120) * 120 / f$x1_fun(2, 120)) < 1e-6)
# the linear form's demand is truncated at zero above its choke price
stopifnot(all(fm$linear$x1_fun(c(30, 50), 100) == 0))
