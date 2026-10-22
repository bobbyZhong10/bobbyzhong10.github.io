# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Marigold data and truth13.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("qr_fns.R")

## ---- chapter card (verbatim) ---------------------------------------------------
tr <- readRDS("truth13.rds")
stopifnot(max(abs(tr$Q1r$qte - tr$Q1$qte)) < 1e-12, max(abs(tr$Q1r$ind_effect_quantiles - tr$Q1r$qte)) > 0.2)

## ---- truth13 fields named on the chapter card --------------------------------------------
P <- tr$params; taus <- P$taus
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))
stopifnot(identical(taus, c(0.1, 0.25, 0.5, 0.75, 0.9)),
          fin(tr$Q1$qte, 5L), fin(tr$Q1$ate), fin(tr$Q1r$qte, 5L), fin(tr$Q1r$ind_effect_quantiles, 5L),
          fin(tr$Q2$uqe), fin(tr$Q3$lqte_median), fin(tr$Q4$beta_tau), fin(tr$Q5$qtt, 5L), fin(tr$Q5$lqte))
sd <- unlist(P$seed)
stopifnot(length(sd) >= 12, !anyDuplicated(sd), all(sd == round(sd)))                 # every block its own seed

## ---- the generators reproduce the saved data from the stored seeds ---------------------------
Q1 <- readRDS("marigold_Q1.rds"); Q1r <- readRDS("marigold_Q1r.rds"); Q2 <- readRDS("marigold_Q2.rds")
Q3 <- readRDS("marigold_Q3.rds"); Q4 <- readRDS("marigold_Q4.rds"); Q5 <- readRDS("marigold_Q5.rds")
Q5f <- readRDS("marigold_Q5f.rds")
same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 0, check.attributes = TRUE))
set.seed(P$seed$Q1); stopifnot(same(gen_Q1(P$Q1), Q1))
set.seed(P$seed$Q1r); stopifnot(same(gen_Q1r(P$Q1, Q1), Q1r))
set.seed(P$seed$Q2); stopifnot(same(gen_Q2(P$Q2), Q2))
set.seed(P$seed$Q3); stopifnot(same(gen_Q3(P$Q3), Q3))
set.seed(P$seed$Q4_T5); r5 <- gen_Q4(P$Q4, P$Q4$T_short)
set.seed(P$seed$Q4_T20); r20 <- gen_Q4(P$Q4, P$Q4$T_long)
stopifnot(same(r5, Q4$T5), same(r20, Q4$T20))
set.seed(P$seed$Q5); stopifnot(same(gen_Q5(P$Q5), Q5))
set.seed(P$seed$Q5f); stopifnot(same(gen_Q5f(P$Q5f), Q5f))

## ---- Q1: AutoPrice experiment (revenue in USD thousands) --------------------------------------
p1 <- P$Q1
stopifnot(nrow(Q1) == p1$n, sum(Q1$D) == p1$n_treat, setequal(unique(Q1$D), 0:1),          # complete randomisation
          identical(Q1$y, ifelse(Q1$D == 1, Q1$y1_true, Q1$y0_true)), all(Q1$y > 0),
          max(abs(Q1$y1_true - exp(p1$logc + p1$k * log(Q1$y0_true)))) < 1e-9)               # rank invariance: y1 = c y0^k
stopifnot(median(Q1$y0_true) > 7, median(Q1$y0_true) < 14, abs(cor(Q1$x, log(Q1$y0_true)) - p1$b / sqrt(p1$v)) < 0.03)
# truth: E[Y(1)] = E[Y(0)] exactly (c chosen for that); QTE(tau) = c q0^k - q0 with q0 = exp(m + sqrt(v) z_tau)
ey0 <- exp(p1$m + p1$v / 2); ey1 <- exp(p1$logc + p1$k * p1$m + p1$k^2 * p1$v / 2)
q0 <- exp(p1$m + sqrt(p1$v) * qnorm(taus))
stopifnot(abs(ey1 - ey0) < 1e-12, tr$Q1$ate == 0, max(abs(tr$Q1$qte - (exp(p1$logc) * q0^p1$k - q0))) < 1e-12,
          tr$Q1$qte[1] < -1, tr$Q1$qte[5] > 2)                                                 # bottom loses, top gains
# population check (1e6 merchants drawn in 01_dgp.R with their own seed)
stopifnot(max(abs(tr$Q1$pop$qte - tr$Q1$qte)) < 0.03, abs(tr$Q1$pop$ate) < 0.03,
          max(abs(tr$Q1r$pop$ind_effect_quantiles - tr$Q1r$ind_effect_quantiles)) < 0.05,
          max(abs(tr$Q1$pop$ind_effect_quantiles - tr$Q1$ind_effect_quantiles)) < 0.03)
# log scale: marginal QTE in log points and the conditional (given X = 0) QTE in log points
stopifnot(max(abs(tr$Q1$qte_log - (p1$logc + (p1$k - 1) * (p1$m + sqrt(p1$v) * qnorm(taus))))) < 1e-12,
          max(abs(tr$Q1$cqte_log - (p1$logc + (p1$k - 1) * (p1$m + p1$s * qnorm(taus))))) < 1e-12)

## ---- Q1r: the same marginals, individual effects re-paired -----------------------------------
stopifnot(nrow(Q1r) == p1$n, identical(Q1r$D, Q1$D), identical(Q1r$x, Q1$x), identical(Q1r$y0_true, Q1$y0_true),
          identical(Q1r$y, ifelse(Q1r$D == 1, Q1r$y1_true, Q1r$y0_true)),
          identical(Q1r$y[Q1r$D == 0], Q1$y[Q1$D == 0]))                                       # control arm identical
z1r <- (log(Q1r$y1_true) - p1$logc) / p1$k; z1r <- (z1r - p1$m) / sqrt(p1$v)
stopifnot(abs(mean(z1r)) < 0.05, abs(sd(z1r) - 1) < 0.05,                                     # Y1r has the marginal of Y1
          abs(cor(z1r, (log(Q1$y0_true) - p1$m) / sqrt(p1$v)) - p1$rho_r) < 0.05)              # Gaussian copula, rho_r

## ---- Q2: Fulfilled-by-Marigold share (log revenue) ------------------------------------------
p2 <- P$Q2
stopifnot(nrow(Q2) == p2$n, setequal(unique(Q2$fbm), 0:1), abs(mean(Q2$fbm) - p2$p) < 0.02,
          max(abs(Q2$logy - log(Q2$y))) < 1e-12, abs(cor(Q2$fbm, Q2$x)) < 0.05)
# UQPE(tau) = (F0(q) - F1(q)) / f(q) at q = F^{-1}(tau): equals the derivative of the exact mixture quantile in p
mq <- function(p, t) uniroot(function(q) mix_cdf_Q2(q, p2, p) - t, c(0, 20), tol = 1e-14)$root
h <- 1e-5
fd <- sapply(taus, function(t) (mq(p2$p + h, t) - mq(p2$p - h, t)) / (2 * h))
stopifnot(max(abs(fd - tr$Q2$uqpe)) < 1e-5, tr$Q2$uqe == tr$Q2$uqpe[taus == p2$tau_uqe],
          max(abs(tr$Q2$cqr - (p2$g + (p2$s1 - p2$s0) * qnorm(taus)))) < 1e-14,
          abs(tr$Q2$uqe - tr$Q2$cqr[taus == p2$tau_uqe]) > 0.1)                                # conditional and unconditional part
stopifnot(max(abs(tr$Q2$q - sapply(taus, function(t) mq(p2$p, t)))) < 1e-9,
          max(abs(tr$Q2$cf - sapply(taus, function(t) mq(p2$p + p2$dp, t)))) < 1e-9)

## ---- Q3: AdBoost adoption with a randomised ad-credit offer (log revenue) ----------------------
p3 <- P$Q3
stopifnot(nrow(Q3) == p3$n, sum(Q3$Z) == p3$n_z1, all(Q3$D == ifelse(Q3$Z == 1, Q3$D1_true, Q3$D0_true)),
          all(Q3$D1_true >= Q3$D0_true),                                                         # monotonicity: no defiers
          max(abs(Q3$logy - ifelse(Q3$D == 1, Q3$logy1_true, Q3$logy0_true))) < 1e-12)
sh <- prop.table(table(factor(Q3$type_true, c("always", "complier", "never"))))
stopifnot(max(abs(sh - c(p3$pi0, p3$pi1 - p3$pi0, 1 - p3$pi1))) < 0.03)
stopifnot(abs(tr$Q3$lqte_median - p3$dmu) < 1e-8,                                               # symmetric complier band
          max(abs(tr$Q3$struct - (p3$dmu + (p3$s1 - p3$s0) * qnorm(taus)))) < 1e-14,
          max(abs(tr$Q3$pop$lqte - tr$Q3$lqte)) < 0.01, abs(tr$Q3$lqte[taus == 0.5] - tr$Q3$lqte_median) < 1e-8)

## ---- Q4: merchant panel with fixed effects (log revenue) --------------------------------------
p4 <- P$Q4
stopifnot(nrow(Q4$T5) == p4$N * p4$T_short, nrow(Q4$T20) == p4$N * p4$T_long,
          all(Q4$T5$x > 0 & Q4$T5$x < 1), all(table(Q4$T20$id) == p4$T_long),
          abs(tr$Q4$beta_tau - (p4$beta + p4$sig * p4$gam * qnorm(p4$tau))) < 1e-14,
          cor(tapply(Q4$T20$x, Q4$T20$id, mean), tapply(Q4$T20$alpha_true, Q4$T20$id, mean)) > 0.2)   # FE correlated with x

## ---- Q5: sharp changes-in-changes (revenue in USD thousands) ---------------------------------
p5 <- P$Q5
stopifnot(nrow(Q5) == 4 * p5$nc, all(table(Q5$S, Q5$T) == p5$nc), identical(Q5$D, Q5$S * Q5$T),
          max(abs(Q5$y - ifelse(Q5$D == 1, Q5$y1_true, Q5$y0_true))) < 1e-12,
          max(abs(Q5$y0_true - h_Q5(Q5$z_true, Q5$T, p5))) < 1e-12)
z1 <- p5$mu1 + p5$sd1 * qnorm(taus)
stopifnot(max(abs(tr$Q5$qtt - p5$delta * h_Q5(z1, 1, p5))) < 1e-12,
          max(abs(tr$Q5$pop$qtt - tr$Q5$qtt)) < 0.05,
          min(abs(tr$Q5$qdid_pop - tr$Q5$qtt)) > 0.3)                                            # QDID is biased here

## ---- Q5f: fuzzy changes-in-changes -----------------------------------------------------
pf <- P$Q5f
stopifnot(nrow(Q5f) == 4 * pf$nc, all(table(Q5f$S, Q5f$T) == pf$nc),
          max(abs(Q5f$y - h_Q5(Q5f$z_true, Q5f$T, pf) * (1 + pf$delta)^Q5f$D)) < 1e-12)
dr <- with(Q5f, tapply(D, list(S, T), mean))
stopifnot(abs(dr["0", "0"] - dr["0", "1"]) < 0.04,                                               # stable control-group rate
          dr["1", "1"] - dr["1", "0"] > 0.35,
          all(Q5f$D[Q5f$type_true == "always"] == 1), all(Q5f$D[Q5f$type_true == "never"] == 0),
          all(Q5f$D[Q5f$type_true == "switcher"] == Q5f$T[Q5f$type_true == "switcher"]),
          all(Q5f$S[Q5f$type_true == "switcher"] == 1))
stopifnot(abs(tr$Q5$lqte - tr$Q5$lqte_grid[taus == pf$tau]) < 1e-12, max(abs(tr$Q5$pop$lqte_grid - tr$Q5$lqte_grid)) < 0.05,
          abs(tr$Q5$pop$late_sw - tr$Q5$late_sw) < 0.05)
