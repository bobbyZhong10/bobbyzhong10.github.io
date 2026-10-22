# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Tamarind data and truth07.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("gmm_fns.R")
tr <- readRDS("truth07.rds")
G1 <- readRDS("tamarind_G1.rds"); G1h <- readRDS("tamarind_G1h.rds"); G2 <- readRDS("tamarind_G2.rds"); G3 <- readRDS("tamarind_G3.rds")
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))
P <- tr$params

## ---- truth07 fields named on the chapter card --------------------------------------------------
stopifnot(fin(tr$G1$theta), fin(tr$G2$theta), fin(tr$G3$theta),
          tr$G1$theta == -1.6, tr$G2$theta == -1.6, tr$G3$theta == -1.6)          # one elasticity, three designs
sd <- unlist(P$seed)
stopifnot(length(sd) >= 8, !anyDuplicated(sd), all(sd == round(sd)))          # every block its own seed

## ---- the generators reproduce the saved data from the stored seeds ---------------------------------
same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 0, check.attributes = TRUE))
set.seed(P$seed$G1);  stopifnot(same(gen_G1(P$G1), G1))
set.seed(P$seed$G1h); stopifnot(same(gen_G1(P$G1, tr$G1h$theta_m), G1h))
set.seed(P$seed$G2);  stopifnot(same(gen_G2(P$G2), G2))
set.seed(P$seed$G3);  stopifnot(same(gen_G3(P$G3), G3))

## ---- G1 and G1h: at least 12 markets, one A/B test each, equal arms ----------------------------------
P1 <- P$G1; r <- length(P1$markets)
stopifnot(r >= 12, tr$G1$r == r, all(P1$n %% 2 == 0), nrow(G1) == sum(P1$n), nrow(G1h) == sum(P1$n),
          identical(levels(G1$market), P1$markets), all(table(G1$market) == P1$n),
          all(tapply(G1$B, G1$market, mean) == 0.5), all(tapply(G1h$B, G1h$market, mean) == 0.5))
for (d in list(G1, G1h)) {                                                       # two prices per market: base and +delta
  lp <- split(d$log_price, d$market)
  stopifnot(all(vapply(seq_len(r), function(m) isTRUE(all.equal(sort(unique(lp[[m]])),
                       log(P1$p0[m]) + c(0, log(1 + P1$delta[m])), tolerance = 1e-12)), logical(1))))
}
stopifnot(all(P1$p0 >= 8.99 & P1$p0 <= 11.49), all(P1$delta >= 0.03 & P1$delta <= 0.12))
med <- tapply(exp(G1$y[G1$B == 0]), G1$market[G1$B == 0], median)               # weekly lessons at the base price
stopifnot(all(abs(med / P1$q0 - 1) < 0.1), all(exp(G1$y) > 0))
# truth: the variance of pi_hat_m and the optimal MD weights
stopifnot(max(abs(tr$G1$avar - P1$sigma^2 * (1 + P1$het^2) * 2 / P1$n / log(1 + P1$delta)^2)) < 1e-15,
          abs(sum(tr$G1$weights) - 1) < 1e-12, abs(tr$G1$se_md - sqrt(1 / sum(1 / tr$G1$avar))) < 1e-15,
          max(tr$G1$weights) / min(tr$G1$weights) > 20)                         # tests differ a lot in precision
# realised per-market slopes are close to their truth at the scale of their standard errors
ms <- market_slopes(G1)
stopifnot(all(abs(ms$pi_hat - tr$G1$theta) / sqrt(tr$G1$avar) < 4), all(abs(ms$se / sqrt(tr$G1$avar) - 1) < 0.1))
# G1h: heterogeneous elasticities, less elastic where income is higher; theta_star is the precision-weighted mean
th <- tr$G1h$theta_m
stopifnot(fin(th, r), abs(mean(th) - P1$theta) < 1e-12, cor(th, P1$income) > 0.99, diff(range(th)) > 1,
          abs(tr$G1h$theta_star - sum(tr$G1$weights * th)) < 1e-12, abs(tr$G1h$theta_star - tr$G1h$theta_mean) > 0.1,
          tr$G1h$power_J_01 > 0.99)

## ---- G2: regions, 10 experiments in natural units (one unbalanced), endogenous price -----------------------
P2 <- P$G2; K2 <- length(P2$step)
z2 <- as.matrix(G2[, paste0("z", seq_len(K2))])
stopifnot(nrow(G2) == P2$n, K2 == 10, all(abs(colSums(z2)) < 1e-9),           # centred assignments: mean exactly 0
          all(vapply(seq_len(K2), function(j) isTRUE(all.equal(sort(unique(z2[, j])), P2$step[j] * (c(0, 1) - P2$share[j]))), logical(1))),
          all(colSums(z2 > 0) == P2$n * P2$share), P2$share[2] == 0.2, all(P2$share[-2] == 0.5),
          max(abs(G2$log_price - log(P2$p_base) - drop(z2 %*% P2$gamma) - G2$v_true)) < 1e-12,
          max(abs(G2$y - log(P2$q_base) - P2$theta * (G2$log_price - log(P2$p_base)) - G2$eps_true)) < 1e-12,
          max(tr$G2$eff_sd) < 0.03, min(tr$G2$eff_sd) > 0.004)                    # each experiment: sd of its price effect 0.4% to 3%
stopifnot(cor(G2$v_true, G2$eps_true) > 0.5,                                    # discretionary promotions follow demand
          coef(lm(y ~ log_price, data = G2))[2] > P2$theta + 0.5,                # so OLS is biased towards zero
          tr$G2$mu2 > 100, tr$G2$ols_plim > P2$theta + 0.5)                      # strong instruments as a set
# heteroskedasticity: the deep-discount regions are noisier; S_0 is not proportional to E[zz'], so efficient GMM beats
# 2SLS asymptotically, and identity weighting (units choose the weights) is far worse
stopifnot(sd(G2$eps_true[z2[, 2] > 0]) > 2 * sd(G2$eps_true[z2[, 2] < 0]),
          tr$G2$se_eff < 0.95 * tr$G2$se_2sls, tr$G2$se_identity > 1.5 * tr$G2$se_eff)
# pop_moments (closed-form D_0 and S_0) agrees with a large sample from the generator
set.seed(1)                                                                      # (check only; not part of the data)
big <- iv_data(gen_G2(P2, n = 400000))
stopifnot(max(abs(-crossprod(big$Z, big$X) / 400000 - tr$G2$D0)) < 0.01,
          max(abs(moment_var(g_iv(tr$G2$par, big)) - tr$G2$S0)) / max(abs(tr$G2$S0)) < 0.02)
stopifnot(P2$R_mc >= 1000, all(P2$R_grid >= 1000), P2$n %in% P2$grid_n, all((P2$grid_n %o% P2$share) %% 1 == 0))

## ---- G3: penny-sized perturbations, weak moments --------------------------------------------------------
P3 <- P$G3; K3 <- length(P3$step)
z3 <- as.matrix(G3[, paste0("z", seq_len(K3))])
stopifnot(nrow(G3) == P3$n, all(abs(colSums(z3)) < 1e-12), all(P3$share == 0.5),
          max(P3$step) / 2 <= 0.02,                                              # at most +/- 2 cents per lesson credit
          max(tr$G3$eff_sd) < 0.002, tr$G3$mu2 < 15, tr$G3$F_pop < 5,
          max(abs(G3$log_price - log(P3$p_base) - drop(z3 %*% P3$gamma) - G3$v_true)) < 1e-12)
fs <- summary(lm(G3$log_price ~ z3))$fstatistic                                  # realised first-stage F is small
stopifnot(fs[["value"]] < 10, P3$R_mc >= 1000)
