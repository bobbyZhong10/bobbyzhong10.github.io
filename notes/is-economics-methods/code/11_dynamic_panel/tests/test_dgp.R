# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Nimbus panels and truth11.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("dyn_fns.R")
tr <- readRDS("truth11.rds"); P <- tr$params
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))
md5_lines <- function(lines) { f <- tempfile(); on.exit(unlink(f)); writeLines(lines, f); unname(tools::md5sum(f)) }

## ---- truth11 fields named on the chapter card ------------------------------------------------
stopifnot(fin(tr$rho), fin(tr$rho_persist), fin(tr$D3$beta), fin(tr$D4$rho),
          tr$rho == 0.70, tr$rho_persist == 0.90, tr$D5$rho == tr$rho)
# every block has its own seed, all different; the inherited blocks keep 808, 7372 (inside gen_dyn) and 303
sd <- unlist(P$seed)
stopifnot(length(sd) >= 15, !anyDuplicated(sd), all(sd == round(sd)), P$seed$D2_mis == 1132, P$seed$D1p_ah_mc == 1161,
          P$seed$D1 == 808, P$seed$D1p == 7372, P$seed$D1T == 303)

## ---- inherited blocks: legacy code byte-identical and first; data regenerate exactly -------------
src <- readLines("01_dgp.R")
stopifnot(md5_lines(src[1:42]) == "af46138f5f415d37cce47517db31c858",          # = the legacy 01_dgp.R
          md5_lines(readLines("03_dynamic_nickell.R")) == "8162c1bf0b557f991df9b6f3dd14f346",
          grep("^gen_dyn <- function", src)[1] == 21, grep("set.seed", src)[1] == 22,
          grep("seed = 808", src)[1] == 32, grep("seed = 7372", src)[1] == 33,
          grep("set.seed", readLines("03_dynamic_nickell.R")) == 43)
e1 <- new.env(); eval(parse(text = src[21:31]), e1)                              # the legacy gen_dyn definition
nb <- readRDS("nimbus.rds")
stopifnot(identical(e1$gen_dyn(N = 500, T = 8, rho = 0.70, seed = 808), nb$dyn_main),
          identical(e1$gen_dyn(N = 500, T = 8, rho = 0.90, seed = 7372), nb$dyn_persist),
          identical(nb$truth, list(rho = 0.70, rho_persist = 0.90)),
          tr$rho == nb$truth$rho, tr$rho_persist == nb$truth$rho_persist,
          identical(deparse(nb$gen_dyn), deparse(e1$gen_dyn)))
# D1 and D1p truth: the exact Nickell probability limits at T = 7 regression periods (8 months)
stopifnot(abs(tr$D1$fe_plim - (0.70 + nickell_bias(0.70, 7))) < 1e-15,
          abs(tr$D1p$fe_plim - (0.90 + nickell_bias(0.90, 7))) < 1e-15,
          identical(tr$D1T$T_obs, c(4, 8, 16, 32)),
          max(abs(tr$D1T$fe_plim - (0.70 + nickell_bias(0.70, c(4, 8, 16, 32) - 1)))) < 1e-15,
          # the legacy row nickell_approx = rho - (1 + rho)/(T_obs - 1) is the order-1/T bias -(1 + rho)/T
          max(abs(tr$D1T$legacy_approx - (0.70 - 1.7 / (c(4, 8, 16, 32) - 1)))) < 1e-15,
          max(abs(tr$D1T$first_factor_plim - (0.70 + nickell_first_factor(0.70, c(4, 8, 16, 32) - 1)))) < 1e-15)

## ---- the generators reproduce the saved data from the stored seeds -------------------------------
blocks <- c(D2 = "nimbus_D2.rds", D3 = "nimbus_D3.rds", D4 = "nimbus_D4.rds")
D <- lapply(blocks, readRDS)
for (b in names(blocks)) {
  set.seed(P$seed[[b]])
  stopifnot(isTRUE(all.equal(get(paste0("gen_", b))(P[[b]]), D[[b]], tolerance = 0, check.attributes = TRUE)),
            is.data.frame(D[[b]]), nrow(D[[b]]) == P[[b]]$N * P[[b]]$T,
            identical(D[[b]]$t, rep(seq_len(P[[b]]$T), P[[b]]$N)), all(D[[b]]$y > 0))
}

## ---- D2: 70 key accounts over 16 months, stationary AR(1); instrument counts exceed N -----------
d <- D$D2; P2 <- P$D2; T2 <- P2$T
stopifnot(P2$N == 70, T2 == 16, P2$rho == 0.70, P2$R_mc >= 500, tr$D2$rho == P2$rho,
          # the misspecified MA(1) twin: the D2 design plus an MA coefficient, 200 samples, its own seed
          identical(P$D2_ma[setdiff(names(P2), "R_mc")], P2[setdiff(names(P2), "R_mc")]), P$D2_ma$ma_coef == 0.6,
          P$D2_ma$R_mc >= 200, P2$mid_lag == 4,
          abs(mean(d$y) / P2$m_mean - 1) < 0.05, max(abs(d$alpha_true - (1 - P2$rho) * d$m_true)) < 1e-12)
# the full system-GMM instrument set has more columns than there are merchants (r > N); collapsed has fewer
stopifnot(tr$D2$r_full == (T2 - 1) * (T2 - 2) / 2 + (T2 - 2) + (T2 - 1), tr$D2$r_full > P2$N,
          tr$D2$r_collapsed == (T2 - 2) + 1 + (T2 - 1), tr$D2$r_collapsed < P2$N)
# closed-form probability limits (N -> infinity) and a large-sample check (seed D2_large)
sw2 <- P2$sigma_e^2 / (1 - P2$rho^2)
stopifnot(abs(tr$D2$pooled_plim - (P2$rho * sw2 + P2$m_sd^2) / (sw2 + P2$m_sd^2)) < 1e-15,
          abs(tr$D2$fe_plim - (P2$rho + nickell_bias(P2$rho, T2 - 1))) < 1e-15,
          abs(tr$D2$ols_diff_plim + (1 - P2$rho) / 2) < 1e-15,
          tr$D2$large$N >= 5e4, abs(tr$D2$large$pooled - tr$D2$pooled_plim) < 0.005,
          abs(tr$D2$large$fe - tr$D2$fe_plim) < 0.005, abs(tr$D2$large$ols_diff - tr$D2$ols_diff_plim) < 0.005)

## ---- D3: lagged GMV and a predetermined promotion budget with feedback -------------------------------
d <- D$D3; P3 <- P$D3
stopifnot(P3$feedback < 0, tr$D3$beta == P3$beta, tr$D3$rho == P3$rho, all(d$x > 0),
          all(c("y", "x", "eps_true", "alpha_true") %in% names(d)))
L <- function(v, k = 1L) lag_within(v, d$id, d$t, k)
r3 <- d$y - P3$rho * L(d$y) - P3$beta * d$x - d$alpha_true - d$eps_true
stopifnot(max(abs(r3), na.rm = TRUE) < 1e-10)
# sequential, not strict, exogeneity: X_it is uncorrelated with eps_it but X_{i,t+1} responds to eps_it
lead_x <- lag_within(d$x, d$id, d$t, -1L)
stopifnot(abs(cor(d$x, d$eps_true)) < 0.03, cor(lead_x, d$eps_true, use = "complete.obs") < -0.3)
# large-sample FE (seed D3_large): biased for beta with feedback; the same draws without feedback
stopifnot(tr$D3$large$N >= 1e5, abs(tr$D3$large$fe_beta - P3$beta) > 0.05, abs(tr$D3$large$fe_rho - P3$rho) > 0.1,
          fin(tr$D3$large$fe_beta_nofeedback), fin(tr$D3$large$fe_rho_nofeedback))

## ---- D4: an onboarding cohort that starts below its long-run level -----------------------------------
d <- D$D4; P4 <- P$D4
stopifnot(P4$start_share < 1, tr$D4$rho == P4$rho, all(d$y_stat_true > 0),
          max(abs(d$y_stat_true - d$y - P4$rho^(d$t - 1) * (1 - P4$start_share) * d$m_true)) < 1e-10)
f1 <- d$t == 1
stopifnot(cor(d$y[f1] - d$m_true[f1], d$alpha_true[f1]) < -0.3,                 # mean stationarity fails
          abs(cor(d$y_stat_true[f1] - d$m_true[f1], d$alpha_true[f1])) < 0.1)     # the twin starts at its level
# large-sample GMM (seed D4_large, N = 20,000): AB near rho, BB far from it, BB fine on the stationary twin
stopifnot(tr$D4$large$N >= 2e4, abs(tr$D4$large$ab - P4$rho) < 0.05, tr$D4$large$bb - P4$rho > 0.1,
          abs(tr$D4$large$bb_stat - P4$rho) < 0.03)
# the same large sample, ONE-step (fix round 2, 11.3b C1): with invalid level moments the one-step and the two-step
# system GMM have different pseudo-true values (about 0.817 against 0.806), so the one-step limit is its own number;
# one-step AB is consistent for rho like the two-step AB (both near 0.60)
stopifnot(fin(tr$D4$large$ab_onestep), fin(tr$D4$large$bb_onestep), abs(tr$D4$large$ab_onestep - P4$rho) < 0.05,
          tr$D4$large$bb_onestep - P4$rho > 0.1, abs(tr$D4$large$bb_onestep - tr$D4$large$bb) > 0.005,
          abs(tr$D4$large$bb_onestep - tr$D4$large$bb) < 0.05)

## ---- D5: the D1 design (legacy gen_dyn) at T_obs = 8, 16, 32 months; exact limits ----------------------
P5 <- P$D5
stopifnot(identical(P5$T_obs, c(8, 16, 32)), P5$R_mc >= 1000, P5$rho == tr$rho, P5$N == 500,
          identical(tr$D5$T, P5$T_obs - 1),
          max(abs(tr$D5$fe_plim - (P5$rho + nickell_bias(P5$rho, tr$D5$T)))) < 1e-15,
          max(abs(tr$D5$hk_plim - hk_correct(tr$D5$fe_plim, tr$D5$T))) < 1e-15,
          # Dhaene and Jochmans (2015) eqs. (2.2), (2.4): each piece of a cut weighted by |S|/T (fix round 2); the two cuts of
          # an odd T have the same limit (floor(T/2) b(floor(T/2)) + ceiling(T/2) b(ceiling(T/2)))/T
          max(abs(tr$D5$spj_plim - sapply(tr$D5$T, function(T) {
            b <- function(k) nickell_bias(P5$rho, k)
            lim_cut <- function(k) (k * b(k) + (T - k) * b(T - k)) / T
            P5$rho + 2 * b(T) - mean(lim_cut(unique(c(floor(T / 2), ceiling(T / 2)))))
          }))) < 1e-15,
          max(abs(round(tr$D5$spj_plim, 4) - c(0.6747, 0.7092, 0.7059))) < 1e-12,               # 11.3a E1: 0.6845, 0.7104, 0.7060 before
          all(abs(tr$D5$spj_plim - P5$rho) < abs(tr$D5$fe_plim - P5$rho)),
          all(abs(tr$D5$hk_plim - P5$rho) < abs(tr$D5$fe_plim - P5$rho)))

## ---- D1p: Anderson-Hsiao with the differenced instrument, Monte Carlo (fix round 2, 11.3b I2) ----------------
# The D1p design (legacy gen_dyn, N = 500, T_obs = 8, rho = 0.90), seed D1p_ah_mc, R = 1,000 samples; the stored
# (estimate, SE) of every sample regenerate from the seed (the first ten are re-drawn here).
Pa <- P$D1p_ah; ah <- readRDS("res_d1.rds")$D1p$ah_diff_mc
stopifnot(is.list(Pa), is.list(ah), Pa$N == 500, Pa$T_obs == 8, Pa$rho == tr$rho_persist, Pa$R_mc >= 1000, ah$R == Pa$R_mc, nrow(ah$draws) == ah$R,
          identical(colnames(ah$draws), c("est", "se")), all(is.finite(ah$draws)))
set.seed(P$seed$D1p_ah_mc)
chk <- t(vapply(1:10, function(r) {
  dd <- e1$gen_dyn(N = Pa$N, T = Pa$T_obs, rho = Pa$rho, burn = Pa$burn, sigma_a = Pa$sigma_a, sigma_e = Pa$sigma_e)
  a <- anderson_hsiao(dd$y, dd$id, dd$t, "diff"); c(est = a$est, se = a$se)
}, numeric(2)))
stopifnot(max(abs(chk - ah$draws[1:10, ])) == 0)
