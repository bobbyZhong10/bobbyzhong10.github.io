# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of dyn_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this). Uses its own simulated data
# (seeds 9101 to 9108), not the chapter's rds files. Every estimator is checked against a package
# (plm, AER::ivreg, sandwich, stats::lm) or against a closed form.
# ------------------------------------------------------------------
source("dyn_fns.R")
suppressPackageStartupMessages(library(plm))
pfe <- function(dd, f = y ~ lag(y, 1)) plm(f, pdata.frame(dd, index = c("id", "t")), model = "within")

## ---- hk_correct: Hahn and Kuersteiner (2001 working paper) eq. (6), exactly --------------------
r <- c(0.3, 0.4246, 0.8); TT <- c(7, 15, 31)
stopifnot(max(abs(hk_correct(r, TT) - ((TT + 1) / TT * r + 1 / TT))) < 1e-15,
          max(abs(hk_correct(r, TT) - (r + (1 + r) / TT))) < 1e-15)

## ---- nickell_bias (T = regression periods) ------------------------------------------------------
# lit_software 13.3: rho = 0.5, T = 10: exact -0.162, first factor -(1 + rho)/(T - 1) = -0.167, -(1 + rho)/T = -0.150
stopifnot(abs(nickell_bias(0.5, 10) - (-0.162)) < 5e-4, abs(nickell_first_factor(0.5, 10) - (-1.5 / 9)) < 1e-15)
# T times the bias tends to -(1 + rho): the order-1/T term -(1 + rho)/T that hk_correct adds back
stopifnot(abs(5000 * nickell_bias(0.7, 5000) + 1.7) < 0.005, abs(5000 * nickell_first_factor(0.7, 5000) + 1.7) < 0.001)
# hk_correct applied to the exact probability limit leaves only an O(1/T^2) bias
for (rho in c(0.3, 0.7, 0.9)) {
  Tr <- c(7, 15, 31, 63)
  left <- hk_correct(rho + nickell_bias(rho, Tr), Tr) - rho
  stopifnot(all(abs(left) < abs(nickell_bias(rho, Tr))), all(diff(abs(left * Tr)) < 0))
}
# the formula against large simulated panels (stationary start, matrix-based FE written out here)
sim_fe <- function(N, Tobs, rho) {
  a <- rnorm(N); y <- matrix(0, N, Tobs + 50)
  y[, 1] <- a / (1 - rho) + rnorm(N, 0, 1 / sqrt(1 - rho^2))
  for (s in 2:ncol(y)) y[, s] <- rho * y[, s - 1] + a + rnorm(N)
  y <- y[, 51:(Tobs + 50)]; Y <- y[, -1]; L <- y[, -Tobs]
  Y <- Y - rowMeans(Y); L <- L - rowMeans(L)
  sum(L * Y) / sum(L^2)
}
set.seed(9101)
for (cfg in list(c(0.5, 5), c(0.7, 8), c(0.9, 11))) {
  stopifnot(abs(sim_fe(1e5, cfg[2], cfg[1]) - (cfg[1] + nickell_bias(cfg[1], cfg[2] - 1))) < 0.005)
}

## ---- fe_ar1 = plm within (coefficients and CR0 SEs), with and without X; any row order -----------
set.seed(9102)
P2 <- list(N = 300, T = 9, rho = 0.7, m_mean = 40, m_sd = 5, sigma_e = 2, burn = 50)
d2 <- gen_D2(P2)
f <- fe_ar1(d2$y, d2$id, d2$t); fp <- pfe(d2)
stopifnot(abs(f$est[["rho"]] - coef(fp)[[1]]) < 1e-12,
          abs(f$se[["rho"]] - sqrt(vcovHC(fp, method = "arellano", type = "HC0")[1, 1])) < 1e-12,
          f$T == P2$T - 1, f$G == P2$N, f$n == P2$N * (P2$T - 1))
o <- sample(nrow(d2))
f_o <- fe_ar1(d2$y[o], d2$id[o], d2$t[o])
stopifnot(max(abs(f_o$est - f$est)) < 1e-12, max(abs(f_o$se - f$se)) < 1e-12)
set.seed(9103)
P3 <- list(N = 400, T = 8, rho = 0.6, beta = 0.5, a_mean = 3, a_sd = 1, x_mean = 6, x_alpha = 0.4, feedback = -0.6,
           eps_x_sd = 0.5, sigma_e = 1, burn = 50)
d3 <- gen_D3(P3)
f3 <- fe_ar1(d3$y, d3$id, d3$t, X = cbind(x = d3$x)); fp3 <- pfe(d3, y ~ lag(y, 1) + x)
stopifnot(max(abs(f3$est - unname(coef(fp3)))) < 1e-12, identical(names(f3$est), c("rho", "x")),
          max(abs(f3$se - sqrt(diag(vcovHC(fp3, method = "arellano", type = "HC0"))))) < 1e-12)

## ---- pooled_ar1 and ols_diff_ar1 = lm with sandwich::vcovCL (HC0, no cluster adjustment) ---------
dd <- d2; dd$yl <- lag_within(dd$y, dd$id, dd$t, 1L); dd <- dd[!is.na(dd$yl), ]
for (tfe in c(FALSE, TRUE)) {
  lmf <- if (tfe) lm(y ~ yl + factor(t), dd) else lm(y ~ yl, dd)
  po <- pooled_ar1(d2$y, d2$id, d2$t, time_fe = tfe)
  stopifnot(abs(po$est - coef(lmf)[["yl"]]) < 1e-10, po$n == nrow(dd),
            abs(po$se - sqrt(sandwich::vcovCL(lmf, cluster = ~id, type = "HC0", cadjust = FALSE)["yl", "yl"])) < 1e-10)
}
dd <- d2; dd$dy <- dd$y - lag_within(dd$y, dd$id, dd$t, 1L); dd$dyl <- lag_within(dd$dy, dd$id, dd$t, 1L)
dd <- dd[!is.na(dd$dy) & !is.na(dd$dyl), ]
lmd <- lm(dy ~ dyl - 1, dd); od <- ols_diff_ar1(d2$y, d2$id, d2$t)
stopifnot(abs(od$est - coef(lmd)[[1]]) < 1e-12,
          abs(od$se - sqrt(sandwich::vcovCL(lmd, cluster = ~id, type = "HC0", cadjust = FALSE)[1, 1])) < 1e-12)

## ---- anderson_hsiao = AER::ivreg (CR0) and, for inst = "level", the just-identified pgmm --------------
dd <- d2; dd$dy <- dd$y - lag_within(dd$y, dd$id, dd$t, 1L); dd$dyl <- lag_within(dd$dy, dd$id, dd$t, 1L)
dd$z_level <- lag_within(dd$y, dd$id, dd$t, 2L); dd$z_diff <- lag_within(dd$dy, dd$id, dd$t, 2L)
for (k in c("level", "diff")) {
  ah <- anderson_hsiao(d2$y, d2$id, d2$t, inst = k)
  sub <- dd[stats::complete.cases(dd[, c("dy", "dyl", paste0("z_", k))]), ]
  iv <- AER::ivreg(stats::as.formula(paste0("dy ~ dyl - 1 | z_", k, " - 1")), data = sub)
  stopifnot(identical(names(ah), c("est", "se", "ci", "n", "inst")), ah$inst == k, ah$n == nrow(sub),
            abs(ah$est - coef(iv)[[1]]) < 1e-10,
            abs(ah$se - sqrt(sandwich::vcovCL(iv, cluster = ~id, type = "HC0", cadjust = FALSE)[1, 1])) < 1e-10,
            max(abs(ah$ci - ci95(ah$est, ah$se))) < 1e-15)
}
stopifnot(anderson_hsiao(d2$y, d2$id, d2$t)$n == P2$N * (P2$T - 2),                  # default inst = "level"
          anderson_hsiao(d2$y, d2$id, d2$t, "diff")$n == P2$N * (P2$T - 3))
pd2 <- pdata.frame(d2, index = c("id", "t"))
g_ah <- pgmm(y ~ lag(y, 1) | lag(y, 2:2), pd2, effect = "individual", model = "onestep", transformation = "d", collapse = TRUE)
ah <- anderson_hsiao(d2$y, d2$id, d2$t, "level")
stopifnot(n_instruments(g_ah) == 1L, abs(ah$est - coef(g_ah)[[1]]) < 1e-8, abs(ah$se - sqrt(vcovHC(g_ah)[1, 1])) < 1e-8)

## ---- spj_fe: Dhaene and Jochmans (2015) eqs. (2.2) and (2.4), the halves weighted by |S|/T -------------
# rho_hat_SPJ = 2 rho_hat_FE - bar_rho_{1/2}; bar_rho_S = sum over the pieces S of the cut (|S|/T) rho_hat_FE,S (2.2);
# bar_rho_{1/2} = the average of the cuts' bar_rho_S (2.4). Even T: one cut, two halves of T/2 periods, weights 1/2
# and 1/2 (the equal average). Odd T: two cuts, (floor(T/2), ceiling(T/2)) and (ceiling(T/2), floor(T/2)) periods,
# each cut's halves weighted by their shares of the T regression periods (not equally).
# (a) T_obs = 9: regression periods t = 2..9 (T = 8, even): halves t = 2..5 (months 1..5) and t = 6..9 (months 5..9)
full <- coef(pfe(d2))[[1]]; h1 <- coef(pfe(d2[d2$t <= 5, ]))[[1]]; h2 <- coef(pfe(d2[d2$t >= 5, ]))[[1]]
sp <- spj_parts(d2$y, d2$id, d2$t)
stopifnot(abs(spj_fe(d2$y, d2$id, d2$t) - (2 * full - (h1 + h2) / 2)) < 1e-12, sp$T == 8, nrow(sp$halves) == 1,
          abs(sp$fe - full) < 1e-12, max(abs(sp$halves - c(h1, h2))) < 1e-12, is.numeric(spj_fe(d2$y, d2$id, d2$t)),
          length(spj_fe(d2$y, d2$id, d2$t)) == 1,
          max(abs(sp$weights - 0.5)) < 1e-15, identical(dim(sp$weights), dim(sp$halves)))      # even T: equal weights
# (b) T_obs = 8: T = 7 (odd): the cuts (3, 4) and (4, 3) regression periods, weights 3/7, 4/7 and 4/7, 3/7
d8 <- d2[d2$t <= 8, ]
f8 <- coef(pfe(d8))[[1]]
hh <- c(coef(pfe(d8[d8$t <= 4, ]))[[1]], coef(pfe(d8[d8$t >= 4, ]))[[1]], coef(pfe(d8[d8$t <= 5, ]))[[1]],
        coef(pfe(d8[d8$t >= 5, ]))[[1]])                         # cut 1: S1, S2 (3 and 4 periods); cut 2: S1, S2 (4 and 3)
w8 <- c(3, 4, 4, 3) / 7
dj8 <- 2 * f8 - ((w8[1] * hh[1] + w8[2] * hh[2]) + (w8[3] * hh[3] + w8[4] * hh[4])) / 2
sp8 <- spj_parts(d8$y, d8$id, d8$t)
stopifnot(abs(spj_fe(d8$y, d8$id, d8$t) - dj8) < 1e-12, nrow(sp8$halves) == 2, sp8$T == 7,
          max(abs(t(sp8$halves) - hh)) < 1e-12, max(abs(t(sp8$weights) - w8)) < 1e-15,
          max(abs(rowSums(sp8$weights) - 1)) < 1e-15, abs(sp8$half_mean - (2 * f8 - dj8)) < 1e-12)
# (c) a small deterministic panel with odd T, against exact arithmetic done by hand (rational numbers).
# 3 merchants, 6 months, T = 5 regression periods (t = 2..6); y_i1..y_i6 =
#     i = 1: 1 3 2 5 4 6;   i = 2: 4 2 3 1 4 2;   i = 3: 0 2 5 3 6 5.
# rho_hat_FE on a piece = sum_i sum_t (L - Lbar_i)(Y - Ybar_i) / sum_i sum_t (L - Lbar_i)^2, L = lag, means over the
# piece's regression periods. Exact fractions (rational arithmetic outside R, numerator/denominator in brackets):
#   all five periods         7/66     (21/5 over 198/5)
#   cut 1, S1 = {2,3}        1/6      (1 over 6)            S2 = {4,5,6}   -1/2     (-7 over 14)
#   cut 2, S1 = {2,3,4}     -1/50     (-1/3 over 50/3)      S2 = {5,6}    -11/19    (-11/2 over 19/2)
# DJ (2.2): bar_1 = (2/5)(1/6) + (3/5)(-1/2) = -7/30;  bar_2 = (3/5)(-1/50) + (2/5)(-11/19) = -1157/4750;
# DJ (2.4): bar_{1/2} = (bar_1 + bar_2)/2 = -1699/7125;  rho_hat_SPJ = 2 (7/66) - bar_{1/2} = 7/33 + 1699/7125 = 35314/78375.
# (The equal-weight average of the four halves would give -2657/11400 and 0.44519: not the same estimator.)
yy <- rbind(c(1, 3, 2, 5, 4, 6), c(4, 2, 3, 1, 4, 2), c(0, 2, 5, 3, 6, 5))
hp <- data.frame(id = rep(1:3, each = 6), t = rep(1:6, times = 3), y = as.vector(t(yy)))
sph <- spj_parts(hp$y, hp$id, hp$t)
fe_hand <- c(7 / 66, 1 / 6, -1 / 2, -1 / 50, -11 / 19)
stopifnot(sph$T == 5, nrow(sph$halves) == 2,
          abs(sph$fe - fe_hand[1]) < 1e-14, max(abs(sph$halves - rbind(fe_hand[2:3], fe_hand[4:5]))) < 1e-14,
          # the same five estimates from plm's within estimator on the sub-panels
          abs(coef(pfe(hp))[[1]] - fe_hand[1]) < 1e-12, abs(coef(pfe(hp[hp$t <= 3, ]))[[1]] - fe_hand[2]) < 1e-12,
          abs(coef(pfe(hp[hp$t >= 3, ]))[[1]] - fe_hand[3]) < 1e-12, abs(coef(pfe(hp[hp$t <= 4, ]))[[1]] - fe_hand[4]) < 1e-12,
          abs(coef(pfe(hp[hp$t >= 4, ]))[[1]] - fe_hand[5]) < 1e-12,
          # weights |S|/T: cut 1 has 2 and 3 periods, cut 2 has 3 and 2
          max(abs(sph$weights - rbind(c(2 / 5, 3 / 5), c(3 / 5, 2 / 5)))) < 1e-15,
          abs(sph$half_mean - (-1699 / 7125)) < 1e-14, abs(sph$spj - 35314 / 78375) < 1e-14,
          abs(spj_fe(hp$y, hp$id, hp$t) - 35314 / 78375) < 1e-14,
          abs(spj_fe(hp$y, hp$id, hp$t) - (2 * fe_hand[1] - mean(fe_hand[2:5]))) > 1e-3)      # not the equal-weight average
# spj removes the order-1/T bias in a large panel (rho = 0.7, T = 8): |bias| below a fifth of FE's
set.seed(9104)
big <- gen_D2(modifyList(P2, list(N = 2e4)))
bs <- spj_parts(big$y, big$id, big$t)
stopifnot(abs(bs$spj - P2$rho) < abs(bs$fe - P2$rho) / 5, abs(bs$fe - (P2$rho + nickell_bias(P2$rho, 8))) < 0.01)
# and for odd T (T_obs = 8, T = 7): the DJ-weighted jackknife is centred on its exact limit
# rho + 2 b(T) - (floor(T/2) b(floor(T/2)) + ceiling(T/2) b(ceiling(T/2)))/T, b = nickell_bias (each half starts stationary)
set.seed(9109)
big8 <- gen_D2(modifyList(P2, list(N = 4e4, T = 8)))
bs8 <- spj_parts(big8$y, big8$id, big8$t)
b7 <- function(k) nickell_bias(P2$rho, k)
lim7 <- P2$rho + 2 * b7(7) - (3 * b7(3) + 4 * b7(4)) / 7
lim7_eq <- P2$rho + 2 * b7(7) - (b7(3) + b7(4)) / 2            # the equal-weight average of the halves (not DJ's)
stopifnot(bs8$T == 7, abs(bs8$spj - lim7) < 0.015, abs(lim7 - lim7_eq) > 5e-3, abs(bs8$spj - lim7) < abs(bs8$spj - lim7_eq))

## ---- n_instruments: the instrument counts r of difference and system GMM ------------------------------
Tobs <- P2$T; dg <- d2[d2$id <= 150, ]; pdg <- pdata.frame(dg, index = c("id", "t"))
cnt <- function(tf, ef, cl) n_instruments(suppressWarnings(pgmm(y ~ lag(y, 1) | lag(y, 2:99), pdg, effect = ef,
                                                                 model = "onestep", transformation = tf, collapse = cl)))
stopifnot(cnt("d", "individual", FALSE) == (Tobs - 1) * (Tobs - 2) / 2,                  # legacy 4.4: (T-1)(T-2)/2
          cnt("d", "individual", TRUE) == Tobs - 2,
          cnt("ld", "individual", FALSE) == (Tobs - 1) * (Tobs - 2) / 2 + (Tobs - 2),     # + one Delta Y per level equation
          cnt("ld", "twoways", FALSE) == (Tobs - 1) * (Tobs - 2) / 2 + (Tobs - 2) + (Tobs - 1),   # + intercept and month dummies
          cnt("ld", "twoways", TRUE) == (Tobs - 2) + 1 + (Tobs - 1),
          is.integer(cnt("d", "individual", TRUE)))
stopifnot(inherits(tryCatch(n_instruments(lm(y ~ t, d2)), error = function(e) e), "error"))

## ---- gmm_summary: robust (Windmeijer) SE, conventional SE, Hansen J, m1 and m2 = plm's own -------------
g2 <- pgmm(y ~ lag(y, 1) | lag(y, 2:99), pd2, effect = "individual", model = "twosteps", transformation = "d")
s2 <- gmm_summary(g2)
sm <- summary(g2, robust = TRUE)
stopifnot(abs(s2$est - coef(g2)[[1]]) < 1e-15, abs(s2$se - sqrt(vcovHC(g2)[1, 1])) < 1e-15,
          abs(s2$se - sm$coefficients[1, 2]) < 1e-12, abs(s2$se_conv - sqrt(vcov(g2)[1, 1])) < 1e-15,
          abs(s2$J - sargan(g2)$statistic) < 1e-10, abs(s2$J_p - sm$sargan$p.value) < 1e-12,
          s2$J_df == s2$r - 1, abs(s2$m2_p - sm$m2$p.value) < 1e-12, abs(s2$m1_p - sm$m1$p.value) < 1e-12,
          s2$r == n_instruments(g2), s2$N == P2$N, max(abs(s2$ci - ci95(s2$est, s2$se))) < 1e-15)
g1 <- pgmm(y ~ lag(y, 1) | lag(y, 2:99), pd2, effect = "individual", model = "onestep", transformation = "d")
s1 <- gmm_summary(g1)
stopifnot(is.na(s1$se_conv), abs(s1$se - sqrt(vcovHC(g1)[1, 1])) < 1e-15)
# plm's one-step Sargan statistic is unscaled (z' A1 z): it is not reported (fix round 1)
stopifnot(!("sargan_onestep_p" %in% names(s1)), !("sargan_onestep_p" %in% names(s2)))

## ---- generators: the stated equations hold exactly; the stationary start is stationary --------------
set.seed(9105)
b2 <- gen_D2(modifyList(P2, list(N = 2e4)))
dev <- b2$y - b2$m_true; v0 <- P2$sigma_e^2 / (1 - P2$rho^2)
stopifnot(abs(mean(dev)) < 0.05, abs(var(dev[b2$t == 1]) / v0 - 1) < 0.05, abs(var(dev[b2$t == P2$T]) / v0 - 1) < 0.05,
          max(abs(b2$alpha_true - (1 - P2$rho) * b2$m_true)) < 1e-12)
e2 <- b2$y - P2$rho * lag_within(b2$y, b2$id, b2$t) - b2$alpha_true
stopifnot(abs(sd(e2, na.rm = TRUE) - P2$sigma_e) < 0.02)
L <- function(v, k = 1L) lag_within(v, d3$id, d3$t, k)
r3 <- d3$y - P3$rho * L(d3$y) - P3$beta * d3$x - d3$alpha_true - d3$eps_true
ex3 <- d3$x - P3$x_mean - P3$x_alpha * d3$alpha_true - P3$feedback * L(d3$eps_true)     # the promotion shock eps_x
stopifnot(max(abs(r3), na.rm = TRUE) < 1e-10, abs(sd(ex3, na.rm = TRUE) - P3$eps_x_sd) < 0.02,
          abs(cor(ex3, L(d3$eps_true), use = "complete.obs")) < 0.05)
# the misspecified twin of D2: eps_it = e_it + ma_coef e_{i,t-1} has variance sigma_e^2 (1 + ma_coef^2), first
# autocovariance ma_coef sigma_e^2 and no second autocovariance (so Y_{i,t-2} is invalid, Y_{i,t-3} is not)
set.seed(9107)
Pma <- modifyList(P2, list(N = 2e4, ma_coef = 0.6))
bm <- gen_D2_ma(Pma)
em <- bm$y - Pma$rho * lag_within(bm$y, bm$id, bm$t) - bm$alpha_true
acv <- function(k) mean(em * lag_within(em, bm$id, bm$t, k), na.rm = TRUE)
stopifnot(abs(acv(0) / (Pma$sigma_e^2 * (1 + Pma$ma_coef^2)) - 1) < 0.02,
          abs(acv(1) / (Pma$ma_coef * Pma$sigma_e^2) - 1) < 0.03, abs(acv(2)) < 0.03 * Pma$sigma_e^2,
          max(abs(bm$alpha_true - (1 - Pma$rho) * bm$m_true)) < 1e-12,
          nrow(bm) == Pma$N * Pma$T, identical(names(bm), names(gen_D2(modifyList(P2, list(N = 5))))))
# with ma_coef = 0 the twin is an AR(1) panel again: the same second moments as gen_D2
set.seed(9108); b0 <- gen_D2_ma(modifyList(Pma, list(ma_coef = 0)))
e0 <- b0$y - Pma$rho * lag_within(b0$y, b0$id, b0$t) - b0$alpha_true
stopifnot(abs(sd(e0, na.rm = TRUE) / Pma$sigma_e - 1) < 0.02,
          abs(mean(e0 * lag_within(e0, b0$id, b0$t, 1L), na.rm = TRUE)) < 0.03 * Pma$sigma_e^2)
set.seed(9106)
P4 <- list(N = 2e4, T = 8, rho = 0.6, m_mean = 15, m_sd = 3, start_share = 0.8, sigma_e = 1)
d4 <- gen_D4(P4)
stopifnot(max(abs(d4$y_stat_true - d4$y - P4$rho^(d4$t - 1) * (1 - P4$start_share) * d4$m_true)) < 1e-10,
          max(abs((d4$y - P4$rho * lag_within(d4$y, d4$id, d4$t)) -
                  (d4$y_stat_true - P4$rho * lag_within(d4$y_stat_true, d4$id, d4$t))), na.rm = TRUE) < 1e-10)
m1st <- d4$t == 1
stopifnot(cor(d4$y[m1st] - d4$m_true[m1st], d4$alpha_true[m1st]) < -0.35,                        # mean stationarity fails
          abs(cor(d4$y_stat_true[m1st] - d4$m_true[m1st], d4$alpha_true[m1st])) < 0.03)            # and holds for the twin

## ---- ci95 and lag_within --------------------------------------------------------------------------------
stopifnot(isTRUE(all.equal(ci95(1, 0.5), c(1 - qnorm(0.975) * 0.5, 1 + qnorm(0.975) * 0.5))))
stopifnot(identical(lag_within(c(1, 2, 3, 10, 20), c(1, 1, 1, 2, 2), c(1, 2, 4, 1, 2)), c(NA, 1, NA, NA, 10)))
