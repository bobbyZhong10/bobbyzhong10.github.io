# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 11.
# Run from the chapter code directory after 01 to 08 (tests/run_tests.R does this).
# The first block is the chapter card's assertions, verbatim except the D2 line, which follows the
# controller's ruling 15 as corrected by ruling 41: the card's "pulled towards FE" is not what plm's system
# GMM does. With the standard one-step weight H = D D', GLS on the differenced equation equals within-groups
# exactly, so DIFFERENCE GMM with many instruments is pulled towards FE; SYSTEM GMM (block D2) is pulled
# mostly towards OLS in levels (a mix with about 0.13 weight on within-groups), so D2's non-instrumented
# counterpart ols_transformed is OLS on the level equation. The second block adds the claims the text makes
# beyond the card (and, in the D2 part, why J_p_full > 0.99 is no evidence: J = N when r >= N).
# ------------------------------------------------------------------
source("../../tools/R/chapter_kit.R")
source("dyn_fns.R")

## ---- chapter card (verbatim; D2 line per rulings 15 and 41) ------------------------------------------------
tr <- readRDS("truth11.rds"); e <- readRDS("res_extra.rds")
stopifnot(covers(e$ah$ci, tr$rho), e$ah$se > e$ab_onestep$se)
stopifnot(e$D2$J_p_full > 0.99, abs(e$D2$est_full - e$D2$ols_transformed) < abs(e$D2$est_collapsed - e$D2$ols_transformed),
          e$D2$n_inst_full > 3 * e$D2$n_inst_collapsed)                                     # 工具扩散把 GMM 拉向不用工具的估计量
stopifnot(all(abs(e$D5$spj - tr$rho) < abs(e$D5$fe - tr$rho)), all(diff(abs(e$D5$fe - tr$rho)) < 0))
stopifnot(!covers(e$D4$bb$ci, tr$D4$rho), covers(e$D4$ab$ci, tr$D4$rho))                    # 初始条件错时 BB 偏
stopifnot(!covers(e$D3$fe$ci, tr$D3$beta), covers(e$D3$ab$ci, tr$D3$beta))

## ---- D1: the inherited sample ------------------------------------------------------------------------
d <- readRDS("res_dynamic.rds")
# the refits in 02 reproduce the legacy two-step estimates; Windmeijer's correction raises the two-step SE
stopifnot(abs(e$ab_twostep$est - d$ab) < 1e-12, abs(e$bb_twostep$est - d$bb) < 1e-12,
          e$ab_twostep$se > e$ab_twostep$se_conv,
          abs(e$D1p$ab$est - d$ab9) < 1e-12, abs(e$D1p$ab$se_conv - d$se_ab9) < 1e-12, abs(e$D1p$bb$se_conv - d$se_bb9) < 1e-12)
# pooled OLS above rho, FE below; both bias corrections move FE towards rho
stopifnot(e$D1$pool$est > tr$rho, e$D1$fe$est < tr$rho, abs(e$D1$fe$est - d$fe) < 1e-12,
          e$D1$hk > e$D1$fe$est, e$D1$hk < tr$rho, abs(e$D1$spj - tr$rho) < abs(e$D1$fe$est - tr$rho))
# Anderson-Hsiao: the differenced instrument is weaker than the level instrument; one-step AB uses all lags
stopifnot(e$ah_diff$se > e$ah$se, e$ab_onestep$r == 21, covers(e$ab_onestep$ci, tr$rho), covers(e$ab_twostep$ci, tr$rho))
# rho near 1 (D1p): difference GMM is far less precise than system GMM (Windmeijer SEs too)
stopifnot(e$D1p$ab$se > 4 * e$D1p$bb$se, covers(e$D1p$bb$ci, tr$rho_persist))
# the legacy Nickell Monte Carlo (40 samples of 300 merchants) agrees with the exact limit
stopifnot(max(abs(d$nick["fe", ] - tr$D1T$fe_plim)) < 0.02)

# Anderson-Hsiao with the differenced instrument at rho = 0.90 (D1p), Monte Carlo of exactly that estimator and design
# (R = 1,000 samples, seed D1p_ah_mc; fix round 2, 11.3b I2): the stored summaries are the ones the draws give. The
# estimator is a just-identified IV estimator, which has no finite moments: its Monte Carlo MEAN is dominated by a few
# extreme draws (nothing is asserted about the mean or the SD); the median is near the truth, the interquartile range is
# wide (a weak instrument), and the 95% interval's coverage is a Monte Carlo proportion with standard error
# sqrt(p (1 - p)/R), about 0.008.
am <- e$D1p$ah_diff_mc; zc <- qnorm(0.975); rho9 <- tr$rho_persist; ad <- am$draws
lo <- ad[, "est"] - zc * ad[, "se"]; hi <- ad[, "est"] + zc * ad[, "se"]
stopifnot(is.list(am), am$R >= 1000, nrow(ad) == am$R, all(is.finite(ad)), all(ad[, "se"] > 0),
          abs(am$cover - mean(lo <= rho9 & rho9 <= hi)) < 1e-15, abs(am$cover_mcse - sqrt(am$cover * (1 - am$cover) / am$R)) < 1e-15,
          abs(am$share_ci_below_truth - mean(hi < rho9)) < 1e-15, abs(am$share_ci_above_truth - mean(lo > rho9)) < 1e-15,
          abs(am$cover + am$share_ci_below_truth + am$share_ci_above_truth - 1) < 1e-15,
          abs(am$median - median(ad[, "est"])) < 1e-15, abs(am$median_bias - (am$median - rho9)) < 1e-15,
          abs(am$mean - mean(ad[, "est"])) < 1e-15, abs(am$iqr - IQR(ad[, "est"])) < 1e-15,
          max(abs(c(am$q025, am$q975) - quantile(ad[, "est"], c(0.025, 0.975)))) < 1e-15,
          am$cover > 0.88, am$cover < 0.99, abs(am$median - rho9) < 0.1, am$iqr > 0.3,         # weak instrument: a wide spread
          am$q025 < rho9, am$q975 > rho9)

## ---- D2: J is mechanically near N when r > N; the collapsed set stays near rho -----------------------
stopifnot(e$D2$n_inst_full > e$D2$N, e$D2$n_inst_collapsed < e$D2$N,
          e$D2$ols_transformed > e$D2$est_full, e$D2$est_full > e$D2$est_collapsed, e$D2$fe < tr$D2$rho)
# J = N exactly when r >= N. The one-step weight matrix is the (generalised) inverse of sum_i g_i g_i', g_i the
# r-vector of unit i's moments Z_i' e_i; it has rank at most N, so J = 1_N' P 1_N, with P the projection on the
# span of the N moment vectors, and P is the identity as soon as the r x N matrix [g_1 ... g_N] has rank N
# (r >= N). So J = N whatever the data: the card line J_p_full > 0.99 is P(chi2_df > N) with df = r - 1 > N,
# a statement about the instrument count and not about the validity of the moments (chapter 7's many-moment
# problem). The power evidence is the misspecified twin below.
stopifnot(abs(e$D2$J_full - e$D2$N) < 1e-6)
pa <- e$D2$path
stopifnot(all(abs(pa$J[pa$r >= e$D2$N] - e$D2$N) < 1e-6), all(pa$J[pa$r < e$D2$N] < e$D2$N - 1e-3), any(pa$r >= e$D2$N))
stopifnot(e$D2$mc$J_full_maxdev_from_N < 1e-6, e$D2$mc$share_J_full_below_005 == 0)
# the Monte Carlo (R_mc samples of the same design): on average full is pulled towards OLS in levels
mc <- e$D2$mc
stopifnot(mc$R >= 500, mc$share_J_full_above_099 > 0.99, mc$mean_full - mc$mean_collapsed > 0.05,
          abs(mc$mean_collapsed - tr$D2$rho) < 0.02, abs(mc$mean_ols - tr$D2$pooled_plim) < 0.01, mc$share_closer > 0.95)
# the lag-depth path (figure): instrument counts rise with the lag depth; the last rows are the two fits above
stopifnot(all(diff(pa$r[!pa$collapse]) > 0), all(diff(pa$r[pa$collapse]) > 0),
          abs(pa$est[!pa$collapse][sum(!pa$collapse)] - e$D2$est_full) < 1e-12,
          abs(pa$est[pa$collapse][sum(pa$collapse)] - e$D2$est_collapsed) < 1e-12)

## ---- D2 misspecified twin: proliferation destroys the power of the Hansen J test -----------------------------
# MA(1) errors (gen_D2_ma, seed D2_mis, 200 samples): Y_{i,t-2} and the level moments are invalid. The full set
# (r > N) cannot reject (J = N); the collapsed set (r = 30) rejects almost always; the middle set (lags 2:4, r = 68 < N
# but close to it) is stored, with its size under the correct model (mc$share_J_mid_below_005).
ms <- e$D2$mis
stopifnot(ms$R >= 200, ms$ma_coef != 0, ms$r_full > e$D2$N, ms$r_mid < e$D2$N, ms$r_collapsed < ms$r_mid,
          ms$rej_full == 0, ms$J_full_maxdev_from_N < 1e-6, ms$rej_collapsed > 0.8,
          is.finite(ms$rej_mid), ms$rej_mid >= 0, ms$rej_mid <= 1, ms$rej_collapsed > ms$rej_mid,
          ms$rej_mid > e$D2$mc$share_J_mid_below_005)                          # the middle set does have some power
# Phase 4 fix (audit ch11-1): the rejection rates above use the J of the one-step fit (chapter Theorem 6). The
# Hansen J of Theorem 4 (iv) is J2, the J of the two-step fit: the same weight at the two-step estimate, which
# minimises the quadratic form, so J2 <= J in every sample and J2 <= N for the full set. One-step system GMM does
# not use the efficient weight, so its J over-rejects; J2 has about the nominal size for the collapsed set (within
# 3 Monte Carlo SEs of 5%) and still rejects the MA(1) twin almost always; the full and the middle set (r near or
# above N) never reject with J2.
stopifnot(mc$J2_minus_J_max <= 1e-8, mc$J2_full_max <= e$D2$N + 1e-6, mc$share_J2_full_below_005 == 0,
          mc$share_J2_full_above_099 > 0.99, abs(mc$share_J2_collapsed_below_005 - 0.05) < 3 * sqrt(0.05 * 0.95 / mc$R),
          mc$share_J_collapsed_below_005 > mc$share_J2_collapsed_below_005,
          ms$J2_minus_J_max <= 1e-8, ms$rej2_full == 0, ms$rej2_collapsed > 0.8, ms$rej2_mid <= ms$rej_mid,
          mc$share_J2_mid_below_005 <= mc$share_J_mid_below_005)

## ---- D3: FE is biased for rho too; its beta is near the large-sample FE limit ----------------------------
stopifnot(!covers(e$D3$fe_rho$ci, tr$D3$rho), abs(e$D3$fe$est - tr$D3$large$fe_beta) < 4 * e$D3$fe$se)

## ---- D4: here the Hansen test does detect the bad level moments (N = 2,000, r < N) -----------------------
stopifnot(e$D4$bb$J_p < 0.05, e$D4$bb$est > tr$D4$rho)

## ---- D5: the Monte Carlo agrees with the exact limits; HK also beats FE; FE's interval undercovers -------
m <- e$D5
mcse <- function(s) 4 * s / sqrt(m$R)
stopifnot(m$R >= 1000, all(abs(m$fe - tr$D5$fe_plim) < mcse(m$sd_fe)), all(abs(m$spj - tr$D5$spj_plim) < mcse(m$sd_spj) + 0.005),
          all(abs(m$hk - tr$rho) < abs(m$fe - tr$rho)), all(m$cover_fe < m$cover_spj), m$cover_fe[1] < 0.05)
