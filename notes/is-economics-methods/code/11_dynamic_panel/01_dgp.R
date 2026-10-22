# ============================================================================
# 01_dgp.R  --  Nimbus marketplace merchant panel, DYNAMIC facet (chapter 11)
# ----------------------------------------------------------------------------
# "Nimbus" is a fictional marketplace. We observe merchants i over months t.
#
#  DYNAMIC facet -- Nickell bias and Arellano-Bond / Blundell-Bond.
#      y_it = rho * y_i,t-1 + alpha_i + eps_it   (persistent GMV).
#      Pooled OLS of y on its lag is biased UP (alpha_i sits in both y and its
#      lag); the within/FE estimator is biased DOWN by the Nickell O(1/T) bias
#      (demeaning correlates the transformed lag with the transformed error);
#      the truth is bracketed between them, and difference/system GMM recover it.
#
# All TRUE parameters known. Seeds fixed (808 and 7372, set inside gen_dyn).
# Scripts use relative paths; run with the working directory set to this
# folder. (The static facet, FE vs RE and the Hausman test, lives in
# code/10_linear_panel_fe_re/.)
# ============================================================================

## ---- dynamic panel: AR(1) with merchant effects ----------------------------
# stationary start; burn-in to erase initial-condition transient
gen_dyn <- function(N, T, rho, burn = 50, sigma_a = 1, sigma_e = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha <- rnorm(N, 0, sigma_a)
  Ttot  <- T + burn
  y <- matrix(0, N, Ttot)
  y[,1] <- alpha/(1-rho) + rnorm(N, 0, sigma_e/sqrt(1-rho^2))  # stationary draw
  for (s in 2:Ttot) y[,s] <- rho*y[,s-1] + alpha + rnorm(N, 0, sigma_e)
  y <- y[, (burn+1):Ttot]                              # keep last T
  data.frame(id = rep(1:N, each = T), t = rep(1:T, times = N),
             y = as.vector(t(y)))
}
dyn_main <- gen_dyn(N = 500, T = 8, rho = 0.70, seed = 808)    # main
dyn_persist <- gen_dyn(N = 500, T = 8, rho = 0.90, seed = 7372) # rho near 1

saveRDS(list(dyn_main = dyn_main, dyn_persist = dyn_persist,
             gen_dyn = gen_dyn,
             truth = list(rho = 0.70, rho_persist = 0.90)),
        "nimbus.rds")

cat(sprintf("dynamic: N=500 T=8, rho=0.70 (main), rho=0.90 (persistent)\n"))
cat(sprintf("  Nickell leading bias -(1+rho)/(T-1) at rho=.7,T=8 = %.3f\n",
            -(1+0.7)/(8-1)))

# ============================================================================
# Lines 1 to 42 above are the inherited legacy blocks D1 (dyn_main, gen_dyn(seed = 808)) and D1p
# (dyn_persist, seed 7372), byte-identical (Phase 3: inherited blocks first and unchanged;
# tests/test_dgp.R checks the md5 of those lines). The inherited Nickell Monte Carlo D1T
# (set.seed(303)) is in 03_dynamic_nickell.R, also unchanged. Everything below is new in v3.
# ----------------------------------------------------------------------------
# 01_dgp.R (v3) -- the Nimbus simulation world for chapter 11 (dynamic panels).
#
# Nimbus is the fictional marketplace of chapter 10 (one world); unit i is a merchant, period t a
# month. The inherited blocks keep the legacy standardised units (alpha_i ~ N(0, 1), eps ~ N(0, 1));
# the new blocks measure GMV Y_it in thousands of USD a month and promotion spend X_it in hundreds
# of USD a month. Blocks:
#   D1   inherited: 500 merchants, 8 months, rho = 0.70, stationary start             nimbus.rds
#   D1p  inherited: the same design with rho = 0.90                                    nimbus.rds
#   D1T  inherited: Nickell Monte Carlo over T_obs = 4, 8, 16, 32 months (in 03)
#   D2   instrument proliferation: 70 key accounts over 16 months, rho = 0.70; system GMM
#        with every available lag (134 instruments, more than merchants) against the
#        collapsed set (30 instruments)                                                 nimbus_D2.rds
#        D2 also has a misspecified twin (MA(1) errors, so Y_{i,t-2} is an invalid instrument), a
#        200-sample Monte Carlo in 04: the evidence that proliferation destroys the power of the
#        Hansen J test (D2 itself: J = N exactly whenever r >= N)
#   D3   predetermined regressor with feedback: after a weak month a merchant raises
#        next month's promotion spend; Y_it = rho Y_{i,t-1} + beta X_it + alpha_i + eps_it   nimbus_D3.rds
#   D4   non-stationary initial condition: an onboarding cohort that starts at 80% of its
#        long-run GMV (mean stationarity fails), and its stationary twin                nimbus_D4.rds
#   D5   bias-corrected FE: the D1 design (the legacy gen_dyn) at T_obs = 8, 16, 32 months,
#        a 1,000-sample Monte Carlo per length in 07_bias_correction_mc.R (data not saved)
#   D1p Anderson-Hsiao Monte Carlo (fix round 2): the D1p design (gen_dyn, N = 500, T_obs = 8, rho = 0.90), 1,000
#        samples, Anderson-Hsiao with the differenced instrument, in 02_d1_more.R (data not saved)
#   truth11.rds   list(rho, rho_persist, D1, D1p, D1T, D2, D3, D4, D5, params): truth, exact and
#                 closed-form probability limits, large-sample checks, parameters and seeds.
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 1.8, 1.9): Y_it, X_it, alpha_i, eps_it; rho the AR coefficient, beta the
# coefficient on X_it; T_obs months observed per merchant and T = T_obs - 1 regression periods
# (see dyn_fns.R: nickell_bias and hk_correct take T).
#
# Seeds (every block its own set.seed; D1 and D1p keep the legacy 808 and 7372 set inside gen_dyn,
# D1T keeps set.seed(303) in 03): D2 1102, D3 1103, D4 1104; Monte Carlo D2 1112 (04), misspecified
# MA(1) twin of D2 1132 (04, fix round 1), D5 1151, 1152, 1153 (07, one per T_obs); Anderson-Hsiao Monte
# Carlo at rho = 0.90 D1p_ah_mc 1161 (02, fix round 2; written down before the Monte Carlo was first run);
# large-sample checks D2 1122, D3 1113, D4 1114.
# Seeds were fixed before any result was seen and were not tuned afterwards.
# Run with the working directory set to this folder (00_run_all.R does this), like the legacy block.
# ============================================================================
source("dyn_fns.R")
suppressPackageStartupMessages(library(plm))

seed <- list(D1 = 808, D1p = 7372, D1T = 303, D2 = 1102, D3 = 1103, D4 = 1104, D2_mc = 1112,
             D5_T8 = 1151, D5_T16 = 1152, D5_T32 = 1153, D2_large = 1122, D3_large = 1113, D4_large = 1114,
             D2_mis = 1132, D1p_ah_mc = 1161)

## =================================================================
## D1, D1p, D1T: truth of the inherited blocks (closed form)
## =================================================================
## gen_dyn: alpha_i ~ N(0, 1), long-run level alpha_i / (1 - rho), eps ~ N(0, 1), stationary start and 50
## burn-in months; T_obs = 8 months, so T = 7 regression periods. Exact Nickell limit rho + nickell_bias(rho, 7);
## its first factor -(1 + rho)/(T - 1); the order-1/T bias -(1 + rho)/T, which is what the legacy
## "nickell_approx" row prints (it divides by T_obs - 1 = T). Pooled OLS of Y_it on (1, Y_{i,t-1}):
## (rho s2w + s2m)/(s2w + s2m), s2w = 1/(1 - rho^2) the within variance, s2m = 1/(1 - rho)^2 that of the level.
pool_plim <- function(rho, s2m, s2w) (rho * s2w + s2m) / (s2w + s2m)
d1_truth <- function(rho, Tobs) {
  TT <- Tobs - 1
  list(rho = rho, T_obs = Tobs, T = TT, fe_plim = rho + nickell_bias(rho, TT),
       first_factor_plim = rho + nickell_first_factor(rho, TT), order1_plim = rho - (1 + rho) / TT,
       pooled_plim = pool_plim(rho, 1 / (1 - rho)^2, 1 / (1 - rho^2)))
}
truth_D1 <- d1_truth(0.70, 8)
truth_D1p <- d1_truth(0.90, 8)
Ts_leg <- c(4, 8, 16, 32)
truth_D1T <- list(rho = 0.70, T_obs = Ts_leg, T = Ts_leg - 1, fe_plim = 0.70 + nickell_bias(0.70, Ts_leg - 1),
                  first_factor_plim = 0.70 + nickell_first_factor(0.70, Ts_leg - 1),
                  legacy_approx = 0.70 - (1 + 0.70) / (Ts_leg - 1))

## =================================================================
## D2: instrument proliferation (70 key accounts, 16 months)
## =================================================================
## Nimbus's 70 largest merchants (key accounts), 16 months. Long-run monthly GMV m_i ~ N(40, 5^2) thousand
## USD, alpha_i = (1 - rho) m_i, rho = 0.70, shocks eps ~ N(0, 2^2) (2,000 USD), stationary start and 50
## burn-in months. System GMM (Blundell-Bond, month effects) with every available lag has
## (T-1)(T-2)/2 + (T-2) + (T-1) = 105 + 14 + 15 = 134 instruments for 70 merchants; the collapsed set
## has 14 + 1 + 15 = 30. Probability limits (N -> infinity): pooled OLS in levels (with or without month
## dummies, the panel is stationary) (rho s2w + s2m)/(s2w + s2m) with s2w = sigma_e^2/(1 - rho^2) and
## s2m = m_sd^2; FE rho + nickell_bias(rho, 15); OLS in first differences -(1 - rho)/2.
## lags: the GMM lag ranges 2:L of the lag-depth path in 04 (99 = every available lag).
## mid_lag: the middle instrument set uses the lags 2:mid_lag (r = 68 < N = 70 with month effects).
P_D2 <- list(N = 70, T = 16, rho = 0.70, m_mean = 40, m_sd = 5, sigma_e = 2, burn = 50, R_mc = 500,
             lags = c(2, 3, 4, 6, 9, 99), mid_lag = 4)
## Misspecified twin of D2 (fix round 1): the same design with MA(1) errors eps_it = e_it + 0.6 e_{i,t-1}, so
## that Y_{i,t-2} is an invalid instrument and a Hansen J test with power must reject (generator gen_D2_ma; the
## 200-sample Monte Carlo is in 04, seed D2_mis). ma_coef 0.6 was chosen from an exploratory run with another
## seed (4242, 60 samples), before the block seed was fixed; the seed was not changed afterwards (ruling 25).
P_D2_ma <- modifyList(P_D2, list(ma_coef = 0.6, R_mc = 200))
## Anderson-Hsiao Monte Carlo at rho = 0.90 (fix round 2, 11.3b I2): exactly the D1p design (the legacy gen_dyn with
## its arguments N = 500, T = T_obs = 8, rho = 0.90, sigma_a = sigma_e = 1, 50 burn-in months) and exactly one
## estimator, anderson_hsiao(inst = "diff"); R_mc = 1,000 samples, seed D1p_ah_mc (02_d1_more.R).
P_D1p_ah <- list(N = 500, T_obs = 8, rho = 0.90, sigma_a = 1, sigma_e = 1, burn = 50, R_mc = 1000)
set.seed(seed$D2)
dD2 <- gen_D2(P_D2)
sw2 <- P_D2$sigma_e^2 / (1 - P_D2$rho^2)
T2 <- P_D2$T
truth_D2 <- list(rho = P_D2$rho, pooled_plim = pool_plim(P_D2$rho, P_D2$m_sd^2, sw2),
                 fe_plim = P_D2$rho + nickell_bias(P_D2$rho, T2 - 1), ols_diff_plim = -(1 - P_D2$rho) / 2,
                 r_full = (T2 - 1) * (T2 - 2) / 2 + (T2 - 2) + (T2 - 1), r_collapsed = (T2 - 2) + 1 + (T2 - 1))
set.seed(seed$D2_large)
big <- gen_D2(modifyList(P_D2, list(N = 5e4)))
truth_D2$large <- list(N = 5e4, pooled = pooled_ar1(big$y, big$id, big$t, time_fe = TRUE)$est,
                       fe = unname(fe_ar1(big$y, big$id, big$t)$est["rho"]), ols_diff = ols_diff_ar1(big$y, big$id, big$t)$est)
rm(big)

## =================================================================
## D3: lagged GMV and a predetermined promotion budget with feedback
## =================================================================
## 4,000 merchants, 8 months. alpha_i ~ N(3, 1). Promotion spend X_it = 6 + 0.4 alpha_i - 0.6 eps_{i,t-1} + eps_x_it,
## eps_x ~ N(0, 0.5^2) (the promotion shock, a code symbol only; about 720 USD a month; after a month 1,000 USD below expectation the merchant spends 60 USD
## more on promotion). GMV Y_it = 0.6 Y_{i,t-1} + 0.5 X_it + alpha_i + eps_it, eps ~ N(0, 1): an extra 100 USD of
## promotion adds 500 USD of GMV this month and 0.5/(1 - 0.6) = 1,250 USD in the long run; long-run mean GMV
## (3 + 0.5 * 7.2)/0.4 = 16,500 USD. X_it is sequentially exogenous: E[eps_it | Y_{i,t-1}, ..., X_it, ..., alpha_i] = 0,
## but X_{i,t+1} depends on eps_it, so FE is inconsistent; difference GMM with Y_{i,t-2}, ... and X_{i,t-1}, ...
## as instruments is consistent. Large-sample FE limits (N = 1e5, seed D3_large) with the feedback and with
## the same draws and feedback = 0 (the FE bias that the lagged dependent variable alone produces).
P_D3 <- list(N = 4000, T = 8, rho = 0.6, beta = 0.5, a_mean = 3, a_sd = 1, x_mean = 6, x_alpha = 0.4, feedback = -0.6,
             eps_x_sd = 0.5, sigma_e = 1, burn = 50)
set.seed(seed$D3)
dD3 <- gen_D3(P_D3)
truth_D3 <- list(beta = P_D3$beta, rho = P_D3$rho, long_run = P_D3$beta / (1 - P_D3$rho),
                 mean_x = P_D3$x_mean + P_D3$x_alpha * P_D3$a_mean,
                 mean_y = (P_D3$a_mean + P_D3$beta * (P_D3$x_mean + P_D3$x_alpha * P_D3$a_mean)) / (1 - P_D3$rho))
fe_large <- function(fb) {
  set.seed(seed$D3_large)
  b <- gen_D3(modifyList(P_D3, list(N = 1e5, feedback = fb)))
  unname(fe_ar1(b$y, b$id, b$t, X = cbind(x = b$x))$est)
}
f_fb <- fe_large(P_D3$feedback); f_0 <- fe_large(0)
truth_D3$large <- list(N = 1e5, fe_rho = f_fb[1], fe_beta = f_fb[2], fe_rho_nofeedback = f_0[1], fe_beta_nofeedback = f_0[2])

## =================================================================
## D4: an onboarding cohort that starts below its long-run level
## =================================================================
## 2,000 merchants who joined Nimbus in the same month, 8 months. Long-run monthly GMV m_i ~ N(15, 3^2) thousand
## USD, alpha_i = (1 - rho) m_i, rho = 0.6, eps ~ N(0, 1). In the first month a merchant reaches a share start_share = 80% of its
## long-run level: Y_i1 = 0.8 m_i + w_i1, w_i1 ~ N(0, 1/(1 - rho^2)). The deviation Y_i1 - m_i = -0.2 m_i + w_i1 is
## correlated with alpha_i (corr_dev_alpha below), so the Blundell-Bond level moments E[Delta Y_{i,t-1}(alpha_i + eps_it)]
## fail; the Arellano-Bond moments do not need the initial condition. y_stat_true: the same draws with start_share = 1.
## Large-sample check (N = 20,000, seed D4_large): two-step AB and BB with month effects, BB on the twin; fix round 2
## (11.3b C1): the ONE-step AB and BB on the same sample (no new draw), because with invalid level moments the one-step
## and the two-step system GMM estimators have different pseudo-true values (the one-step 0.817, the two-step 0.806).
P_D4 <- list(N = 2000, T = 8, rho = 0.6, m_mean = 15, m_sd = 3, start_share = 0.8, sigma_e = 1)
set.seed(seed$D4)
dD4 <- gen_D4(P_D4)
sw4 <- P_D4$sigma_e^2 / (1 - P_D4$rho^2)
truth_D4 <- list(rho = P_D4$rho, start_share = P_D4$start_share,
                 corr_dev_alpha = -(1 - P_D4$start_share) * P_D4$m_sd / sqrt((1 - P_D4$start_share)^2 * P_D4$m_sd^2 + sw4))
set.seed(seed$D4_large)
big <- pdata.frame(gen_D4(modifyList(P_D4, list(N = 2e4))), index = c("id", "t"))
g4 <- function(f, tf, model = "twosteps") coef(pgmm(f, big, effect = "twoways", model = model, transformation = tf))[[1]]
truth_D4$large <- list(N = 2e4, ab = g4(y ~ lag(y, 1) | lag(y, 2:99), "d"), bb = g4(y ~ lag(y, 1) | lag(y, 2:99), "ld"),
                       bb_stat = g4(y_stat_true ~ lag(y_stat_true, 1) | lag(y_stat_true, 2:99), "ld"),
                       ab_onestep = g4(y ~ lag(y, 1) | lag(y, 2:99), "d", "onestep"),
                       bb_onestep = g4(y ~ lag(y, 1) | lag(y, 2:99), "ld", "onestep"))
rm(big)

## =================================================================
## D5: bias-corrected FE at T_obs = 8, 16, 32 (the D1 design)
## =================================================================
## gen_dyn(N = 500, T = T_obs, rho = 0.70) of the legacy block, R_mc = 1000 samples per length (07). Exact limits:
## FE rho + b(T); HK hk_correct(FE limit, T) (linear in rho_hat_FE); half-panel jackknife with the halves weighted by
## their shares |S|/T of the regression periods (Dhaene and Jochmans 2015 eqs. 2.2 and 2.4; fix round 2): a cut (k, T - k)
## has the limit (k b(k) + (T - k) b(T - k))/T, b = nickell_bias (each half starts from a stationary month); for even T
## k = T/2 gives b(T/2); for odd T the two cuts (floor, ceiling) and (ceiling, floor) have the same limit
## (floor(T/2) b(floor(T/2)) + ceiling(T/2) b(ceiling(T/2)))/T, so the SPJ limit is rho + 2 b(T) minus that. (Before the fix
## the four halves had equal weight, rho + 2 b(T) - (b(floor) + b(ceiling))/2: 0.6845, 0.7104, 0.7060 for T = 7, 15, 31.)
P_D5 <- list(N = 500, T_obs = c(8, 16, 32), rho = 0.70, sigma_a = 1, sigma_e = 1, burn = 50, R_mc = 1000)
TT5 <- P_D5$T_obs - 1
b5 <- function(TT) nickell_bias(P_D5$rho, TT)
fe5 <- P_D5$rho + b5(TT5)
truth_D5 <- list(rho = P_D5$rho, T_obs = P_D5$T_obs, T = TT5, fe_plim = fe5, hk_plim = hk_correct(fe5, TT5),
                 spj_plim = P_D5$rho + 2 * b5(TT5) - (floor(TT5 / 2) * b5(floor(TT5 / 2)) +
                                                        ceiling(TT5 / 2) * b5(ceiling(TT5 / 2))) / TT5)

## =================================================================
## Save
## =================================================================
saveRDS(dD2, "nimbus_D2.rds"); saveRDS(dD3, "nimbus_D3.rds"); saveRDS(dD4, "nimbus_D4.rds")
tr <- list(rho = 0.70, rho_persist = 0.90, D1 = truth_D1, D1p = truth_D1p, D1T = truth_D1T, D2 = truth_D2,
           D3 = truth_D3, D4 = truth_D4, D5 = truth_D5,
           params = list(seed = seed, D2 = P_D2, D2_ma = P_D2_ma, D1p_ah = P_D1p_ah, D3 = P_D3, D4 = P_D4, D5 = P_D5))
saveRDS(tr, "truth11.rds")

cat(sprintf("D1: exact Nickell limit %.4f (first factor %.4f, order 1/T %.4f), pooled limit %.4f\n",
            tr$D1$fe_plim, tr$D1$first_factor_plim, tr$D1$order1_plim, tr$D1$pooled_plim))
cat(sprintf("D2: pooled limit %.4f (large %.4f), FE limit %.4f (large %.4f), OLS-diff limit %.4f; r = %d / %d for N = %d\n",
            tr$D2$pooled_plim, tr$D2$large$pooled, tr$D2$fe_plim, tr$D2$large$fe, tr$D2$ols_diff_plim,
            tr$D2$r_full, tr$D2$r_collapsed, P_D2$N))
cat(sprintf("D3: large-sample FE rho %.4f beta %.4f (no feedback: %.4f %.4f); truth rho %.2f beta %.2f\n",
            tr$D3$large$fe_rho, tr$D3$large$fe_beta, tr$D3$large$fe_rho_nofeedback, tr$D3$large$fe_beta_nofeedback,
            tr$D3$rho, tr$D3$beta))
cat(sprintf("D4: corr(Y_i1 - m_i, alpha_i) = %.3f; large-sample two-step AB %.4f BB %.4f BB(twin) %.4f; truth %.2f\n",
            tr$D4$corr_dev_alpha, tr$D4$large$ab, tr$D4$large$bb, tr$D4$large$bb_stat, tr$D4$rho))
cat(sprintf("D4: large-sample one-step AB %.4f BB %.4f (same N = 20,000 sample)\n", tr$D4$large$ab_onestep, tr$D4$large$bb_onestep))
cat(sprintf("D5: T = %s; FE limit %s; HK %s; SPJ %s\n", paste(TT5, collapse = "/"),
            paste(sprintf("%.4f", fe5), collapse = "/"), paste(sprintf("%.4f", tr$D5$hk_plim), collapse = "/"),
            paste(sprintf("%.4f", tr$D5$spj_plim), collapse = "/")))
cat("Saved nimbus_D2.rds, nimbus_D3.rds, nimbus_D4.rds, truth11.rds\n")
