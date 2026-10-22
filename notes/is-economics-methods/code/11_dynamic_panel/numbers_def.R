# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers11(dir = "."): every quotable number of chapter 11 as a flat named list
# of finite numeric scalars (10_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth11.rds, res_dynamic.rds (the legacy 03) and res_extra.rds (08) in `dir`.
#
# Name scheme: <prefix>_<object>_<field>. Prefixes: dgp (01_dgp.R: seeds, parameters, truth, exact and
# closed-form limits, large-sample checks), d1 and d1p (blocks D1, D1p: 02 and the legacy 03), d1t (the
# legacy Nickell Monte Carlo of 03), d2 (04), d3 (05), d4 (06), d5 (07). Conventions:
#   T_obs, T          months observed per merchant; T = T_obs - 1 regression periods
#   *_est, *_se       estimate and its SE: CR0 clustered by merchant (pooled OLS, FE, Anderson-Hsiao), the
#                     robust one-step SE (GMM "ab1", D2 to D4 one-step), the Windmeijer-corrected two-step SE
#                     ("ab2", "bb2", "*_twostep"); *_se_conv is the uncorrected two-step SE
#   *_lb, *_ub        95% interval est -/+ qnorm(0.975) se of the same row
#   *_covers_truth    1 when that interval contains the truth, 0 otherwise
#   *_r               number of instruments; *_J, *_J_df, *_J_p Hansen J (sargan(fit, "twosteps")) and its
#                     p-value; *_m1 (z statistic), *_m2 and *_m2_p the Arellano-Bond serial-correlation tests
#                     (robust variance)
#   *_plim            probability limit as N -> infinity (exact Nickell formula for FE; for the half-panel jackknife the
#                     limit of the Dhaene-Jochmans (2015) estimator with the halves weighted by |S|/T, eqs. 2.2 and 2.4)
#   *_spj_*           half-panel jackknife: *_spj_split<c>_half<h> are the FE estimates on the two pieces of cut c (odd T: cut 1
#                     has floor(T/2) and ceiling(T/2) periods, cut 2 the reverse), *_spj_half_mean = bar_rho_{1/2}, the cuts'
#                     |S|/T-weighted piece averages averaged over the cuts (eq. 2.4), *_spj_est = 2 FE - bar_rho_{1/2}
#   *_usd             a parameter converted to US dollars (the blocks record thousands or hundreds of USD)
#   d1_legacy_*       as the legacy code computes them (Sargan with two-step weights, mtest with the
#                     conventional variance)
#   d1p_ah_diff_mc_*  Monte Carlo (seed D1p_ah_mc, 1,000 samples of the D1p design) of Anderson-Hsiao with the differenced
#                     instrument: cover = share of nominal 95% intervals containing 0.90, median, mean (no finite mean: see
#                     below), iqr and quantiles of the estimates
#   dgp_D4_large_*    large-sample (N = 20,000, one draw) estimates: two-step ab, bb, bb_stat; one-step ab_onestep, bb_onestep
#   d2_ma_*           the Monte Carlo of the misspecified (MA(1)) twin of D2: rej_* = rejection rate of the
#                     Hansen J test at 5% for the full (r > N), middle (lags 2:mid_lag, r < N) and collapsed
#                     set; mean_J_*; mean_est_*
#   d2_mc_*J_*, d2_ma_rej_*, d2_ma_mean_J_*   in the two D2 Monte Carlos: the J of the one-step fit (moments and
#                     weight from the one-step residuals; chapter Theorem 6, equal to N when r >= N)
#   d2_mc_*J2_*, d2_ma_rej2_*, d2_ma_mean_J2_*  the two-step Hansen J of chapter Theorem 4 (iv) (sargan() of the
#                     two-step fit; J2 <= J in every sample: flags d2_mc_J2_le_J, d2_ma_J2_le_J; d2_mc_J2_full_max
#                     is the largest J2 of the full set, at most N: flag d2_mc_J2_full_le_N)
# P-values (every key that ends in _p, or in _p_full or _p_collapsed, apart from means and shares): one rule.
# A p-value of 1e-4 or more is stored as is. A p-value below 1e-4 is stored as the flag <key>_lt_0001 = 1
# (<key>_p_lt_0001 when the key does not end in _p) and, when it is positive, as log10 of the p-value in
# <key>_log10 (<key>_p_log10): 1e-33 does not round-trip through JSON text bit for bit, 0 would be a false
# value, and the text says "p < 0.0001" with the log10 where it needs the size. No p-value is dropped or
# stored as 0. put() applies the rule to every such key. The one-step Sargan p-value of plm (always 0) is not
# a test and is not stored (dyn_fns.R, gmm_summary).
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json.
# ------------------------------------------------------------------
collect_numbers11 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 08 first")
    readRDS(p)
  }
  tr <- rd("truth11.rds"); dy <- rd("res_dynamic.rds"); x <- rd("res_extra.rds")
  P <- tr$params

  nums <- list()
  put_raw <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers11: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers11: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers11: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  # a p-value key: ends in _p, _p_full or _p_collapsed, and is not a mean or a share of p-values
  is_pkey <- function(name) grepl("(^|_)p$|_p_(full|collapsed)$", name) && !grepl("mean|share", name)
  put <- function(name, value) {
    if (!is_pkey(name)) return(put_raw(name, value))
    p <- as.numeric(value)
    if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) stop("collect_numbers11: ", name, " is not a p-value (", p, ")")
    if (p >= 1e-4) return(put_raw(name, p))
    base <- if (grepl("_p$", name)) name else paste0(name, "_p")
    put_raw(paste0(sub("_p$", "", base), "_p_lt_0001"), 1)                    # flag: p < 1e-4
    if (p > 0) put_raw(paste0(sub("_p$", "", base), "_p_log10"), log10(p))   # log10 of the p-value, where finite
    invisible(NULL)
  }
  put_est <- function(prefix, s, truth = NULL) {
    put(paste0(prefix, "_est"), s$est); put(paste0(prefix, "_se"), s$se)
    put(paste0(prefix, "_lb"), s$ci[1]); put(paste0(prefix, "_ub"), s$ci[2])
    if (!is.null(truth)) put(paste0(prefix, "_covers_truth"), as.numeric(s$ci[1] <= truth && truth <= s$ci[2]))
  }
  put_gmm <- function(prefix, s, truth = NULL) {
    put_est(prefix, s, truth)
    if (!is.na(s$se_conv)) { put(paste0(prefix, "_se_conv"), s$se_conv); put(paste0(prefix, "_se_ratio"), s$se / s$se_conv) }
    # m1_p is as small as 1e-68 here: put() stores it as the flag m1_p_lt_0001 and log10 of the p-value (rule above).
    # There is no one-step Sargan key: plm's sargan(fit, "onestep") is an unscaled statistic (dyn_fns.R, gmm_summary).
    for (k in c("r", "J", "J_df", "J_p", "m1", "m1_p", "m2", "m2_p")) put(paste0(prefix, "_", k), s[[k]])
  }
  put_list <- function(prefix, l, keys = names(l)) for (k in keys) put(paste0(prefix, "_", k), l[[k]])

  ## ---- dgp: seeds, parameters, truth ---------------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  put("dgp_rho", tr$rho); put("dgp_rho_persist", tr$rho_persist)
  for (b in c("D1", "D1p")) {
    z <- tr[[b]]
    put_list(paste0("dgp_", b), z, c("T_obs", "T", "fe_plim", "first_factor_plim", "order1_plim", "pooled_plim"))
    put(paste0("dgp_", b, "_fe_bias"), z$fe_plim - z$rho)
    put(paste0("dgp_", b, "_first_factor"), z$first_factor_plim - z$rho)
    put(paste0("dgp_", b, "_order1_bias"), z$order1_plim - z$rho)
  }
  put("dgp_D1_N", 500)
  z <- tr$D1T
  for (j in seq_along(z$T_obs)) {
    nm <- paste0("dgp_D1T_T", z$T_obs[j])
    put(paste0(nm, "_fe_plim"), z$fe_plim[j]); put(paste0(nm, "_first_factor_plim"), z$first_factor_plim[j])
    put(paste0(nm, "_legacy_approx"), z$legacy_approx[j]); put(paste0(nm, "_fe_bias"), z$fe_plim[j] - z$rho)
  }
  P2 <- P$D2
  put_list("dgp_D2", P2, c("N", "T", "rho", "m_mean", "m_sd", "sigma_e", "burn", "R_mc", "mid_lag"))
  put_list("dgp_D2_ma", P$D2_ma, c("ma_coef", "R_mc"))
  put("dgp_D2_m_mean_usd", 1000 * P2$m_mean); put("dgp_D2_m_sd_usd", 1000 * P2$m_sd); put("dgp_D2_sigma_e_usd", 1000 * P2$sigma_e)
  put_list("dgp_D2", tr$D2, c("pooled_plim", "fe_plim", "ols_diff_plim", "r_full", "r_collapsed"))
  put_list("dgp_D2_large", tr$D2$large)
  put("dgp_D2_r_full_diff", (P2$T - 1) * (P2$T - 2) / 2); put("dgp_D2_r_full_level", P2$T - 2)
  put("dgp_D2_r_time", P2$T - 1)
  P3 <- P$D3
  put_list("dgp_D3", P3, c("N", "T", "rho", "beta", "a_mean", "a_sd", "x_mean", "x_alpha", "feedback", "eps_x_sd", "sigma_e", "burn"))
  put_list("dgp_D3", tr$D3, c("long_run", "mean_x", "mean_y"))
  put_list("dgp_D3_large", tr$D3$large)
  put("dgp_D3_large_fe_beta_bias", tr$D3$large$fe_beta - P3$beta); put("dgp_D3_large_fe_rho_bias", tr$D3$large$fe_rho - P3$rho)
  put("dgp_D3_beta_usd_per_100", 1000 * P3$beta)                 # USD of GMV per 100 USD of promotion, this month
  put("dgp_D3_long_run_usd_per_100", 1000 * tr$D3$long_run)      # ... in the long run
  put("dgp_D3_mean_y_usd", 1000 * tr$D3$mean_y); put("dgp_D3_mean_x_usd", 100 * tr$D3$mean_x)
  put("dgp_D3_feedback_usd", -100 * P3$feedback)                 # extra promotion after a month 1,000 USD below expectation
  put("dgp_D3_eps_x_sd_usd", 100 * P3$eps_x_sd)                   # sd of the promotion shock eps_x, in USD
  P4 <- P$D4
  put_list("dgp_D4", P4, c("N", "T", "rho", "m_mean", "m_sd", "start_share", "sigma_e"))
  put("dgp_D4_start_share_pct", 100 * P4$start_share); put("dgp_D4_m_mean_usd", 1000 * P4$m_mean); put("dgp_D4_m_sd_usd", 1000 * P4$m_sd)
  put("dgp_D4_corr_dev_alpha", tr$D4$corr_dev_alpha)
  put_list("dgp_D4_large", tr$D4$large)
  put_list("dgp_D1p_ah", P$D1p_ah, c("N", "T_obs", "rho", "R_mc"))   # the Anderson-Hsiao Monte Carlo at rho = 0.90 (02)
  P5 <- P$D5
  put_list("dgp_D5", P5, c("N", "rho", "R_mc"))
  for (j in seq_along(P5$T_obs)) {
    nm <- paste0("dgp_D5_T", P5$T_obs[j])
    put(paste0(nm, "_T"), tr$D5$T[j])
    put(paste0(nm, "_fe_plim"), tr$D5$fe_plim[j]); put(paste0(nm, "_hk_plim"), tr$D5$hk_plim[j])
    put(paste0(nm, "_spj_plim"), tr$D5$spj_plim[j])
  }

  ## ---- d1, d1p: the inherited samples (02; legacy 03) ------------------------------------------------
  for (b in c("D1", "D1p")) {
    z <- x[[b]]; pre <- tolower(b); truth <- if (b == "D1") tr$rho else tr$rho_persist
    put_est(paste0(pre, "_pool"), z$pool, truth); put_est(paste0(pre, "_fe"), z$fe, truth)
    put(paste0(pre, "_fe_bias"), z$fe$est - truth)
    put_est(paste0(pre, "_hk"), list(est = z$hk, se = z$hk_se, ci = z$hk_ci), truth)
    put(paste0(pre, "_spj_est"), z$spj); put(paste0(pre, "_spj_half_mean"), z$spj_half_mean)
    for (i in seq_len(nrow(z$spj_halves))) for (h in 1:2) put(sprintf("%s_spj_split%d_half%d", pre, i, h), z$spj_halves[i, h])
    put_est(paste0(pre, "_ah"), z$ah, truth); put(paste0(pre, "_ah_n"), z$ah$n)
    put_est(paste0(pre, "_ah_diff"), z$ah_diff, truth); put(paste0(pre, "_ah_diff_n"), z$ah_diff$n)
    put_gmm(paste0(pre, "_ab1"), z$ab_onestep, truth); put_gmm(paste0(pre, "_ab2"), z$ab, truth); put_gmm(paste0(pre, "_bb2"), z$bb, truth)
    put(paste0(pre, "_T"), z$T)
  }
  # Monte Carlo of the Anderson-Hsiao estimator with the differenced instrument in the D1p design (rho = 0.90): R samples,
  # coverage of the nominal 95% interval (and its Monte Carlo SE; the shares of intervals entirely below and entirely above
  # 0.90), median and mean of the estimates, the median's bias, the interquartile range and the 2.5% and 97.5% quantiles.
  # A just-identified IV estimator has no finite moments: the mean is dominated by a few extreme draws (no SD is stored).
  put_list("d1p_ah_diff_mc", x$D1p$ah_diff_mc, c("R", "cover", "cover_mcse", "share_ci_below_truth", "share_ci_above_truth",
                                                  "median", "median_bias", "iqr", "q025", "q975"))   # no mean or SD: a just-identified IV has no finite moments
  put("d1_ah_se_over_ab1_se", x$ah$se / x$ab_onestep$se)
  put("d1p_ab2_se_over_bb2_se", x$D1p$ab$se / x$D1p$bb$se)
  # The legacy "Sargan p 0.290" is the two-step Hansen J p-value (plm::sargan(fit, "twosteps") of the two-step
  # fit, the "Sargan test" line of summary.pgmm); it equals d1_ab2_J_p. Not a one-step Sargan statistic.
  put("d1_legacy_sargan_p", dy$sargan_p); put("d1_legacy_ar1_p", dy$ar1_p); put("d1_legacy_ar2_p", dy$ar2_p)

  ## ---- d1t: the legacy Nickell Monte Carlo (03, set.seed(303)) ---------------------------------------
  put("d1t_R", 40); put("d1t_N", 300)
  for (j in seq_len(ncol(dy$nick))) {
    nm <- paste0("d1t_T", colnames(dy$nick)[j])
    put(paste0(nm, "_fe_mc"), dy$nick["fe", j]); put(paste0(nm, "_ab_mc"), dy$nick["ab", j])
    put(paste0(nm, "_fe_mc_minus_plim"), dy$nick["fe", j] - tr$D1T$fe_plim[j])
  }

  ## ---- d2: instrument proliferation (04) --------------------------------------------------------------
  D2 <- x$D2
  put_est("d2_full", list(est = D2$est_full, se = D2$se_full, ci = D2$ci_full), tr$D2$rho)
  put_est("d2_collapsed", list(est = D2$est_collapsed, se = D2$se_collapsed, ci = D2$ci_collapsed), tr$D2$rho)
  put_list("d2", D2, c("J_full", "J_df_full", "J_p_full", "J_collapsed", "J_df_collapsed", "J_p_collapsed", "m2_p_full",
                       "m2_p_collapsed", "ols_transformed", "ols_transformed_se", "ols_levels_no_month", "fe", "fe_se",
                       "ols_diff", "n_inst_full", "n_inst_collapsed", "N", "T"))
  put("d2_full_minus_ols", D2$est_full - D2$ols_transformed); put("d2_collapsed_minus_ols", D2$est_collapsed - D2$ols_transformed)
  put("d2_full_minus_collapsed", D2$est_full - D2$est_collapsed)
  pa <- D2$path
  for (i in seq_len(nrow(pa))) {
    nm <- sprintf("d2_path_%s_lag%d", if (pa$collapse[i]) "collapsed" else "full", pa$max_lag[i])
    put_list(nm, as.list(pa[i, ]), c("r", "est", "se", "lo", "hi", "J", "J_df", "J_p"))
  }
  put_list("d2_mc", D2$mc, c("R", "mean_full", "mean_collapsed", "mean_ols", "mean_fe", "sd_full", "sd_collapsed",
                             "share_J_full_above_099", "mean_J_p_full", "mean_J_p_collapsed", "share_J_collapsed_below_005",
                             "share_J_full_below_005", "share_J_mid_below_005", "mean_J_p_mid", "share_closer",
                             "share_J2_full_above_099", "share_J2_full_below_005", "share_J2_mid_below_005",
                             "share_J2_collapsed_below_005", "mean_J2_p_full", "mean_J2_p_mid", "mean_J2_p_collapsed",
                             "J2_full_max"))
  # the two-step Hansen J2 never exceeds the J of the one-step fit (same weight, the two-step estimate minimises
  # the quadratic form), so J2 <= N for the full set (r >= N); flags: 1 = holds in every sample (1e-8 tolerance)
  put("d2_mc_J2_le_J", as.numeric(D2$mc$J2_minus_J_max <= 1e-8))
  put("d2_mc_J2_full_le_N", as.numeric(D2$mc$J2_full_max <= D2$N + 1e-6))
  put("d2_ma_J2_le_J", as.numeric(D2$mis$J2_minus_J_max <= 1e-8))
  put("d2_mc_bias_full", D2$mc$mean_full - tr$D2$rho); put("d2_mc_bias_collapsed", D2$mc$mean_collapsed - tr$D2$rho)
  # J is N exactly for the full set (r >= N): the p-value depends on the instrument count alone. Flags (1 = the
  # largest |J - N| is below 1e-6; the deviations themselves are rounding noise, about 1e-10): the D2 sample, its
  # 500 Monte Carlo samples, the 200 samples of the misspecified twin
  put("d2_J_full_equals_N", as.numeric(abs(D2$J_full - D2$N) < 1e-6))
  put("d2_mc_J_full_equals_N", as.numeric(D2$mc$J_full_maxdev_from_N < 1e-6))
  put("d2_ma_J_full_equals_N", as.numeric(D2$mis$J_full_maxdev_from_N < 1e-6))
  # misspecified MA(1) twin: Hansen J rejection rates at 5% (the power of the test), 200 samples
  put_list("d2_ma", D2$mis, c("R", "ma_coef", "mid_lag", "r_full", "r_mid", "r_collapsed", "rej_full", "rej_mid", "rej_collapsed",
                              "mean_J_full", "mean_J_mid", "mean_J_collapsed", "mean_est_full", "mean_est_mid",
                              "mean_est_collapsed", "rej2_full", "rej2_mid", "rej2_collapsed", "mean_J2_full",
                              "mean_J2_mid", "mean_J2_collapsed"))

  ## ---- d3: predetermined regressor with feedback (05) -----------------------------------------------
  D3 <- x$D3
  put_est("d3_fe_beta", D3$fe, tr$D3$beta); put_est("d3_fe_rho", D3$fe_rho, tr$D3$rho)
  for (k in c("ab", "ab_twostep", "ab_endog", "ab_exog")) {
    put_gmm(paste0("d3_", k, "_beta"), D3[[k]], tr$D3$beta)
    put_est(paste0("d3_", k, "_rho"), D3[[paste0(k, "_rho")]], tr$D3$rho)
  }
  put("d3_fe_beta_bias", D3$fe$est - tr$D3$beta); put("d3_fe_rho_bias", D3$fe_rho$est - tr$D3$rho)

  ## ---- d4: non-stationary initial condition (06) ------------------------------------------------------
  D4 <- x$D4
  for (k in c("ab", "bb", "ab_twostep", "bb_twostep", "ab_stat", "bb_stat")) put_gmm(paste0("d4_", k), D4[[k]], tr$D4$rho)
  put_est("d4_fe", D4$fe, tr$D4$rho)
  put("d4_bb_bias", D4$bb$est - tr$D4$rho)

  ## ---- d5: bias-corrected FE Monte Carlo (07) ---------------------------------------------------------
  D5 <- x$D5
  put("d5_R", D5$R); put("d5_N", D5$N)
  for (j in seq_along(D5$T_obs)) {
    nm <- paste0("d5_T", D5$T_obs[j])
    for (k in c("fe", "hk", "spj", "sd_fe", "sd_hk", "sd_spj", "rmse_fe", "rmse_hk", "rmse_spj", "mean_se_fe",
                "cover_fe", "cover_hk", "cover_spj")) put(paste0(nm, "_", k), D5[[k]][j])
    for (k in c("fe", "hk", "spj")) put(paste0(nm, "_bias_", k), D5[[k]][j] - tr$rho)
  }
  put("d5_mcse_coverage_R1000", sqrt(0.95 * 0.05 / D5$R))

  nums
}
