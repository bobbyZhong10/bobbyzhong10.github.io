# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers09(dir = "."): every quotable number of chapter 9 as a flat named
# list of finite numeric scalars (07_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth09.rds and res_system.rds, res_contagion_r0.rds, res_sur.rds, res_weak.rds, res_overid.rds in `dir`.
#
# Name scheme: <prefix>_<object>_<field>. Prefixes: dgp (01_dgp.R: parameters, truth, population and
# large-sample quantities), sys (02, block S1), mis (02, block S1m), sur (03, block S2), weak (04, block S3),
# overid (05, blocks S4 and S4v). Conventions:
#   *_est, *_se             point estimate and standard error: 2SLS and OLS HC0 (eq. 5.34) unless named
#                           *_se_homo (eq. 5.27, n - k); 3SLS eq. 8.45 unless named *_se_robust
#   *_lb, *_ub              95% interval est -/+ qnorm(0.975) se of the same row
#   *_covers_truth          1 when that interval contains the true value, 0 otherwise
#   *_p                     p-value of the statistic named just before it
#   cover_*                 Monte Carlo coverage rate of nominal 95% intervals or sets
#   mc_reject_*             Monte Carlo rejection rate of a 5% test (_01: of a 1% test)
#   Pi_<z>_<y>              reduced-form coefficient of exogenous variable z in the equation for y
#   mu2_<v>                 weak-IV grid point with concentration parameter v (decimal point written p)
#   r_minus0p7              population SE ratio when corr(coupon, fee) = -0.7
#   ar_type                 Anderson-Rubin set shape: 1 bounded interval, 2 two rays, 3 the whole line
#   ar_bounded, ar_n_pieces always stored (the set's shape must not change the key set): ar_bounded = 1 when the
#                           set is a bounded interval, 0 otherwise; ar_n_pieces = number of disjoint pieces
#                           (bounded interval 1, two rays 2, the whole line 1, empty 0). The end-point keys
#                           weak_sample_ar_lo / _hi exist ONLY when finite, i.e. for a bounded interval (both are
#                           its end points) and for two rays (lo closes the left ray (-Inf, lo], hi opens the
#                           right ray [hi, Inf)); for the whole line neither key exists. Read ar_type, ar_bounded
#                           or ar_n_pieces before quoting lo or hi.
#   *_n_denominator         2SLS homoskedastic SE with sigma_hat^2 = sum eps_hat^2 / n (not n - k), the same
#                           convention as the 3SLS SE (Omega_hat / n); the ratio with the same suffix divides the
#                           3SLS SE by it
#   dgp_S1m_r0_*, dgp_S1m_r_minus0p7_*   large-sample (n = 1e6, 02b_contagion_r0.R, own seeds) limits of the
#                           S1m buyer coefficients n_sellers and coupon under 2SLS and 3SLS, with
#                           corr(coupon, fee) = 0 and -0.7; *_bias = limit - truth (n_sellers 1.2, coupon 0.4);
#                           *_three_minus_tsls = 3SLS limit - 2SLS limit
# Units: buyers in thousands, n_sellers in hundreds, ad_budget in USD thousands, coupon, fee,
# free_ship_min, basket and referral bonus in USD. A non-finite value is an error (put() stops).
# Rounding-level diagnostics (differences that are zero in exact arithmetic, *_max_abs*, sur_pkg_*) are
# stored to 3 significant digits (put_tiny): only their order of magnitude means anything.
# ------------------------------------------------------------------
collect_numbers09 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 05 first")
    readRDS(p)
  }
  tr <- rd("truth09.rds"); s1 <- rd("res_system.rds"); s2 <- rd("res_sur.rds"); s3 <- rd("res_weak.rds"); s4 <- rd("res_overid.rds")
  s1c <- rd("res_contagion_r0.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers09: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers09: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers09: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  put_ci <- function(prefix, est, se, ci, truth = NULL) {
    put(paste0(prefix, "_est"), est); put(paste0(prefix, "_se"), se)
    put(paste0(prefix, "_lb"), ci[1]); put(paste0(prefix, "_ub"), ci[2])
    if (!is.null(truth)) put(paste0(prefix, "_covers_truth"), covers01(ci, truth))
  }
  put_list <- function(prefix, l, fields = names(l)) for (k in fields) put(paste0(prefix, "_", k), l[[k]])
  put_tiny <- function(name, value) put(name, signif(value, 3))
  pid <- function(x) gsub("-", "minus", gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)))
  P <- tr$params
  b <- tr$S1$buyer[["n_sellers"]]

  ## ---- dgp: seeds and parameters ---------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  put("dgp_N_large", P$N_large)
  P1 <- P$S1
  put("dgp_S1_n", P1$n)
  put("dgp_S1_ad_mean", P1$ad0); put("dgp_S1_ad_sd", P1$ad_sd); put("dgp_S1_coupon_mean", P1$cp0); put("dgp_S1_coupon_sd", P1$cp_sd)
  put("dgp_S1_fee_mean", P1$fee0); put("dgp_S1_fee_sd", P1$fee_sd); put("dgp_S1_corr_coupon_fee", P1$r)
  put("dgp_S1_mean_buyers", P1$EB); put("dgp_S1_mean_sellers", P1$ES)
  put("dgp_S1_sd_eps_buyer", P1$sB); put("dgp_S1_sd_eps_seller", P1$sS); put("dgp_S1_corr_eps", P1$rho)
  for (k in names(tr$S1$buyer)) put(paste0("dgp_S1_buyer_", k), tr$S1$buyer[[k]])
  for (k in names(tr$S1$seller)) put(paste0("dgp_S1_seller_", k), tr$S1$seller[[k]])
  put("dgp_S1_feedback_product", P1$bBS * P1$bSB); put("dgp_S1_multiplier", tr$S1$multiplier)
  put("dgp_S1_usd_per_buyer_from_ads", 1 / P1$gad)                                   # USD 1,000 of ads brings 1,000 gad buyers
  put("dgp_S1_buyers_per_100_sellers", 1000 * P1$bBS); put("dgp_S1_sellers_per_1000_buyers", 100 * P1$bSB)
  Om <- tr$S1$Omega; Ov <- tr$S1$Omega_v
  put("dgp_S1_Omega_11", Om[1, 1]); put("dgp_S1_Omega_12", Om[1, 2]); put("dgp_S1_Omega_22", Om[2, 2])
  put("dgp_S1_Omega_v_11", Ov[1, 1]); put("dgp_S1_Omega_v_12", Ov[1, 2]); put("dgp_S1_Omega_v_22", Ov[2, 2])
  put("dgp_S1_corr_eB_v", tr$S1$corr_eB_v)
  Pi <- tr$S1$Pi
  for (zn in rownames(Pi)) for (yn in colnames(Pi)) put(sprintf("dgp_S1_Pi_%s_%s", zn, yn), Pi[zn, yn])
  for (m in c("buyer", "seller")) {
    id <- tr$S1$ident[[m]]
    for (k in c("n_excluded", "n_endog", "order", "rank", "rank_ok", "overid")) put(sprintf("dgp_S1_ident_%s_%s", m, k), id[[k]])
    put(sprintf("dgp_S1_rank_general_%s", m), tr$S1$rank_general[[m]])
  }
  put("dgp_S1_no_fee_order", tr$S1$ident_no_fee$order); put("dgp_S1_no_fee_rank", tr$S1$ident_no_fee$rank)
  put("dgp_S1_no_fee_rank_ok", tr$S1$ident_no_fee$ok)
  pa <- tr$S1$pop_avar
  put("dgp_S1_pop_se_ratio", pa$se_ratio); put("dgp_S1_pop_se_gain_pct", 100 * (1 - pa$se_ratio))
  put("dgp_S1_pop_se_ratio_independent", pa$se_ratio_independent)
  put("dgp_S1_pop_se_2sls_n", pa$se_2sls_n); put("dgp_S1_pop_se_3sls_n", pa$se_3sls_n)
  for (i in which(pa$by_r$r <= 0)) put(paste0("dgp_S1_pop_se_ratio_r_", pid(pa$by_r$r[i])), pa$by_r$se_ratio[i])
  pl <- tr$S1$plim
  put_list("dgp_S1_plim", pl)
  put("dgp_S1_plim_ols_bias", pl$ols_buyer - b); put("dgp_S1_plim_ols_bias_pct", 100 * (pl$ols_buyer / b - 1))
  # S1m
  put("dgp_S1m_kappa", P$S1m$kappa)
  put("dgp_S1m_kappa_times_mean_buyers", P$S1m$kappa * P1$EB)                      # implied direct coupon effect on sellers (hundreds per USD)
  put("dgp_S1m_slope_coupon_plus1usd_pct", 100 * P$S1m$kappa / P1$bSB)               # change of the seller slope per USD of coupon (%)
  put("dgp_S1m_slope_coupon_plus1sd", P1$bSB + P$S1m$kappa * P1$cp_sd)
  put("dgp_S1m_slope_coupon_plus1sd_pct", 100 * P$S1m$kappa * P1$cp_sd / P1$bSB)
  plm <- tr$S1m$plim
  put_list("dgp_S1m_plim", plm)
  put("dgp_S1m_plim_three_sls_bias", plm$three_sls_buyer - b); put("dgp_S1m_plim_three_sls_bias_pct", 100 * (plm$three_sls_buyer / b - 1))
  put("dgp_S1m_plim_tsls_bias", plm$tsls_buyer - b)
  # S1m at corr(coupon, fee) = 0 (and, for comparison, -0.7 on a different draw): which buyer coefficient is contaminated
  for (d in c("r0", "r_minus0p7")) {
    cr <- s1c[[d]]; nm <- paste0("dgp_S1m_", d)
    put(paste0("dgp_seed_S1m_", d, "_large"), cr$seed); put(paste0(nm, "_corr_coupon_fee_sample"), cr$corr_coupon_fee)
    for (v in c("n_sellers", "coupon")) {
      put(paste0(nm, "_", v, "_tsls"), cr[[paste0(v, "_tsls")]]); put(paste0(nm, "_", v, "_three_sls"), cr[[paste0(v, "_three_sls")]])
      put(paste0(nm, "_", v, "_tsls_bias"), cr[[paste0(v, "_tsls")]] - s1c$truth[[v]])
      put(paste0(nm, "_", v, "_three_sls_bias"), cr[[paste0(v, "_three_sls")]] - s1c$truth[[v]])
      put(paste0(nm, "_", v, "_three_minus_tsls"), cr[[paste0(v, "_three_sls")]] - cr[[paste0(v, "_tsls")]])
    }
    put(paste0(nm, "_seller_buyers_tsls"), cr$seller_buyers_tsls); put(paste0(nm, "_J_seller_per_n"), cr$J_seller_per_n)
  }
  # S2
  P2 <- P$S2
  for (k in c("n", "ship0", "ship_sd", "prem0", "prem_ship", "prem_sd", "loy_ship", "loy_cut", "s1", "s2", "rho")) put(paste0("dgp_S2_", k), P2[[k]])
  for (k in names(tr$S2$beta_same_X)) put(paste0("dgp_S2_beta_", k), tr$S2$beta_same_X[[k]])
  for (k in c("p_loyalty", "E_basket", "E_orders", "revenue_effect")) put(paste0("dgp_S2_", k), tr$S2[[k]])
  # S3
  P3 <- P$S3
  for (k in c("mu2", "fee_sd", "pi_f", "sigma_v", "corr_eB_v", "plim_ols")) put(paste0("dgp_S3_", k), tr$S3[[k]])
  put("dgp_S3_n", P3$n); put("dgp_S3_R_mc", P3$R_mc); put("dgp_S3_R_grid", P3$R_grid); put("dgp_S3_corr_coupon_fee", P3$r)
  put("dgp_S3_fee_sd_cents", 100 * tr$S3$fee_sd); put("dgp_S3_plim_ols_bias", tr$S3$plim_ols - b)
  for (i in seq_len(nrow(tr$S3$grid))) {
    g <- tr$S3$grid[i, ]
    put(paste0("dgp_S3_grid_mu2_", pid(g$mu2), "_fee_sd"), g$fee_sd); put(paste0("dgp_S3_grid_mu2_", pid(g$mu2), "_plim_ols"), g$plim_ols)
  }
  # S4
  P4 <- P$S4
  for (k in c("n", "fee0", "fee_sd", "pay_lo", "pay_hi", "ref_lo", "ref_hi", "ref_mid", "dpay", "dref", "theta", "eta", "R_mc", "R_power"))
    put(paste0("dgp_S4_", k), P4[[k]])
  put("dgp_S4_aS", P4$aS); put("dgp_S4_overid_df", tr$S4$overid_df)
  put("dgp_S4_theta_buyers_per_10usd", 10 * 1000 * P4$theta)                         # buyers per USD 10 of referral bonus
  t2 <- 2 * P4$eta * sqrt(3); hmax <- exp(P4$eta * sqrt(3)) / sqrt(sinh(t2) / t2); hmin <- exp(-P4$eta * sqrt(3)) / sqrt(sinh(t2) / t2)
  put("dgp_S4_h_max", hmax); put("dgp_S4_h_min", hmin); put("dgp_S4_var_ratio_max_min", (hmax / hmin)^2)

  ## ---- sys: block S1 ---------------------------------------------------------------------------
  put_ci("sys_ols_buyer", s1$ols$buyer_est, s1$ols$buyer_se, s1$ols$buyer_ci, b)
  put("sys_ols_buyer_se_homo", s1$ols$buyer_se_homo); put("sys_ols_buyer_minus_truth", s1$ols$buyer_est - b)
  put_ci("sys_ols_seller", s1$ols$seller_est, s1$ols$seller_se, s1$ols$seller_ci, tr$S1$seller[["buyers"]])
  put_ci("sys_tsls_buyer", s1$tsls$buyer_est, s1$tsls$buyer_se, s1$tsls$buyer_ci, b)
  put("sys_tsls_buyer_se_homo", s1$tsls$buyer_se_homo); put("sys_tsls_buyer_minus_truth", s1$tsls$buyer_est - b)
  put_ci("sys_tsls_seller", s1$tsls$seller_est, s1$tsls$seller_se, s1$tsls$seller_ci, tr$S1$seller[["buyers"]])
  put("sys_tsls_seller_se_homo", s1$tsls$seller_se_homo)
  put_ci("sys_three_sls_buyer", s1$three_sls$buyer_est, s1$three_sls$buyer_se, s1$three_sls$buyer_ci, b)
  put("sys_three_sls_buyer_se_robust", s1$three_sls$buyer_se_robust); put("sys_three_sls_buyer_minus_truth", s1$three_sls$buyer_est - b)
  put_ci("sys_three_sls_seller", s1$three_sls$seller_est, s1$three_sls$seller_se, s1$three_sls$seller_ci, tr$S1$seller[["buyers"]])
  put("sys_three_minus_tsls_buyer_n_sellers", s1$diff$buyer_n_sellers)
  put("sys_three_minus_tsls_buyer_max_abs", max(abs(s1$diff$buyer))); put_tiny("sys_three_minus_tsls_seller_max_abs", max(abs(s1$diff$seller)))
  put("sys_se_tsls_homo", s1$se$tsls); put("sys_se_three_sls", s1$se$three_sls)
  put("sys_se_tsls_robust", s1$se$tsls_robust); put("sys_se_three_sls_robust", s1$se$three_sls_robust)
  put("sys_se_ratio", s1$se$ratio); put("sys_se_ratio_robust", s1$se$ratio_robust); put("sys_se_gain_pct", 100 * (1 - s1$se$ratio))
  put("sys_tsls_se_n_denominator", s1$se$tsls_n_denominator); put("sys_se_ratio_n_denominator", s1$se$ratio_n_denominator)
  put("sys_se_gain_pct_n_denominator", 100 * (1 - s1$se$ratio_n_denominator))
  put("sys_rho_hat", s1$three_sls$rho_hat)
  put("sys_Omega_hat_11", s1$three_sls$Sigma[1, 1]); put("sys_Omega_hat_12", s1$three_sls$Sigma[1, 2]); put("sys_Omega_hat_22", s1$three_sls$Sigma[2, 2])
  put_list("sys_fs", s1$first_stage)
  put("sys_cf_rho", s1$endog$rho); put("sys_cf_se", s1$endog$se); put("sys_cf_t", s1$endog$t); put("sys_cf_p", s1$endog$p)
  put("sys_J_seller", s1$overid$J_seller); put("sys_J_seller_p", s1$overid$J_seller_p); put("sys_J_seller_df", s1$overid$J_seller_df)
  put("sys_sargan_seller", s1$overid$sargan_seller); put("sys_sargan_seller_p", s1$overid$sargan_seller_p)
  for (m in c("buyer", "seller")) {
    dg <- s1$diagnostics[[m]]
    put(sprintf("sys_diag_%s_weak_F", m), dg["Weak instruments", "statistic"])
    put(sprintf("sys_diag_%s_wu_hausman", m), dg["Wu-Hausman", "statistic"]); put(sprintf("sys_diag_%s_wu_hausman_p", m), dg["Wu-Hausman", "p-value"])
  }
  put("sys_diag_seller_sargan", s1$diagnostics$seller["Sargan", "statistic"])
  Ph <- s1$identification$Pi_hat
  for (zn in rownames(Ph)) for (yn in colnames(Ph)) put(sprintf("sys_Pi_hat_%s_%s", zn, yn), Ph[zn, yn])
  put_tiny("sys_pkg_tsls_max_abs_diff", max(abs(s1$tsls_hand$b - s1$tsls_pkg$b)))
  put_tiny("sys_pkg_three_sls_max_abs_diff", max(abs(s1$three_sls_hand$b - s1$three_sls_pkg$b)))
  put("sys_n", s1$n)
  put_list("sys_mc", s1$mc$S1)
  put("sys_mc_sd_gain_pct", 100 * (1 - s1$mc$S1$sd_ratio))
  put("sys_mc_bias_tsls", s1$mc$S1$mean_tsls - b); put("sys_mc_bias_three_sls", s1$mc$S1$mean_three_sls - b)
  put("mc_mcse_coverage_R1000", sqrt(0.95 * 0.05 / 1000)); put("mc_mcse_coverage_R2000", sqrt(0.95 * 0.05 / 2000))
  put("mc_mcse_rate05_R2000", sqrt(0.05 * 0.95 / 2000))

  ## ---- mis: block S1m --------------------------------------------------------------------------
  ms <- s1$misspec
  put_ci("mis_tsls_buyer", ms$tsls_est, ms$tsls_se, ms$tsls_ci, b); put("mis_tsls_bias", ms$tsls_bias)
  put_ci("mis_three_sls_buyer", ms$three_sls_est, ms$three_sls_se, ms$three_sls_ci, b); put("mis_three_sls_bias", ms$three_sls_bias)
  put("mis_three_sls_bias_pct", 100 * ms$three_sls_bias / b); put("mis_three_minus_tsls", ms$three_sls_minus_tsls)
  put("mis_bias_ratio", abs(ms$three_sls_bias) / abs(ms$tsls_bias))
  put_ci("mis_tsls_seller", ms$seller_tsls_est, ms$seller_tsls_se, ms$seller_tsls_ci, tr$S1$seller[["buyers"]])
  put("mis_three_sls_seller_est", ms$seller_three_sls_est); put_tiny("mis_three_minus_tsls_seller_max_abs", max(abs(ms$seller_diff)))
  put("mis_ols_buyer_est", ms$ols_est); put("mis_rho_hat", ms$rho_hat)
  put("mis_J_seller", ms$J_seller); put("mis_J_seller_p", ms$J_seller_p)
  put("mis_sargan_seller", ms$sargan_seller); put("mis_sargan_seller_p", ms$sargan_seller_p)
  put("mis_se_three_sls", ms$se_three_sls); put("mis_se_tsls_homo", ms$se_tsls_homo)
  put_list("mis_mc", s1$mc$S1m)
  put("mis_mc_bias_tsls", s1$mc$S1m$mean_tsls - b); put("mis_mc_bias_three_sls", s1$mc$S1m$mean_three_sls - b)
  put("mis_mc_bias_three_sls_in_sd", (s1$mc$S1m$mean_three_sls - b) / s1$mc$S1m$sd_three_sls)

  ## ---- sur: block S2 ---------------------------------------------------------------------------
  cname <- function(x) gsub(":", "_", x)
  sx <- s2$same_X
  put_tiny("sur_same_max_abs_sur_minus_ols", sx$max_abs_sur_minus_ols); put("sur_same_rho_hat", sx$rho_hat)
  for (j in seq_along(sx$sur)) {
    nm <- paste0("sur_same_", cname(names(sx$sur)[j]))
    put(paste0(nm, "_est"), sx$sur[[j]]); put(paste0(nm, "_se"), sx$se_sur[[j]]); put(paste0(nm, "_se_lm"), sx$se_lm[[j]])
  }
  dx <- s2$diff_X
  put("sur_diff_rho_hat", dx$rho_hat)
  put("sur_diff_Omega_hat_11", dx$Sigma[1, 1]); put("sur_diff_Omega_hat_12", dx$Sigma[1, 2]); put("sur_diff_Omega_hat_22", dx$Sigma[2, 2])
  for (j in seq_along(dx$sur)) {
    nm <- paste0("sur_diff_", cname(names(dx$sur)[j]))
    put(paste0(nm, "_sur_est"), dx$sur[[j]]); put(paste0(nm, "_ols_est"), dx$ols[[j]])
    put(paste0(nm, "_se_sur"), dx$se_sur[[j]]); put(paste0(nm, "_se_ols"), dx$se_ols[[j]]); put(paste0(nm, "_se_lm"), dx$se_lm[[j]])
    put(paste0(nm, "_se_ratio"), dx$se_ratio[[j]]); put(paste0(nm, "_se_gain_pct"), 100 * (1 - dx$se_ratio[[j]]))
    put(paste0(nm, "_covers_truth"), s2$covers_diff[[j]])
  }
  for (sp in c("same_X", "diff_X")) {
    rv <- s2$revenue[[sp]]
    nm <- paste0("sur_revenue_", sub("_X", "", sp))
    put_ci(nm, rv$est, rv$se_full, rv$ci, tr$S2$revenue_effect)
    put(paste0(nm, "_se_blockdiag"), rv$se_blockdiag); put(paste0(nm, "_se_ratio"), rv$se_ratio)
    put(paste0(nm, "_se_understated_pct"), 100 * (1 - 1 / rv$se_ratio))
    put(paste0(nm, "_corr_ship"), rv$corr_ship); put(paste0(nm, "_w_orders"), rv$w_orders); put(paste0(nm, "_w_basket"), rv$w_basket)
  }
  for (k in names(s2$pkg)) put_tiny(paste0("sur_pkg_", k), s2$pkg[[k]])
  put_list("sur_sample", s2$sample)

  ## ---- weak: block S3 --------------------------------------------------------------------------
  sm <- s3$sample
  put_ci("weak_sample_tsls", sm$tsls_est, sm$tsls_se, sm$wald_ci, b)
  put("weak_sample_tsls_se_homo", sm$tsls_se_homo)
  put("weak_sample_F_fs", sm$F_fs); put("weak_sample_F_eff", sm$F_eff); put("weak_sample_pi_fee", sm$pi_fee)
  put("weak_sample_ar_type", match(sm$ar$type, c("bounded", "two_rays", "real_line", "empty")))
  put("weak_sample_ar_bounded", sm$ar$type == "bounded")
  put("weak_sample_ar_n_pieces", c(bounded = 1, two_rays = 2, real_line = 1, empty = 0)[[sm$ar$type]])
  # lo / hi exist only when finite (bounded interval or two rays; none for the whole line): see the header
  if (is.finite(sm$ar$lo)) put("weak_sample_ar_lo", sm$ar$lo)
  if (is.finite(sm$ar$hi)) put("weak_sample_ar_hi", sm$ar$hi)
  put("weak_sample_ar_covers_truth", sm$ar_covers)
  put("weak_sample_ar_stat_at_truth", sm$ar_at_truth$stat); put("weak_sample_ar_stat_at_truth_p", sm$ar_at_truth$p)
  put("weak_sample_ar_stat_at_zero", sm$ar_at_zero$stat); put("weak_sample_ar_stat_at_zero_p", sm$ar_at_zero$p)
  put("weak_sample_ar_df2", sm$ar_at_truth$df2); put("weak_sample_ar_crit", qf(0.95, 1, sm$ar_at_truth$df2))
  put("weak_sample_liml_est", sm$liml_est); put("weak_sample_fuller_est", sm$fuller_est); put("weak_sample_fuller_se", sm$fuller_se)
  put("weak_sample_ols_est", sm$ols_est); put("weak_sample_ols_se", sm$ols_se)
  put("weak_cover_wald", s3$cover$wald); put("weak_cover_wald_homo", s3$cover$wald_homo); put("weak_cover_ar", s3$cover$ar)
  put("weak_cover_R", s3$cover$R); put("weak_cover_mu2", s3$cover$mu2)
  mfields <- c("R", "mean_F", "median_F", "share_F_below_10", "cover_wald", "cover_wald_homo", "cover_ar", "share_ar_bounded",
               "share_ar_two_rays", "share_ar_real_line", "median_est", "q25_est", "q75_est", "median_bias", "median_se_r",
               "share_wald_reject_zero")
  put_list("weak_mc", s3$main, mfields)
  for (i in seq_len(nrow(s3$grid))) put_list(paste0("weak_grid_mu2_", pid(s3$grid$mu2[i])), s3$grid[i, ], mfields)

  ## ---- overid: blocks S4 and S4v ------------------------------------------------------------------
  for (k in c("valid", "invalid")) {
    o <- s4[[k]]; nm <- paste0("overid_", k)
    put_ci(paste0(nm, "_tsls"), o$tsls_est, o$tsls_se, o$tsls_ci, tr$S4$beta[["n_sellers"]])
    put_ci(paste0(nm, "_tsls_two_valid"), o$tsls_valid2_est, o$tsls_valid2_se, o$tsls_valid2_ci, tr$S4$beta[["n_sellers"]])
    for (f in c("J", "J_p", "J_df", "gmm_est", "sargan", "sargan_p", "J_valid2", "J_valid2_p", "J_valid2_df", "F_fs", "F_eff"))
      put(paste0(nm, "_", f), o[[f]])
  }
  v <- s4$valid; iv <- s4$invalid
  for (f in c("R", "mc_reject_J", "mc_reject_sargan", "mc_reject_J_01", "mc_reject_sargan_01", "mc_mean_est", "mc_sd_est"))
    put(paste0("overid_valid_", f), v[[f]])
  for (f in c("R_power", "mc_reject_J", "mc_reject_J_01", "mc_reject_sargan", "mc_mean_est", "mc_sd_est", "mc_bias"))
    put(paste0("overid_invalid_", f), iv[[f]])
  put("overid_invalid_tsls_bias", iv$tsls_est - tr$S4$beta[["n_sellers"]])
  put("overid_invalid_diag_sargan", iv$diagnostics["Sargan", "statistic"])
  put("overid_crit_chi2_2_05", s4$crit$chi2_2_05); put("overid_crit_chi2_2_01", s4$crit$chi2_2_01)

  nums
}

covers01 <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])
