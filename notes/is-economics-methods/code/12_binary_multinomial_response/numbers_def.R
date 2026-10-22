# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers12(dir = "."): every quotable number of chapter 12 as a flat named
# list of finite numeric scalars (10_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth12.rds, res_helix.rds and res_extra.rds in `dir`.
#
# Name scheme: <prefix>_<block>_<field>. Prefixes: dgp (01_dgp.R: seeds, parameters, truth, large-sample
# checks), helix (04_helix.R, inherited block B1), cs (05_cross_section.R, blocks B2 to B4), panel
# (06_panel.R, block B5), mult (07_multinomial.R, blocks B6 and B7). Conventions:
#   *_est, *_se          estimate and standard error (APEs: delta method plus the sampling variance of X,
#                        see avg_effect_inference in bin_fns.R; *_se_delta the delta-method part alone)
#   *_lb, *_ub           95% interval of the same row (est -/+ 1.96 se unless named otherwise:
#                        *_pct_lb, *_pct_ub are bootstrap percentiles)
#   *_covers_*           1 when that interval contains the named value, 0 otherwise
#   *_log10p             log10 of a p-value (ruling 30: tiny p-values are stored as log10)
#   ame                  the legacy field name of the APE (notation 1.10 rule 6: APE in the text)
#   coefficient names    int = intercept; B2 x1, x2, x1x2 (beta_jj'); B3 X (log cloud spend), D, vhat (rho);
#                        B7 pro_*, ent_* (beta_Pro, beta_Ent, base Free) with int, x (engagement), s (log seats)
#   mc_*, *_mc_*         Monte Carlo summaries; cover_* coverage rates of nominal 95% intervals
#   money                *_usd in US dollars, *_kusd in thousands of US dollars
# Block B5 has two "T = 4" designs that are never mixed (controller ruling 42):
#   panel_T4_*, panel_T8_*      the ONE sample (helix_B5.rds, T = 8); panel_T4_* is its months 1 to 4. Truths:
#                               dgp_B5_T8first4_* (months 1 to 4 of the T = 8 panel, alpha_i built on the 8-month
#                               mean) and dgp_B5_T8_* / dgp_B5_ape; panel_T*_truth_* repeat the benchmark used and
#                               panel_T*_*_minus_truth, *_covers_truth, *_minus_plim, *_covers_plim compare with it
#   panel_mc_T4_*, panel_mc_T8_* the Monte Carlo; panel_mc_T4_* draws a 4-month DGP (alpha_i built on the 4-month
#                               mean). Truths: dgp_B5_T4_*, dgp_B5_ape_T4
# B5 truths are EXACT 2-D Gauss-Hermite quadrature values (controller ruling 44): because (X_it, alpha_i) has the same
# law in every month, dgp_B5_T8first4_{ape,rate,pooled_plim,pooled_ape_plim} equal the dgp_B5_T8_* / dgp_B5_ape values;
# dgp_B5_*_draw are the large-draw CHECKS of the same quantities (Monte Carlo error: about 0.003 for the pooled limit,
# 0.00013 for the APE); the CRE-probit limits and the share of movers are large-draw values.
# B5 key meanings: panel_T*_ape_fe and panel_T*_ape_bc are the all-firm APEs (over all firm-months; alpaca::getAPEs' delta
# already divides by nobs.full, the never/always-adopting firms count as zero effects); *_ape_fe_kept is the movers-only
# APE of the FE MLE (delta / share_kept, compared with ape_true_movers; there is no movers-only bias-corrected APE);
# panel_T*_{fe,bc,clogit,pooled}_covers_truth test coverage of beta = 0.8; *_covers_plim use the design's pooled-logit
# limit (T8 sample: dgp_B5_T8_pooled_plim; T4 sample: dgp_B5_T8first4_pooled_plim); panel_T*_ape_*_covers_truth and
# *_pooled_ape_covers_truth use the design's true APE (panel_T*_truth_ape); ape_*_se = population-APE standard error
# (getAPEs' conditional variance ape_*_se_cond plus the firm-clustered sampling variance ape_se_sample, squared).
# Decomposition of the all-firm FE APE error (panel_mc_T*_ape_fe_kept_part, _drop_part, _sampling_part; the one-sample
# panel_T*_ape_fe_{kept,drop}_part split the error against that panel's true APE ape_true_sample): kept_part = the kept
# firms' overstated effects (positive), drop_part = the dropped firms' zero estimated effect (negative), sampling_part
# = mean per-panel true APE minus the population truth; they add up to panel_mc_T*_ape_fe_bias exactly.
# B6 Bond and Lang: *_bl_ratio_h = exp(*_bl_log_ratio_h) is E[h | x1, bill, G = 1] / E[h | x1, bill, G = 0] under
# h(Y*) = -exp(-Y*), CONDITIONAL on the covariates (the unconditional E[h | G] is minus infinity: lognormal bill).
# B7: mult_B7_iia_* / mult_B7_nested_* with hm = the Hausman-McFadden test (hm_test = mlogit::hmftest, drop Free);
# *_suest_* = the stacked-score (suest-type, PSD) contrast on the same samples and Monte Carlo draws.
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json.
# ------------------------------------------------------------------
collect_numbers12 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 08 first")
    readRDS(p)
  }
  tr <- rd("truth12.rds"); h <- rd("res_helix.rds"); x <- rd("res_extra.rds")
  P <- tr$params

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers12: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers12: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers12: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  put_vec <- function(prefix, v, nms) for (k in seq_along(nms)) put(paste0(prefix, "_", nms[k]), v[[k]])
  put_ci <- function(prefix, ci) { put(paste0(prefix, "_lb"), ci[[1]]); put(paste0(prefix, "_ub"), ci[[2]]) }
  put_est <- function(prefix, est, se, ci = NULL) {
    put(paste0(prefix, "_est"), est); put(paste0(prefix, "_se"), se)
    if (is.null(ci)) ci <- c(est - qnorm(0.975) * se, est + qnorm(0.975) * se)
    put_ci(prefix, ci)
  }
  covers_num <- function(ci, v) as.numeric(ci[[1]] <= v && v <= ci[[2]])
  log10p <- function(p) { if (p <= 0) stop("p-value underflow"); log10(p) }

  ## ---- dgp: seeds ---------------------------------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])

  ## ---- dgp: B1 (inherited) -----------------------------------------------------------------------
  T1 <- tr$B1
  put("dgp_B1_n", P$B1$n); put("dgp_B1_beta0", T1$beta[1]); put("dgp_B1_beta1", T1$beta[2])
  put("dgp_B1_oracle_N", 5e6)
  for (k in c("ame", "ame_exact", "mem", "base_rate", "base_rate_exact", "lpm_slope_pop", "lpm_intercept_pop",
              "share_above_pop", "share_below_pop", "probit_ape_pseudo", "logit_over_probit_pseudo"))
    put(paste0("dgp_B1_", k), T1[[k]])
  put("dgp_B1_share_out_pop", T1$share_above_pop + T1$share_below_pop)
  put("dgp_B1_probit_pseudo_beta0", T1$probit_pseudo[1]); put("dgp_B1_probit_pseudo_beta1", T1$probit_pseudo[2])
  put("dgp_B1_ame_over_mem", T1$ame / T1$mem)
  put("dgp_logistic_var", pi^2 / 3); put("dgp_logistic_sd", pi / sqrt(3))
  put("dgp_B1_consistency_R", 400)
  for (n in c(100, 400, 1600, 6400, 25600, 102400)) put(paste0("dgp_B1_consistency_n", n), n)

  ## ---- dgp: B2 -------------------------------------------------------------------------------------
  P2 <- P$B2; T2 <- tr$B2
  put("dgp_B2_n", P2$n); put("dgp_B2_rho_x", P2$rho_x)
  put_vec("dgp_B2_beta", P2$beta, c("int", "x1", "x2", "x1x2"))
  for (k in c("avg_cross", "share_pos", "rate", "ape1", "ape2", "naive_int_me")) put(paste0("dgp_B2_", k), T2[[k]])
  for (k in c("N", "avg_cross", "share_pos", "rate")) put(paste0("dgp_B2_large_", k), T2$large[[k]])

  ## ---- dgp: B3 -------------------------------------------------------------------------------------
  P3 <- P$B3; T3 <- tr$B3
  put("dgp_B3_n", P3$n); put("dgp_B3_spend_median_usd", 1000 * P3$spend_median_kusd); put("dgp_B3_sd_logspend", P3$sd_logspend)
  put("dgp_B3_p_reminder", P3$p_reminder); put_vec("dgp_B3_pi", P3$pi, c("int", "X", "Z"))
  put("dgp_B3_sigma_v", P3$sigma_v); put_vec("dgp_B3_beta", P3$beta, c("int", "X")); put("dgp_B3_beta_D", P3$beta_D)
  put("dgp_B3_corr_eps_v", P3$r); put("dgp_B3_rho", T3$rho); put("dgp_B3_sigma_e", T3$sigma_e)
  put_vec("dgp_B3_scaled", T3$scaled, c("int", "X", "D", "vhat"))
  for (k in c("ape", "rate", "partial_r2_Z")) put(paste0("dgp_B3_", k), T3[[k]])
  put("dgp_B3_large_N", T3$large$N); put("dgp_B3_large_ape", T3$large$ape); put("dgp_B3_large_rate", T3$large$rate)
  put_vec("dgp_B3_naive_plim", T3$large$naive_coef, c("int", "X", "D"))
  put("dgp_B3_naive_ape_plim", T3$large$naive_ape); put("dgp_B3_naive_ape_plim_over_ape", T3$large$naive_ape / T3$ape)

  ## ---- dgp: B4 -------------------------------------------------------------------------------------
  P4 <- P$B4; T4 <- tr$B4
  put("dgp_B4_n", P4$n); put_vec("dgp_B4_beta", P4$beta, c("int", "X")); put("dgp_B4_sigma_a", P4$sigma_a)
  for (k in c("scale", "beta_att", "ape", "rate", "lpm_slope_pop")) put(paste0("dgp_B4_", k), T4[[k]])
  put("dgp_B4_beta_att_int", T4$beta_att_vec[1]); put("dgp_B4_att_factor", 1 / T4$scale)
  put("dgp_B4_large_N", T4$large$N); put("dgp_B4_large_ape", T4$large$ape); put("dgp_B4_large_rate", T4$large$rate)
  put_vec("dgp_B4_large_short", T4$large$short_coef, c("int", "X"))

  ## ---- dgp: B5 -------------------------------------------------------------------------------------
  P5 <- P$B5; T5 <- tr$B5
  for (k in c("N", "T", "T_short", "sigma_c", "sigma_a", "xi0", "xi", "beta", "R_mc", "N_large")) put(paste0("dgp_B5_", k), P5[[k]])
  put("dgp_B5_ape", T5$ape); put("dgp_B5_ape_T4", T5$ape_T4)
  put("dgp_B5_T8first4_ape", T5$T8first4$ape)                       # months 1 to 4 of the T = 8 panel (not the 4-month DGP)
  put("dgp_B5_quad_nodes", T5$T8$quad_nodes)
  for (tt in c("T8", "T4", "T8first4")) {
    L <- T5[[tt]]
    for (k in c("ape", "rate", "pooled_plim", "pooled_ape_plim")) put(paste0("dgp_B5_", tt, "_", k, "_draw"), L$draw[[k]])   # large-draw checks
    for (k in c("rate", "pooled_plim", "pooled_ape_plim", "cre_probit_ape_plim", "share_movers")) put(paste0("dgp_B5_", tt, "_", k), L[[k]])
    put_vec(paste0("dgp_B5_", tt, "_cre_probit_plim"), L$cre_probit_plim, c("int", "X", "Xbar"))
    put(paste0("dgp_B5_", tt, "_pooled_bias"), L$pooled_plim - P5$beta)
  }

  ## ---- dgp: B6 -------------------------------------------------------------------------------------
  P6 <- P$B6; T6 <- tr$B6
  put("dgp_B6_n", P6$n); put_vec("dgp_B6_beta", P6$beta, c("x1", "bill_kusd")); put("dgp_B6_delta_G", P6$delta_G)
  put("dgp_B6_sigma_self_serve", P6$sigma[1]); put("dgp_B6_sigma_managed", P6$sigma[2])
  put_vec("dgp_B6_mu", P6$cuts, 1:4)
  put("dgp_B6_bill_median_usd", 1000 * P6$bill_median_kusd); put("dgp_B6_sd_logbill", P6$sd_logbill)
  put("dgp_B6_mean_bill_usd", 1000 * T6$mean_bill_kusd); put("dgp_B6_p_managed", P6$p_managed); put("dgp_B6_R_mc", P6$R_mc)
  put("dgp_B6_ape_top", T6$ape_top)
  put_vec("dgp_B6_share", T6$shares, paste0("r", 1:5))
  put_vec("dgp_B6_share_self_serve", T6$shares_G["self_serve", ], paste0("r", 1:5))
  put_vec("dgp_B6_share_managed", T6$shares_G["managed", ], paste0("r", 1:5))
  put("dgp_B6_mean_rating_self_serve", T6$mean_rating_G[1]); put("dgp_B6_mean_rating_managed", T6$mean_rating_G[2])
  put("dgp_B6_bl_latent_diff", T6$bl_latent_diff); put("dgp_B6_bl_log_ratio_h", T6$bl_log_ratio_h)
  put("dgp_B6_bl_ratio_h", exp(T6$bl_log_ratio_h))
  put("dgp_B6_large_N", T6$large$N); put("dgp_B6_large_ape_top", T6$large$ape_top)

  ## ---- dgp: B7 -------------------------------------------------------------------------------------
  P7 <- P$B7; T7 <- tr$B7
  put("dgp_B7_n", P7$n); put("dgp_B7_lambda", P7$lambda_nested); put("dgp_B7_R_mc", P7$R_mc)
  for (m in c("Pro", "Enterprise")) put_vec(paste0("dgp_B7_", if (m == "Pro") "pro" else "ent"), P7$coef[m, ], c("int", "x", "s"))
  lv <- c("free", "pro", "ent")
  put_vec("dgp_B7_share", T7$shares, lv); put_vec("dgp_B7_share_nested", T7$shares_nested, lv)
  put_vec("dgp_B7_ape_x", T7$ape_x, lv); put_vec("dgp_B7_ape_s", T7$ape_s, lv)
  put_vec("dgp_B7_within_paid_nested", T7$within_paid_nested, c("int", "x", "s"))
  put_vec("dgp_B7_ent_minus_pro", P7$coef["Enterprise", ] - P7$coef["Pro", ], c("int", "x", "s"))
  put("dgp_B7_large_N", T7$large$N)

  ## ---- helix: block B1 results (04_helix.R; legacy 02 and 03) ---------------------------------------
  L <- h$lpm
  for (k in c("intercept", "slope", "se", "se_intercept", "se_hc1", "se_hc0", "min_phat", "max_phat", "share_out",
              "share_above", "share_below", "n_above", "n_below")) put(paste0("helix_lpm_", k), L[[k]])
  put_ci("helix_lpm_slope_hc1", L$ci_hc1)
  G <- h$logit
  for (k in c("beta0", "beta1", "se0_fisher", "se1_fisher", "se0_sand", "se1_sand", "sand_over_fisher", "ame",
              "ame_se", "ame_se_delta", "pea", "mean_x", "iter", "se1_glm", "se1_hc0", "loglik"))
    put(paste0("helix_logit_", k), G[[k]])
  put("helix_logit_glm_maxdiff_log10", log10(G$glm_maxdiff))           # ruling 30: tiny numbers as log10
  put_ci("helix_logit_ame", G$ame_ci); put_ci("helix_logit_beta1", G$ci1)
  put("helix_logit_beta1_covers_truth", covers_num(G$ci1, T1$beta[2]))
  put("helix_logit_ame_covers_truth", covers_num(G$ame_ci, T1$ame_exact))
  put("helix_logit_ame_minus_lpm_slope", G$ame - L$slope)
  pr <- h$probit
  for (k in c("beta0", "beta1", "se1", "ape", "ape_se", "pea", "ratio_logit_probit", "loglik")) put(paste0("helix_probit_", k), pr[[k]])
  put_ci("helix_probit_ape", pr$ape_ci)
  cons <- h$consistency
  for (j in seq_len(ncol(cons))) for (k in rownames(cons)) put(paste0("helix_cons_n", colnames(cons)[j], "_", k), cons[k, j])
  put("helix_cons_sd_ratio_100_400", cons["sd", "100"] / cons["sd", "400"])
  put("helix_cons_sd_ratio_25600_102400", cons["sd", "25600"] / cons["sd", "102400"])
  for (k in c("n", "adoption_rate", "mean_x", "sd_x")) put(paste0("helix_sample_", k), h$sample[[k]])

  ## ---- cs: block B2 -------------------------------------------------------------------------------
  B2 <- x$B2; cn2 <- c("int", "x1", "x2", "x1x2")
  for (j in 1:4) put_est(paste0("cs_B2_coef_", cn2[j]), B2$coef[j], B2$se[j])
  put("cs_B2_se_interaction_robust", B2$se_interaction_robust); put("cs_B2_z_interaction", B2$z_interaction)
  put("cs_B2_interaction_log10p", log10p(B2$p_interaction)); put("cs_B2_lr_interaction", B2$lr_interaction)
  put("cs_B2_lr_interaction_log10p", log10p(B2$p_lr_interaction))
  put_est("cs_B2_avg_cross", B2$avg_cross, B2$avg_cross_se, B2$avg_cross_ci)
  put("cs_B2_avg_cross_se_delta", B2$avg_cross_se_delta); put("cs_B2_avg_cross_z", B2$avg_cross_z)
  put("cs_B2_avg_cross_covers_truth", B2$covers_truth); put("cs_B2_avg_cross_minus_truth", B2$avg_cross - T2$avg_cross)
  put("cs_B2_share_pos", B2$share_pos); put("cs_B2_cross_min", B2$cross_min); put("cs_B2_cross_max", B2$cross_max)
  put_est("cs_B2_ape1", B2$ape1, B2$ape1_se, B2$ape1_ci); put_est("cs_B2_ape2", B2$ape2, B2$ape2_se, B2$ape2_ci)
  put("cs_B2_ape1_se_marginaleffects", B2$ape1_se_me)
  put("cs_B2_ape1_covers_truth", covers_num(B2$ape1_ci, T2$ape1)); put("cs_B2_ape2_covers_truth", covers_num(B2$ape2_ci, T2$ape2))
  put("cs_B2_naive_int_me_est", B2$naive_int_me); put_ci("cs_B2_naive_int_me", B2$naive_int_me_ci)
  put("cs_B2_rate", B2$rate); put("cs_B2_n", B2$n)

  ## ---- cs: block B3 -------------------------------------------------------------------------------
  B3 <- x$B3; cn3 <- c("int", "X", "D")
  for (j in 1:3) put_est(paste0("cs_B3_naive_coef_", cn3[j]), B3$naive$coef[j], B3$naive$se[j])
  put_est("cs_B3_naive_ape", B3$naive$ape, B3$naive$ape_se, B3$naive$ape_ci)
  put("cs_B3_naive_ape_covers_truth", B3$naive$covers_truth); put("cs_B3_naive_ape_over_truth", B3$naive$ape / T3$ape)
  cf <- B3$cf; cn3b <- c("int", "X", "D", "vhat")
  for (j in 1:4) {
    put_est(paste0("cs_B3_cf_coef_", cn3b[j]), cf$beta[j], cf$beta_se[j])
    put(paste0("cs_B3_cf_coef_", cn3b[j], "_covers_scaled_truth"), covers_num(cf$beta_ci[[j]], T3$scaled[j]))
  }
  put_est("cs_B3_cf_ape", cf$ape, cf$ape_se, cf$ape_ci); put_ci("cs_B3_cf_ape_pct", cf$ape_ci_pct)
  put("cs_B3_cf_ape_covers_truth", cf$covers_truth); put("cs_B3_cf_ape_pct_covers_truth", cf$covers_truth_pct)
  put("cs_B3_cf_ape_minus_truth", cf$ape - T3$ape)
  put_vec("cs_B3_cf_pi", cf$pi, c("int", "X", "Z")); put("cs_B3_cf_rho_t", cf$rho_t); put("cs_B3_cf_rho_log10p", log10p(cf$rho_p))
  put("cs_B3_first_stage_F", cf$first_stage_F); put("cs_B3_cf_B", cf$B)
  iv <- B3$ivmle
  for (j in 1:3) put_est(paste0("cs_B3_ivmle_coef_", cn3[j]), iv$beta[j], iv$se_beta[j])
  put("cs_B3_ivmle_sigma_v", iv$sigma_v); put("cs_B3_ivmle_corr", iv$r); put("cs_B3_ivmle_loglik", iv$loglik)
  put_est("cs_B3_ivmle_ape", iv$ape, iv$ape_se, iv$ape_ci); put("cs_B3_ivmle_ape_covers_truth", iv$covers_truth)
  put("cs_B3_ivmle_beta_D_covers_truth", covers_num(iv$ci_beta_D, T3$beta_D[1]))
  put("cs_B3_ivmle_minus_cf_ape", iv$ape - cf$ape)
  put_est("cs_B3_tsls", B3$tsls$est, B3$tsls$se, B3$tsls$ci); put("cs_B3_tsls_covers_ape", B3$tsls$covers_ape)
  put_est("cs_B3_ols_lpm", B3$ols_lpm$est, B3$ols_lpm$se)
  put("cs_B3_rate", B3$rate); put("cs_B3_n", B3$n)
  put("cs_B3_mean_spend_usd", 1000 * B3$mean_spend_kusd); put("cs_B3_median_spend_usd", 1000 * B3$median_spend_kusd)

  ## ---- cs: block B4 -------------------------------------------------------------------------------
  B4 <- x$B4
  put_est("cs_B4_beta", B4$beta, B4$beta_se, B4$beta_ci); put("cs_B4_coef_int", B4$coef[1])
  put("cs_B4_beta_covers_att", B4$covers_att); put("cs_B4_beta_covers_beta", covers_num(B4$beta_ci, T4$beta))
  put("cs_B4_beta_over_truth", B4$beta / T4$beta)
  put_est("cs_B4_ape", B4$ape, B4$ape_se, B4$ape_ci); put("cs_B4_ape_covers_truth", B4$covers_ape)
  for (j in 1:3) put_est(paste0("cs_B4_long_coef_", c("int", "X", "a")[j]), B4$long$coef[j], B4$long$se[j])
  put_est("cs_B4_long_ape", B4$long$ape, B4$long$ape_se, B4$long$ape_ci); put("cs_B4_long_ape_covers_truth", B4$long$covers_ape)
  put("cs_B4_ratio_short_long", B4$ratio_short_long)
  put_est("cs_B4_lpm_slope", B4$lpm_slope, B4$lpm_se, B4$lpm_ci)
  put("cs_B4_lpm_covers_ape", covers_num(B4$lpm_ci, T4$ape))
  put("cs_B4_rate", B4$rate); put("cs_B4_n", B4$n)

  ## ---- panel: block B5 ------------------------------------------------------------------------------
  B5 <- x$B5
  for (tt in c("T8", "T4")) {
    S <- B5[[tt]]; pf <- paste0("panel_", tt, "_")
    for (k in c("fe", "bc", "clogit")) {
      put_est(paste0(pf, k), S[[k]]$est, S[[k]]$se, S[[k]]$ci)
      put(paste0(pf, k, "_covers_truth"), covers_num(S[[k]]$ci, P5$beta))
    }
    put_est(paste0(pf, "pooled"), S$pooled$est, S$pooled$se, S$pooled$ci); put(paste0(pf, "pooled_se_model"), S$pooled$se_model)
    put(paste0(pf, "pooled_covers_truth"), covers_num(S$pooled$ci, P5$beta))
    put_est(paste0(pf, "pooled_ape"), S$pooled$ape, S$pooled$ape_se, S$pooled$ape_ci)
    put(paste0(pf, "cmle_est"), S$cmle$est); put(paste0(pf, "cmle_se_model"), S$cmle$se_model)
    put(paste0(pf, "cmle_se_robust"), S$cmle$se_robust); put_ci(paste0(pf, "cmle_robust"), S$cmle$ci_robust)
    put(paste0(pf, "cmle_n_units_used"), S$cmle$n_units_used)
    put(paste0(pf, "cmle_diff_clogit_log10"), log10(max(S$cmle$max_diff_clogit, 1e-300)))
    for (k in c("share_kept", "n_firms_kept", "ape_fe", "ape_bc", "ape_fe_se", "ape_bc_se", "ape_fe_se_cond", "ape_bc_se_cond",
                "ape_se_sample", "ape_fe_kept", "ape_fe_kept_part", "ape_fe_drop_part", "ape_true_sample", "ape_true_movers")) put(paste0(pf, k), S[[k]])
    put_ci(paste0(pf, "ape_fe"), S$fe$ape_ci); put_ci(paste0(pf, "ape_bc"), S$bc$ape_ci)
    # benchmark truth of this design and the comparisons (T8: the T = 8 truths; T4: months 1 to 4 of the T = 8 panel)
    put(paste0(pf, "truth_ape"), S$truth$ape); put(paste0(pf, "truth_pooled_plim"), S$truth$pooled_plim)
    for (k in names(S$vs)) put(paste0(pf, k), S$vs[[k]])
  }
  cre <- B5$cre
  for (j in 1:3) put_est(paste0("panel_cre_coef_", c("int", "X", "Xbar")[j]), cre$beta[j], cre$se[j])
  put_est("panel_cre_ape", cre$ape[[1]], cre$ape_se[[1]], cre$ape_ci[[1]]); put("panel_cre_ape_covers_truth", covers_num(cre$ape_ci[[1]], T5$ape))
  rr <- B5$cre_re
  for (j in 1:3) {
    nm <- c("int", "X", "Xbar")[j]
    put_est(paste0("panel_cre_re_coef_", nm), rr$coef[j], rr$se[j]); put(paste0("panel_cre_re_scaled_", nm), rr$scaled_coef[j])
  }
  put("panel_cre_re_sigma_a", rr$sigma_a); put("panel_cre_re_ape", rr$ape)
  for (tt in c("T4", "T8")) {
    m <- B5$mc[[tt]]; pf <- paste0("panel_mc_", tt, "_")
    for (k in c("R", "fe_mean", "fe_sd", "fe_bias", "fe_rel_bias", "bc_mean", "bc_sd", "bc_bias", "bc_rel_bias",
                "cl_mean", "cl_sd", "cl_bias", "pooled_mean", "pooled_bias", "cover_fe", "cover_bc", "cover_cl",
                "mean_se_cl", "mean_se_bc", "ape_true", "ape_fe_bias", "ape_bc_bias", "ape_fe_rel_bias", "ape_bc_rel_bias",
                "ape_movers_true", "ape_fe_kept_bias", "ape_fe_kept_rel_bias",
                "ape_fe_kept_part", "ape_fe_drop_part", "ape_fe_sampling_part",
                "share_kept", "ape_fe_sd", "ape_bc_sd", "ape_fe_mean_se", "ape_bc_mean_se", "ape_fe_cover", "ape_bc_cover",
                "ape_true_sample_mean")) put(paste0(pf, k), m[[k]])
  }
  put("panel_fe_bias_T4", B5$fe$bias_T4); put("panel_fe_bias_T8", B5$fe$bias_T8)
  put("panel_bc_bias_T4", B5$bc$bias_T4); put("panel_bc_bias_T8", B5$bc$bias_T8)
  put("panel_fe_bias_ratio_T4_T8", B5$fe$bias_T4 / B5$fe$bias_T8)
  for (k in c("N", "T", "rate", "share_movers", "share_movers_T4")) put(paste0("panel_sample_", k), B5$sample[[k]])

  ## ---- mult: block B6 -------------------------------------------------------------------------------
  B6 <- x$B6; cn6 <- c("mu1", "mu2", "mu3", "mu4", "x1", "bill_kusd", "G", "log_sigma_G")
  for (j in seq_along(cn6)) put_est(paste0("mult_B6_coef_", cn6[j]), B6$coef[j], B6$se[j])
  put("mult_B6_sigma_managed", B6$sigma_managed)
  put("mult_B6_x1_minus_truth_in_se", (B6$coef[["x1"]] - P6$beta[1]) / B6$se[["x1"]])
  put_est("mult_B6_ape_top", B6$ape_top, B6$ape_top_se, B6$ape_top_ci); put("mult_B6_ape_top_se_delta", B6$ape_top_se_delta)
  put("mult_B6_ape_top_covers_truth", B6$covers_truth)                       # 0: the ruling-29 miss
  put("mult_B6_ape_top_minus_truth_in_se", (B6$ape_top - T6$ape_top) / B6$ape_top_se)
  put_vec("mult_B6_ape", B6$ape_by_rating, paste0("r", 1:5))
  for (j in 1:7) put_est(paste0("mult_B6_homo_coef_", cn6[j]), B6$homo$coef[j], B6$homo$se[j])
  put("mult_B6_homo_ape_top_est", B6$homo$ape_top); put_ci("mult_B6_homo_ape_top", B6$homo$ape_top_ci)
  put("mult_B6_lr_scale", B6$lr_scale); put("mult_B6_lr_scale_log10p", log10p(B6$p_lr_scale))
  put("mult_B6_bl_latent_diff", B6$bl_latent_diff); put("mult_B6_bl_log_ratio_h", B6$bl_log_ratio_h)
  put("mult_B6_bl_ratio_h", exp(B6$bl_log_ratio_h))
  put_vec("mult_B6_share", B6$shares, paste0("r", 1:5))
  put_vec("mult_B6_share_self_serve", B6$shares_G["self_serve", ], paste0("r", 1:5))
  put_vec("mult_B6_share_managed", B6$shares_G["managed", ], paste0("r", 1:5))
  put("mult_B6_mean_rating_self_serve", B6$mean_rating_G[["self_serve"]]); put("mult_B6_mean_rating_managed", B6$mean_rating_G[["managed"]])
  put("mult_B6_n", B6$n); put("mult_B6_share_managed", B6$share_managed)
  put("mult_B6_mean_bill_usd", 1000 * B6$mean_bill_kusd); put("mult_B6_median_bill_usd", 1000 * B6$median_bill_kusd)
  for (k in c("R", "mean_est", "sd_est", "mean_se", "cover", "bias", "mean_b1", "share_est_above_sample"))
    put(paste0("mult_B6_mc_", k), B6$mc[[k]])

  ## ---- mult: block B7 -------------------------------------------------------------------------------
  I <- x$B7$iia; Nn <- x$B7$nested
  for (m in c("Pro", "Enterprise")) for (j in 1:3) {
    pf <- paste0("mult_B7_", if (m == "Pro") "pro" else "ent", "_", c("int", "x", "s")[j])
    put_est(pf, I$coef[m, j], I$se[m, j], c(I$ci_lb[m, j], I$ci_ub[m, j]))
    put(paste0(pf, "_covers_truth"), as.numeric(I$covers[m, j]))
    put(paste0(pf, "_nested_est"), Nn$coef[m, j])
  }
  put_vec("mult_B7_ape_x", I$ape_x, lv); put_vec("mult_B7_ape_x_se", I$ape_x_se, lv)
  put_vec("mult_B7_ape_s", I$ape_s, lv); put_vec("mult_B7_ape_s_se", I$ape_s_se, lv)
  for (k in seq_along(lv)) put(paste0("mult_B7_ape_x_", lv[k], "_covers_truth"),
                               covers_num(I$ape_x[k] + c(-1, 1) * qnorm(0.975) * I$ape_x_se[k], T7$ape_x[k]))
  put_vec("mult_B7_share", I$shares, lv)
  put_vec("mult_B7_full_ent_minus_pro", I$full_ent_minus_pro, c("int", "x", "s"))
  put("mult_B7_iia_hm_stat", I$hm_stat); put("mult_B7_iia_hm_df", I$hm_df); put("mult_B7_iia_hm_p", I$hm_p)
  put("mult_B7_iia_hm_stat_mlogit", I$hm_stat_mlogit)
  put("mult_B7_iia_hm_ent_stat", I$hm_ent_stat); put("mult_B7_iia_hm_ent_p", I$hm_ent_p)
  for (k in c("mc_R", "mc_reject", "mc_reject01", "mc_reject10", "mc_neg_share", "mc_mean_stat", "mc_median_stat",
              "mc_mean_stat_pos", "mc_mcse", "mc_nonpsd_share", "mc_median_eigen_ratio")) put(paste0("mult_B7_iia_", k), I[[k]])
  # stacked-score (suest-type) contrast: single IIA sample and the same Monte Carlo draws
  su <- I$suest
  put("mult_B7_iia_suest_stat", su$hm_stat); put("mult_B7_iia_suest_df", su$hm_df); put("mult_B7_iia_suest_p", su$hm_p)
  for (k in c("mc_R", "mc_reject01", "mc_reject", "mc_reject10", "mc_neg_share", "mc_mean_stat", "mc_median_stat", "mc_min_eigen_min"))
    put(paste0("mult_B7_iia_suest_", k), su[[k]])
  put("mult_B7_iia_suest_mc_reject_minus_005_in_mcse", (su$mc_reject - 0.05) / I$mc_mcse)
  put("mult_B7_iia_mc_reject_minus_005_in_mcse", (I$mc_reject - 0.05) / I$mc_mcse)
  put("mult_B7_nested_hm_stat", Nn$hm_stat); put("mult_B7_nested_hm_df", Nn$hm_df); put("mult_B7_nested_hm_log10p", Nn$hm_log10p)
  put("mult_B7_nested_hm_stat_mlogit", Nn$hm_stat_mlogit)
  put("mult_B7_nested_suest_stat", Nn$suest$hm_stat); put("mult_B7_nested_suest_df", Nn$suest$hm_df)
  put("mult_B7_nested_suest_log10p", Nn$suest$hm_log10p)
  put("mult_B7_nested_hm_ent_stat", Nn$hm_ent_stat); put("mult_B7_nested_hm_ent_log10p", Nn$hm_ent_log10p)
  put_vec("mult_B7_nested_full_ent_minus_pro", Nn$full_ent_minus_pro, c("int", "x", "s"))
  put_vec("mult_B7_nested_paid_logit", Nn$paid_logit, c("int", "x", "s"))
  put_vec("mult_B7_nested_share", Nn$shares, lv)

  nums
}
