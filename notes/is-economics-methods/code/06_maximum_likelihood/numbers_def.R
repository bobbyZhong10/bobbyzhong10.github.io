# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers06(dir = "."): every quotable number of chapter 6 as a
# flat named list of finite numeric scalars (07_summary.R writes it to numbers.json with
# chapter_kit.R's write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth06.rds, the four sorrel_C*.rds data sets and res_cens.rds, res_qmle.rds,
# res_grouped.rds, res_twostep.rds in `dir`.
#
# Name scheme: <prefix>_<block or object>_<field>. Prefixes: dgp (01_dgp.R: parameters, truth,
# realised sample facts), cens (02, block C1), qmle (03, block C2), grouped (04, block C3),
# twostep (05, block C4). Conventions:
#   *_est, *_se          point estimate and its standard error (the SE named in the key when there are
#                        several: *_se_model = Hessian A^-1 / n, *_se_opg = B^-1 / n, *_se_sandwich,
#                        *_se_naive, *_se_corrected, *_se_glm, *_se_cluster)
#   *_lb, *_ub           95% interval est -/+ qnorm(0.975) se of the same row (the SE named in the key,
#                        otherwise the row's *_se)
#   *_covers_*           1 when that interval contains the named value, 0 otherwise
#   *_cover_*            Monte Carlo coverage rate of nominal 95% intervals
#   b0, b_fee, b_ref, b_u, log_kappa, kappa   the Weibull PH parameters (beta_0, beta_fee, beta_ref, beta_u,
#                        log kappa, kappa); fee in USD per month; times in months
#   m01 ... m24          grouped-duration interval m (block C3)
# Scales of the block C1 naive-regression keys (OLS of log tenure on x, ols_weibull):
#   cens_ols_*_fee_*     the PH (hazard) scale of truth06$C1$beta: beta_fee = -kappa_hat b_fee with
#                        kappa_hat = pi / sqrt(6 s2) from the same regression; on censored spells this mixes
#                        two distortions, the attenuated slope b_fee and the shrunken residual variance
#                        (cens_ols_observed_kappa = 1.014 against the true 0.9)
#   cens_ols_*_slope_*,  the raw log-tenure (AFT) slope b_fee, no kappa involved; its truth is
#   cens_aft_*           cens_aft_truth = -beta_fee / kappa
# cens_mle_fee_est (with _se, _lb, _ub) and cens_mle_b_fee_est are the same number (the first is the
# interval row the card's test uses, the second part of the full parameter listing); both are kept.
# ------------------------------------------------------------------
collect_numbers06 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 05 first")
    readRDS(p)
  }
  tr <- rd("truth06.rds"); C1 <- rd("sorrel_C1.rds"); C2 <- rd("sorrel_C2.rds"); C4 <- rd("sorrel_C4.rds")
  c1 <- rd("res_cens.rds"); c2 <- rd("res_qmle.rds"); c3 <- rd("res_grouped.rds"); c4 <- rd("res_twostep.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers06: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers06: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers06: ", name, " is not finite")
    nums[[name]] <<- value
    invisible(NULL)
  }
  put_ci <- function(prefix, est, se, ci = NULL) {
    put(paste0(prefix, "_est"), est); put(paste0(prefix, "_se"), se)
    if (is.null(ci)) ci <- c(est - qnorm(0.975) * se, est + qnorm(0.975) * se)
    put(paste0(prefix, "_lb"), ci[1]); put(paste0(prefix, "_ub"), ci[2])
  }
  covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])
  mm <- function(m) sprintf("m%02d", m)
  pn <- c("b0", "b_fee", "b_ref", "log_kappa")
  P <- tr$params

  ## ---- dgp: seeds, design, parameters, truth, realised samples --------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  P1 <- P$C1
  put("dgp_months", P1$months); put("dgp_growth", P1$growth); put("dgp_growth_pct", P1$growth - 1)
  for (i in seq_along(P1$fees)) {
    put(sprintf("dgp_fee_%d_usd", P1$fees[i]), P1$fees[i]); put(sprintf("dgp_fee_%d_share", P1$fees[i]), P1$p_fee[i])
  }
  put("dgp_C1_n", P1$n); put("dgp_C1_p_ref", P1$p_ref)
  for (i in 1:3) put(paste0("dgp_C1_", pn[i]), P1$beta[i])
  put("dgp_C1_kappa", P1$kappa); put("dgp_C1_beta_aft", tr$C1$beta_aft)
  put("dgp_C1_hr_fee_per_usd", exp(P1$beta[2])); put("dgp_C1_pct_hazard_per_usd", exp(P1$beta[2]) - 1)
  put("dgp_C1_hr_fee_per_15usd", exp(15 * P1$beta[2])); put("dgp_C1_pct_hazard_per_15usd", exp(15 * P1$beta[2]) - 1)
  put("dgp_C1_hr_ref", exp(P1$beta[3])); put("dgp_C1_pct_hazard_ref", exp(P1$beta[3]) - 1)
  put("dgp_C1_median_fee60", tr$C1$median_fee60); put("dgp_C1_S12_fee60", tr$C1$S12_fee60)
  put("dgp_C1_share_cens_pop", tr$C1$share_cens_pop)
  put("dgp_C1_sample_share_cens", 1 - mean(C1$delta)); put("dgp_C1_sample_n_exit", sum(C1$delta))
  put("dgp_C1_sample_share_ref", mean(C1$referred)); put("dgp_C1_sample_mean_fee_usd", mean(C1$fee_usd))

  P2 <- P$C2
  for (k in c("n", "p_ref", "g0", "g_ref", "med_deal", "p_deal", "med_reg", "p_reg", "theta_fee", "fee_ref", "N_large", "R_mc"))
    put(paste0("dgp_C2_", k), P2[[k]])
  put("dgp_C2_deal_share_nonref", tr$C2$deal_share_nonref); put("dgp_C2_deal_share_ref", tr$C2$deal_share_ref)
  put("dgp_C2_pct_time_per_usd", exp(-P2$theta_fee) - 1)
  put("dgp_C2_reg_hazard_peak_month", P2$med_reg * (P2$p_reg - 1)^(1 / P2$p_reg))
  put("dgp_C2_pseudo_kappa", tr$C2$pseudo)
  for (i in 1:4) put(paste0("dgp_C2_pseudo_", pn[i]), tr$C2$pseudo_vec[i])
  for (i in 1:4) put(paste0("dgp_C2_pseudo_", pn[i], "_se_large"), tr$C2$large$se_sand[i])
  put("dgp_C2_quad_kappa", tr$C2$quad$kappa)
  for (i in 1:4) put(paste0("dgp_C2_quad_", pn[i]), tr$C2$quad$par[i])
  put("dgp_C2_pseudo_minus_quad_kappa", tr$C2$pseudo - tr$C2$quad$kappa)
  put("dgp_C2_large_share_cens", tr$C2$large$share_cens)
  put("dgp_C2_median_fee60", tr$C2$median_fee60); put("dgp_C2_share_cens_pop", tr$C2$share_cens_pop)
  hz <- tr$C2$hazard_fee60
  dh <- diff(hz$h)
  i_min <- which(dh[-1] > 0 & dh[-length(dh)] < 0)[1] + 1                 # first local minimum of the true hazard
  i_max <- which(dh[-1] < 0 & dh[-length(dh)] > 0)[1] + 1                 # the local maximum after it
  put("dgp_C2_hazard_fee60_min_month", hz$t[i_min]); put("dgp_C2_hazard_fee60_min", hz$h[i_min])
  put("dgp_C2_hazard_fee60_max_month", hz$t[i_max]); put("dgp_C2_hazard_fee60_max", hz$h[i_max])
  put("dgp_C2_sample_share_cens", 1 - mean(C2$delta)); put("dgp_C2_sample_share_deal", mean(C2$deal_true))

  P3 <- P$C3
  put("dgp_C3_M", P3$M)
  for (m in seq_len(P3$M)) put(paste0("dgp_C3_eta_true_", mm(m)), tr$C3$eta[m])

  P4 <- P$C4
  put("dgp_C4_n", P4$n); put("dgp_C4_R_mc", P4$R_mc); put("dgp_C4_kappa", P4$kappa)
  for (i in 1:3) put(sprintf("dgp_C4_gamma%d", i - 1), P4$gamma[i])
  for (i in 1:3) put(paste0("dgp_C4_", c("b0", "b_fee", "b_u")[i]), P4$beta[i])
  for (i in seq_along(P4$hh_sizes)) put(sprintf("dgp_C4_hh%d_share", P4$hh_sizes[i]), P4$p_hh[i])
  put("dgp_C4_mean_usage_pop", tr$C4$mean_usage); put("dgp_C4_var_u", tr$C4$var_u); put("dgp_C4_sd_u", sqrt(tr$C4$var_u))
  put("dgp_C4_hr_u", exp(P4$beta[3])); put("dgp_C4_pct_hazard_u", exp(P4$beta[3]) - 1)
  put("dgp_C4_sample_mean_usage", mean(C4$usage_wk1)); put("dgp_C4_sample_share_zero_usage", mean(C4$usage_wk1 == 0))
  put("dgp_C4_sample_share_cens", 1 - mean(C4$delta))

  ## ---- cens: block C1 --------------------------------------------------------------------------------
  b <- tr$C1$beta
  put("cens_n", c1$n); put("cens_n_exit", c1$n_exit); put("cens_share", c1$share_cens); put("cens_mean_t_observed", c1$mean_t_observed)
  for (nm in c("ols_observed", "ols_uncensored", "ols_latent")) {
    o <- c1[[nm]]
    pre <- paste0("cens_", nm)
    put_ci(paste0(pre, "_fee"), o$est, o$se, o$ci)
    put(paste0(pre, "_fee_covers_truth"), covers_num(o$ci, b)); put(paste0(pre, "_fee_over_truth"), o$est / b)
    put(paste0(pre, "_fee_minus_truth_in_se"), (o$est - b) / o$se)
    put(paste0(pre, "_kappa"), o$kappa); put(paste0(pre, "_kappa_se"), o$se_kappa)
    put_ci(paste0(pre, "_slope"), o$slope, o$slope_se, o$slope_ci)
    put(paste0(pre, "_s2"), o$s2)
    put(paste0(pre, "_b0"), o$beta[1]); put(paste0(pre, "_b_ref"), o$beta[3]); put(paste0(pre, "_b_ref_se"), o$se_beta[3])
  }
  put("cens_ols_uncensored_n", c1$ols_uncensored$n)
  m1 <- c1$mle
  put_ci("cens_mle_fee", m1$est, m1$se, m1$ci)
  put("cens_mle_fee_covers_truth", covers_num(m1$ci, b)); put("cens_mle_fee_minus_truth_in_se", (m1$est - b) / m1$se)
  put("cens_mle_fee_se_opg", c1$se_opg); put("cens_mle_fee_se_sandwich", c1$se_sand)
  put("cens_mle_fee_opg_lb", m1$ci_opg[1]); put("cens_mle_fee_opg_ub", m1$ci_opg[2])
  put("cens_mle_fee_sandwich_lb", m1$ci_sand[1]); put("cens_mle_fee_sandwich_ub", m1$ci_sand[2])
  for (i in 1:4) {
    put(paste0("cens_mle_", pn[i], "_est"), m1$par[i]); put(paste0("cens_mle_", pn[i], "_se_model"), m1$se_hess_vec[i])
    put(paste0("cens_mle_", pn[i], "_se_opg"), m1$se_opg_vec[i]); put(paste0("cens_mle_", pn[i], "_se_sandwich"), m1$se_sand_vec[i])
  }
  put_ci("cens_mle_ref", m1$ref$est, m1$ref$se, m1$ref$ci)
  put("cens_mle_ref_covers_truth", covers_num(m1$ref$ci, tr$C1$beta_vec[3]))
  put_ci("cens_mle_kappa", m1$kappa, m1$se_kappa, m1$ci_kappa)
  put("cens_mle_kappa_covers_truth", covers_num(m1$ci_kappa, tr$C1$shape))
  put("cens_mle_hr_fee_per_15usd", exp(15 * m1$est)); put("cens_mle_median_fee60", m1$median_fee60)
  put("cens_mle_loglik", -m1$nll); put("cens_mle_iter", m1$iter)
  put("cens_se_hess", c1$se_hess); put("cens_se_opg", c1$se_opg); put("cens_se_sand", c1$se_sand)
  put("cens_ratio_hess_opg", c1$se_hess / c1$se_opg); put("cens_ratio_sand_hess", c1$se_sand / c1$se_hess)
  put("cens_ols_observed_over_mle_fee", c1$ols_observed$est / m1$est)
  sv <- c1$survreg
  for (i in 1:3) {
    put(paste0("cens_survreg_aft_", c("b0", "b_fee", "b_ref")[i]), sv$coef[i])
    put(paste0("cens_survreg_aft_", c("b0", "b_fee", "b_ref")[i], "_se"), sv$se_coef[i])
  }
  put("cens_survreg_scale", sv$scale); put("cens_survreg_kappa", 1 / sv$scale); put("cens_survreg_loglik", sv$loglik)
  put("cens_survreg_ph_b_fee", sv$ph_par[2]); put("cens_survreg_ph_b_fee_se", sv$ph_se[2])
  put("cens_survreg_max_abs_diff", sv$max_abs_diff)
  tri <- c1$trinity
  for (k in c("W", "LR", "LM", "df", "crit_5pct", "LR_survreg")) put(paste0("cens_trinity_", k), tri[[k]])
  # p-values near 1e-15: their log10 (a 16-digit text of 1.29e-15 reads back one ulp off in jsonlite's parser)
  for (k in c("W", "LR", "LM")) put(paste0("cens_trinity_log10_p_", k), pchisq(tri[[k]], tri$df, lower.tail = FALSE, log.p = TRUE) / log(10))
  put("cens_trinity_W_minus_LR_rel", (tri$W - tri$LR) / tri$LR); put("cens_trinity_LM_minus_LR_rel", (tri$LM - tri$LR) / tri$LR)
  put("cens_trinity_max_gap_over_LR", max(abs(tri$W - tri$LR), abs(tri$LM - tri$LR)) / tri$LR)
  for (i in 1:3) put(paste0("cens_trinity_restricted_", pn[i]), tri$par_r[i])
  a <- c1$aft
  put("cens_aft_truth", a$truth)
  put_ci("cens_aft_ols_observed", a$ols_observed, a$ols_observed_se, a$ols_observed_ci)
  put_ci("cens_aft_ols_latent", a$ols_latent, a$ols_latent_se, a$ols_latent_ci)
  put_ci("cens_aft_mle", a$mle, a$mle_se, a$mle_ci)
  put("cens_aft_ols_observed_covers_truth", covers_num(a$ols_observed_ci, a$truth))
  put("cens_aft_mle_covers_truth", covers_num(a$mle_ci, a$truth))
  put("cens_aft_ols_observed_over_truth", a$ols_observed / a$truth)

  ## ---- qmle: block C2 -----------------------------------------------------------------------------------
  ks <- tr$C2$pseudo
  q <- c2$qmle
  put("qmle_n", c2$n); put("qmle_n_exit", c2$n_exit); put("qmle_share_cens", c2$share_cens); put("qmle_share_deal", c2$share_deal)
  put_ci("qmle_kappa", q$est, q$se, q$ci)
  put("qmle_kappa_se_model", c2$se_hess); put("qmle_kappa_se_opg", c2$se_opg); put("qmle_kappa_se_sandwich", c2$se_sand)
  put("qmle_kappa_model_lb", q$ci_model[1]); put("qmle_kappa_model_ub", q$ci_model[2])
  put("qmle_kappa_opg_lb", q$ci_opg[1]); put("qmle_kappa_opg_ub", q$ci_opg[2])
  put("qmle_kappa_covers_pseudo", c2$covers_sandwich); put("qmle_kappa_model_covers_pseudo", c2$covers_model)
  put("qmle_kappa_opg_covers_pseudo", c2$covers_opg)
  put("qmle_kappa_minus_pseudo", q$est - ks); put("qmle_kappa_minus_pseudo_in_se", (q$est - ks) / q$se)
  put("qmle_ratio_hess_opg", c2$ratio_hess_opg); put("qmle_ratio_sand_hess", c2$ratio_sand_hess)
  put("qmle_ratio_sand_opg", c2$se_sand / c2$se_opg)
  for (i in 1:4) {
    put(paste0("qmle_", pn[i], "_est"), q$par[i]); put(paste0("qmle_", pn[i], "_se_model"), q$se_model_vec[i])
    put(paste0("qmle_", pn[i], "_se_opg"), q$se_opg_vec[i]); put(paste0("qmle_", pn[i], "_se_sandwich"), q$se_sand_vec[i])
  }
  put("qmle_loglik", -q$nll)
  fe <- c2$fee
  put("qmle_fee_est", fe$est); put("qmle_fee_se_model", fe$se_model); put("qmle_fee_se_opg", fe$se_opg)
  put("qmle_fee_se_sandwich", fe$se_sandwich); put("qmle_fee_lb", fe$ci[1]); put("qmle_fee_ub", fe$ci[2])
  put("qmle_fee_model_lb", fe$ci_model[1]); put("qmle_fee_model_ub", fe$ci_model[2])
  put("qmle_fee_covers_pseudo", fe$covers_sandwich); put("qmle_fee_model_covers_pseudo", fe$covers_model)
  vu <- c2$vuong
  # the two-sided p-value is about 1e-317 (a subnormal double that does not survive JSON): its log10 instead
  put("qmle_vuong_stat", vu$stat); put("qmle_vuong_log10_p", (log(2) + pnorm(-abs(vu$stat), log.p = TRUE)) / log(10))
  put("qmle_vuong_mean", vu$mean); put("qmle_vuong_sd", vu$sd)
  put("qmle_vuong_sum", vu$sum); put("qmle_vuong_loglik_weibull", vu$loglik_weibull); put("qmle_vuong_loglik_loglogistic", vu$loglik_loglogistic)
  put("qmle_loglogistic_scale", vu$loglogistic_scale); put("qmle_loglogistic_shape", 1 / vu$loglogistic_scale)
  for (i in 1:3) put(paste0("qmle_loglogistic_aft_", c("b0", "b_fee", "b_ref")[i]), vu$loglogistic_coef[i])
  mc <- c2$mc
  for (k in c("R", "n", "cover_model", "cover_opg", "cover_sandwich", "mean_est", "sd_est", "mean_se_model", "mean_se_opg",
              "mean_se_sandwich", "share_ratio_gt_0p3", "mean_ratio_model_opg"))
    put(paste0("qmle_mc_", k), mc[[k]])
  put("qmle_mc_bias", mc$mean_est - ks); put("qmle_mc_sd_over_mean_se_model", mc$sd_est / mc$mean_se_model)
  put("qmle_mc_sd_over_mean_se_sandwich", mc$sd_est / mc$mean_se_sandwich)
  for (k in names(mc$fee)) put(paste0("qmle_mc_fee_", k), mc$fee[[k]])
  put("qmle_mcse_coverage", sqrt(0.95 * 0.05 / mc$R))

  ## ---- grouped: block C3 ------------------------------------------------------------------------------
  put("grouped_max_abs_diff", c3$max_abs_diff); put("grouped_nll_gap", c3$nll_gap)
  put("grouped_iter_newton", c3$iter_grouped); put("grouped_iter_glm", c3$iter_glm)
  put("grouped_n_rows", c3$n_rows); put("grouped_n_exits", c3$n_exits); put("grouped_n_subscribers", c3$n_subscribers)
  for (nm in c("fee", "ref")) {
    g <- c3[[nm]]
    put_ci(paste0("grouped_", nm), g$est, g$se, g$ci)
    put(paste0("grouped_", nm, "_se_glm"), g$se_glm); put(paste0("grouped_", nm, "_se_cluster"), g$se_cluster)
    put(paste0("grouped_", nm, "_covers_truth"), g$covers_truth)
  }
  put("grouped_ratio_cluster_glm_fee", c3$ratio_cluster_glm_fee); put("grouped_ratio_cluster_glm_ref", c3$ratio_cluster_glm_ref)
  put("grouped_fee_minus_continuous", c3$fee_minus_continuous)
  for (m in seq_along(c3$eta)) {
    put(paste0("grouped_eta_", mm(m), "_est"), c3$eta[m]); put(paste0("grouped_eta_", mm(m), "_se"), c3$eta_se[m])
    put(paste0("grouped_exits_", mm(m)), c3$exits_by_month[m]); put(paste0("grouped_at_risk_", mm(m)), c3$at_risk_by_month[m])
  }

  ## ---- twostep: block C4 ----------------------------------------------------------------------------------
  bu <- tr$C4$beta
  put("twostep_n", c4$n); put("twostep_share_cens", c4$share_cens); put("twostep_mean_usage", c4$mean_usage)
  put("twostep_share_zero_usage", c4$share_zero_usage)
  for (i in 1:3) { put(sprintf("twostep_gamma%d_est", i - 1), c4$gamma[i]); put(sprintf("twostep_gamma%d_se", i - 1), c4$gamma_se[i]) }
  for (i in 1:4) put(paste0("twostep_", c("b0", "b_fee", "b_u", "log_kappa")[i], "_est"), c4$par[i])
  put("twostep_kappa_est", exp(c4$par[4]))
  put("twostep_est", c4$est); put("twostep_se_naive", c4$se_naive); put("twostep_se_naive_sandwich", c4$se_naive_sand)
  put("twostep_se_corrected", c4$se_corrected); put("twostep_se_ratio", c4$se_ratio)
  put("twostep_naive_lb", c4$naive$ci[1]); put("twostep_naive_ub", c4$naive$ci[2]); put("twostep_naive_covers_truth", c4$naive$covers)
  put("twostep_corrected_lb", c4$corrected$ci[1]); put("twostep_corrected_ub", c4$corrected$ci[2])
  put("twostep_corrected_covers_truth", c4$corrected$covers)
  put_ci("twostep_oracle", c4$oracle$est, c4$oracle$se, c4$oracle$ci); put("twostep_oracle_covers_truth", c4$oracle$covers)
  put("twostep_hr_u_est", exp(c4$est)); put("twostep_pct_hazard_u_est", exp(c4$est) - 1)
  for (r in 1:4) for (s in 1:3) put(sprintf("twostep_F2_%s_gamma%d", c("b0", "b_fee", "b_u", "log_kappa")[r], s - 1), c4$F2[r, s])
  mc4 <- c4$mc
  for (k in c("R", "n", "cover_naive", "cover_naive_sand", "cover_corrected", "cover_oracle", "mean_est", "sd_est", "mean_se_naive",
              "mean_se_naive_sand", "mean_se_corrected", "mean_oracle", "sd_oracle", "mean_se_oracle", "share_corrected_ge_naive"))
    put(paste0("twostep_mc_", k), mc4[[k]])
  put("twostep_mc_bias", mc4$mean_est - bu); put("twostep_mc_sd_over_mean_se_naive", mc4$sd_est / mc4$mean_se_naive)
  put("twostep_mc_sd_over_mean_se_corrected", mc4$sd_est / mc4$mean_se_corrected)
  put("twostep_mcse_coverage", sqrt(0.95 * 0.05 / mc4$R))

  nums
}
