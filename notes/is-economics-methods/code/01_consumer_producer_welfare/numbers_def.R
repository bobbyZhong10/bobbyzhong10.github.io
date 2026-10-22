# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers01(dir = "."): every quotable number of chapter 1 as a flat named list
# of finite numeric scalars (06_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth01.rds, plover_W1.rds and res_welfare.rds, res_rp.rds, res_risk.rds in `dir`.
#
# Name scheme: <prefix>_<block or object>_<field>. Prefixes: dgp (01_dgp.R: seeds, parameters, realised
# sample facts, simulation truth), w0 (block W0, the hand block), welfare (02, blocks W1 and W1g), rp (03,
# block W2), risk (04, block W3). Conventions:
#   loss-positive welfare   cv, ev, dcs are losses (positive for a price rise), in USD per subscriber per month
#   *_mean                  average over the W1 (or W1g) subscribers; truth fields are dgp_W1_*, analyst
#                           fields welfare_*
#   cv12, *_out             the loss from raising p1 from 2 to 12 USD (Ruling P7; out of the 2 to 4 support)
#   eta                     income elasticity of premium demand at p1 = 2 USD (E_1w in the text)
#   rep                     representative consumer: dgp_W1_cv_rep_pop is the population LES with b = mean beta,
#                           welfare_cv_rep the analyst's (income-term form at mean income)
#   forms                   linear, loglinear, income (fit_demand_forms in welfare_fns.R); the headline
#                           welfare_dcs, welfare_willig, welfare_vartia, welfare_cv_rep are the income-term form
#   check                   Ruling P13: the Willig formula against exact CVs (welfare_check_*)
#   p2p5, rho0p5            decimal point written p (p1 = 2.5 USD, rho = 0.5)
#   flags                   0/1 numbers (w0_garp_pass, risk_rabin_bet_rejected_at_A_lo, ...)
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json.
# ------------------------------------------------------------------
collect_numbers01 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 04 first")
    readRDS(p)
  }
  tr <- rd("truth01.rds"); W1 <- rd("plover_W1.rds")
  wf <- rd("res_welfare.rds"); rp <- rd("res_rp.rds"); rk <- rd("res_risk.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers01: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers01: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers01: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  pid <- function(x) gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))
  P <- tr$params

  ## ---- dgp: seeds and parameters ----------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  P1 <- P$W1
  put("dgp_W1_N", P1$N); put("dgp_W1_n_per_arm", P1$N / length(P1$p_grid)); put("dgp_W1_n_arms", length(P1$p_grid))
  put("dgp_W1_p_lo", min(P1$p_grid)); put("dgp_W1_p_hi", max(P1$p_grid)); put("dgp_W1_p_step", diff(P1$p_grid)[1])
  for (k in c("p0", "p1", "p_out", "p2", "w_median", "w_sdlog", "beta_mean", "beta_phi")) put(paste0("dgp_W1_", k), P1[[k]])
  put("dgp_W1_p_out_over_p_lo", P1$p_out / min(P1$p_grid))
  put("dgp_W1_beta_shape1", P1$beta_phi * P1$beta_mean); put("dgp_W1_beta_shape2", P1$beta_phi * (1 - P1$beta_mean))
  put("dgp_W1_beta_sd", sqrt(P1$beta_mean * (1 - P1$beta_mean) / (P1$beta_phi + 1)))
  put("dgp_W1_w_mean_pop", P1$w_median * exp(P1$w_sdlog^2 / 2))
  put("dgp_W1_g1_lo", P1$g1_range[1]); put("dgp_W1_g1_hi", P1$g1_range[2])
  put("dgp_W1_g2_lo", P1$g2_range[1]); put("dgp_W1_g2_hi", P1$g2_range[2])
  # realised sample (W1; W1g has the same incomes, subsistence quantities and prices)
  q <- function(x, p) unname(stats::quantile(x, p, type = 7))
  put("dgp_W1_sample_w_mean", mean(W1$w)); put("dgp_W1_sample_w_median", stats::median(W1$w))
  put("dgp_W1_sample_w_min", min(W1$w)); put("dgp_W1_sample_w_max", max(W1$w))
  put("dgp_W1_sample_w_p10", q(W1$w, 0.1)); put("dgp_W1_sample_w_p90", q(W1$w, 0.9))
  put("dgp_W1_sample_beta_mean", mean(W1$beta_true)); put("dgp_W1_sample_beta_sd", stats::sd(W1$beta_true))
  put("dgp_W1_sample_beta_min", min(W1$beta_true)); put("dgp_W1_sample_beta_max", max(W1$beta_true))
  put("dgp_W1_sample_cor_w_beta", stats::cor(W1$w, W1$beta_true))
  put("dgp_W1_sample_x1_mean", mean(W1$x1)); put("dgp_W1_sample_spend_share_mean", mean(W1$p1 * W1$x1 / W1$w))
  put("dgp_W1_sample_min_supernumerary_at_p_out", min(W1$w - P1$p_out * W1$g1_true - P1$p2 * W1$g2_true))
  # truth, W1 and W1g
  for (v in c("W1", "W1g")) {
    t <- tr[[v]]; pre <- paste0("dgp_", v, "_")
    for (k in c("cv_mean", "ev_mean", "dcs_mean", "eta_mean", "willig_mean", "cv12_mean", "ev12_mean", "dcs12_mean",
                "cv_rep_pop", "ev_rep_pop", "cv12_rep_pop", "cov_x0_beta", "cov_spend0_beta", "cov_w_beta", "gap_second_order"))
      put(paste0(pre, k), t[[k]])
    for (k in c("b", "g1", "g2", "wbar")) put(paste0(pre, "rep_", k), t$rep[[k]])
    put(paste0(pre, "cv_mean_minus_rep_pop"), t$cv_mean - t$cv_rep_pop)
    put(paste0(pre, "cv_minus_dcs_mean"), t$cv_mean - t$dcs_mean)
    put(paste0(pre, "dcs_minus_ev_mean"), t$dcs_mean - t$ev_mean)
    put(paste0(pre, "cv_minus_dcs_rel"), (t$cv_mean - t$dcs_mean) / t$cv_mean)
    put(paste0(pre, "dcs_over_w_mean"), mean(t$dcs_i / W1$w))
    put(paste0(pre, "share_ev_lt_dcs_lt_cv"), mean(t$ev_i < t$dcs_i & t$dcs_i < t$cv_i))
  }
  P2 <- P$W2
  for (k in c("n_users", "n_budgets", "share_rational", "m", "noise_sd", "n_bronars_per_user")) put(paste0("dgp_W2_", k), P2[[k]])
  put("dgp_W2_n_rational", sum(tr$W2$type == "rational")); put("dgp_W2_n_noisy", sum(tr$W2$type == "noisy"))
  put("dgp_W2_price_lo", P2$price_range[1]); put("dgp_W2_price_hi", P2$price_range[2])
  put("dgp_W2_alpha_lo", P2$alpha_range[1]); put("dgp_W2_alpha_hi", P2$alpha_range[2])
  put("dgp_W2_rho_lo", P2$rho_range[1]); put("dgp_W2_rho_hi", P2$rho_range[2])
  put("dgp_W2_sigma_lo", 1 / (1 - P2$rho_range[1])); put("dgp_W2_sigma_hi", 1 / (1 - P2$rho_range[2]))   # elasticity of substitution
  P3 <- P$W3
  for (k in c("loss", "prob", "buy_at", "decline_at", "rabin_loss", "rabin_gain", "big_loss")) put(paste0("dgp_W3_", k), P3[[k]])
  for (pr in P3$prices) put(paste0("dgp_W3_price_", pid(pr)), pr)

  ## ---- w0: the hand block -------------------------------------------------------------------
  W0 <- tr$W0
  scal <- names(W0)[vapply(W0, function(z) is.numeric(z) && length(z) == 1L, logical(1))]
  for (k in scal) put(paste0("w0_", k), W0[[k]])
  for (i in 1:2) for (j in 1:2) put(sprintf("w0_S%d%d", i, j), W0$S[i, j])
  put("w0_S_eig1", W0$S_eig[1]); put("w0_S_eig2", W0$S_eig[2])
  for (j in seq_along(W0$crra_rho)) put(paste0("w0_crra_rp_local_rho", pid(W0$crra_rho[j])), W0$crra_rp_local[j])
  put("w0_garp_pass", as.numeric(W0$garp_pass))
  put("w0_cv_minus_dcs", W0$cv - W0$dcs); put("w0_dcs_minus_ev", W0$dcs - W0$ev)
  put("w0_cv_minus_ev", W0$cv - W0$ev); put("w0_cv_minus_ev_over_dcs", (W0$cv - W0$ev) / W0$dcs)   # section 1: CV - EV and its base
  put("w0_h1_drop_subst", W0$x1_0 - W0$h1_p1_u0); put("w0_h1_drop_income", W0$h1_p1_u0 - W0$x1_1)

  ## ---- welfare: blocks W1 and W1g (02_welfare.R) -------------------------------------------------
  cv_true <- tr$W1$cv_mean; cvg_true <- tr$W1g$cv_mean
  for (k in c("dcs", "willig", "vartia", "ev", "willig_ev", "cv_rep", "cv_rep_g", "spread_in_range", "spread_outside"))
    put(paste0("welfare_", k), wf[[k]])
  put("welfare_dcs_minus_cv_true", wf$dcs - cv_true)
  put("welfare_willig_minus_cv_true", wf$willig - cv_true)
  put("welfare_vartia_minus_cv_true", wf$vartia - cv_true)
  put("welfare_willig_ratio", abs(wf$willig - cv_true) / abs(wf$dcs - cv_true))
  put("welfare_share_cs_error_removed", 1 - abs(wf$willig - cv_true) / abs(wf$dcs - cv_true))
  put("welfare_cv_rep_minus_cv_true", wf$cv_rep - cv_true)
  put("welfare_cv_rep_rel_err", (wf$cv_rep - cv_true) / cv_true)
  put("welfare_cv_rep_g_minus_cv_true", wf$cv_rep_g - cvg_true)
  put("welfare_cv_rep_g_rel_err", (wf$cv_rep_g - cvg_true) / cvg_true)
  put("welfare_rep_bias_ratio", abs(wf$cv_rep - cv_true) / abs(wf$cv_rep_g - cvg_true))
  put("welfare_cv_rep_minus_rep_pop", wf$cv_rep - tr$W1$cv_rep_pop)
  for (v in c("W1", "W1g")) {
    r <- wf[[v]]; t <- tr[[v]]; pre <- paste0("welfare_", v, "_")
    put(paste0(pre, "spread_in_range"), r$spread_in_range); put(paste0(pre, "spread_outside"), r$spread_outside)
    put(paste0(pre, "wbar"), r$wbar)
    for (fn in names(r$forms)) {
      f <- r$forms[[fn]]; fp <- paste0(pre, fn, "_")
      for (k in c("dcs", "willig", "willig_ev", "cv", "ev", "cv_rep", "cv_out", "dcs_out", "eta", "r2",
                  "x1_at_p0_mean", "x1_at_pout_mean"))
        put(paste0(fp, k), f[[k]])
      for (cn in names(f$coef)) put(paste0(fp, "coef_", cn), f$coef[[cn]])
      put(paste0(fp, "cv_rel_err"), (f$cv - t$cv_mean) / t$cv_mean)
      put(paste0(fp, "cv_out_rel_err"), (f$cv_out - t$cv12_mean) / t$cv12_mean)
      put(paste0(fp, "willig_ratio"), abs(f$willig - t$cv_mean) / abs(f$dcs - t$cv_mean))
    }
    put(paste0(pre, "loglinear_smear"), r$forms$loglinear$smear)
    put(paste0(pre, "linear_choke_at_wbar"), r$forms$linear$choke_at_wbar)
    for (i in seq_len(nrow(r$arms))) {
      ap <- paste0(pre, "arm_p", pid(r$arms$p1[i]), "_")
      put(paste0(ap, "x1_mean"), r$arms$x1_mean[i]); put(paste0(ap, "n"), r$arms$n[i]); put(paste0(ap, "w_mean"), r$arms$w_mean[i])
    }
  }
  # Ruling P13: the Willig formula against exact CVs
  c0 <- wf$willig_check$W0
  for (k in c("willig", "willig_minus", "willig_coef1", "share_removed", "share_removed_ev")) put(paste0("welfare_check_W0_", k), c0[[k]])
  for (v in c("W1", "W1g")) {
    cc <- wf$willig_check[[v]]
    for (k in c("willig_mean", "willig_coef1_mean", "willig_minus_mean", "share_removed", "share_removed_ev",
                "mean_abs_err_dcs", "mean_abs_err_willig", "max_rel_err_willig", "min_share_removed_i"))
      put(paste0("welfare_check_", v, "_", k), cc[[k]])
    put(paste0("welfare_check_", v, "_all_closer"), as.numeric(cc$all_closer))
  }

  ## ---- rp: block W2 (03_revealed_pref.R) ---------------------------------------------------------
  for (k in c("n_users", "n_rational", "n_noisy", "n_bronars", "n_budgets", "pass_rational", "ccei_rational_min",
              "pass_noisy", "ccei_noisy_mean", "ccei_noisy_median", "ccei_noisy_min", "ccei_noisy_p10",
              "share_noisy_ccei_ge_095", "pass_bronars", "power", "ccei_bronars_mean", "ccei_bronars_median",
              "ccei_bronars_min", "ccei_bronars_p10", "share_bronars_ccei_ge_095", "pass_all_users",
              "n_viol_noisy_mean_failers"))
    put(paste0("rp_", k), rp[[k]])

  ## ---- risk: block W3 (04_risk.R) ----------------------------------------------------------------
  put("risk_expected_loss", rk$expected_loss)
  for (j in seq_along(rk$prices)) {
    put(paste0("risk_A_star_price", pid(rk$prices[j])), rk$A_star[[j]])
    put(paste0("risk_markup_price", pid(rk$prices[j])), rk$markup_over_expected_loss[j])
  }
  put("risk_A_lo", rk$A_interval[1]); put("risk_A_hi", rk$A_interval[2])
  for (k in c("A_rabin", "A_big", "loss_never_accepted_at_A_rabin", "loss_never_accepted_at_A_lo",
              "loss_never_accepted_at_A_hi", "rabin_bet_rejected_at_A_lo", "rabin_bet_rejected_at_A_hi"))
    put(paste0("risk_", k), rk[[k]])

  nums
}
