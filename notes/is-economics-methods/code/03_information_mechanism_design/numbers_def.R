# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers03(dir = "."): every quotable number of chapter 3 as a flat named list
# of finite numeric scalars (07_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth03.rds, vireo_V1.rds and res_auction.rds, res_menu.rds, res_len.rds, res_match.rds in `dir`.
#
# Name scheme: <prefix>_<block or object>_<field>. Prefixes: dgp (01_dgp.R: seeds, parameters, realised
# sample facts, simulation truth), v0 (block V0, the hand block), auc (02, block V1), menu (03, block V2),
# len (04, block V3), match (05, block V4). Conventions:
#   units               auction values, bids and revenues in USD per click (revenue: per auction); menu types
#                       theta in USD per thousand calls, q in thousand calls a month, t, U, profit in USD a month
#                       per customer; LEN amounts in thousand USD a month per moderator, *_usd in USD
#   fp, sp, *_res       first price, second price; *_res with the optimal reserve r* (auc_rev_reserve is the
#                       second price with r*, the card's name)
#   n2 ... n6           bidder count; tau0p5 = quantile level 0.5
#   mc_*                Monte Carlo over 1,000 replications (Ruling P23): *_mean, *_sd, and *_lo, *_hi for the
#                       2.5% and 97.5% quantiles of the replications
#   reg, irr            the regular and the bimodal (ironing) type density of block V2; two = the two-type menu
#   without, report,    the LEN contract on output alone, and with the user-report rate, a peer team's output,
#   peer, noise         a pure-noise metric
#   P, R (match_)       workers propose (P) or firms propose (R); ranks are positions on the agent's own true list
#   receiver_gain_*     match_receiver_gain_any (the card's field) counts profitable firm REPORTS, random plus
#                       targeted truncations, so a firm can count twice; match_receiver_gain_n_firms counts the
#                       distinct firms that gain, match_receiver_gain_random_n_firms those among the random ones
#   p                   decimal point written p (b = 0.05 -> b0p05)
#   flags               0/1 numbers (v0_gsp_truth_lef, match_stable_both, ...)
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json; fields that are NA by
# construction (the regular density has no bunch) are left out.
# Sources mech_fns.R (run from the chapter code directory, as 07_summary.R and tests/run_tests.R do) for value_dist.
# ------------------------------------------------------------------
source("mech_fns.R")                                         # value_dist: the V1 value quantiles
collect_numbers03 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 05 first")
    readRDS(p)
  }
  tr <- rd("truth03.rds"); V1 <- rd("vireo_V1.rds")
  a1 <- rd("res_auction.rds"); s2 <- rd("res_menu.rds"); k3 <- rd("res_len.rds"); m4 <- rd("res_match.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers03: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers03: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers03: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  pid <- function(x) gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))
  P <- tr$params
  q <- function(x, p) unname(stats::quantile(x, p, type = 7))

  ## ---- dgp: seeds and parameters ----------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  P1 <- P$V1; T1 <- tr$V1
  put("dgp_V1_F_mean", P1$F$mean); put("dgp_V1_F_sd", P1$F$sd); put("dgp_V1_F_lower", P1$F$lower); put("dgp_V1_F_upper", P1$F$upper)
  put("dgp_V1_n_min", min(P1$n_set)); put("dgp_V1_n_max", max(P1$n_set)); put("dgp_V1_n_auc_per_n", P1$n_auc_per_n)
  put("dgp_V1_n_auc", max(V1$auction)); put("dgp_V1_n_bids", nrow(V1))
  put("dgp_V1_h", P1$h); put("dgp_V1_tau_lo", min(P1$taus)); put("dgp_V1_tau_hi", max(P1$taus)); put("dgp_V1_tau_step", diff(P1$taus)[1])
  put("dgp_V1_n_tau", length(P1$taus)); put("dgp_V1_R_mc", P1$R_mc)
  put("dgp_V1_r_star", T1$r_star); put("dgp_V1_mean_value", T1$mean_value); put("dgp_V1_median_value", T1$median_value)
  put("dgp_V1_pdf_lower", T1$pdf_lower); put("dgp_V1_pdf_upper", T1$pdf_upper)
  for (i in seq_along(T1$rev_true$n)) {
    k <- T1$rev_true$n[i]
    put(sprintf("dgp_V1_rev_true_n%d", k), T1$rev_true$none[i]); put(sprintf("dgp_V1_rev_true_res_n%d", k), T1$rev_true$reserve[i])
    put(sprintf("dgp_V1_gain_true_n%d", k), T1$rev_true$reserve[i] - T1$rev_true$none[i])
    put(sprintf("dgp_V1_p_unsold_n%d", k), T1$p_unsold[i])
  }
  put("dgp_V1_rev_true_pooled", T1$rev_true_pooled); put("dgp_V1_rev_true_res_pooled", T1$rev_reserve_true_pooled)
  put("dgp_V1_gain_true_pooled", T1$rev_reserve_true_pooled - T1$rev_true_pooled)
  D1 <- value_dist(T1$F)
  for (tau in c(0.05, 0.25, 0.5, 0.75, 0.95)) put(paste0("dgp_V1_v_quantile_tau", pid(tau)), D1$quantile(tau))
  put("dgp_V1_sample_v_mean", mean(V1$v_true)); put("dgp_V1_sample_v_min", min(V1$v_true)); put("dgp_V1_sample_v_max", max(V1$v_true))
  put("dgp_V1_sample_share_below_r", mean(V1$v_true < T1$r_star))
  P2 <- P$V2
  for (k in c("theta_lo", "theta_hi", "n_grid", "cost2")) put(paste0("dgp_V2_", k), P2[[k]])
  put("dgp_V2_reg_mean", P2$regular$mean); put("dgp_V2_reg_sd", P2$regular$sd)
  for (j in 1:2) {
    put(sprintf("dgp_V2_irr_weight%d", j), P2$irregular$weights[j]); put(sprintf("dgp_V2_irr_mean%d", j), P2$irregular$means[j])
    put(sprintf("dgp_V2_irr_sd%d", j), P2$irregular$sds[j])
  }
  for (k in names(P2$two_type)) put(paste0("dgp_V2_two_", k), P2$two_type[[k]])
  P3 <- P$V3
  put("dgp_V3_A", P3$A); put("dgp_V3_c2", P3$c2); put("dgp_V3_sigma2", P3$sigma2); put("dgp_V3_sigma", sqrt(P3$sigma2))
  for (nm in names(P3$signals)) for (k in c("k", "var", "rho")) put(sprintf("dgp_V3_%s_%s", nm, k), P3$signals[[nm]][[k]])
  for (nm in names(P3$signals)) put(sprintf("dgp_V3_%s_sd", nm), sqrt(P3$signals[[nm]]$var))
  put("dgp_V3_N_sim", P3$N_sim); put("dgp_V3_R_sim", P3$R_sim)
  P4 <- P$V4
  for (k in c("n_workers", "n_firms", "sd_common", "sd_idio", "n_trunc_prop", "n_swap_prop", "n_trunc_recv")) put(paste0("dgp_V4_", k), P4[[k]])
  put("dgp_V4_common_share", P4$sd_common^2 / (P4$sd_common^2 + P4$sd_idio^2))

  ## ---- v0: the hand block -------------------------------------------------------------------
  V0 <- tr$V0
  scal <- names(V0)[vapply(V0, function(z) (is.numeric(z) || is.logical(z)) && length(z) == 1L, logical(1))]
  for (k in scal) put(paste0("v0_", k), V0[[k]])
  for (i in 1:3) put(sprintf("v0_v%d", i), V0$v[i])
  for (i in 1:2) put(sprintf("v0_ctr%d", i), V0$ctr[i])
  for (i in seq_along(V0$re_n)) {
    put(sprintf("v0_re_spa_n%d", V0$re_n[i]), V0$re_spa[i]); put(sprintf("v0_re_fpa_n%d", V0$re_n[i]), V0$re_fpa[i])
  }
  for (i in seq_along(V0$shading_n)) put(sprintf("v0_shading_n%d", V0$shading_n[i]), V0$shading[i])
  for (i in seq_along(V0$akerlof_p)) {
    put(paste0("v0_akerlof_p_", pid(V0$akerlof_p[i])), V0$akerlof_p[i])
    put(paste0("v0_akerlof_wtp_at_", pid(V0$akerlof_p[i])), V0$akerlof_wtp[i])
  }
  for (i in seq_along(V0$ct_b)) {
    put(paste0("v0_ct_b_", pid(V0$ct_b[i])), V0$ct_b[i]); put(paste0("v0_ct_N_b", pid(V0$ct_b[i])), V0$ct_N[i])
  }
  for (i in seq_along(V0$ct_partition_b0p05)) put(sprintf("v0_ct_partition_b0p05_a%d", i - 1L), V0$ct_partition_b0p05[i])
  for (i in 1:3) {
    put(sprintf("v0_da_worker_opt_w%d", i), V0$da_worker_opt[i]); put(sprintf("v0_da_firm_opt_w%d", i), V0$da_firm_opt[i])
  }
  wr <- partner_rank_list(V0$da_worker_opt, V0$da_worker_pref); fr <- partner_rank_list(V0$da_firm_opt, V0$da_worker_pref)
  for (i in 1:3) { put(sprintf("v0_da_worker_opt_rank_w%d", i), wr[i]); put(sprintf("v0_da_firm_opt_rank_w%d", i), fr[i]) }

  ## ---- auc: block V1 (02_auction.R) -------------------------------------------------------------
  for (k in c("N_auc", "r_star", "rev_fp", "rev_sp", "rev_fp_reserve", "rev_reserve", "diff", "se_diff", "diff_res",
              "se_diff_res", "gain", "se_gain", "se_rev_sp", "se_rev_fp", "share_unsold",
              "gpv_med_relerr", "gpv_max_relerr", "gpv_mean_relerr"))
    put(paste0("auc_", k), a1[[k]])
  put("auc_diff_over_se", a1$diff / a1$se_diff); put("auc_diff_res_over_se", a1$diff_res / a1$se_diff_res)
  put("auc_gain_over_se", a1$gain / a1$se_gain)
  put("auc_rev_sp_minus_truth", a1$rev_sp - T1$rev_true_pooled); put("auc_rev_fp_minus_truth", a1$rev_fp - T1$rev_true_pooled)
  put("auc_gain_pct_of_rev_sp", 100 * a1$gain / a1$rev_sp)
  bn <- as.data.frame(a1$by_n)       # data.frames: the same indexing whether or not data.table is attached
  for (i in seq_len(nrow(bn))) {
    pre <- sprintf("auc_n%d_", bn$n[i])
    for (k in c("rev_fp", "rev_sp", "rev_fp_res", "rev_sp_res", "se_diff", "se_diff_res", "gain_reserve", "se_gain",
                "share_unsold", "shading", "bid_over_value"))
      put(paste0(pre, k), bn[[k]][i])
    put(paste0(pre, "diff"), bn$rev_fp[i] - bn$rev_sp[i])
    put(paste0(pre, "diff_over_se"), (bn$rev_fp[i] - bn$rev_sp[i]) / bn$se_diff[i])
    put(paste0(pre, "gain_pct_of_rev_sp"), 100 * bn$gain_reserve[i] / bn$rev_sp[i])
  }
  gb <- as.data.frame(a1$gpv_by_n)
  for (i in seq_len(nrow(gb))) {
    pre <- sprintf("auc_gpv_n%d_", gb$n[i])
    put(paste0(pre, "med_relerr"), gb$med_relerr[i]); put(paste0(pre, "max_relerr"), gb$max_relerr[i])
    put(paste0(pre, "mean_relerr"), gb$mean_relerr[i])
  }
  gp <- as.data.frame(a1$gpv)
  for (k in c(2L, 4L, 6L)) for (tau in c(0.25, 0.5, 0.75)) {
    r <- gp[gp$n == k & abs(gp$tau - tau) < 1e-9, ]
    pre <- sprintf("auc_gpv_n%d_tau%s_", k, pid(tau))
    put(paste0(pre, "b"), r$b); put(paste0(pre, "v_hat"), r$v_hat); put(paste0(pre, "v_true"), r$v_true)
    put(paste0(pre, "relerr"), r$relerr)
  }
  mc <- a1$mc
  put("auc_mc_R", mc$R)
  for (k in c("diff", "diff_res", "gain")) {
    put(paste0("auc_mc_", k, "_mean"), mc[[paste0(k, "_mean")]]); put(paste0("auc_mc_", k, "_sd"), mc[[paste0(k, "_sd")]])
    put(paste0("auc_mc_", k, "_lo"), mc[[paste0(k, "_range")]][1]); put(paste0("auc_mc_", k, "_hi"), mc[[paste0(k, "_range")]][2])
  }
  for (k in c("rev_fp", "rev_sp", "rev_fp_res", "rev_sp_res")) {
    put(paste0("auc_mc_", k, "_mean"), mc[[paste0(k, "_mean")]])
    put(paste0("auc_mc_", k, "_lo"), mc[[paste0(k, "_range")]][1]); put(paste0("auc_mc_", k, "_hi"), mc[[paste0(k, "_range")]][2])
  }
  for (k in c("se_diff_mean", "share_within_3se", "share_within_3se_res", "share_gain_positive", "gpv_relerr_mean"))
    put(paste0("auc_mc_", k), mc[[k]])
  put("auc_mc_gpv_relerr_lo", mc$gpv_relerr_range[1]); put("auc_mc_gpv_relerr_hi", mc$gpv_relerr_range[2])
  put("auc_mc_sd_over_se", mc$diff_sd / a1$se_diff)
  gbn <- as.data.frame(mc$gain_by_n)
  for (i in seq_len(nrow(gbn))) {
    g <- gbn[i, ]; pre <- sprintf("auc_mc_n%d_", g$n)
    put(paste0(pre, "gain_lo"), g$gain_lo); put(paste0(pre, "gain_hi"), g$gain_hi)
    put(paste0(pre, "diff_lo"), g$diff_lo); put(paste0(pre, "diff_hi"), g$diff_hi)
  }
  band <- as.data.frame(mc$gpv_band)
  put("auc_mc_gpv_band_cover_share", mean(band$lo <= band$v_true & band$v_true <= band$hi))
  put("auc_mc_gpv_max_rel_bias", max(abs(band$mean - band$v_true) / band$v_true))
  for (k in c(2L, 4L, 6L)) {
    r <- band[band$n == k & abs(band$tau - 0.5) < 1e-9, ]
    put(sprintf("auc_mc_gpv_n%d_tau0p5_lo", k), r$lo); put(sprintf("auc_mc_gpv_n%d_tau0p5_hi", k), r$hi)
  }

  ## ---- menu: block V2 (03_menu.R) ----------------------------------------------------------------
  for (nm in c("regular", "irregular")) {
    s <- s2[[nm]]; pre <- paste0("menu_", c(regular = "reg", irregular = "irr")[[nm]], "_")
    for (k in c("profit", "surplus", "rent", "surplus_efficient", "dwl", "theta_exclusion", "theta_psi0", "share_served",
                "q_top", "q_eff_top", "ic_gain", "n_bunches", "bunch_width", "bunch_q", "bunch_mass", "mean_theta",
                "profit_naive_if_truthful", "ic_gain_naive", "rent_share", "profit_share_eff"))
      put(paste0(pre, k), s[[k]])
    for (k in c("bunch_lo", "bunch_hi", "bunch_psi_bar", "psi_fall_lo", "psi_fall_hi")) if (is.finite(s[[k]])) put(paste0(pre, k), s[[k]])
    put(paste0(pre, "psi_monotone"), s$psi_monotone)
    for (tt in names(s$q_at)) {
      put(paste0(pre, "q_at_", tt), s$q_at[[tt]]); put(paste0(pre, "q_eff_at_", tt), s$q_eff_at[[tt]])
      put(paste0(pre, "t_at_", tt), s$t_at[[tt]]); put(paste0(pre, "U_at_", tt), s$U_at[[tt]])
      put(paste0(pre, "distortion_at_", tt), s$q_eff_at[[tt]] - s$q_at[[tt]])
    }
  }
  tw <- s2$two_type
  for (k in c("q_L", "q_H", "t_L", "t_H", "rent_H", "profit_L", "profit_H", "profit", "q_L_fb", "q_H_fb", "t_L_fb", "t_H_fb", "profit_fb"))
    put(paste0("menu_two_", k), tw[[k]])
  put("menu_two_profit_loss", tw$profit_fb - tw$profit)

  ## ---- len: block V3 (04_contract.R) ---------------------------------------------------------------
  cs <- as.data.frame(k3$cases)
  for (i in seq_len(nrow(cs))) for (k in c("sigma_s2", "beta", "b_x", "b_signal", "effort", "effort_cost", "risk_premium",
                                            "profit", "pay_sd", "w0")) {
    put(sprintf("len_%s_%s", cs$case[i], k), cs[[k]][i])
    if (k %in% c("profit", "risk_premium", "effort_cost", "w0")) put(sprintf("len_%s_%s_usd", cs$case[i], k), 1000 * cs[[k]][i])
  }
  put("len_b_without", k3$b_without); put("len_b_with_signal", k3$b_with_signal)
  put("len_profit_without", k3$profit_without); put("len_profit_with_signal", k3$profit_with_signal)
  put("len_profit_gain", k3$profit_with_signal - k3$profit_without)
  put("len_profit_gain_usd", 1000 * (k3$profit_with_signal - k3$profit_without))
  put("len_profit_gain_pct", 100 * (k3$profit_with_signal / k3$profit_without - 1))
  put("len_peer_profit_gain_usd", 1000 * (k3$profit_peer - k3$profit_without))
  sm <- k3$sim
  for (k in c("R", "N", "gain_mean", "gain_sd", "share_gain_positive", "profit_without_mean", "profit_with_mean",
              "ce_without_mean", "ce_with_mean"))
    put(paste0("len_sim_", k), sm[[k]])
  for (k in c("gain", "profit_without", "profit_with", "ce_without", "ce_with")) {
    put(paste0("len_sim_", k, "_lo"), sm[[paste0(k, "_range")]][1]); put(paste0("len_sim_", k, "_hi"), sm[[paste0(k, "_range")]][2])
  }
  put("len_sim_gain_mean_usd", 1000 * sm$gain_mean)
  put("len_sim_gain_lo_usd", 1000 * sm$gain_range[1]); put("len_sim_gain_hi_usd", 1000 * sm$gain_range[2])

  ## ---- match: block V4 (05_matching.R) -------------------------------------------------------------
  for (k in c("n_workers", "n_firms", "stable_P", "stable_R", "stable_both", "n_diff_partner",
              "share_worker_top1_P", "share_worker_top1_R", "share_firm_top1_P", "share_firm_top1_R",
              "n_prop_trunc", "n_prop_swap", "n_prop_misreports", "proposer_gain_any", "n_prop_changed", "n_prop_worse",
              "prop_max_gain", "n_recv_random", "receiver_gain_random", "receiver_gain_targeted",
              "targeted_get_receiver_optimal", "targeted_mean_gain", "targeted_max_gain", "random_recv_worse",
              "receiver_gain_any", "receiver_gain_n_firms", "receiver_gain_random_n_firms"))
    put(paste0("match_", k), m4[[k]])
  put("match_share_diff_partner", m4$n_diff_partner / m4$n_workers)
  for (side in c("worker", "firm")) for (d in c("P", "R")) {
    r <- m4[[sprintf("%s_rank_%s", side, d)]]
    put(sprintf("match_%s_rank_%s_mean", side, d), mean(r)); put(sprintf("match_%s_rank_%s_median", side, d), stats::median(r))
  }
  put("match_worker_rank_gain_P", mean(m4$worker_rank_R) - mean(m4$worker_rank_P))
  put("match_firm_rank_gain_R", mean(m4$firm_rank_P) - mean(m4$firm_rank_R))

  nums
}

# helper used above
partner_rank_list <- function(match, prefs) vapply(seq_along(prefs), function(i) as.numeric(match(match[i], prefs[[i]])), numeric(1))
