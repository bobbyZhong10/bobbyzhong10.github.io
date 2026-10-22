# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers02(dir = "."): every quotable number of chapter 2 as a flat named list
# of finite numeric scalars (07_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth02.rds, osprey_G1.rds and res_qre.rds, res_gp.rds, res_learn.rds, res_gg.rds in `dir`.
#
# Name scheme: <prefix>_<block or object>_<field>. Prefixes: dgp (01_dgp.R: seeds, parameters, realised
# sample facts, simulation truth), g0 (block G0, the hand block), qre (02, block G1), gp (03, block G2),
# learn (04, block G3), gg (05, block G4). Conventions:
#   payoffs              G0 and G3 in million USD a month; G1 in thousand USD a month per local market;
#                        G2 prices in USD per unit, quantities in thousand units, profits in thousand USD a month
#   u11, u22             G1's payoff parameters: Osprey's gain from S1 over S2 when Boreal is on S1 (truth 3)
#                        and from S2 over S1 when Boreal is on S2 (truth 1); chi is the QRE precision
#   unique, pooled,      the three G1 inversions (02_qre.R header): draw_* plain OLS on the displayed draw,
#   by_market            mc_* over 1,000 replications (mean, 2.5% to 97.5% range, sd); no standard errors (P23)
#   low_, fixed_strategy G2 at delta_low = 0.5: low_* the re-optimised cartel, fixed_strategy_* the delta = 0.95
#                        strategies held fixed (03_green_porter.R header)
#   s<sigma>             a noise level, decimal point written p (s0p1 = 0.1, s0p05 = 0.05); x<value> a G1 shifter,
#                        minus sign written m (xm0p2 = -0.2)
#   value_per_period     (1 - delta) Vc, the collusive value as a per-month annuity
#   cell11 .. cell22     cells of a 2x2 game, [Osprey's action, Boreal's action]
#   flags                0/1 numbers (g0_pd_high_dominated, gp_osd_ok, ...)
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json.
# ------------------------------------------------------------------
collect_numbers02 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 05 first")
    readRDS(p)
  }
  tr <- rd("truth02.rds"); G1 <- rd("osprey_G1.rds")
  qr <- rd("res_qre.rds"); gp <- rd("res_gp.rds"); ln <- rd("res_learn.rds"); gg <- rd("res_gg.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers02: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers02: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers02: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  pid <- function(x) gsub("-", "m", gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)))
  cells <- function(pre, M) for (i in 1:2) for (j in 1:2) put(sprintf("%s_cell%d%d", pre, i, j), M[i, j])
  P <- tr$params

  ## ---- dgp: seeds ------------------------------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])

  ## ---- dgp: G1 ----------------------------------------------------------------------------------
  P1 <- P$G1
  put("dgp_G1_n_markets", nrow(G1))
  put("dgp_G1_n_markets_unique", P1$n_markets_unique); put("dgp_G1_n_markets_high", P1$n_markets_high)
  put("dgp_G1_n_dec", P1$n_dec); put("dgp_G1_kappa", P1$kappa)
  put("dgp_G1_chi_unique", P1$lambda[["unique"]]); put("dgp_G1_chi_high", P1$lambda[["high"]])
  put("dgp_G1_n_x_unique", length(P1$x_unique)); put("dgp_G1_n_x_high", length(P1$x_high))
  put("dgp_G1_x_unique_lo", min(P1$x_unique)); put("dgp_G1_x_unique_hi", max(P1$x_unique))
  put("dgp_G1_x_high_lo", min(P1$x_high)); put("dgp_G1_x_high_hi", max(P1$x_high))
  put("dgp_G1_x_step", diff(P1$x_unique)[1])
  put("dgp_G1_markets_per_x_unique", P1$n_markets_unique / length(P1$x_unique))
  put("dgp_G1_markets_per_x_high", P1$n_markets_high / length(P1$x_high))
  put("dgp_G1_decisions_unique", P1$n_markets_unique * P1$n_dec); put("dgp_G1_decisions_high", P1$n_markets_high * P1$n_dec)
  cells("dgp_G1_A0", P1$A0); cells("dgp_G1_B0", P1$B0)
  put("dgp_G1_payoff_u11", tr$G1$payoff); put("dgp_G1_payoff_u22", tr$G1$payoff_u22)
  put("dgp_G1_payoff_slope", tr$G1$payoff + tr$G1$payoff_u22)               # how Osprey's payoff difference moves with q
  put("dgp_G1_bifurcation_chi", tr$G1$bifurcation[["lambda"]])
  put("dgp_G1_bifurcation_p", tr$G1$bifurcation[["p"]]); put("dgp_G1_bifurcation_q", tr$G1$bifurcation[["q"]])
  qx <- tr$G1$qre_x0                     # the x = 0 QRE at chi = 2 are the dgp_G1_high_x0_* keys below
  put("dgp_G1_x0_unique_p", qx$unique[1, "p"]); put("dgp_G1_x0_unique_q", qx$unique[1, "q"])
  eq <- tr$G1$eq
  eu <- eq[eq$block == "unique", ]
  put("dgp_G1_unique_p_min", min(eu$p)); put("dgp_G1_unique_p_max", max(eu$p))
  put("dgp_G1_unique_q_min", min(eu$q)); put("dgp_G1_unique_q_max", max(eu$q))
  for (x in P1$x_high) {
    e <- eq[eq$block == "high" & abs(eq$x - x) < 1e-12, ]
    for (role in c("S1", "S2", "middle")) {
      put(sprintf("dgp_G1_high_x%s_%s_p", pid(x), role), e$p[e$role == role])
      put(sprintf("dgp_G1_high_x%s_%s_q", pid(x), role), e$q[e$role == role])
    }
    put(sprintf("dgp_G1_high_x%s_prob_S1_first", pid(x)), stats::plogis(P1$kappa * x))
  }
  put("dgp_G1_share_S1_first", tr$G1$share_S1_first)
  gu <- G1[G1$block == "unique", ]; gh <- G1[G1$block == "high", ]
  put("dgp_G1_freq_osprey_unique", sum(gu$k_osprey) / sum(gu$n_dec)); put("dgp_G1_freq_boreal_unique", sum(gu$k_boreal) / sum(gu$n_dec))
  put("dgp_G1_freq_osprey_high", sum(gh$k_osprey) / sum(gh$n_dec)); put("dgp_G1_freq_boreal_high", sum(gh$k_boreal) / sum(gh$n_dec))

  ## ---- dgp: G2 ----------------------------------------------------------------------------------
  P2 <- P$G2
  for (k in c("a", "cost", "delta", "n_pun", "K", "T", "delta_low")) put(paste0("dgp_G2_", k), P2[[k]])
  put("dgp_G2_n_market_months", P2$K * P2$T); put("dgp_G2_years", P2$T / 12); put("dgp_G2_n_sigma", length(P2$sigma))
  for (s in P2$sigma) put(paste0("dgp_G2_sigma_s", pid(s)), s)
  for (k in c("q", "price", "profit")) {
    put(paste0("dgp_G2_monopoly_", k), tr$G2$monopoly[[k]]); put(paste0("dgp_G2_cournot_", k), tr$G2$cournot[[k]])
  }
  put("dgp_G2_monopoly_total_q", 2 * tr$G2$monopoly$q); put("dgp_G2_cournot_total_q", 2 * tr$G2$cournot$q)
  put("dgp_G2_deviation_q", tr$G2$deviation$q); put("dgp_G2_deviation_profit", tr$G2$deviation$profit)
  put("dgp_G2_deviation_price", P2$a - tr$G2$monopoly$q - tr$G2$deviation$q)
  put("dgp_G2_grim_delta_star", (tr$G2$deviation$profit - tr$G2$monopoly$profit) / (tr$G2$deviation$profit - tr$G2$cournot$profit))

  ## ---- dgp: G3 ----------------------------------------------------------------------------------
  P3 <- P$G3
  put("dgp_G3_T", P3$T); put("dgp_G3_R", P3$R); put("dgp_G3_n_checkpoints", length(P3$checkpoints))
  cells("dgp_G3_chicken_A", P3$chicken_A); cells("dgp_G3_chicken_B", P3$chicken_B); cells("dgp_G3_mp_A", P3$mp_A)
  t3 <- tr$G3
  put("dgp_G3_chicken_mixed_p", t3$chicken_mixed[1]); put("dgp_G3_chicken_mixed_q", t3$chicken_mixed[2])
  put("dgp_G3_chicken_mixed_payoff", t3$chicken_mixed_payoff)
  put("dgp_G3_chicken_ce3_payoff", t3$chicken_ce3_payoff); put("dgp_G3_chicken_ce3_violation", t3$chicken_ce3_violation)
  put("dgp_G3_chicken_ce3_total", 2 * t3$chicken_ce3_payoff)
  put("dgp_G3_chicken_ne_total", sum(t3$chicken_pure_payoffs[1, ]))
  put("dgp_G3_mp_mixed_p", t3$mp_mixed[1]); put("dgp_G3_mp_mixed_q", t3$mp_mixed[2])

  ## ---- dgp: G4 ----------------------------------------------------------------------------------
  P4 <- P$G4
  put("dgp_G4_prior_mean", P4$prior_mean); put("dgp_G4_prior_sd", P4$prior_sd); put("dgp_G4_prior_precision", 1 / P4$prior_sd^2)
  put("dgp_G4_n_pop", P4$n_pop); put("dgp_G4_n_mc", P4$n_mc)
  put("dgp_G4_theta_lo", min(P4$theta_grid)); put("dgp_G4_theta_hi", max(P4$theta_grid)); put("dgp_G4_theta_step", diff(P4$theta_grid)[1])
  put("dgp_G4_theta_star_limit", tr$G4$theta_star_limit)
  for (k in seq_along(P4$sigma)) put(paste0("dgp_G4_sigma_s", pid(P4$sigma[k])), P4$sigma[k])

  ## ---- g0: the hand block -----------------------------------------------------------------------
  G0 <- tr$G0
  scal <- names(G0)[vapply(G0, function(z) (is.numeric(z) || is.logical(z)) && length(z) == 1L && is.null(dim(z)), logical(1))]
  for (k in scal) put(paste0("g0_", k), as.numeric(G0[[k]]))
  cells("g0_pd_A", G0$pd_A); cells("g0_compat_A", G0$compat_A); cells("g0_compat_B", G0$compat_B); cells("g0_bne_A", G0$bne_A)
  for (j in seq_along(G0$grim_deltas)) put(paste0("g0_grim_slack_d", pid(G0$grim_deltas[j])), G0$grim_slack[j])
  for (j in seq_along(G0$rub_delta)) {
    put(paste0("g0_rub_mover_d", pid(G0$rub_delta[j])), G0$rub_mover[j])
    put(paste0("g0_rub_responder_d", pid(G0$rub_delta[j])), G0$rub_responder[j])
    put(paste0("g0_rub_mover_pct_d", pid(G0$rub_delta[j])), 100 * G0$rub_mover[j])
    put(paste0("g0_rub_responder_pct_d", pid(G0$rub_delta[j])), 100 * G0$rub_responder[j])
  }
  put("g0_brw_r_osprey", G0$brw_r[["osprey"]]); put("g0_brw_r_boreal", G0$brw_r[["boreal"]])
  dt_lab <- c("year", "month", "week", "day")
  for (j in seq_along(G0$brw_dt)) {
    put(paste0("g0_brw_dt_", dt_lab[j]), G0$brw_dt[j])
    put(paste0("g0_brw_share_osprey_first_", dt_lab[j]), G0$brw_share_osprey_first[j])
    put(paste0("g0_brw_share_osprey_second_", dt_lab[j]), G0$brw_share_osprey_second[j])
  }
  put("g0_pd_ne_osprey_action", G0$pd_ne[[1]][1]); put("g0_pd_ne_boreal_action", G0$pd_ne[[1]][2])
  put("g0_seq_boreal_follow_S1", G0$seq_boreal_follow[1]); put("g0_seq_boreal_follow_S2", G0$seq_boreal_follow[2])
  put("g0_compat_n_pure", length(G0$compat_pure))
  put("g0_folk_minmax_osprey", G0$folk_minmax[1]); put("g0_folk_minmax_boreal", G0$folk_minmax[2])
  for (j in seq_len(nrow(G0$folk_ir))) {                 # counter-clockwise from the minmax corner
    put(sprintf("g0_folk_ir_v%d_osprey", j), G0$folk_ir[j, 1]); put(sprintf("g0_folk_ir_v%d_boreal", j), G0$folk_ir[j, 2])
  }

  ## ---- qre: block G1 (02_qre.R) ----------------------------------------------------------------------
  # draw_*: plain OLS on the displayed draw; mc_*: over the Monte Carlo replications, the mean estimate, the
  # 2.5% to 97.5% range of the estimates and their sd (Ruling P23: no standard-error formulas)
  for (k in c("unique", "pooled", "by_market")) {
    z <- qr[[k]]; pre <- paste0("qre_", k, "_")
    put(paste0(pre, "draw_u11"), z$draw[["u11"]]); put(paste0(pre, "draw_u22"), z$draw[["u22"]])
    put(paste0(pre, "draw_u11_minus_truth"), z$draw[["u11"]] - qr$payoff)
    put(paste0(pre, "n_groups"), z$n_groups)
    put(paste0(pre, "mc_mean_u11"), z$coef[["u11"]]); put(paste0(pre, "mc_mean_u22"), z$coef[["u22"]])
    put(paste0(pre, "mc_sd_u11"), z$sd[["u11"]]); put(paste0(pre, "mc_sd_u22"), z$sd[["u22"]])
    put(paste0(pre, "mc_range_lo_u11"), z$ci[["lo"]]); put(paste0(pre, "mc_range_hi_u11"), z$ci[["hi"]])
    put(paste0(pre, "mc_range_lo_u22"), z$ci_u22[["lo"]]); put(paste0(pre, "mc_range_hi_u22"), z$ci_u22[["hi"]])
    put(paste0(pre, "mc_range_covers_u11"), as.numeric(z$ci[["lo"]] <= qr$payoff && qr$payoff <= z$ci[["hi"]]))
    put(paste0(pre, "mc_range_covers_u22"), as.numeric(z$ci_u22[["lo"]] <= qr$payoff_u22 && qr$payoff_u22 <= z$ci_u22[["hi"]]))
    put(paste0(pre, "mc_mean_u11_minus_truth"), z$coef[["u11"]] - qr$payoff)
    put(paste0(pre, "mc_mean_u11_rel_bias"), z$coef[["u11"]] / qr$payoff - 1)
  }
  put("qre_mc_R", qr$R)
  mx <- qr$pooled_mixture
  for (i in seq_len(nrow(mx))) {
    pre <- sprintf("qre_pooled_x%s_", pid(mx$x[i]))
    put(paste0(pre, "share_S1"), mx$share_S1[i]); put(paste0(pre, "p"), mx$p[i]); put(paste0(pre, "q"), mx$q[i])
  }

  ## ---- gp: block G2 (03_green_porter.R) ----------------------------------------------------------------
  for (l in seq_along(gp$sigma)) {
    pre <- sprintf("gp_s%s_", pid(gp$sigma[l]))
    for (k in c("trigger", "q_c", "price_collusive", "beta", "value_per_period", "value_war_start_per_period",
                "war_share_theory", "war_freq", "war_freq_se", "wars_per_market_year", "mean_price_collusion",
                "share_markets_with_war"))
      put(paste0(pre, k), gp[[k]][l])
    if (is.finite(gp$mean_price_war[l])) put(paste0(pre, "mean_price_war"), gp$mean_price_war[l])
    put(paste0(pre, "osd"), as.numeric(gp$osd_each[l]))
    put(paste0(pre, "q_c_total"), 2 * gp$q_c[l])
    put(paste0(pre, "price_minus_trigger"), gp$price_collusive[l] - gp$trigger[l])
    if (gp$sigma[l] > 0) put(paste0(pre, "trigger_gap_in_sd"), (gp$trigger[l] - gp$price_collusive[l]) / gp$sigma[l])
    put(paste0(pre, "value_gain_over_cournot"), gp$value_per_period[l] - gp$cournot_profit)
    put(paste0(pre, "value_share_of_collusive_gain"), (gp$value_per_period[l] - gp$cournot_profit) / (gp$monopoly_profit - gp$cournot_profit))
    if (gp$beta[l] > 0) put(paste0(pre, "expected_months_between_wars"), 1 / gp$beta[l])
    # patience at delta_low (03_green_porter.R header): (a) the delta = 0.95 strategies held fixed;
    # (b) the best collusion at delta_low with the trigger re-optimised
    put(paste0(pre, "fixed_strategy_pass_at_delta_low"), as.numeric(gp$fixed_strategy_pass_low[l]))
    put(paste0(pre, "fixed_strategy_gain_at_delta_low"), gp$fixed_strategy_gain_low[l])
    put(paste0(pre, "fixed_strategy_q_dev_at_delta_low"), gp$fixed_strategy_q_dev_low[l])
    put(paste0(pre, "low_collusive"), as.numeric(gp$collusive_low[l]))
    put(paste0(pre, "low_value_per_period"), gp$value_low_per_period[l])
    put(paste0(pre, "low_q_c"), gp$q_c_low[l])
    if (gp$collusive_low[l]) {                            # (no gain key when it collapses to Cournot: it would be rounding noise)
      put(paste0(pre, "low_value_gain_over_cournot"), gp$value_low_per_period[l] - gp$cournot_profit)
      put(paste0(pre, "low_trigger"), gp$trigger_low[l]); put(paste0(pre, "low_beta"), gp$beta_low[l])
      put(paste0(pre, "low_war_share_theory"), gp$war_share_theory_low[l])
      put(paste0(pre, "low_price_collusive"), P2$a - 2 * gp$q_c_low[l])
    }
  }
  put("gp_osd_ok", as.numeric(gp$osd_ok))
  put("gp_fixed_strategies_pass_at_delta_low", as.numeric(gp$osd_low_delta))     # the card's g2$osd_low_delta
  put("gp_n_collusive_at_delta_low", sum(gp$collusive_low))
  put("gp_cournot_profit", gp$cournot_profit); put("gp_cournot_price", gp$cournot_price); put("gp_cournot_q", gp$cournot_q)
  put("gp_monopoly_profit", gp$monopoly_profit); put("gp_monopoly_price", gp$monopoly_price)

  ## ---- learn: block G3 (04_learning.R) ------------------------------------------------------------------
  for (g in c("compat", "chicken")) {
    pre <- paste0("learn_", g, "_")
    put(paste0(pre, "ce_violation_T"), ln$ce_violation_T[[g]])
    cells(paste0(pre, "mu_T"), ln$mu_T[[g]])
    put(paste0(pre, "final_cell"), ln$final_cell[[g]]); put(paste0(pre, "lock_in_period"), ln$lock_in_period[[g]])
    put(paste0(pre, "reps_max_violation_T"), ln$reps_max_violation[[g]]); put(paste0(pre, "reps_mean_violation_T"), ln$reps_mean_violation[[g]])
    put(paste0(pre, "reps_share_zero_T"), ln$reps_share_zero_T[[g]])
    sc <- ln$reps_share_cell[[g]]
    for (j in 1:4) put(sprintf("%sreps_share_cell%s", pre, c("11", "12", "21", "22")[j]), sc[[j]])
    pa <- ln$paths[ln$paths$game == g, ]
    for (t in c(10, 100, 1000, 10000, 50000)) {
      r <- pa[pa$t == t, ]
      if (nrow(r) != 1L) stop("collect_numbers02: checkpoint ", t, " missing")
      put(sprintf("%sreps_mean_violation_t%d", pre, t), r$viol_reps_mean); put(sprintf("%sreps_max_violation_t%d", pre, t), r$viol_reps_max)
    }
  }
  put("learn_fp_p", ln$fp_freq[["p"]]); put("learn_fp_q", ln$fp_freq[["q"]])
  cells("learn_fp_joint", ln$fp_joint)
  put("learn_fp_ce_violation", ln$fp_ce_violation); put("learn_fp_max_dev_last_half", ln$fp_max_dev_last_half)

  ## ---- gg: block G4 (05_global_game.R) -------------------------------------------------------------------
  for (k in seq_along(gg$sigma)) {
    pre <- sprintf("gg_s%s_", pid(gg$sigma[k]))
    for (f in c("theta_star", "dist_to_limit", "theta_star_sim", "m_post", "El_exact", "El_sim", "mc_indiff", "min_slope",
                "share_at_half", "share_at_055", "flat_prior_threshold"))
      put(paste0(pre, f), gg[[f]][k])
  }
  put("gg_sigma_max", gg$sigma_max)
  me <- gg$multiple_example
  put("gg_multi_sigma", me$sigma); put("gg_multi_prior_mean", me$prior_mean); put("gg_multi_prior_sd", me$prior_sd)
  put("gg_multi_sigma_max", me$sigma_max); put("gg_multi_n_roots", length(me$roots))
  for (j in seq_along(me$roots)) put(paste0("gg_multi_root", j), me$roots[j])
  nums
}
