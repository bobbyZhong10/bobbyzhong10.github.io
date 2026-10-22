# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers04(dir = "."): every quotable number of chapter 4 as a flat named list of finite
# numeric scalars (07_summary.R writes it to numbers.json with chapter_kit.R's write_numbers_json;
# tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth04.rds, cirrus_C1.rds, cirrus_C2.rds, cirrus_C3.rds and res_logit.rds, res_pd.rds, res_bundle.rds,
# res_leader.rds in `dir`.
#
# Name scheme: <prefix>_<block or object>_<field>. Prefixes: dgp (01_dgp.R: seeds, parameters, realised sample
# facts, simulation truth), c0 (block C0, the hand block), logit (02, block C1), pd (03, block C2), bundle (04,
# block C3), leader (05, block C4). Conventions:
#   cirrus, garnet       the two firms (C1, C4); in pass-through keys <price>_<cost>: logit_ols_bertrand_cirrus_garnet
#                        is dp_Cirrus / dc_Garnet under Bertrand (the text's PT_{Garnet, Cirrus})
#   bertrand, joint      conduct kappa = 0 and kappa = 1
#   uniform, coarse, fine   the C2 regimes (K = 1, 4, 50 segments); K<k> a rung of the C2 ladder
#   J<k>                 the C3 bundle of the first k goods
#   _1, _2, ...          elements of a vector, in order (c0_seg_p_1 is segment 1's price)
#   p                    decimal point written p in a name (n100, e2, rate0p05)
#   pct                  a ratio minus 1 (0.19 = 19 percent); pp: a difference of shares
#   flags                0/1 numbers
# All money is USD (C0: per unit of operating activity; C1: per unit, market size 1; C2: per merchant per month;
# C3: per consumer, values in [0, 1] per good). A non-finite value is an error (put() stops).
# ------------------------------------------------------------------
collect_numbers04 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 05 first")
    readRDS(p)
  }
  tr <- rd("truth04.rds"); C1 <- rd("cirrus_C1.rds"); C2 <- rd("cirrus_C2.rds"); C3 <- rd("cirrus_C3.rds")
  lg <- rd("res_logit.rds"); pd <- rd("res_pd.rds"); bu <- rd("res_bundle.rds"); ld <- rd("res_leader.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers04: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers04: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers04: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  putv <- function(name, v) for (i in seq_along(v)) put(paste0(name, "_", i), v[[i]])
  pid <- function(x) gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))
  firm <- c("cirrus", "garnet")
  P <- tr$params

  ## ---- dgp: seeds, parameters, sample facts, truth -----------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  P1 <- P$C1
  put("dgp_C1_M", P1$M); put("dgp_C1_alpha", P1$alpha)
  put("dgp_C1_delta_cirrus", P1$delta[["Cirrus"]]); put("dgp_C1_delta_garnet", P1$delta[["Garnet"]])
  put("dgp_C1_c_base_cirrus", P1$c_base[1]); put("dgp_C1_c_base_garnet", P1$c_base[2])
  put("dgp_C1_cost_halfwidth", P1$cost_halfwidth); put("dgp_C1_xi_sd", P1$xi_sd)
  put("dgp_C1_hom_a", P1$hom_a); put("dgp_C1_hom_b", P1$hom_b); put("dgp_C1_hom_eta_sd", P1$hom_eta_sd)
  put("dgp_C1_sample_c_cirrus_mean", mean(C1$c1)); put("dgp_C1_sample_c_garnet_mean", mean(C1$c2))
  put("dgp_C1_sample_c_cirrus_min", min(C1$c1)); put("dgp_C1_sample_c_cirrus_max", max(C1$c1))
  put("dgp_C1_sample_cor_c_cirrus_c_garnet", stats::cor(C1$c1, C1$c2))
  T1 <- tr$C1
  for (j in 1:2) {
    put(paste0("dgp_C1_p_bertrand_base_", firm[j]), T1$p_bertrand_base[j]); put(paste0("dgp_C1_s_bertrand_base_", firm[j]), T1$s_bertrand_base[j])
    put(paste0("dgp_C1_p_joint_base_", firm[j]), T1$p_joint_base[j]); put(paste0("dgp_C1_s_joint_base_", firm[j]), T1$s_joint_base[j])
    put(paste0("dgp_C1_markup_bertrand_base_", firm[j]), T1$p_bertrand_base[j] - T1$c_base[j])
  }
  put("dgp_C1_s0_bertrand_base", 1 - sum(T1$s_bertrand_base)); put("dgp_C1_s0_joint_base", 1 - sum(T1$s_joint_base))
  put("dgp_C1_markup_joint_base", T1$p_joint_base[1] - T1$c_base[1])
  for (cond in c("bertrand", "joint")) for (what in c("base", "mean")) {
    M <- T1[[paste0("pt_", cond, "_", what)]]
    for (k in 1:2) for (j in 1:2) put(sprintf("dgp_C1_pt_%s_%s_%s_%s", cond, what, firm[k], firm[j]), M[k, j])
  }
  put("dgp_C1_pt_own_cirrus_bertrand_min", T1$pt_own1_bertrand_range[1]); put("dgp_C1_pt_own_cirrus_bertrand_max", T1$pt_own1_bertrand_range[2])
  put("dgp_C1_pt_own_cirrus_joint_min", T1$pt_own1_joint_range[1]); put("dgp_C1_pt_own_cirrus_joint_max", T1$pt_own1_joint_range[2])
  put("dgp_C1_pt_cournot_hom", T1$pt_cournot_hom)
  P2 <- P$C2
  put("dgp_C2_N", P2$N); put("dgp_C2_sigma", P2$sigma); put("dgp_C2_cost", P2$cost); put("dgp_C2_R_mc", P2$R_mc)
  for (k in names(P2$beta)) put(paste0("dgp_C2_beta_", k), P2$beta[[k]])
  put("dgp_C2_median_wtp_at_const", exp(P2$beta[["const"]]))
  for (i in seq_along(P2$cat_levels)) put(paste0("dgp_C2_cat_prob_", P2$cat_levels[i]), P2$cat_probs[i])
  for (k in names(P2$K)) put(paste0("dgp_C2_K_", k), P2$K[[k]])
  put("dgp_C2_n_ladder", length(P2$K_ladder)); put("dgp_C2_K_ladder_max", max(P2$K_ladder))
  q <- function(x, p) unname(stats::quantile(x, p, type = 7))
  w <- C2$wtp
  put("dgp_C2_sample_n_train", sum(C2$train)); put("dgp_C2_sample_n_eval", sum(!C2$train))
  put("dgp_C2_sample_wtp_mean", mean(w)); put("dgp_C2_sample_wtp_median", stats::median(w))
  put("dgp_C2_sample_wtp_min", min(w)); put("dgp_C2_sample_wtp_max", max(w))
  put("dgp_C2_sample_wtp_p10", q(w, 0.1)); put("dgp_C2_sample_wtp_p90", q(w, 0.9))
  put("dgp_C2_sample_wtp_mean_train", mean(w[C2$train])); put("dgp_C2_sample_wtp_mean_eval", mean(w[!C2$train]))
  put("dgp_C2_sample_n_wtp_values", length(unique(w)))
  for (lv in P2$cat_levels) put(paste0("dgp_C2_sample_share_", lv), mean(C2$category == lv))
  put("dgp_C2_r2_latent", tr$C2$r2_latent); put("dgp_C2_efficient", tr$C2$efficient); put("dgp_C2_efficient_eval", tr$C2$efficient_eval)
  P3 <- P$C3
  put("dgp_C3_N", P3$N); put("dgp_C3_n_bundles", length(P3$J)); put("dgp_C3_J_max", max(P3$J))
  for (k in seq_along(P3$J)) put(paste0("dgp_C3_J", P3$J[k], "_mean_value"), mean(C3[, k]))
  put("dgp_C3_mean_value_per_good", tr$C3$mean_value)
  put("dgp_C3_rev_per_good_separate", tr$C3$rev_per_good_separate); put("dgp_C3_dwl_share_separate", tr$C3$dwl_share_separate)

  ## ---- c0: the hand block --------------------------------------------------------------------------------
  C0 <- tr$C0
  # vectors indexed by period t (Hotelling's Rule), by the number of firms n (Cournot path, pass-through) or by
  # the elasticity e get that index in the key; other vectors (segments, merchant types) are numbered _1, _2
  by_t <- c("hot_rent", "hot_price", "hot_pv"); by_n_path <- c("path_p", "path_Q", "path_profit_total", "path_dwl")
  skip <- c("hot_t", "pt_n", "path_n", "ptm_ce_e")
  for (k in setdiff(names(C0), skip)) {
    v <- C0[[k]]
    if (k %in% by_t) for (i in seq_along(v)) put(sprintf("c0_%s_t%d", k, C0$hot_t[i]), v[i])
    else if (k %in% by_n_path) for (i in seq_along(v)) put(sprintf("c0_%s_n%d", k, C0$path_n[i]), v[i])
    else if (k == "pt_cournot") for (i in seq_along(v)) put(sprintf("c0_pt_cournot_n%d", C0$pt_n[i]), v[i])
    else if (k == "ptm_ce") for (i in seq_along(v)) put(sprintf("c0_ptm_ce_e%s", pid(C0$ptm_ce_e[i])), v[i])
    else if (length(v) == 1L) put(paste0("c0_", k), v) else putv(paste0("c0_", k), v)
  }
  put("c0_dwl_share_m", C0$dwl_m / C0$ts_competitive)
  put("c0_cournot_profit_total_over_monopoly", C0$cournot_profit_total / C0$profit_m)
  put("c0_vert_pd_minus_pm", C0$vert_pd - C0$vert_pm)
  put("c0_menu_gain_over_single", C0$profit_menu - C0$single_bench)

  ## ---- logit: block C1 (02_logit.R) ---------------------------------------------------------------------------
  o <- lg$ols
  for (i in seq_len(nrow(o))) {
    pre <- sprintf("logit_ols_%s_%s_%s_", o$conduct[i], tolower(o$price[i]), tolower(o$cost[i]))
    put(paste0(pre, "est"), o$estimate[i]); put(paste0(pre, "se"), o$se[i]); put(paste0(pre, "truth_mean"), o$truth_mean[i])
    put(paste0(pre, "est_minus_truth"), o$estimate[i] - o$truth_mean[i])
    put(paste0(pre, "ci_lo"), o$estimate[i] - 1.96 * o$se[i]); put(paste0(pre, "ci_hi"), o$estimate[i] + 1.96 * o$se[i])
  }
  for (cond in c("bertrand", "joint")) for (k in 1:2) {
    i <- which(o$conduct == cond & tolower(o$price) == firm[k])[1]
    put(sprintf("logit_ols_%s_%s_r2", cond, firm[k]), o$r2[i]); put(sprintf("logit_ols_%s_%s_sigma", cond, firm[k]), o$sigma[i])
  }
  for (k in c("pt_own_bertrand", "pt_cross_bertrand", "pt_own_joint", "pt_cross_joint",
              "se_own_bertrand", "se_cross_bertrand", "se_own_joint", "se_cross_joint")) put(paste0("logit_", k), lg[[k]])
  h <- lg$hom
  for (j in 1:2) {
    put(paste0("logit_hom_cournot_coef_", firm[j]), h$cournot_coef[j]); put(paste0("logit_hom_cournot_se_", firm[j]), h$cournot_se[j])
    put(paste0("logit_hom_bertrand_coef_", firm[j]), h$bertrand_coef[j]); put(paste0("logit_hom_bertrand_se_", firm[j]), h$bertrand_se[j])
  }
  put("logit_hom_cournot_coef_sum", sum(h$cournot_coef)); put("logit_hom_bertrand_coef_sum", sum(h$bertrand_coef))
  for (k in c("p_cournot_mean", "p_bertrand_mean", "markup_cournot_mean", "markup_bertrand_mean")) put(paste0("logit_hom_", k), h[[k]])
  put("logit_mono_p", lg$mono_example$p); put("logit_mono_s", lg$mono_example$s); put("logit_mono_pt", lg$mono_example$pt)
  pc <- lg$pt_curves
  for (i in seq_len(nrow(pc))) put(sprintf("logit_ptcurve_%s_n%d", pc$model[i], pc$n[i]), pc$pt[i])
  put("logit_ptcurve_delta", lg$pt_curves_par$delta); put("logit_ptcurve_ce_e", lg$pt_curves_par$ce_e)
  for (j in 1:2) {
    for (k in c("lerner_bertrand", "elas_abs_bertrand", "lerner_joint", "elas_abs_joint", "mc_if_bertrand_at_joint"))
      put(sprintf("logit_%s_%s", k, firm[j]), lg$lerner[[k]][j])
    put(paste0("logit_inv_elas_bertrand_", firm[j]), 1 / lg$lerner$elas_abs_bertrand[j])
    put(paste0("logit_mc_error_if_bertrand_", firm[j]), lg$lerner$mc_if_bertrand_at_joint[j] - T1$c_base[j])
    put(paste0("logit_mc_error_pct_if_bertrand_", firm[j]), lg$lerner$mc_if_bertrand_at_joint[j] / T1$c_base[j] - 1)
  }
  for (k in names(lg$sample)) put(paste0("logit_sample_", k), lg$sample[[k]])

  ## ---- pd: block C2 (03_personalized.R) ---------------------------------------------------------------------
  reg <- names(pd$K)
  for (r in reg) {
    put(paste0("pd_K_", r), pd$K[[r]])
    for (k in c("profit", "cs", "ts", "served", "mean_price_paid")) put(sprintf("pd_%s_%s", k, r), pd[[k]][[r]])
  }
  for (r in c("coarse", "fine")) {
    put(sprintf("pd_profit_gain_pct_%s", r), pd$profit[[r]] / pd$profit[["uniform"]] - 1)
    put(sprintf("pd_cs_change_pct_%s", r), pd$cs[[r]] / pd$cs[["uniform"]] - 1)
    put(sprintf("pd_ts_change_pct_%s", r), pd$ts[[r]] / pd$ts[["uniform"]] - 1)
    put(sprintf("pd_served_change_pp_%s", r), pd$served[[r]] - pd$served[["uniform"]])
  }
  put("pd_profit_fine_over_coarse_pct", pd$profit[["fine"]] / pd$profit[["coarse"]] - 1)
  for (i in seq_len(nrow(pd$ladder))) for (k in c("profit", "cs", "ts", "served", "mean_price_paid"))
    put(sprintf("pd_ladder_K%d_%s", pd$ladder$K[i], k), pd$ladder[[k]][i])
  put("pd_ladder_best_K", pd$ladder$K[which.max(pd$ladder$profit)]); put("pd_ladder_best_profit", max(pd$ladder$profit))
  for (i in seq_len(nrow(pd$ols))) {
    tm <- c("(Intercept)" = "const", sales_z = "sales", tenure_z = "tenure", categoryelectronics = "electronics",
            categoryhome = "home")[[pd$ols$term[i]]]
    put(paste0("pd_ols_", tm, "_est"), pd$ols$estimate[i]); put(paste0("pd_ols_", tm, "_se"), pd$ols$se[i])
    put(paste0("pd_ols_", tm, "_truth"), pd$ols$truth[i])
  }
  put("pd_ols_r2", pd$r2); put("pd_n_train", pd$n_train); put("pd_n_eval", pd$n_eval)
  put("pd_price_uniform", pd$prices$uniform)
  putv("pd_price_coarse", pd$prices$coarse)
  put("pd_price_fine_min", min(pd$prices$fine)); put("pd_price_fine_max", max(pd$prices$fine))
  put("pd_price_fine_median", stats::median(pd$prices$fine))
  put("pd_n_train_per_coarse_segment", pd$prices$coarse_n_train[1]); put("pd_n_train_per_fine_segment", pd$prices$fine_n_train[1])
  for (k in names(pd$by_merchant)) put(paste0("pd_fine_vs_uniform_", k), pd$by_merchant[[k]])
  put("pd_efficient", pd$efficient)
  for (k in c("uniform_profit", "uniform_price", "efficient", "profit", "total", "cs_max", "n_segments", "n_values",
              "largest_segment_weight", "min_price", "max_price", "uniform_gap")) put(paste0("pd_bbm_", k), pd$bbm[[k]])
  put("pd_bbm_cs_max_over_cs_uniform", pd$bbm$cs_max / pd$cs[["uniform"]])
  put("pd_first_degree_profit", pd$bbm$efficient); put("pd_first_degree_cs", 0)
  for (nm in c("linear", "linear_open")) {
    L <- pd[[nm]]
    for (k in c("p_uniform", "Q_uniform", "profit_uniform", "CS_uniform", "W_uniform", "n_served_uniform",
                "Q_pd", "profit_pd", "CS_pd", "W_pd", "n_served_pd")) put(sprintf("pd_%s_%s", nm, k), L[[k]])
    putv(sprintf("pd_%s_p_pd", nm), L$p_pd); putv(sprintf("pd_%s_q_pd", nm), L$q_pd); putv(sprintf("pd_%s_q_uniform", nm), L$q_uniform)
  }
  putv("pd_linear_a", pd$linear_a); putv("pd_linear_open_a", pd$linear_open_a)
  put("pd_mc_R", pd$mc$R)
  ms <- pd$mc$summary
  for (i in seq_len(nrow(ms))) for (k in c("mean", "lo", "hi", "share_positive")) put(sprintf("pd_mc_%s_%s", ms$stat[i], k), ms[[k]][i])

  ## ---- bundle: block C3 (04_bundle.R) -------------------------------------------------------------------------
  for (i in seq_along(bu$J)) for (k in c("price", "rev", "dwl", "share_buy", "cs", "efficient", "price_per_good",
                                         "rev_per_good", "dwl_share", "cs_per_good", "value_per_good", "sd_bundle_value"))
    put(sprintf("bundle_J%d_%s", bu$J[i], k), bu[[k]][i])
  nJ <- length(bu$J)
  put("bundle_rev_per_good_gain_pct", bu$rev_per_good[nJ] / bu$rev_per_good[1] - 1)
  put("bundle_rev_per_good_share_of_mean_J100", bu$rev_per_good[nJ] / tr$C3$mean_value)

  ## ---- leader: block C4 (05_leader.R) -------------------------------------------------------------------------
  # scenarios: sim (simultaneous), cirrus_leads (res_leader.rds field suffix seq), garnet_leads (seq_garnet), joint
  scen <- c(sim = "sim", cirrus_leads = "seq", garnet_leads = "seq_garnet", joint = "joint")
  for (sc in names(scen)) {
    f <- scen[[sc]]
    for (j in 1:2) {
      put(sprintf("leader_p_%s_%s", sc, firm[j]), ld[[paste0("p_", f)]][j])
      put(sprintf("leader_profit_%s_%s", sc, firm[j]), ld[[paste0("profit_", f)]][j])
      put(sprintf("leader_s_%s_%s", sc, firm[j]), ld[[paste0("s_", f)]][j])
    }
    put(sprintf("leader_cs_%s", sc), ld[[paste0("cs_", f)]])
  }
  for (sc in c("cirrus_leads", "garnet_leads")) for (j in 1:2) {
    f <- scen[[sc]]
    put(sprintf("leader_dp_%s_%s", sc, firm[j]), ld[[paste0("p_", f)]][j] - ld$p_sim[j])
    put(sprintf("leader_dprofit_pct_%s_%s", sc, firm[j]), ld[[paste0("profit_", f)]][j] / ld$profit_sim[j] - 1)
  }
  put("leader_flag_both_prices_above_sim", as.numeric(all(ld$p_seq > ld$p_sim) && all(ld$p_seq_garnet > ld$p_sim)))

  nums
}
