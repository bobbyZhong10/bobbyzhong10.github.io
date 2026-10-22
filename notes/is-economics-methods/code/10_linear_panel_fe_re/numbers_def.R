# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers10(dir = "."): every quotable number of chapter 10 as a flat named
# list of finite numeric scalars (10_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth10.rds, res_static.rds (legacy 02) and res_extra.rds (08) in `dir`.
#
# Name scheme: <prefix>_<object>_<field>. Prefixes: dgp (01_dgp.R: parameters, seeds, truth,
# closed-form probability limits, large-sample checks), static (the legacy res_static.rds, block P1),
# p1 to p7 and p5b (blocks P1 to P7: 03 to 07), mundlak (the identities collected by 08). Conventions:
#   *_est, *_se, *_lb, *_ub   estimate, cluster-robust (by merchant, CR0) standard error, 95% interval
#                             est -/+ qnorm(0.975) se; *_se_classic is the non-robust textbook formula
#                             (the between estimator and Hausman-Taylor have only *_se_classic or *_se
#                             as stated in the script)
#   *_covers_truth            1 when that interval contains the true value, 0 otherwise
#   *_stat, *_df, *_p         test statistic, degrees of freedom, p value
#   sd_*, mean_*, size_*      Monte Carlo sampling SD, mean, rejection rate
#   0p05                      a decimal point written p (nominal level 0.05)
#   usd                       an amount in US dollars (Y is GMV in thousands of USD, X promotion spend in
#                             hundreds of USD; P1 keeps the legacy standardised units)
# ------------------------------------------------------------------
collect_numbers10 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 08 first")
    readRDS(p)
  }
  tr <- rd("truth10.rds"); st <- rd("res_static.rds"); x <- rd("res_extra.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers10: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L || !is.finite(value)) stop("collect_numbers10: ", name, " is not a finite scalar")
    nums[[name]] <<- value
    invisible(NULL)
  }
  put_list <- function(prefix, L, keys = names(L)) for (k in keys) put(paste0(prefix, "_", k), L[[k]])
  put_row <- function(prefix, r) {                 # est, se, ci (and se_classic, covers_truth when present)
    put(paste0(prefix, "_est"), r$est)
    if (!is.null(r$se)) put(paste0(prefix, "_se"), r$se)
    if (!is.null(r$se_classic)) put(paste0(prefix, "_se_classic"), r$se_classic)
    put(paste0(prefix, "_lb"), r$ci[1]); put(paste0(prefix, "_ub"), r$ci[2])
    if (!is.null(r$covers_truth)) put(paste0(prefix, "_covers_truth"), r$covers_truth)
  }
  put_test <- function(prefix, w) { put(paste0(prefix, "_stat"), w$stat); put(paste0(prefix, "_df"), w$df); put(paste0(prefix, "_p"), w$p) }
  pid <- function(v) gsub("\\.", "p", format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))
  P <- tr$params

  ## ---- dgp: seeds, parameters, truth --------------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  put_list("dgp_P1", P$P1)
  put_list("dgp_P1", tr$P1, c("plim_pols", "plim_be", "plim_re", "varphi_plim", "w_plim", "varphi_components"))
  put("dgp_P1_plim_pols_minus_beta", tr$P1$plim_pols - tr$P1$beta); put("dgp_P1_plim_re_minus_beta", tr$P1$plim_re - tr$P1$beta)
  put("dgp_P1_plim_be_minus_beta", tr$P1$plim_be - tr$P1$beta)
  put_list("dgp_P1_large", tr$P1$large)

  put_list("dgp_P2", P$P2)
  put("dgp_P2_beta_usd_per_100usd", 1000 * P$P2$beta)                       # GMV (USD) per 100 USD of promotion
  put("dgp_P2_x_mean_usd", 100 * P$P2$x_mean)

  put_list("dgp_P3", P$P3)
  put_list("dgp_P3", tr$P3, c("fe_plim", "fd_plim"))
  put("dgp_P3_fe_plim_minus_beta", tr$P3$fe_plim - tr$P3$beta); put("dgp_P3_fd_plim_minus_beta", tr$P3$fd_plim - tr$P3$beta)
  put("dgp_P3_within_var_x", P$P3$feedback^2 * P$P3$sigma_e^2 + P$P3$x_sd^2)
  put_list("dgp_P3_large", tr$P3$large)

  put_list("dgp_P4", P$P4)
  put("dgp_P4_gamma_usd", 1000 * P$P4$gamma)

  P5 <- P$P5
  put_list("dgp_P5", P5[setdiff(names(P5), c("s_types", "s_probs", "design"))])
  for (j in seq_along(P5$s_types)) {
    sk <- pid(P5$s_types[j])
    put(paste0("dgp_P5_type_s", sk, "_prob"), P5$s_probs[j]); put(paste0("dgp_P5_type_s", sk, "_sd_usd"), 100 * P5$s_types[j])
    put(paste0("dgp_P5_type_s", sk, "_slope"), tr$P5$type_slope[j])
    put(paste0("dgp_P5_type_s", sk, "_weight_fe"), tr$P5$type_weight_fe[j])
  }
  put_list("dgp_P5", tr$P5, c("fe_plim", "Es2", "Vs2"))                      # avg_slope is in P5 (parameters)
  put("dgp_P5_fe_plim_minus_avg", tr$P5$fe_plim - tr$P5$avg_slope)
  put("dgp_P5_avg_slope_usd_per_usd", 10 * tr$P5$avg_slope); put("dgp_P5_fe_plim_usd_per_usd", 10 * tr$P5$fe_plim)
  put_list("dgp_P5_large", tr$P5$large)
  put("dgp_P5b_kappa", P$P5b$kappa); put("dgp_P5b_avg_slope", tr$P5b$avg_slope); put("dgp_P5b_fe_plim", tr$P5b$fe_plim)
  put_list("dgp_P5b_large", tr$P5b$large)

  put_list("dgp_P6", P$P6)
  put("dgp_P6_floor_usd", 1000 * P$P6$floor)
  put("dgp_P6_hazard_alpha_m1p5", pnorm(P$P6$haz0 + P$P6$haz_alpha * (-1.5)))
  put("dgp_P6_hazard_alpha_0", pnorm(P$P6$haz0)); put("dgp_P6_hazard_alpha_1p5", pnorm(P$P6$haz0 + P$P6$haz_alpha * 1.5))
  put_list("dgp_P6_large", tr$P6$large)
  put("dgp_P6_large_fe_eps_minus_beta", tr$P6$large$fe_eps - tr$P6$beta)

  put_list("dgp_P7", P$P7)

  ## ---- static: the legacy res_static.rds (block P1) ---------------------------------------------
  for (k in names(st)) put(paste0("static_", k), st[[k]])

  ## ---- p1: block P1 beyond the legacy script ------------------------------------------------------
  X1 <- x$P1
  for (k in c("pool", "re", "fe", "fd")) put_row(paste0("p1_", k), X1[[k]])
  put_row("p1_be", X1$be)
  put("p1_twfe_est", x$twfe)
  put_list("p1", X1, c("varphi", "sigma2_eps", "sigma2_alpha", "share_alpha", "lsdv", "lsdv_se_classic",
                       "mundlak_se_x", "mundlak_se_xbar", "lead_p", "lead_stat", "lead_delta", "lead_se", "lead_n"))
  put("p1_re_w_between", X1$re_decomp$w); put("p1_re_S_W", X1$re_decomp$S_W); put("p1_re_S_B", X1$re_decomp$S_B)
  put("p1_re_decomp_check", X1$re_decomp$check)
  put("p1_mundlak_b_x", x$mundlak$b_x); put("p1_mundlak_b_xbar", x$mundlak$b_xbar)
  put("p1_mundlak_b_x_minus_fe", x$mundlak$b_x - st$fe); put("p1_be_minus_fe", X1$be$est - st$fe)
  put("p1_mundlak_tw_b_x", x$mundlak_tw$b_x); put("p1_mundlak_tw_b_x_minus_twfe", x$mundlak_tw$b_x - x$twfe)
  put_test("p1_haus_classic", X1$haus_classic); put_test("p1_haus_robust", X1$haus_robust)
  put_test("p1_haus_robust_tw", X1$haus_robust_tw); put_test("p1_haus_aux", X1$haus_aux)
  put("p1_haus_classic_minus_static", X1$haus_classic$stat - st$haus_stat)
  put_list("p1_fd_ar1", X1$fd_ar1)
  ch <- X1$chamberlain
  put("p1_chamberlain_est", ch$beta); put("p1_chamberlain_se", ch$se); put("p1_chamberlain_lb", ch$ci[1])
  put("p1_chamberlain_ub", ch$ci[2]); put_test("p1_chamberlain_overid", ch)
  for (j in seq_along(ch$xi)) put(sprintf("p1_chamberlain_xi_t%d", j), ch$xi[j])
  put("p1_chamberlain_xi_sum", sum(ch$xi))
  put("p1_pool_minus_beta", X1$pool$est - tr$P1$beta); put("p1_re_minus_beta", X1$re$est - tr$P1$beta)

  ## ---- p2: FE vs FD -----------------------------------------------------------------------------
  X2 <- x$P2
  put_list("p2_mc", X2, setdiff(names(X2), "sample"))
  for (k in c("iid", "rw")) {
    S <- X2$sample[[k]]
    put_row(paste0("p2_", k, "_fe"), S$fe); put_row(paste0("p2_", k, "_fd"), S$fd)
    put_list(paste0("p2_", k, "_fd_ar1"), S$fd_ar1)
    put_list(paste0("p2_", k), S, c("pwartest_p", "pwfdtest_fd_p", "pwfdtest_fe_p"))
  }

  ## ---- p3: feedback ------------------------------------------------------------------------------
  X3 <- x$P3
  put_row("p3_fe", X3$fe); put("p3_fe_covers_plim", X3$fe$covers_plim)
  put_row("p3_fd", X3$fd); put("p3_fd_covers_plim", X3$fd$covers_plim)
  put_list("p3", X3, c("lead_p", "lead_stat", "lead_delta", "lead_se", "lead_n"))
  put("p3_chamberlain_est", X3$chamberlain$beta); put("p3_chamberlain_se", X3$chamberlain$se)
  put("p3_chamberlain_lb", X3$chamberlain$ci[1]); put("p3_chamberlain_ub", X3$chamberlain$ci[2])
  put_test("p3_chamberlain_overid", X3$chamberlain)
  put("p3_fe_minus_beta", X3$fe$est - tr$P3$beta); put("p3_fd_minus_beta", X3$fd$est - tr$P3$beta)

  ## ---- p4: Hausman-Taylor -------------------------------------------------------------------------
  X4 <- x$P4
  for (k in c("ht", "re", "pool", "re_cl", "ht_beta1", "ht_beta2", "re_beta2", "fe_beta1", "fe_beta2")) put_row(paste0("p4_", k), X4[[k]])
  put_list("p4", X4, c("varphi_ht", "varphi_re", "ht_iv_maxdiff", "share_certified", "corr_w_x1bar", "corr_w_alpha",
                       "first_stage_coef", "first_stage_F", "N", "T"))
  put_test("p4_haus_robust", X4$haus_robust)
  put("p4_re_minus_gamma", X4$re$est - tr$P4$gamma); put("p4_ht_minus_gamma", X4$ht$est - tr$P4$gamma)
  put("p4_ht_est_usd", 1000 * X4$ht$est); put("p4_re_est_usd", 1000 * X4$re$est)

  ## ---- p5 and p5b: random slopes ---------------------------------------------------------------------
  for (b in c("P5", "P5b")) {
    X5 <- x[[b]]; pf <- tolower(b)
    put_row(paste0(pf, "_fe"), X5$fe); put_row(paste0(pf, "_mg"), X5$mg)
    put(paste0(pf, "_fe_covers_plim"), X5$fe_covers_plim)
    put_list(pf, X5, c("pooled", "sample_avg_b", "sample_weighted_b"))
    put_test(paste0(pf, "_test_xbar"), X5$test_xbar); put_test(paste0(pf, "_test_s2"), X5$test_s2)
    bt <- X5$by_type
    for (j in seq_len(nrow(bt))) {
      sk <- pid(bt$s[j])
      put(paste0(pf, "_type_s", sk, "_n"), bt$n[j]); put(paste0(pf, "_type_s", sk, "_mean_b"), bt$mean_b[j])
      put(paste0(pf, "_type_s", sk, "_weight_fe"), bt$weight_fe[j])
    }
    put(paste0(pf, "_fe_minus_avg"), X5$fe$est - tr[[b]]$avg_slope)
  }

  ## ---- p6: unbalanced panel and attrition ----------------------------------------------------------
  X6 <- x$P6
  for (k in c("full", "attr_alpha", "attr_eps")) {
    put_row(paste0("p6_", k), X6[[k]])
    put_list(paste0("p6_", k), X6[[k]], c("n", "G", "mean_Ti", "share_Ti_full", "share_Ti_1", "keep"))
  }
  put("p6_attr_eps_bias", X6$attr_eps$bias); put("p6_attr_eps_z_vs_truth", X6$attr_eps$z_vs_truth)
  put_list("p6", X6, c("alpha_leavers_alpha", "alpha_stayers_alpha", "mean_y_full", "share_below_floor", "mundlak_se_x"))
  put_list("p6_re_alpha", X6$re_alpha)
  put_test("p6_haus_robust_alpha", X6$haus_robust_alpha)

  ## ---- mundlak identities (08) -----------------------------------------------------------------------
  put("mundlak_unbal_b_x", x$mundlak_unbal$b_x); put("mundlak_unbal_b_xbar", x$mundlak_unbal$b_xbar)
  put("mundlak_fe_unbal", x$fe_unbal); put("mundlak_unbal_b_x_minus_fe", x$mundlak_unbal$b_x - x$fe_unbal)
  put("mundlak_tw_unbal_b_x", x$mundlak_tw_unbal$b_x); put("mundlak_twfe_unbal", x$twfe_unbal)
  put("mundlak_tw_unbal_b_x_minus_twfe", x$mundlak_tw_unbal$b_x - x$twfe_unbal)
  put("mundlak_tw_naive_unbal", x$mundlak_tw_naive_unbal)
  put("mundlak_tw_naive_unbal_minus_twfe", x$mundlak_tw_naive_unbal - x$twfe_unbal)

  ## ---- p7: Hausman size Monte Carlo -------------------------------------------------------------------
  X7 <- x$P7
  put_list("p7", X7, setdiff(names(X7), c("curve", "sample")))
  for (j in seq_len(nrow(X7$curve))) {
    a <- pid(X7$curve$nominal[j])
    put(paste0("p7_curve_", a, "_classic"), X7$curve$classic[j]); put(paste0("p7_curve_", a, "_robust"), X7$curve$robust[j])
  }
  put_list("p7_sample", X7$sample)

  nums
}
