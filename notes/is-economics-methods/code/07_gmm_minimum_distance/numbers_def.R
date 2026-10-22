# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers07(dir = "."): every quotable number of chapter 7 as a flat named
# list of finite numeric scalars (06_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth07.rds and res_md.rds, res_gmm_mc.rds, res_weak.rds in `dir`, and sources gmm_fns.R from `dir` (pop_moments).
#
# Name scheme: <prefix>_<block or object>_<field>. Prefixes: dgp (01_dgp.R: parameters, seeds, truth,
# population moments), md (02, blocks G1 and G1h, and pooled GMM on the G1 data), gmm (03, block G2),
# weak (04, block G3). Conventions:
#   *_est, *_se          point estimate and standard error (the SE named in the estimator's description:
#                        robust sandwich for 2SLS and identity weights, conventional (D'WD)^-1 / n for two-step,
#                        iterated and CUE, Windmeijer for *_windmeijer)
#   *_lb, *_ub           95% interval est -/+ qnorm(0.975) se of the same row
#   *_covers*            1 when that interval contains the named value (the truth unless named), 0 otherwise
#   cover_*              Monte Carlo coverage rate of nominal 95% intervals
#   <market>             lower-case market name (atlanta ... seattle); z1 ... z10 the G2 experiments, z1 ... z4
#                        the G3 perturbations, intercept the constant moment
#   grid_n<nn>           G2 Monte Carlo at n = nn regions
# A non-finite value (Inf, NA, NaN) is an error (put() stops), as in chapter 5: every number listed here exists
# in the saved results (the G2 and G3 S-set bounds are finite). A value whose shortest decimal
# text jsonlite does not read back exactly (jsonlite misreads a few per cent of 17-digit numbers) is stored
# as signif(value, d) for the largest d <= 16 that round-trips (a change of at most 1e-15 relative).
# ------------------------------------------------------------------
json_safe <- function(x) {
  txt <- function(y) {
    for (d in 15:17) { s <- sprintf(paste0("%.", d, "g"), y); if (as.numeric(s) == y) return(s) }
    s
  }
  rt <- function(y) isTRUE(as.numeric(jsonlite::parse_json(paste0("[", txt(y), "]"))[[1]]) == y)
  if (rt(x)) return(x)
  for (d in 16:10) { y <- signif(x, d); if (rt(y)) return(y) }
  stop("json_safe: no round-tripping text for ", format(x, digits = 17))
}

collect_numbers07 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 04 first")
    readRDS(p)
  }
  tr <- rd("truth07.rds"); md <- rd("res_md.rds"); g2 <- rd("res_gmm_mc.rds"); g3 <- rd("res_weak.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers07: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers07: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers07: ", name, " is not finite (", value, ")")
    nums[[name]] <<- json_safe(value)
    invisible(NULL)
  }
  put_ci <- function(prefix, est, se, ci, with_est = TRUE) {
    if (with_est) put(paste0(prefix, "_est"), est)
    put(paste0(prefix, "_se"), se)
    put(paste0(prefix, "_lb"), ci[[1]]); put(paste0(prefix, "_ub"), ci[[2]])
  }
  P <- tr$params
  th0 <- tr$G1$theta

  ## ---- dgp: seeds, parameters, truth ------------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  put("dgp_theta", th0)
  P1 <- P$G1
  mk <- tolower(P1$markets)
  put("dgp_G1_r", length(mk)); put("dgp_G1_n_total", sum(P1$n)); put("dgp_G1_het", P1$het); put("dgp_G1_R_mc", P1$R_mc)
  put("dgp_G1_n_min", min(P1$n)); put("dgp_G1_n_max", max(P1$n))
  put("dgp_G1_delta_min", min(P1$delta)); put("dgp_G1_delta_max", max(P1$delta))
  put("dgp_G1_p0_min", min(P1$p0)); put("dgp_G1_p0_max", max(P1$p0))
  put("dgp_G1_q0_min", min(P1$q0)); put("dgp_G1_q0_max", max(P1$q0))
  put("dgp_G1_sigma_min", min(P1$sigma)); put("dgp_G1_sigma_max", max(P1$sigma))
  put("dgp_G1_se_md", tr$G1$se_md)
  put("dgp_G1_se_market_min", min(tr$G1$se_market)); put("dgp_G1_se_market_max", max(tr$G1$se_market))
  put("dgp_G1_weight_min", min(tr$G1$weights)); put("dgp_G1_weight_max", max(tr$G1$weights))
  put("dgp_G1_weight_ratio_max_min", max(tr$G1$weights) / min(tr$G1$weights))
  for (j in seq_along(mk)) {
    m <- mk[j]
    put(paste0("dgp_G1_", m, "_n"), P1$n[j]); put(paste0("dgp_G1_", m, "_delta"), P1$delta[j])
    put(paste0("dgp_G1_", m, "_p0"), P1$p0[j]); put(paste0("dgp_G1_", m, "_p1"), P1$p0[j] * (1 + P1$delta[j]))
    put(paste0("dgp_G1_", m, "_q0"), P1$q0[j]); put(paste0("dgp_G1_", m, "_sigma"), P1$sigma[j])
    put(paste0("dgp_G1_", m, "_income"), P1$income[j]); put(paste0("dgp_G1_", m, "_households"), P1$households[j])
    put(paste0("dgp_G1_", m, "_se"), tr$G1$se_market[j]); put(paste0("dgp_G1_", m, "_weight"), tr$G1$weights[j])
    put(paste0("dgp_G1h_", m, "_theta_m"), tr$G1h$theta_m[j])
  }
  put("dgp_G1h_kappa", P1$kappa)
  for (k in c("theta_star", "theta_mean", "theta_hh", "lambda_J", "power_J_01")) put(paste0("dgp_G1h_", k), tr$G1h[[k]])
  put("dgp_G1h_theta_m_min", min(tr$G1h$theta_m)); put("dgp_G1h_theta_m_max", max(tr$G1h$theta_m))
  put("dgp_G1h_theta_star_minus_mean", tr$G1h$theta_star - tr$G1h$theta_mean)

  P2 <- P$G2
  K2 <- length(P2$step)
  for (k in c("n", "sigma_v", "rho", "sigma_eps", "het", "het_on", "p_base", "q_base", "R_mc", "R_grid")) put(paste0("dgp_G2_", k), P2[[k]])
  put("dgp_G2_K", K2); put("dgp_G2_r", K2 + 1); put("dgp_G2_k", 2); put("dgp_G2_df", K2 - 1)
  put("dgp_G2_sigma_eps_deep", P2$sigma_eps * (1 + P2$het)); put("dgp_G2_pack_price", 8 * P2$p_base)
  put("dgp_G2_beta0", tr$G2$beta0)
  for (j in seq_along(P2$grid_n)) put(sprintf("dgp_G2_grid_n%d", j), P2$grid_n[j])
  for (j in seq_len(K2)) {
    pz <- sprintf("dgp_G2_z%d", j)
    put(paste0(pz, "_step"), P2$step[j]); put(paste0(pz, "_half_step"), P2$step[j] / 2); put(paste0(pz, "_share"), P2$share[j])
    put(paste0(pz, "_gamma"), P2$gamma[j]); put(paste0(pz, "_arm_effect"), P2$gamma[j] * P2$step[j])
    put(paste0(pz, "_eff_sd"), tr$G2$eff_sd[j])
  }
  for (k in c("mu2", "r2_first_stage", "ols_plim", "se_2sls", "se_eff", "se_identity")) put(paste0("dgp_G2_", k), tr$G2[[k]])
  put("dgp_G2_ols_plim_minus_theta", tr$G2$ols_plim - th0)
  put("dgp_G2_se_eff_over_2sls", tr$G2$se_eff / tr$G2$se_2sls); put("dgp_G2_se_identity_over_eff", tr$G2$se_identity / tr$G2$se_eff)
  # Common asymptotic SE of efficient GMM and 2SLS if every experiment were balanced (share 0.5), by pop_moments
  # (gmm_fns.R, read into a private environment). G2 is not balanced because its deep-discount test (z2) puts only 20% of
  # regions in the upper arm and is four times as noisy there, so S_0 is not proportional to E[zz'] and efficient GMM beats 2SLS.
  fns <- new.env(); sys.source(file.path(dir, "gmm_fns.R"), envir = fns)
  bal <- fns$pop_moments(modifyList(P2, list(share = rep(0.5, K2))))
  stopifnot(abs(bal$se_eff / bal$se_2sls - 1) < 1e-8)
  put("dgp_G2_se_balanced", bal$se_eff)
  mom <- c("intercept", paste0("z", seq_len(K2)))
  for (j in seq_along(mom)) put(paste0("dgp_G2_lambda0_", mom[j]), tr$G2$sensitivity[j])

  P3 <- P$G3
  K3 <- length(P3$step)
  for (k in c("n", "sigma_v", "rho", "sigma_eps", "R_mc")) put(paste0("dgp_G3_", k), P3[[k]])
  put("dgp_G3_K", K3); put("dgp_G3_beta0", tr$G3$beta0)
  for (j in seq_len(K3)) {
    pz <- sprintf("dgp_G3_z%d", j)
    put(paste0(pz, "_step"), P3$step[j]); put(paste0(pz, "_half_step"), P3$step[j] / 2)
    put(paste0(pz, "_half_step_cents"), 100 * P3$step[j] / 2); put(paste0(pz, "_gamma"), P3$gamma[j])
    put(paste0(pz, "_eff_sd"), tr$G3$eff_sd[j])
  }
  for (k in c("mu2", "F_pop", "ols_plim", "r2_first_stage")) put(paste0("dgp_G3_", k), tr$G3[[k]])
  put("dgp_G3_chi2_crit", qchisq(0.95, K3)); put("dgp_G2_chi2_crit_S", qchisq(0.95, K2)); put("dgp_G2_chi2_crit_J", qchisq(0.95, K2 - 1))
  put("dgp_G1_chi2_crit_J", qchisq(0.95, length(mk) - 1))

  ## ---- md: blocks G1, G1h, pooled GMM -------------------------------------------------------------
  for (b in c("G1", "G1h")) {
    x <- md[[b]]
    pre <- paste0("md_", b)
    put_ci(pre, x$theta, x$se, x$ci)
    put(paste0(pre, "_J"), x$J); put(paste0(pre, "_J_p"), x$J_p); put(paste0(pre, "_df"), x$df)
    put(paste0(pre, "_ivw_mean"), x$ivw_mean); put(paste0(pre, "_n_market_ci_cover"), x$n_market_ci_cover)
    put_ci(paste0(pre, "_equal"), x$equal$theta, x$equal$se, x$equal$ci)
    put(paste0(pre, "_se_equal_over_md"), x$equal$se / x$se)
    for (j in seq_along(mk)) {
      s <- x$slopes[j, ]
      stopifnot(tolower(as.character(s$market)) == mk[j])
      pm <- paste0(pre, "_", mk[j])
      put(paste0(pm, "_pi_hat"), s$pi_hat); put(paste0(pm, "_se"), s$se)
      put(paste0(pm, "_lb"), s$lo); put(paste0(pm, "_ub"), s$hi); put(paste0(pm, "_weight"), s$weight)
    }
  }
  put("md_G1_covers_truth", md$G1$covers_truth); put("md_G1_minus_truth", md$G1$theta - th0)
  for (k in c("mc_R", "mc_reject_J", "mc_reject_J_01", "mc_reject_J_10", "mc_mean_J", "mc_cover", "mc_mean_theta",
              "mc_sd_theta", "mc_mean_se")) put(paste0("md_G1_", k), md$G1[[k]])
  put("md_G1_mc_bias", md$G1$mc_mean_theta - th0)
  put("md_G1_mc_mcse_reject_J", sqrt(0.05 * 0.95 / md$G1$mc_R))
  for (k in c("covers_star", "covers_mean", "covers_hh", "equal_covers_mean")) put(paste0("md_G1h_", k), md$G1h[[k]])
  put("md_G1h_minus_star", md$G1h$theta - tr$G1h$theta_star); put("md_G1h_minus_mean", md$G1h$theta - tr$G1h$theta_mean)
  g <- md$gmm$two_step
  put("md_gmm_two_step_est", g$theta); put("md_gmm_two_step_se", g$se); put("md_gmm_two_step_J", g$J)
  put("md_gmm_two_step_J_p", g$J_p); put("md_gmm_two_step_df", g$df); put("md_gmm_first_step_est", g$theta1)
  put("md_gmm_n_obs", g$n); put("md_gmm_n_moments", g$n_moments); put("md_gmm_n_par", g$n_par)
  u <- md$gmm$unrestricted
  put("md_gmm_unrestricted_est", u$theta); put("md_gmm_unrestricted_se", u$se); put("md_gmm_unrestricted_J", u$J)
  put("md_vs_gmm", md$md_vs_gmm); put("md_vs_gmm_abs", abs(md$md_vs_gmm)); put("md_vs_gmm_over_se", abs(md$md_vs_gmm) / md$G1$se)
  put("md_vs_gmm_unrestricted", md$md_vs_gmm_unrestricted); put("md_vs_gmm_unrestricted_abs", abs(md$md_vs_gmm_unrestricted))
  put("md_J_vs_gmm_unrestricted", md$J_md_vs_gmm_unrestricted); put("md_se_vs_gmm_unrestricted", md$se_md_vs_gmm_unrestricted)
  put("md_vs_gmm_first_step", md$G1$theta - g$theta1)

  ## ---- gmm: block G2 ------------------------------------------------------------------------------
  S <- g2$sample
  for (k in names(S$est)) {
    e <- S$est[[k]]
    pre <- paste0("gmm_G2_", k)
    put_ci(pre, e$est, e$se, e$ci, with_est = k != "windmeijer")
    put(paste0(pre, "_covers"), e$covers)
  }
  put("gmm_G2_beta0_est", S$beta0_two_step)
  put("gmm_G2_J", S$J); put("gmm_G2_J_p", S$J_p); put("gmm_G2_J_df", S$J_df)
  put("gmm_G2_J_iterated", S$J_iterated); put("gmm_G2_J_cue", S$J_cue); put("gmm_G2_iter_iterated", S$iter_iterated)
  put("gmm_G2_ols_est", S$ols$est); put("gmm_G2_ols_se", S$ols$se)
  put("gmm_G2_first_stage_F", S$first_stage_F); put("gmm_G2_first_stage_r2", S$first_stage_r2)
  put("gmm_G2_se_ratio_windmeijer", S$se_ratio_windmeijer)
  put("gmm_G2_se_identity_over_eff", g2$se_identity / g2$se_eff); put("gmm_G2_se_one_step_over_eff", g2$se_one_step / g2$se_eff)
  for (j in seq_along(mom)) {
    put(paste0("gmm_G2_lambda_", mom[j]), S$lambda_theta[j]); put(paste0("gmm_G2_lambda_std_", mom[j]), S$lambda_theta_std[j])
  }
  put("gmm_G2_s_set_lb", S$s_set[["lo"]]); put("gmm_G2_s_set_ub", S$s_set[["hi"]]); put("gmm_G2_s_min", S$s_min)
  put("gmm_G2_s_argmin", S$s_argmin); put("gmm_G2_s_df", S$s_df)
  est_names <- c("one_step", "identity", "two_step_naive", "windmeijer", "iterated", "cue")
  put_mc <- function(pre, s, full = TRUE) {
    put(paste0(pre, "_R"), s$R)
    for (k in est_names) {
      put(paste0(pre, "_cover_", k), s[[paste0("cover_", k)]])
      put(paste0(pre, "_sd_", k), s[[paste0("sd_", k)]]); put(paste0(pre, "_mean_se_", k), s[[paste0("mean_se_", k)]])
      if (full) {
        put(paste0(pre, "_mean_", k), s[[paste0("mean_", k)]]); put(paste0(pre, "_median_", k), s[[paste0("median_", k)]])
        put(paste0(pre, "_bias_", k), s[[paste0("mean_", k)]] - th0)
        put(paste0(pre, "_mean_se_over_sd_", k), s[[paste0("mean_se_", k)]] / s[[paste0("sd_", k)]])
      }
    }
    for (k in c("reject_J_two_step", "reject_J_iterated", "reject_J_cue")) put(paste0(pre, "_", k), s[[k]])
  }
  put_mc("gmm_mc", g2$mc)
  put("gmm_mc_n", g2$mc$n); put("gmm_mc_mean_iter_iterated", g2$mc$mean_iter_iterated)
  put("gmm_mc_sd_identity_over_two_step", g2$mc$sd_identity / g2$mc$sd_two_step_naive)
  put("gmm_mc_sd_one_step_over_two_step", g2$mc$sd_one_step / g2$mc$sd_two_step_naive)
  for (i in seq_len(nrow(g2$grid))) put_mc(sprintf("gmm_grid_n%d", g2$grid$n[i]), g2$grid[i, ], full = FALSE)
  put("gmm_mcse_coverage_R1000", sqrt(0.95 * 0.05 / 1000)); put("gmm_mcse_coverage_R2000", sqrt(0.95 * 0.05 / 2000))

  ## ---- weak: block G3 ------------------------------------------------------------------------------
  W <- g3$sample
  put_ci("weak_G3_two_step", W$two_step$est, W$two_step$se, W$two_step$ci)
  put_ci("weak_G3_windmeijer", W$windmeijer$est, W$windmeijer$se, W$windmeijer$ci, with_est = FALSE)
  put_ci("weak_G3_2sls", W$tsls$est, W$tsls$se, W$tsls$ci)
  put_ci("weak_G3_cue", W$cue$est, W$cue$se, W$cue$ci)
  put("weak_G3_ols_est", W$ols$est)
  put("weak_G3_wald_covers", W$wald_covers); put("weak_G3_s_covers", W$s_covers)
  for (k in c("S0", "S0_p", "S_inf", "unbounded", "first_stage_F", "J", "crit", "df", "s_set_pieces", "s_min", "s_argmin"))
    put(paste0("weak_G3_", k), W[[k]])
  put("weak_G3_s_set_lb", W$s_set[["lo"]]); put("weak_G3_s_set_ub", W$s_set[["hi"]])
  put("weak_G3_s_set_length", W$s_set[["hi"]] - W$s_set[["lo"]])
  put("weak_G3_wald_length", W$two_step$ci[["hi"]] - W$two_step$ci[["lo"]])
  put("weak_G3_grid_lo", W$grid[1]); put("weak_G3_grid_hi", W$grid[2])
  for (k in names(g3$cover)) put(paste0("weak_mc_cover_", k), g3$cover[[k]])
  for (k in names(g3$mc)) put(paste0("weak_mc_", k), g3$mc[[k]])
  put("weak_mc_bias_median_two_step", g3$mc$median_two_step - th0)

  nums
}
