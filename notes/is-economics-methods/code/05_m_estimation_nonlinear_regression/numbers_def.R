# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers05(dir = "."): every quotable number of chapter 5 as a
# flat named list of finite numeric scalars (07_summary.R writes it to numbers.json with
# chapter_kit.R's write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth05.rds and res_logs.rds, res_nls.rds, res_twostep.rds, res_mc.rds in `dir`.
#
# Name scheme: <prefix>_<block or object>_<field>. Prefixes: dgp (01_dgp.R: parameters, truth,
# realised sample facts), logs (02, block M1), nls (03, blocks M2 and M2m), twostep (04, block M3),
# mc (05, block M4). Conventions:
#   *_est, *_se          point estimate and HC0 sandwich standard error (unless named otherwise:
#                        *_se_model, *_se_naive, *_se_semi, *_se_quasi, *_se_hetero)
#   *_lb, *_ub           95% interval est -/+ qnorm(0.975) se of the same row
#   *_covers_*           1 when that interval contains the named value, 0 otherwise
#   cover_*              Monte Carlo coverage rate of nominal 95% intervals
#   usd, cents, kusd     sales in USD, cents, thousands of USD; a<k> is the unit multiplier a = k
#                        (decimal point written p: a0p001 = thousands of USD)
#   theta1_0p8, size_3   Monte Carlo grid points (size_poisson: the Poisson case, size = Inf)
#   M4 overdispersion    avg_cond_var_to_mean = E[Var(y | x)] / E[y]; pearson_dispersion = the Pearson
#                        statistic and its target E[Var(y | x) / E[y | x]]; marginal_var_to_mean = Var(y) / E[y]
#                        (includes Var(mu) from x). Three different quantities: never compare across them.
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json.
# ------------------------------------------------------------------
collect_numbers05 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 05 first")
    readRDS(p)
  }
  tr <- rd("truth05.rds"); lg <- rd("res_logs.rds"); nl <- rd("res_nls.rds")
  ts <- rd("res_twostep.rds"); mc <- rd("res_mc.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers05: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers05: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers05: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  put_est <- function(prefix, est, se, ci = NULL) {
    put(paste0(prefix, "_est"), est); put(paste0(prefix, "_se"), se)
    if (is.null(ci)) ci <- c(est - qnorm(0.975) * se, est + qnorm(0.975) * se)
    put(paste0(prefix, "_lb"), ci[1]); put(paste0(prefix, "_ub"), ci[2])
  }
  pid <- function(x) gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))
  P <- tr$params

  ## ---- dgp: seeds and parameters ---------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  P1 <- P$M1
  for (k in c("n", "k", "sigma", "beta_tenure", "delta")) put(paste0("dgp_M1_", k), P1[[k]])
  for (j in seq_along(P1$categories)) {
    cn <- P1$categories[j]
    put(paste0("dgp_M1_pi_", cn), P1$pi[j]); put(paste0("dgp_M1_p0_", cn), P1$p0[j])
    put(paste0("dgp_M1_p1_", cn), P1$k * P1$p0[j]); put(paste0("dgp_M1_median0_", cn), P1$median0[j])
  }
  put("dgp_M1_log_k", log(P1$k)); put("dgp_M1_pct_intensive", exp(P1$delta) - 1)
  for (k in c("log_ratio", "ratio", "share_sell0", "share_sell1", "ext_effect", "mean_y0_usd", "mean_y1_usd"))
    put(paste0("dgp_M1_", k), tr$M1[[k]])
  put("dgp_M1_pct_effect", tr$M1$ratio - 1)
  put("dgp_M1_share_zero0", 1 - tr$M1$share_sell0); put("dgp_M1_share_zero1", 1 - tr$M1$share_sell1)
  pp <- tr$M1$pop
  for (k in c("n", "log_ratio", "mean_y0_usd", "mean_y1_usd", "share_sell0", "share_sell1", "ate_log1p_kusd",
              "ate_log1p_usd", "ate_log1p_cents", "ate_asinh_usd", "ate_asinh_cents"))
    put(paste0("dgp_M1_pop_", k), pp[[k]])
  put("dgp_M1_pop_ate_log1p_cents_minus_usd", pp$ate_log1p_cents - pp$ate_log1p_usd)
  for (k in names(lg$sample)) put(paste0("dgp_M1_sample_", k), lg$sample[[k]])

  P2 <- P$M2
  for (k in c("n", "listings_median", "sd_log_listings", "p_free_ship", "sigma_eta")) put(paste0("dgp_M2_", k), P2[[k]])
  for (j in 1:3) put(sprintf("dgp_M2_theta%d", j - 1), P2$theta[j])
  put("dgp_M2_mean_at_median_usd", 1000 * exp(P2$theta[1])); put("dgp_M2_pct_free_ship", exp(P2$theta[3]) - 1)
  P2m <- P$M2m
  for (k in c("n", "sd_x", "a0", "a1", "a2", "sigma_eta")) put(paste0("dgp_M2m_", k), P2m[[k]])
  put("dgp_M2m_pseudo_b0", tr$M2m$pseudo_vec[1]); put("dgp_M2m_pseudo_b1", tr$M2m$pseudo_vec[2])
  put("dgp_M2m_pseudo_ppml_b0", tr$M2m$pseudo_ppml_vec[1]); put("dgp_M2m_pseudo_ppml_b1", tr$M2m$pseudo_ppml_vec[2])
  put("dgp_M2m_large_N", tr$M2m$large$N)
  put("dgp_M2m_large_b0", tr$M2m$large$pseudo_vec[1]); put("dgp_M2m_large_b1", tr$M2m$large$pseudo_vec[2])
  put("dgp_M2m_pseudo_ppml_minus_nls_b1", tr$M2m$pseudo_ppml - tr$M2m$pseudo)
  P3 <- P$M3
  for (k in c("n", "sigma_e", "sigma_eta", "R_mc", "R_grid")) put(paste0("dgp_M3_", k), P3[[k]])
  for (j in 1:3) put(sprintf("dgp_M3_gamma%d", j - 1), P3$gamma[j])
  for (j in 1:2) put(sprintf("dgp_M3_theta%d", j - 1), P3$theta[j])
  put("dgp_M3_mean_at_median_usd", 1000 * exp(P3$theta[1])); put("dgp_M3_r2_first_stage", tr$M3$r2_first_stage)
  P4 <- P$M4
  for (k in c("n", "size", "R_mc", "R_grid")) put(paste0("dgp_M4_", k), P4[[k]])
  for (j in 1:2) put(sprintf("dgp_M4_beta%d", j - 1), P4$beta[j])
  put("dgp_M4_mean_at_average_rating", exp(P4$beta[1]))
  # three different overdispersion quantities (population values, see 01_dgp.R M4); never compare across them
  put("dgp_M4_avg_cond_var_to_mean", tr$M4$avg_cond_var_to_mean)     # E[Var(y | x)] / E[y], the Monte Carlo x-axis
  put("dgp_M4_pearson_dispersion", tr$M4$pearson_dispersion)         # E[Var(y | x) / E[y | x]], target of the Pearson statistic
  put("dgp_M4_marginal_var_to_mean", tr$M4$marginal_var_to_mean)     # Var(y) / E[y], includes Var(mu) from x

  ## ---- logs: block M1 ------------------------------------------------------------------------
  L <- lg$log1p
  for (u in c("kusd", "usd", "cents")) put_est(paste0("logs_log1p_", u), L[[u]], L[[paste0("se_", u)]], L[[paste0("ci_", u)]])
  put("logs_log1p_cents_minus_usd", L$cents - L$usd); put("logs_log1p_cents_over_usd", L$cents / L$usd)
  for (u in c("usd", "cents")) put_est(paste0("logs_asinh_", u), lg$asinh[[u]], lg$asinh[[paste0("se_", u)]])
  pu <- lg$ppml$usd
  put_est("logs_ppml_usd", pu$est, pu$se, pu$ci)
  put("logs_ppml_usd_se_model", pu$se_model); put("logs_ppml_usd_pct", pu$pct)
  put("logs_ppml_usd_pct_lb", pu$pct_ci[1]); put("logs_ppml_usd_pct_ub", pu$pct_ci[2])
  put("logs_ppml_usd_covers_truth", covers_num(pu$ci, tr$M1$log_ratio))
  put("logs_ppml_usd_minus_truth", pu$est - tr$M1$log_ratio)
  put_est("logs_ppml_cents", lg$ppml$cents$est, lg$ppml$cents$se, lg$ppml$cents$ci)
  put("logs_ppml_usd_minus_cents_abs", abs(pu$est - lg$ppml$cents$est))
  put_est("logs_ppml_raw", lg$ppml$raw$est, lg$ppml$raw$se, lg$ppml$raw$ci)
  put("logs_ppml_raw_covers_truth", covers_num(lg$ppml$raw$ci, tr$M1$log_ratio))
  put("logs_ppml_fepois_est", lg$ppml$fepois$est); put("logs_ppml_fepois_se_hetero", lg$ppml$fepois$se_hetero)
  put("logs_ppml_fepois_se_ratio", lg$ppml$fepois$se_ratio)
  put_est("logs_ext", lg$ext$est, lg$ext$se, lg$ext$ci)
  put("logs_ext_times_log100", lg$ext$est * log(100))
  put("logs_ext_covers_truth", covers_num(lg$ext$ci, tr$M1$ext_effect))
  put_est("logs_logpos", lg$logpos$est, lg$logpos$se, lg$logpos$ci); put("logs_logpos_n", lg$logpos$n)
  # in this DGP the seller-only comparison identifies delta exactly (no selection bias; see 02_logs.R)
  put("logs_logpos_covers_delta", covers_num(lg$logpos$ci, P1$delta)); put("logs_logpos_minus_delta", lg$logpos$est - P1$delta)
  put("logs_log1p_usd_covers_pop_ate", covers_num(L$ci_usd, tr$M1$pop$ate_log1p_usd))
  put("logs_log1p_cents_covers_pop_ate", covers_num(L$ci_cents, tr$M1$pop$ate_log1p_cents))
  sw <- lg$sweep
  for (k in -3:4) {
    j <- which.min(abs(log10(sw$a) - k))
    stopifnot(abs(log10(sw$a[j]) - k) < 1e-9)
    nm <- paste0("logs_sweep_a", pid(10^k))
    put(paste0(nm, "_est"), sw$est[j]); put(paste0(nm, "_se"), sw$se[j]); put(paste0(nm, "_ate_pop"), sw$ate_pop[j])
  }
  put("logs_sweep_min_est", min(sw$est)); put("logs_sweep_max_est", max(sw$est))

  ## ---- nls: blocks M2 and M2m ----------------------------------------------------------------
  for (j in 1:3) {
    nm <- sprintf("nls_M2_theta%d", j - 1)
    put_est(nm, nl$theta[j], nl$se_robust[j], nl$ci[[j]])
    put(paste0(nm, "_se_model"), nl$se_model[j]); put(paste0(nm, "_se_full"), nl$se_full[j])
    put(paste0(nm, "_nls_est"), nl$nls_coef[j]); put(paste0(nm, "_covers_truth"), nl$covers[j])
    put(paste0(nm, "_se_robust_over_model"), nl$se_robust[j] / nl$se_model[j])
  }
  put("nls_M2_gn_iter", nl$gn$iter); put("nls_M2_nls_iter", nl$nls_iter); put("nls_M2_gn_minus_nls_max", nl$max_abs_diff)
  put("nls_M2_sigma2", nl$sigma2)
  ms <- nl$misspec
  put("nls_M2m_b0_est", ms$theta[1])
  put_est("nls_M2m_b1", ms$est, ms$se, ms$ci)
  put("nls_M2m_b1_se_semi", ms$se_semi); put("nls_M2m_b1_lb_semi", ms$ci_semi[1]); put("nls_M2m_b1_ub_semi", ms$ci_semi[2])
  put("nls_M2m_b1_se_model", ms$se_model); put("nls_M2m_gn_iter", ms$iter)
  put("nls_M2m_b1_covers_pseudo", ms$covers); put("nls_M2m_b1_covers_a1", ms$covers_a1)
  put("nls_M2m_b1_minus_pseudo", ms$est - tr$M2m$pseudo)
  put("nls_M2m_ppml_b0_est", ms$ppml$theta[1])
  put_est("nls_M2m_ppml_b1", ms$ppml$est, ms$ppml$se, ms$ppml$ci)
  put("nls_M2m_ppml_b1_covers_pseudo_ppml", ms$ppml$covers); put("nls_M2m_ppml_b1_covers_nls_pseudo", ms$ppml$covers_nls_pseudo)

  ## ---- twostep: block M3 --------------------------------------------------------------------
  put("twostep_theta0_est", ts$theta[1]); put("twostep_theta1_est", ts$est)
  for (j in 1:3) put(sprintf("twostep_gamma%d_est", j - 1), ts$gamma[j])
  put("twostep_se_naive", ts$se_naive); put("twostep_se_corrected", ts$se_corrected); put("twostep_se_ratio", ts$se_ratio)
  put("twostep_naive_lb", ts$ci_naive[1]); put("twostep_naive_ub", ts$ci_naive[2])
  put("twostep_corrected_lb", ts$ci_corrected[1]); put("twostep_corrected_ub", ts$ci_corrected[2])
  put("twostep_covers_naive", ts$covers_naive); put("twostep_covers_corrected", ts$covers_corrected)
  put("twostep_r2_first_stage", ts$r2_first_stage); put("twostep_gn_iter", ts$iter)
  for (i in 1:2) for (j in 1:3) put(sprintf("twostep_F2_theta%d_gamma%d", i - 1, j - 1), ts$F2[i, j])
  put_est("twostep_direct", ts$direct$est, ts$direct$se, ts$direct$ci)
  put("twostep_direct_over_truth", ts$direct$est / P3$theta[2])
  put_est("twostep_oracle", ts$oracle$est, ts$oracle$se, ts$oracle$ci)
  put_mc3 <- function(prefix, s) {
    for (k in c("R", "cover_naive", "cover_corrected", "mean_est", "sd_est", "mean_se_naive", "mean_se_corrected"))
      put(paste0(prefix, "_", k), s[[k]])
  }
  put_mc3("twostep_mc", ts$mc)
  put("twostep_mc_bias", ts$mc$mean_est - P3$theta[2])
  for (i in seq_len(nrow(ts$grid))) put_mc3(paste0("twostep_grid_theta1_", pid(ts$grid$theta1[i])), ts$grid[i, ])

  ## ---- mc: block M4 -----------------------------------------------------------------------------
  put_mc4 <- function(prefix, s) {
    for (k in c("R", "avg_cond_var_to_mean", "cover_model", "cover_sandwich", "cover_quasi", "mean_est", "sd_est",
                "mean_se_model", "mean_se_sandwich", "mean_se_quasi"))
      put(paste0(prefix, "_", k), s[[k]])
  }
  put_mc4("mc", mc)
  put("mc_bias", mc$mean_est - P4$beta[2]); put("mc_sd_over_mean_se_model", mc$sd_est / mc$mean_se_model)
  for (i in seq_len(nrow(mc$grid))) {
    sz <- mc$grid$size[i]
    put_mc4(paste0("mc_grid_size_", if (is.infinite(sz)) "poisson" else pid(sz)), mc$grid[i, ])
  }
  s <- mc$sample
  put("mc_sample_beta0_est", s$beta[1]); put("mc_sample_beta1_est", s$est)
  put("mc_sample_se_model", s$se_model); put("mc_sample_se_sandwich", s$se_sandwich); put("mc_sample_se_quasi", s$se_quasi)
  put("mc_sample_model_lb", s$ci_model[1]); put("mc_sample_model_ub", s$ci_model[2])
  put("mc_sample_sandwich_lb", s$ci_sandwich[1]); put("mc_sample_sandwich_ub", s$ci_sandwich[2])
  put("mc_sample_model_covers_beta1", covers_num(s$ci_model, P4$beta[2]))
  put("mc_sample_sandwich_covers_beta1", covers_num(s$ci_sandwich, P4$beta[2]))
  put("mc_sample_beta1_minus_truth_in_sandwich_se", (s$est - P4$beta[2]) / s$se_sandwich)
  for (k in c("pearson_dispersion", "se_ratio", "marginal_var_to_mean", "mean_y", "share_zero")) put(paste0("mc_sample_", k), s[[k]])
  put("mc_mcse_coverage_R1000", sqrt(0.95 * 0.05 / 1000)); put("mc_mcse_coverage_R2000", sqrt(0.95 * 0.05 / 2000))

  nums
}

covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])
