# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers08(dir = "."): every quotable number of chapter 8 as a
# flat named list of finite numeric scalars (07_summary.R writes it to numbers.json with
# chapter_kit.R's write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads halcyon.rds (01_dgp.R) and res_delta.rds (02), res_cluster.rds (03), res_hac.rds (04),
# res_test.rds (05) in `dir`.
#
# Name scheme: <prefix>_<block>_<field>. Prefixes: dgp (01_dgp.R: true parameters and sample facts),
# delta (02, WTP ratio: Delta method against the pairs bootstrap; main = strong denominator,
# weak = weak denominator), cluster (03, city rollout: G40 = one sample with G = 40 cities,
# size_G40 / size_G10 = Monte Carlo rejection rates of a true null at nominal 5%, G10 = one G = 10
# sample with a real effect), hac (04, AR(1) daily series), test (05, logit trinity and the Wald
# reparameterization demo), const (fixed constants the text prints). Conventions:
#   *_se_iid, *_se_hc1     default OLS and heteroskedasticity-robust HC1 standard errors
#   *_se_cr1, *_se_cr2     cluster-robust (city) CR1 (vcovCL type HC1) and CR2 (type HC2)
#   *_t_*                  estimate / that standard error
#   *_ci_*_lb, *_ci_*_ub   95% interval bounds (delta: est -/+ 1.96 se; perc: bootstrap 2.5% and
#                          97.5% quantiles; pt: percentile-t, studentized with each resample's Delta SE)
#   *_cover_*              Monte Carlo coverage of nominal 95% intervals (2000 series)
#   cluster_G10_p_wcb      wild cluster bootstrap p (Rademacher, null imposed). With G = 10 clusters there
#                          are 2^10 = 1024 sign patterns, so wild_cluster_p() enumerates all of them (exact,
#                          no seed) and 76 of 1024 count (0.0742). Two patterns (all +1 and all -1) give
#                          exactly |t_obs|; they count as ">=" up to a relative tolerance of 1e-8, so the
#                          last bits of floating point cannot flip them (the nearest non-tie is 0.0078
#                          away in |t|). The earlier version drew B = 1999 random patterns, where such
#                          ties fell either way: 145 of 1999 (0.0725) in the frozen v2 run and 148 of
#                          1999 (0.0740) in R 4.6.1. Every other stored v2 value agrees to a relative 1e-10.
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json.
# ------------------------------------------------------------------
collect_numbers08 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 05 first")
    readRDS(p)
  }
  H <- rd("halcyon.rds"); dl <- rd("res_delta.rds"); cl <- rd("res_cluster.rds")
  hc <- rd("res_hac.rds"); te <- rd("res_test.rds")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers08: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers08: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers08: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }

  ## ---- dgp: true parameters and sample facts (01_dgp.R) -------------------------------------
  tr <- H$truth
  put("dgp_wtp_true", tr$WTP); put("dgp_b_time", tr$b_time); put("dgp_b_fee", tr$b_fee)
  put("dgp_tau", tr$tau); put("dgp_icc_rho", tr$rho); put("dgp_moulton", tr$moulton)
  put("dgp_beta_ts", tr$beta_ts); put("dgp_G", tr$G); put("dgp_n_c", tr$n_c)
  mk <- formals(H$mkclus)
  put("dgp_sigma_a", mk$sigma_a); put("dgp_sigma_e", mk$sigma_e)
  put("dgp_n_wtp", nrow(H$wtp_main)); put("dgp_n_wtp_weak", nrow(H$wtp_weak))
  put("dgp_n_sessions", nrow(H$clus_main)); put("dgp_T", nrow(H$ts_dat))

  ## ---- delta: WTP ratio, Delta method against the bootstrap (02) ----------------------------
  for (blk in c("main", "weak")) {
    b <- dl[[blk]]
    pre <- paste0("delta_", blk, "_")
    put(paste0(pre, "wtp"), b$wtp_hat); put(paste0(pre, "se_delta"), b$se_delta)
    put(paste0(pre, "se_boot"), b$se_boot); put(paste0(pre, "skew"), b$skew)
    put(paste0(pre, "ci_delta_lb"), b$ci_delta[1]); put(paste0(pre, "ci_delta_ub"), b$ci_delta[2])
    put(paste0(pre, "ci_perc_lb"), b$ci_perc[1]); put(paste0(pre, "ci_perc_ub"), b$ci_perc[2])
    put(paste0(pre, "ci_pt_lb"), b$ci_pt[1]); put(paste0(pre, "ci_pt_ub"), b$ci_pt[2])
    put(paste0(pre, "n_boot"), length(b$wstar))
  }
  put("delta_weak_b_fee", dl$base_w["bf"]); put("delta_weak_se_b_fee", dl$base_w["se_bf"])
  put("delta_weak_b_time", dl$base_w["bt"])

  ## ---- cluster: one G = 40 rollout, size study, one G = 10 sample (03) ----------------------
  put("cluster_G40_tau_hat", cl$tau_hat)
  for (k in c("iid", "hc1", "cr1", "cr2")) {
    se <- cl[[paste0("se_", k)]]
    put(paste0("cluster_G40_se_", k), se); put(paste0("cluster_G40_t_", k), cl$tau_hat / se)
  }
  # printed as 4 x 10^-43; kept to 13 significant digits because jsonlite 2.0.0 reads the 16-digit text of
  # this number (near 1e-43) back one ulp off, which write_numbers_json's read-back check rejects
  put("cluster_G40_p_iid", signif(2 * pnorm(-abs(cl$tau_hat / cl$se_iid)), 13))
  put("cluster_G40_p_cr1_t", 2 * pt(-abs(cl$tau_hat / cl$se_cr1), df = tr$G - 1))
  put("cluster_G40_inflation_cr1_iid", cl$se_cr1 / cl$se_iid)
  for (k in c("iid", "cr1_z", "cr1_t")) {
    put(paste0("cluster_size_G40_", k), cl$s40[[k]]); put(paste0("cluster_size_G10_", k), cl$s10[[k]])
  }
  put("cluster_G10_t", cl$t10); put("cluster_G10_p_z", cl$p10_z); put("cluster_G10_p_t", cl$p10_t)
  put("cluster_G10_p_wcb", cl$p10_wcb)

  ## ---- hac: AR(1) daily series (04) ---------------------------------------------------------
  put("hac_b", hc$b); put("hac_se_iid", hc$se_iid); put("hac_se_nw", hc$se_nw); put("hac_lag", hc$L)
  put("hac_t_iid", hc$b / hc$se_iid); put("hac_t_nw", hc$b / hc$se_nw)
  put("hac_ratio_nw_iid", hc$se_nw / hc$se_iid)
  put("hac_cover_iid", hc$cov_iid); put("hac_cover_nw", hc$cov_nw)

  ## ---- test: logit trinity and Wald reparameterization (05) ---------------------------------
  put("test_g1_hat", te$g1_hat)
  put("test_W", te$W); put("test_LR", te$LR); put("test_LM", te$LM)
  put("test_p_W", te$pW); put("test_p_LR", te$pLR); put("test_p_LM", te$pLM)
  put("test_repar_W_direct", te$W_direct); put("test_repar_W_log", te$W_repar)
  put("test_repar_p_direct", te$pWd); put("test_repar_p_log", te$pWr)
  put("test_repar_null_g1", 0.5)

  ## ---- const: fixed constants printed in the text --------------------------------------------
  put("const_z_0p975", qnorm(0.975))            # 1.96 in theta_k -/+ 1.96 SE
  put("const_ci_lower_prob", 0.025)             # 2.5% quantile of a percentile interval
  put("const_ci_upper_prob", 0.975)             # 97.5% quantile
  nums
}
