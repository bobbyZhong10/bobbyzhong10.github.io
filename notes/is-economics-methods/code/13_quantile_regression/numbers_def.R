# ------------------------------------------------------------------
# numbers_def.R -- collect_numbers13(dir = "."): every quotable number of chapter 13 as a flat
# named list of finite numeric scalars (08_summary.R writes it to numbers.json with chapter_kit.R's
# write_numbers_json; tests/test_numbers_fresh.R checks the file against it).
#
# Reads truth13.rds, res_qte.rds and res_extra.rds in `dir`.
#
# Name scheme: <prefix>_<block>_<object>_<field>. Prefixes: dgp (01_dgp.R: seeds, parameters, truth,
# population checks), qte (02, blocks Q1 and Q1r), rif (03, block Q2: RIF, conditional QR, distribution
# regression), iv (04, block Q3), panel (05, block Q4), cic (06, blocks Q5 and Q5f), cross (06, the
# toy data for the chapter-20 cross-check). Conventions:
#   tNN                  quantile level tau = NN / 100 (t10, t25, t50, t75, t90)
#   *_est, *_se          point estimate and standard error (bootstrap unless named: *_se_ols is the OLS
#                        sandwich that treats the RIF as known; QR rows use summary.rq(se = "nid"))
#   *_lb, *_ub           95% interval of the same row (est -/+ 1.96 se; IV-QR: the inverted-test set)
#   *_covers_*           1 when that interval contains the named true value, 0 otherwise (IV-QR rows: _covers_struct
#                        is coverage of beta_D(tau), the IV-QR target; the LQTE, the complier effect, equals beta_D
#                        only at tau = 0.5, so _covers_lqte exists for iv_ivqr_t50 only)
#   units                Q1, Q5, Q5f: monthly revenue in USD thousands; Q2, Q3, Q4: log monthly revenue
#                        (log USD), effects in log points; *_usd: a level in USD
#   panel_mc_*           Monte Carlo (R = 1000) mean, bias, sd, rmse of estimator (pooled, dummies, canay)
#                        at panel length T5 or T20 (\mathcal{T} = 5 or 20)
# A non-finite value is an error (put() stops), as in chapter_kit.R's write_numbers_json. No p-values are
# stored (ruling 30 would store a tiny one as log10).
# ------------------------------------------------------------------
collect_numbers13 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop(f, " is missing: run 01 to 06 first")
    readRDS(p)
  }
  tr <- rd("truth13.rds"); q <- rd("res_qte.rds"); x <- rd("res_extra.rds")
  for (b in c("Q2", "Q3", "Q4", "Q5", "cross")) if (is.null(x[[b]])) stop("res_extra.rds has no element ", b, ": run 03 to 06")

  nums <- list()
  put <- function(name, value) {
    if (name %in% names(nums)) stop("collect_numbers13: duplicate key ", name)
    value <- as.numeric(value)
    if (length(value) != 1L) stop("collect_numbers13: ", name, " is not a scalar")
    if (!is.finite(value)) stop("collect_numbers13: ", name, " is not finite (", value, ")")
    nums[[name]] <<- value
    invisible(NULL)
  }
  tid <- function(t) sprintf("t%02d", round(100 * t))
  cov01 <- function(ci, v) as.numeric(v >= ci[1] && v <= ci[2])
  put_ci <- function(prefix, est, se, lb, ub) {
    put(paste0(prefix, "_est"), est)
    if (!is.null(se)) put(paste0(prefix, "_se"), se)
    put(paste0(prefix, "_lb"), lb); put(paste0(prefix, "_ub"), ub)
  }
  put_vec <- function(prefix, v, taus, suffix = "") for (j in seq_along(taus)) put(paste0(prefix, "_", tid(taus[j]), suffix), v[j])
  P <- tr$params; taus <- P$taus

  ## ---- dgp: seeds, parameters, truth ----------------------------------------------------------
  for (k in names(P$seed)) put(paste0("dgp_seed_", k), P$seed[[k]])
  put("dgp_n_taus", length(taus))
  P1 <- P$Q1
  for (k in c("n", "n_treat", "m", "b", "s", "v", "k", "logc", "rho_r", "n_pop")) put(paste0("dgp_Q1_", k), P1[[k]])
  put("dgp_Q1_n_control", P1$n - P1$n_treat)
  put("dgp_Q1_median_y0", exp(P1$m)); put("dgp_Q1_sd_log", sqrt(P1$v)); put("dgp_Q1_c", exp(P1$logc))
  put("dgp_Q1_corr_x_logy0", P1$b / sqrt(P1$v))
  T1 <- tr$Q1
  put("dgp_Q1_ate", T1$ate); put("dgp_Q1_ate_check_abs", abs(T1$ate_check)); put("dgp_Q1_mean_y", T1$mean_y)
  put("dgp_Q1_sd_y0", T1$sd_y0); put("dgp_Q1_sd_y1", T1$sd_y1)
  put("dgp_Q1_ind_min", T1$ind_min); put("dgp_Q1_y_at_ind_min", T1$y_at_ind_min); put("dgp_Q1_share_lose", T1$share_lose)
  put_vec("dgp_Q1_qte", T1$qte, taus); put_vec("dgp_Q1_q0", T1$q0, taus); put_vec("dgp_Q1_q1", T1$q1, taus)
  put_vec("dgp_Q1_qte_log", T1$qte_log, taus); put_vec("dgp_Q1_cqte_log", T1$cqte_log, taus)
  put_vec("dgp_Q1_ind", T1$ind_effect_quantiles, taus)
  put_vec("dgp_Q1_pop_qte", T1$pop$qte, taus); put("dgp_Q1_pop_ate", T1$pop$ate)
  put_vec("dgp_Q1_pop_ind", T1$pop$ind_effect_quantiles, taus)
  T1r <- tr$Q1r
  put_vec("dgp_Q1r_qte", T1r$qte, taus); put_vec("dgp_Q1r_ind", T1r$ind_effect_quantiles, taus)
  put_vec("dgp_Q1r_ind_minus_qte", T1r$ind_effect_quantiles - T1r$qte, taus)
  put("dgp_Q1r_share_lose", T1r$share_lose); put("dgp_Q1r_pop_share_lose", T1r$pop$share_lose)
  put_vec("dgp_Q1r_pop_ind", T1r$pop$ind_effect_quantiles, taus)
  put("dgp_Q1r_max_abs_ind_minus_qte", max(abs(T1r$ind_effect_quantiles - T1r$qte)))

  P2 <- P$Q2; T2 <- tr$Q2
  for (k in c("n", "n_fbm", "p", "mu", "g", "bx", "s0", "s1", "tau_uqe", "dp")) put(paste0("dgp_Q2_", k), P2[[k]])
  put("dgp_Q2_median0_usd", exp(P2$mu)); put("dgp_Q2_pct_fbm", exp(P2$g) - 1)
  put("dgp_Q2_sd0_marginal", sqrt(P2$bx^2 + P2$s0^2)); put("dgp_Q2_sd1_marginal", sqrt(P2$bx^2 + P2$s1^2))
  put_vec("dgp_Q2_q", T2$q, taus); put_vec("dgp_Q2_q_usd", exp(T2$q), taus)
  put_vec("dgp_Q2_uqpe", T2$uqpe, taus); put_vec("dgp_Q2_cqr", T2$cqr, taus); put_vec("dgp_Q2_marg_diff", T2$marg_diff, taus)
  put_vec("dgp_Q2_cf", T2$cf, taus); put_vec("dgp_Q2_uqe_discrete", T2$uqe_discrete, taus); put_vec("dgp_Q2_pct_cf", T2$pct_cf, taus)
  put("dgp_Q2_uqe", T2$uqe); put("dgp_Q2_uqe_minus_cqr", T2$uqe - T2$cqr[taus == P2$tau_uqe])

  P3 <- P$Q3; T3 <- tr$Q3
  for (k in c("n", "n_z1", "mu0", "dmu", "s0", "s1", "bx", "rho", "pi0", "pi1", "tau", "grid_step", "n_pop")) put(paste0("dgp_Q3_", k), P3[[k]])
  put("dgp_Q3_median0_usd", exp(P3$mu0)); put("dgp_Q3_pct_dmu", exp(P3$dmu) - 1)
  for (k in c("share_always", "share_complier", "share_never", "late_log", "lqte_median")) put(paste0("dgp_Q3_", k), T3[[k]])
  put_vec("dgp_Q3_lqte", T3$lqte, taus); put_vec("dgp_Q3_struct", T3$struct, taus)
  put_vec("dgp_Q3_pop_lqte", T3$pop$lqte, taus); put("dgp_Q3_pop_late_log", T3$pop$late_log); put("dgp_Q3_pop_n_complier", T3$pop$n_complier)

  P4 <- P$Q4; T4 <- tr$Q4
  for (k in c("N", "T_short", "T_long", "mu", "sa", "cx", "beta", "sig", "gam", "tau", "R")) put(paste0("dgp_Q4_", k), P4[[k]])
  put("dgp_Q4_beta_tau", T4$beta_tau); put_vec("dgp_Q4_beta", T4$beta_grid, taus)

  P5 <- P$Q5; T5 <- tr$Q5
  for (k in c("nc", "A", "g", "s", "kap", "mu1", "sd1", "delta", "nc_pop")) put(paste0("dgp_Q5_", k), P5[[k]])
  put_vec("dgp_Q5_qtt", T5$qtt, taus); put_vec("dgp_Q5_qdid", T5$qdid_pop, taus); put_vec("dgp_Q5_qdid_bias", T5$qdid_pop - T5$qtt, taus)
  put_vec("dgp_Q5_q11", T5$q11, taus); put_vec("dgp_Q5_cf", T5$cf, taus); put("dgp_Q5_att", T5$att)
  put_vec("dgp_Q5_pop_qtt", T5$pop$qtt, taus)
  P5f <- P$Q5f
  for (k in c("nc", "rho", "p_at1", "p_sw1", "p_at0", "tau", "nc_pop")) put(paste0("dgp_Q5f_", k), P5f[[k]])
  put("dgp_Q5f_lqte", T5$lqte); put_vec("dgp_Q5f_lqte", T5$lqte_grid, taus); put("dgp_Q5f_late_sw", T5$late_sw)
  put_vec("dgp_Q5f_pop_lqte", T5$pop$lqte_grid, taus); put("dgp_Q5f_pop_late_sw", T5$pop$late_sw); put("dgp_Q5f_pop_n_switcher", T5$pop$n_switcher)

  ## ---- qte: blocks Q1 and Q1r ------------------------------------------------------------------
  for (w in c("Q1", "Q1r")) {
    r <- if (w == "Q1") q else q$Q1r
    tq <- if (w == "Q1") tr$Q1 else tr$Q1r
    pre <- paste0("qte_", w)
    put_ci(paste0(pre, "_ols"), r$ols$est, r$ols$se, r$ols$ci[1], r$ols$ci[2])
    put(paste0(pre, "_ols_covers_zero"), cov01(r$ols$ci, 0))
    put_ci(paste0(pre, "_olsx"), r$ols_x$est, r$ols_x$se, r$ols_x$ci[1], r$ols_x$ci[2])
    for (j in seq_along(taus)) {
      a <- r$qte[j, ]; stopifnot(abs(a$tau - taus[j]) < 1e-12)
      put_ci(paste0(pre, "_", tid(taus[j])), a$est, a$se, a$lo, a$hi)
      put(paste0(pre, "_", tid(taus[j]), "_covers_truth"), cov01(c(a$lo, a$hi), tq$qte[j]))
      put(paste0(pre, "_log_", tid(taus[j])), r$qte_log$est[j])
      cq <- r$cqr_log[j, ]
      put_ci(paste0(pre, "_cqr_log_", tid(taus[j])), cq$est, cq$se, cq$lo, cq$hi)
    }
    for (k in names(r$sample)) put(paste0(pre, "_sample_", k), r$sample[[k]])
  }
  for (j in seq_along(taus)) put(paste0("qte_Q1_cqr_log_", tid(taus[j]), "_covers_truth"),
                                 cov01(c(q$cqr_log$lo[j], q$cqr_log$hi[j]), tr$Q1$cqte_log[j]))
  put("qte_Q1r_max_abs_qte_diff_vs_Q1", q$Q1r$max_abs_qte_diff_vs_Q1)

  ## figure grid (fig_13_qte): pointwise 95% intervals at tau = 0.05, ..., 0.95. The grid size and the number of
  ## levels whose interval covers the true QTE (the true QTE is the same in Q1 and Q1r), a 0/1 flag per level, and
  ## the one Q1 interval the text may quote as an example of a pointwise miss (tau = 0.60). Cf. rif_n_covers_uqpe_fine.
  tf <- P$taus_fine; tq_fine <- T1$fine$qte
  stopifnot(length(tq_fine) == length(tf), max(abs(T1$fine$tau - tf)) < 1e-12)
  for (w in c("Q1", "Q1r")) {
    f <- (if (w == "Q1") q else q$Q1r)$qte_fine
    stopifnot(nrow(f) == length(tf), max(abs(f$tau - tf)) < 1e-12)
    hit <- f$lo <= tq_fine & tq_fine <= f$hi
    put(paste0("qte_", w, "_n_covers_fine"), sum(hit)); put(paste0("qte_", w, "_n_fine"), nrow(f))
    for (j in seq_along(tf)) put(paste0("qte_", w, "_fine_", tid(tf[j]), "_covers_truth"), as.numeric(hit[j]))
  }
  j60 <- which(abs(tf - 0.6) < 1e-12); f1 <- q$qte_fine[j60, ]
  put_ci("qte_Q1_fine_t60", f1$est, f1$se, f1$lo, f1$hi); put("dgp_Q1_fine_qte_t60", tq_fine[j60])

  ## ---- rif: block Q2 --------------------------------------------------------------------------
  Q2 <- x$Q2; g2 <- Q2$grid
  for (t in taus) {
    a <- g2[abs(g2$tau - t) < 1e-12, ]; stopifnot(nrow(a) == 1)
    put_ci(paste0("rif_", tid(t)), a$rif_est, a$rif_se, a$rif_lo, a$rif_hi); put(paste0("rif_", tid(t), "_se_ols"), a$rif_se_ols)
    put_ci(paste0("rif_cqr_", tid(t)), a$cqr_est, a$cqr_se, a$cqr_lo, a$cqr_hi)
    put(paste0("rif_", tid(t), "_covers_uqpe"), cov01(c(a$rif_lo, a$rif_hi), a$uqpe_true))
    put(paste0("rif_cqr_", tid(t), "_covers_cond"), cov01(c(a$cqr_lo, a$cqr_hi), a$cqr_true))
  }
  put("rif_covers_uqe", cov01(Q2$rif$ci, tr$Q2$uqe)); put("rif_cqr_covers_uqe", cov01(Q2$cqr$ci, tr$Q2$uqe))
  put("rif_n_covers_uqpe_fine", sum(g2$rif_lo <= g2$uqpe_true & g2$uqpe_true <= g2$rif_hi)); put("rif_n_fine", nrow(g2))
  put("rif_cqr_n_covers_cond_fine", sum(g2$cqr_lo <= g2$cqr_true & g2$cqr_true <= g2$cqr_hi))
  put("rif_B", Q2$B); put("rif_dr_thresholds", Q2$dr_thresholds); put("rif_dr_n_extreme", Q2$dr_n_extreme)
  put_vec("rif_dr_q_fit", Q2$dr$q_fit, taus); put_vec("rif_dr_q_cf", Q2$dr$q_cf, taus)
  put_vec("rif_dr_uqe_discrete", Q2$dr$uqe_discrete, taus); put_vec("rif_q_sample", Q2$dr$q_sample, taus)
  for (k in names(Q2$sample)) put(paste0("rif_sample_", k), Q2$sample[[k]])

  ## ---- iv: block Q3 ---------------------------------------------------------------------------
  Q3 <- x$Q3
  put("iv_first_stage", Q3$first_stage); put("iv_B", Q3$B)
  put("iv_grid_min", Q3$ivqr$grid_range[1]); put("iv_grid_max", Q3$ivqr$grid_range[2])
  put("iv_ivqr_t50_tstat_at_truth", Q3$ivqr$tstat_at_truth)
  for (j in seq_along(taus)) {
    tg <- tid(taus[j])
    a <- Q3$qr_tab[j, ]; put_ci(paste0("iv_qr_", tg), a$est, a$se, a$lo, a$hi)
    b <- Q3$iv_tab[j, ]; put_ci(paste0("iv_ivqr_", tg), b$est, NULL, b$lo, b$hi)
    put(paste0("iv_ivqr_", tg, "_connected"), as.numeric(b$connected)); put(paste0("iv_ivqr_", tg, "_n_accepted"), b$n_accepted)
    # IV-QR (Chernozhukov-Hansen, rank similarity) targets the structural effect beta_D(tau) for the whole
    # population, so its coverage is stored against beta_D(tau) at every tau (_covers_struct). The LQTE is the
    # complier effect, a different object; the two coincide only at tau = 0.5, where the symmetric design makes
    # both 0.10. The IV-QR interval is therefore compared with the LQTE at the median only
    # (iv_ivqr_t50_covers_lqte); away from the median, coverage of the LQTE by IV-QR is not a meaningful statistic.
    put(paste0("iv_ivqr_", tg, "_covers_struct"), cov01(c(b$lo, b$hi), b$struct_true))
    if (abs(taus[j] - 0.5) < 1e-12) put("iv_ivqr_t50_covers_lqte", cov01(c(b$lo, b$hi), b$lqte_true))
    l <- Q3$lqte_tab[j, ]; put_ci(paste0("iv_lqte_", tg), l$est, l$se, l$lo, l$hi)
    put(paste0("iv_lqte_", tg, "_covers_truth"), cov01(c(l$lo, l$hi), l$truth))
  }
  put("iv_ivqr_t50_covers_lqte_median", cov01(Q3$ivqr$ci, tr$Q3$lqte_median))
  put("iv_qr_t50_covers_lqte_median", cov01(Q3$qr$ci, tr$Q3$lqte_median))
  put("iv_qr_t50_minus_truth", Q3$qr$est - tr$Q3$lqte_median)
  put_ci("iv_wald", Q3$wald$est, Q3$wald$se, Q3$wald$ci[1], Q3$wald$ci[2])
  put("iv_wald_covers_late", cov01(Q3$wald$ci, tr$Q3$late_log))
  for (k in names(Q3$sample)) put(paste0("iv_sample_", k), Q3$sample[[k]])

  ## ---- panel: block Q4 (Monte Carlo means at every tau; single samples) ----------------------------
  Q4 <- x$Q4
  put("panel_R", Q4$R); put("panel_tau", Q4$tau)
  for (e in c("pooled", "dummies", "canay")) for (Tn in c("T5", "T20")) for (j in seq_along(taus)) {
    pre <- paste0("panel_mc_", e, "_", Tn, "_", tid(taus[j]))
    put(paste0(pre, "_mean"), Q4$mc_mean[e, j, Tn]); put(paste0(pre, "_bias"), Q4$mc_bias[e, j, Tn])
    put(paste0(pre, "_sd"), Q4$mc_sd[e, j, Tn]); put(paste0(pre, "_rmse"), Q4$mc_rmse[e, j, Tn])
    put(paste0(pre, "_mcse"), Q4$mcse_mean[e, j, Tn])
    put(paste0("panel_single_", e, "_", Tn, "_", tid(taus[j])), Q4$single[[Tn]][e, j])
  }
  put("panel_canay_T20_abs_bias_t90", abs(Q4$canay_T20 - tr$Q4$beta_tau))
  put("panel_dummies_T5_abs_bias_t90", abs(Q4$dummies_T5 - tr$Q4$beta_tau))

  ## ---- cic: blocks Q5 and Q5f -----------------------------------------------------------------
  Q5 <- x$Q5
  put("cic_B", Q5$B)
  for (j in seq_along(taus)) {
    tg <- tid(taus[j])
    put_ci(paste0("cic_", tg), Q5$cic$est[j], Q5$cic$se[j], Q5$cic$ci[[j]][1], Q5$cic$ci[[j]][2])
    put(paste0("cic_", tg, "_covers_truth"), cov01(Q5$cic$ci[[j]], tr$Q5$qtt[j]))
    put(paste0("cic_", tg, "_cf"), Q5$cic$cf[j]); put(paste0("cic_", tg, "_q11"), Q5$cic$q11[j])
    put_ci(paste0("cic_qdid_", tg), Q5$qdid$est[j], Q5$qdid$se[j], Q5$qdid$ci[[j]][1], Q5$qdid$ci[[j]][2])
    put(paste0("cic_qdid_", tg, "_covers_truth"), cov01(Q5$qdid$ci[[j]], tr$Q5$qtt[j]))
    put(paste0("cic_lqte_grid_", tg, "_est"), Q5$lqte_grid$est[j])
  }
  g5 <- Q5$grid
  put("cic_n_covers_fine", sum(g5$cic_lo <= g5$qtt_true & g5$qtt_true <= g5$cic_hi)); put("cic_n_fine", nrow(g5))
  # per-level 0/1 flags on the figure grid (fig_13_cic) and the one CIC interval the text may quote as a pointwise
  # miss (tau = 0.20), with its true QTT
  for (j in seq_len(nrow(g5))) put(paste0("cic_fine_", tid(g5$tau[j]), "_covers_truth"),
                                   as.numeric(g5$cic_lo[j] <= g5$qtt_true[j] && g5$qtt_true[j] <= g5$cic_hi[j]))
  j20 <- which(abs(g5$tau - 0.2) < 1e-12)
  put_ci("cic_fine_t20", g5$cic_est[j20], g5$cic_se[j20], g5$cic_lo[j20], g5$cic_hi[j20]); put("dgp_Q5_fine_qtt_t20", g5$qtt_true[j20])
  put("cic_qdid_n_covers_fine", sum(g5$qdid_lo <= g5$qtt_true & g5$qtt_true <= g5$qdid_hi))
  L <- Q5$lqte
  put("cic_lqte_tau", L$tau); put_ci("cic_lqte", L$est, L$se, L$ci[1], L$ci[2]); put("cic_lqte_covers_truth", cov01(L$ci, tr$Q5$lqte))
  put("cic_lqte_num", L$num); put("cic_lqte_first_stage", L$first_stage); put("cic_lqte_B", L$B)
  put_ci("cic_wald", L$wald, L$se_wald, L$ci_wald[1], L$ci_wald[2]); put("cic_wald_covers_truth", cov01(L$ci_wald, tr$Q5$late_sw))
  for (k in names(Q5$sample$means)) put(paste0("cic_sample_mean_", k), Q5$sample$means[[k]])
  for (k in names(Q5$sample$medians)) put(paste0("cic_sample_median_", k), Q5$sample$medians[[k]])
  for (k in names(Q5$sample$d_rate)) put(paste0("cic_sample_f_rate_", k), Q5$sample$d_rate[[k]])

  ## ---- cross: chapter-20 cross-check toys -------------------------------------------------------
  for (k in names(x$cross$num13)) put(paste0("cross_num13_", k), x$cross$num13[[k]])
  put("cross_n_toys", length(x$cross$toys))
  nums
}
