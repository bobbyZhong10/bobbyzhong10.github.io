# ------------------------------------------------------------------
# numbers_def.R -- the number registry of chapter 21.
#
# collect_numbers21(dir = ".") -> named list of finite numeric scalars, read from the
# res_*.rds files (02 to 06), truth21.rds and the block metadata in aster_*.rds, plus four
# constants of the legacy R1 block read from the 01_dgp.R lines that use them (aster_R1.rds
# keeps the legacy meta, which does not store them; the stored scores are checked against
# the score constants). Nothing is estimated here, no random numbers are drawn and no
# design value is typed in: the list is a pure function of those files. A field that is
# not a finite number stops the collection with its key. 08_summary.R writes the list to
# numbers.json; tests/test_numbers_fresh.R recomputes it and compares. Task 8 may add keys
# (ruling R13), never rename them.
#
# Key scheme: <block>_<object>_<field>, block in r1, r1c, r1h, r2, r3, r4, r4g.
#   *_est, *_se          conventional point estimate and standard error
#   *_est_bc, *_se_rob   bias-corrected estimate and robust standard error (rdrobust)
#   *_ci_lo, *_ci_hi     the object's headline 95% interval: robust bias-corrected for the
#                        rdrobust, rdmc and rdms objects; the exceptions say what they are:
#                        r2_conv (conventional), r2_crv (conventional, clustered by support
#                        point), r2_honest_* (honest, Armstrong and Kolesar), r2_locrand
#                        (randomization test inversion), r4_extrap* (conventional local
#                        linear at the fixed h of the key, stacked sandwich se), r2_mc_*lr*
#                        (Neyman interval in the fixed window); *_ci_conv_lo / _hi is the
#                        conventional interval where both exist
#   *_h, *_b             main and pilot bandwidth; *_n_eff_left / _right effective n
#   *_p, *_p_rob         p-values (stored as computed, also when tiny or 0)
#   *_truth, truth_*     simulation truth;  dgp_*  design parameters and sample facts
#   r2_mc_*              Monte Carlo (coverage of the truth, mean length, ...)
#   *_gap_*, *_ratio_*   differences and ratios the text compares
# Integers (sample sizes, seeds) are registered too; the number checker ignores them.
# ------------------------------------------------------------------

collect_numbers21 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop("collect_numbers21: ", f, " is missing (run 01 to 06 first)", call. = FALSE)
    readRDS(p)
  }
  tr <- rd("truth21.rds")
  core <- rd("res_core.rds"); disc <- rd("res_discrete.rds"); kink <- rd("res_kink.rds")
  multi <- rd("res_multi.rds"); extra <- rd("res_extra.rds")
  R1 <- rd("aster_R1.rds"); m_R1 <- R1$meta
  m_R1c <- rd("aster_R1c.rds")$meta; m_R1h <- rd("aster_R1h.rds")$meta
  m_R2 <- rd("aster_R2.rds")$meta; m_R3 <- rd("aster_R3.rds")$meta; m_R4 <- rd("aster_R4.rds")$meta
  m_R4g <- rd("aster_R4g.rds")$meta

  nums <- list()
  put <- function(name, value) {
    if (!is.character(name) || length(name) != 1L || !nzchar(name)) stop("collect_numbers21: bad key", call. = FALSE)
    if (name %in% names(nums)) stop("collect_numbers21: duplicate key ", name, call. = FALSE)
    if (is.logical(value)) value <- as.numeric(value)
    if (!(is.numeric(value) && length(value) == 1L && is.finite(value)))
      stop("collect_numbers21: ", name, " is not a finite numeric scalar", call. = FALSE)
    nums[[name]] <<- as.numeric(value)
    invisible(NULL)
  }
  put_ci <- function(name, ci) { put(paste0(name, "_lo"), ci[[1]]); put(paste0(name, "_hi"), ci[[2]]) }
  # every numeric field of a result list: scalars as is; ci* and window as _lo / _hi;
  # n_eff and n_support_in_h as _left / _right; nested lists recursively; skip the rest
  put_obj <- function(prefix, obj, skip = character(0)) {
    for (nm in setdiff(names(obj), skip)) {
      v <- obj[[nm]]; key <- paste0(prefix, "_", nm)
      if (is.list(v) && !is.data.frame(v)) {
        put_obj(key, v)
      } else if ((is.numeric(v) || is.logical(v)) && length(v) == 1L) {
        put(key, v)                                   # put() stops on NA, NaN, Inf
      } else if (is.numeric(v) && length(v) == 2L && (grepl("(^|_)ci(_|$)", nm) || nm == "window")) {
        put_ci(key, v)
      } else if (is.numeric(v) && length(v) == 2L && nm %in% c("n_eff", "n_support_in_h")) {
        put(paste0(key, "_left"), v[[1]]); put(paste0(key, "_right"), v[[2]])
      } else if (is.numeric(v) && length(v) == 2L) {
        # a pair this function does not know how to name: stop instead of dropping it silently
        # (register it by hand and pass its name in `skip`)
        stop("collect_numbers21: ", key, " is a length-2 numeric field that put_obj cannot name; ",
             "register it explicitly and skip it here", call. = FALSE)
      }
    }
  }
  id <- function(v) gsub("-", "m", gsub("\\.", "p", format(v, scientific = FALSE, drop0trailing = TRUE)))  # 4.25 -> 4p25

  ## ==================================================================
  ## R1 (inherited block): design, truth, sample facts, res_core.rds
  ## ==================================================================
  put("dgp_r1_seed", tr$seeds[["R1"]]); put("dgp_r1_n", m_R1$N); put("dgp_r1_cutoff", m_R1$c_cut)
  put("dgp_r1_sigma_e", m_R1$sigma_e)
  # constants of the legacy R1 block that its meta does not store (aster_R1.rds must stay the
  # legacy file): the score distribution and the fuzzy badge probabilities, read from the
  # 01_dgp.R lines that draw with them (test_reproduce.R and test_dgp.R lock that stream)
  r1_code <- function(pattern) {
    src <- readLines(file.path(dir, "01_dgp.R"), warn = FALSE)
    hit <- Filter(length, regmatches(src, regexec(pattern, src)))
    if (length(hit) != 1L) stop("collect_numbers21: 01_dgp.R has ", length(hit), " lines matching ", pattern, call. = FALSE)
    as.numeric(hit[[1]][-1])
  }
  sc <- r1_code("^raw +<- rnorm\\(N, mean = ([0-9.]+), sd = ([0-9.]+)\\)")
  pb <- r1_code("^p_badge <- ifelse\\(score >= c_cut, ([0-9.]+), ([0-9.]+)\\)")     # above, below
  stopifnot(abs(m_R1$first_stage_jump_true - (pb[1] - pb[2])) < 1e-12,
            abs(m_R1$frac_badge_above - pb[1]) < 0.02, abs(m_R1$frac_badge_below - pb[2]) < 0.02)
  # data-level guard for the score constants: the stored scores are the draw clamped to [45, 100],
  # so their mean is within 0.1 and their sd within 0.25 of the constants (sampling error and
  # clamping at both ends); a wrong source line would put them off by far more
  stopifnot(abs(mean(R1$dat$score) - sc[1]) < 0.3, abs(stats::sd(R1$dat$score) - sc[2]) < 0.5)
  put("dgp_r1_score_mean", sc[1]); put("dgp_r1_score_sd", sc[2])
  put("dgp_r1_p_badge_below", pb[2]); put("dgp_r1_p_badge_above", pb[1])
  put("r1_truth_jump", tr$R1$jump); put("r1_truth_first_stage", tr$R1$first_stage_jump)
  put_obj("r1_sample", core$sample)
  put("r1_naive", core$naive)
  put("r1_ratio_naive_to_truth", core$naive / tr$R1$jump)
  put("r1_gap_naive_minus_truth", core$naive - tr$R1$jump)
  put("r1_share_naive_selection_bias", (core$naive - tr$R1$jump) / core$naive)
  put_obj("r1_sharp", core$sharp)
  put("r1_bw_half_est", core$bw_sens[[1]]); put("r1_bw_dbl_est", core$bw_sens[[2]])
  put_obj("r1_bw_half", core$bw_sens_detail$half, skip = "est"); put_obj("r1_bw_dbl", core$bw_sens_detail$dbl, skip = "est")
  put_obj("r1_cer", core$cer); put_obj("r1_p2", core$p2)
  put("r1_ratio_cer_h_to_mse_h", core$cer$h / core$sharp$h)
  put("r1_poly4_est", core$poly[[1]]); put("r1_poly6_est", core$poly[[2]]); put("r1_poly8_est", core$poly[[3]])
  put_obj("r1_density", core$density)
  put_obj("r1_fuzzy", core$fuzzy)
  put("r1_gap_fuzzy_minus_sharp", core$fuzzy$est - core$sharp$est)
  put("r1_ratio_fuzzy_h_to_sharp_h", core$fuzzy$h / core$sharp$h)
  ## power (06_covs_donut_power.R, on the R1 design)
  put_obj("r1_power", extra$power)
  put("r1_power_pct", 100 * extra$power$power)

  ## ==================================================================
  ## R1c (covariates)
  ## ==================================================================
  put("dgp_r1c_seed", tr$seeds[["R1c"]]); put("dgp_r1c_n", m_R1c$N); put("r1c_truth_jump", tr$R1c$jump)
  for (nm in c("sd_a", "sd_prev", "sd_eps", "tenure_coef")) put(paste0("dgp_r1c_", nm), m_R1c[[nm]])
  put_obj("r1c_nocovs", extra$nocovs)
  put_obj("r1c_covs", extra$covs, skip = "coef_covs")
  stopifnot(setequal(names(extra$covs$coef_covs), c("gmv_prev", "tenure")))   # names attached by matching in 06
  put("r1c_covs_coef_gmv_prev", extra$covs$coef_covs[["gmv_prev"]]); put("r1c_covs_coef_tenure", extra$covs$coef_covs[["tenure"]])
  put_obj("r1c_covs_same_h", extra$covs_same_h)
  put_obj("r1c_ratio", extra$covs_ratio)
  put("r1c_ci_len_reduction_pct", 100 * (1 - extra$covs_ratio$len_rob))
  put("r1c_ci_len_reduction_same_h_pct", 100 * (1 - extra$covs_ratio$len_rob_same_h))
  put("r1c_covs_partial_r2", extra$covs_partial_r2)
  for (nm in names(extra$balance)) put_obj(paste0("r1c_balance_", nm), extra$balance[[nm]])

  ## ==================================================================
  ## R1h (heaping, donut)
  ## ==================================================================
  put("dgp_r1h_seed", tr$seeds[["R1h"]]); put("dgp_r1h_n", m_R1h$N); put("r1h_truth_jump", tr$R1h$jump)
  put("dgp_r1h_push_max", m_R1h$push_max); put("dgp_r1h_push_width", m_R1h$push_width)
  put("dgp_r1h_share_heaped", m_R1h$share_heaped); put("dgp_r1h_share_heaped_pct", 100 * m_R1h$share_heaped)
  put_obj("r1h_heap", extra$heap, skip = c("push_max", "push_width", "share_heaped", "N"))
  put("r1h_ratio_at_80_to_at_79", extra$heap$n_at_80 / extra$heap$n_at_79)
  put("r1h_gap_heaped_minus_true80_mean_y", extra$heap$mean_y_heaped - extra$heap$mean_y_true_80)
  for (nm in c("heap_naive", "donut", "donut_heap")) {
    put_obj(paste0("r1h_", nm), extra[[nm]])
    put(paste0("r1h_", nm, "_abs_error"), abs(extra[[nm]]$est - tr$R1h$jump))
  }

  ## ==================================================================
  ## R2 (discrete running variable, honest inference, local randomization)
  ## ==================================================================
  P <- m_R2$par
  put("dgp_r2_seed", tr$seeds[["R2"]]); put("dgp_r2_n", m_R2$N); put("dgp_r2_cutoff", m_R2$cutoff)
  put("r2_truth_jump", tr$R2$jump); put("r2_truth_m_bound", tr$R2$M_bound)
  for (nm in c("K", "load", "shift", "sd_e", "xmin", "base", "beta_q", "tau0", "page_gain", "page_mid", "page_scale",
               "sd_eps")) put(paste0("dgp_r2_", nm), P[[nm]])
  for (k in seq_along(P$th)) put(sprintf("dgp_r2_threshold_%d_%d", k, k + 1), P$th[k])
  cef <- m_R2$cef
  m_part <- max(abs(diff(cef$m0, differences = 2))) / (1 / P$K)^2
  put("dgp_r2_m_bound_m0_part", m_part); put("dgp_r2_m_bound_tau_part", m_R2$M_bound - m_part)
  # evaluation points: the two support points of the local randomization window, the next
  # two support points right of the cutoff, and the top of the scale
  wl <- disc$locrand$window[1]; wr <- disc$locrand$window[2]
  x_eval <- unique(c(wl, cef$x[cef$x >= m_R2$cutoff][1:3], max(cef$x)))
  stopifnot(any(abs(cef$x - wr) < 1e-9), length(x_eval) == 5L)
  for (xx in x_eval) {
    i <- which(abs(cef$x - xx) < 1e-9); stopifnot(length(i) == 1L)
    put(sprintf("dgp_r2_m0_at_%s", id(xx)), cef$m0[i]); put(sprintf("dgp_r2_tau_at_%s", id(xx)), cef$tau[i])
    put(sprintf("dgp_r2_mass_at_%s", id(xx)), cef$mass[i])
  }
  il <- which(abs(cef$x - wl) < 1e-9); ir <- which(abs(cef$x - wr) < 1e-9)
  stopifnot(length(il) == 1L, length(ir) == 1L, identical(disc$mc$lr_window, disc$locrand$window))
  lr_target <- cef$m0[ir] + cef$tau[ir] - cef$m0[il]
  put("dgp_r2_lr_window_target", lr_target)                     # E[y | wr] - E[y | wl]
  put("dgp_r2_lr_window_bias", lr_target - tr$R2$jump)
  put("dgp_r2_lr_window_m0_gap", cef$m0[ir] - cef$m0[il])            # Y(0) part of the window bias
  put(sprintf("dgp_r2_tau_at_%s_minus_jump", id(wr)), cef$tau[ir] - tr$R2$jump)
  put_obj("r2_sample", disc$sample)
  put_obj("r2_conv", disc$conv); put_obj("r2_robust", disc$robust); put_obj("r2_crv", disc$crv)
  put("r2_ratio_crv_se_to_conv_se", disc$crv$se / disc$conv$se)
  put_obj("r2_honest_mbound", disc$honest); put_obj("r2_honest_mlb", disc$honest_lb)
  put_obj("r2_honest_mest", disc$honest_mest)
  put_obj("r2_smooth_bound", disc$smooth_bound)
  put_obj("r2_locrand", disc$locrand)
  put("r2_winselect_window_lo", disc$winselect$window[1]); put("r2_winselect_window_hi", disc$winselect$window[2])
  for (k in seq_len(min(3L, length(disc$winselect$balance_p)))) {
    put(sprintf("r2_winselect_w%d_right", k), disc$winselect$w_right[k])
    put(sprintf("r2_winselect_w%d_balance_p", k), disc$winselect$balance_p[k])
  }
  S <- disc$sensitivity
  for (k in seq_along(S$w_right)) {
    put(sprintf("r2_sens_w%d_lo", k), S$w_left[k]); put(sprintf("r2_sens_w%d_hi", k), S$w_right[k])
    put(sprintf("r2_sens_w%d_est", k), S$est[k])
    for (j in seq_along(S$tlist)) put(sprintf("r2_sens_w%d_p_tau_%s", k, id(S$tlist[j])), S$p[j, k])
  }
  mc <- disc$mc
  put_obj("r2_mc", mc[setdiff(names(mc), c("draws", "lr_window", "cores", "runtime_sec"))])
  put("r2_mc_lr_window_lo", mc$lr_window[1]); put("r2_mc_lr_window_hi", mc$lr_window[2])
  for (nm in grep("^cover_", names(mc), value = TRUE)) put(paste0("r2_mc_", nm, "_pct"), 100 * mc[[nm]])
  put("r2_mc_bias_conv", mc$mean_est_conv - tr$R2$jump)
  put("r2_mc_ratio_len_honest_to_conv", mc$len_honest / mc$len_conv)
  put("r2_ratio_len_honest_mbound_to_conv", diff(disc$honest$ci) / diff(disc$conv$ci))

  ## ==================================================================
  ## R3 (regression kink)
  ## ==================================================================
  put("dgp_r3_seed", tr$seeds[["R3"]]); put("dgp_r3_n", m_R3$N); put("dgp_r3_cutoff", m_R3$cutoff)
  put("dgp_r3_floor_from", m_R3$floor_from); put("dgp_r3_dslope_b", m_R3$slope_right - m_R3$slope_left)
  put("dgp_r3_dslope_b_per_100", 100 * (m_R3$slope_right - m_R3$slope_left))
  put("r3_truth_kink_effect", tr$R3$kink_effect); put("r3_truth_level_jump", tr$R3$level_jump)
  put("r3_truth_dslope_y", tr$R3$kink_effect * tr$R3$dslope_b)
  # design constants of gen_R3, from the block's meta: elasticity of GMV to last quarter's
  # GMV, intercept, residual sd, lognormal median and sd of log x, sample range
  put("dgp_r3_elasticity_x", m_R3$elasticity_x); put("dgp_r3_intercept", m_R3$intercept)
  put("dgp_r3_sd_eps", m_R3$sd_eps); put("dgp_r3_x_median", m_R3$x_median); put("dgp_r3_sdlog", m_R3$sdlog)
  put("dgp_r3_x_lo", m_R3$x_lo); put("dgp_r3_x_hi", m_R3$x_hi)
  # the schedule: base rate up to the kink, floor from floor_from on, slope change between
  bf <- m_R3$b_fun
  stopifnot(isTRUE(all.equal(bf(c(m_R3$x_lo, m_R3$cutoff, m_R3$floor_from, m_R3$x_hi)),
                             c(m_R3$rate_base, m_R3$rate_base, m_R3$rate_floor, m_R3$rate_floor))),
            isTRUE(all.equal(m_R3$rate_floor - m_R3$rate_base, m_R3$slope_right * (m_R3$floor_from - m_R3$cutoff))))
  put("dgp_r3_rate_base", m_R3$rate_base); put("dgp_r3_rate_floor", m_R3$rate_floor)
  put_obj("r3_sample", kink$sample)
  put_obj("r3_rk", kink$rk); put_obj("r3_rk_p1", kink$rk_p1); put_obj("r3_rk_hand", kink$rk_hand)
  put_obj("r3_rk_fuzzy", kink$rk_fuzzy); put_obj("r3_level", kink$level_jump)
  put("r3_gap_rk_p1_minus_truth", kink$rk_p1$est - tr$R3$kink_effect)
  put("r3_gap_rk_minus_truth", kink$rk$est - tr$R3$kink_effect)
  put("r3_rk_pct_per_pp", 100 * kink$rk$est)

  ## ==================================================================
  ## R4 (multi-cutoff, extrapolation)
  ## ==================================================================
  G <- m_R4$categories
  put("dgp_r4_seed", tr$seeds[["R4"]]); put("dgp_r4_n", m_R4$N)
  for (i in seq_len(nrow(G))) {
    pre <- sprintf("dgp_r4_c%d", G$cutoff[i])
    put(paste0(pre, "_share"), G$share[i]); put(paste0(pre, "_mu"), G$mu[i]); put(paste0(pre, "_sd"), G$sd[i])
    put(paste0(pre, "_alpha"), G$alpha[i])
  }
  # the effect line tau(x) = a + s (x - c_1) of gen_R4, from the true cutoff effects: a is the
  # effect at the lowest cutoff, s the slope from the first to the second cutoff; the third
  # cutoff and the extrapolation point must lie on the same line
  cuts4 <- sort(m_R4$cutoffs); tc <- tr$R4$tau_c[as.character(cuts4)]
  a4 <- tc[[1]]; s4 <- (tc[[2]] - tc[[1]]) / (cuts4[2] - cuts4[1])
  stopifnot(abs(tc[[3]] - (a4 + s4 * (cuts4[3] - cuts4[1]))) < 1e-12,
            abs(tr$R4$tau_extrap - (a4 + s4 * (tr$R4$x_extrap - cuts4[1]))) < 1e-12)
  put("dgp_r4_tau_intercept", a4); put("dgp_r4_tau_slope", s4)
  put("dgp_r4_alpha_gap_vintage_home", G$alpha[G$cutoff == m_R4$extrap_cutoff] - G$alpha[G$cutoff == m_R4$ref_cutoff])
  for (cc in names(tr$R4$tau_c)) {
    put(sprintf("r4_truth_tau_c%s", cc), tr$R4$tau_c[[cc]]); put(sprintf("r4_truth_omega_c%s", cc), tr$R4$omega[[cc]])
  }
  put("r4_truth_tau_extrap", tr$R4$tau_extrap); put("r4_truth_x_extrap", tr$R4$x_extrap)
  put("r4_truth_tau_pooled", tr$R4$tau_pooled); put("r4_truth_simple_mean", mean(tr$R4$tau_c))
  put("r4_truth_gap_pooled_minus_mean", tr$R4$tau_pooled - mean(tr$R4$tau_c))
  bc <- multi$by_cutoff
  for (i in seq_len(nrow(bc))) {
    pre <- sprintf("r4_c%d", bc$cutoff[i])
    for (nm in c("est", "est_bc", "se", "se_rob", "h", "n_left", "n_right", "n_sample")) put(paste0(pre, "_", nm), bc[[nm]][i])
    put_ci(paste0(pre, "_ci"), bc$ci[[i]])
    put(paste0(pre, "_weight"), multi$weights[[as.character(bc$cutoff[i])]])
    put(paste0(pre, "_est_common_h"), multi$by_cutoff_common$est[[as.character(bc$cutoff[i])]])
    put(paste0(pre, "_rdmc_nh_weight"), multi$rdmc_weighted$weights[[as.character(bc$cutoff[i])]])
  }
  put_obj("r4_pooled", multi$pooled)
  put_obj("r4_pooled_norm", multi$pooled_norm)
  put("r4_simple_mean", multi$simple_mean)
  put("r4_gap_pooled_minus_simple_mean", multi$pooled$est - multi$simple_mean)
  put("r4_gap_pooled_minus_truth_mean", multi$pooled$est - mean(tr$R4$tau_c))
  put("r4_gap_pooled_norm_minus_target", multi$pooled_norm$est - multi$pooled_norm$target)
  put("r4_rdmc_weighted_est", multi$rdmc_weighted$est); put_ci("r4_rdmc_weighted_ci", multi$rdmc_weighted$ci)
  put("r4_common_h_pooled", multi$by_cutoff_common$pooled)
  ex <- multi$extrap
  put_obj("r4_extrap", ex[c("est", "se", "ci", "se_indep", "h", "naive_cross", "bias_at_70", "truth")])
  for (nm in names(ex$parts)) { put(paste0("r4_extrap_part_", nm), ex$parts[[nm]]); put(paste0("r4_extrap_n_", nm), ex$n[[nm]]) }
  for (nm in names(ex$sens)) put_obj(paste0("r4_extrap_", nm), ex$sens[[nm]])
  put("r4_extrap_gap_naive_cross_minus_truth", ex$naive_cross - ex$truth)
  put_obj("r4_sample", list(N = multi$sample$N4))

  ## ==================================================================
  ## R4g (geographic RD)
  ## ==================================================================
  put("dgp_r4g_seed", tr$seeds[["R4g"]]); put("dgp_r4g_n", m_R4g$N); put("dgp_r4g_radius", m_R4g$radius)
  # design constants of gen_R4g, from the block's meta
  put("dgp_r4g_share_urban", m_R4g$share_urban); put("dgp_r4g_urban_scale_km", m_R4g$urban_scale)
  put("dgp_r4g_suburb_lo_km", m_R4g$suburb_lo); put("dgp_r4g_suburb_hi_km", m_R4g$suburb_hi)
  put("dgp_r4g_tau_center", m_R4g$tau_center); put("dgp_r4g_tau_east_slope", m_R4g$tau_east)
  put("dgp_r4g_base", m_R4g$base); put("dgp_r4g_demand_amp", m_R4g$demand_amp)
  put("dgp_r4g_demand_scale_km", m_R4g$demand_scale); put("dgp_r4g_north_slope", m_R4g$north)
  put("dgp_r4g_sd_eps", m_R4g$sd_eps)
  stopifnot(isTRUE(all.equal(tr$R4g$tau_b, m_R4g$tau_center + m_R4g$tau_east * m_R4g$boundary$b1 / m_R4g$radius)))
  put("r4g_n_inside", multi$sample$n_inside); put("r4g_n_outside", multi$sample$n_outside)
  geo <- multi$geo
  for (k in seq_len(nrow(geo))) {
    pre <- sprintf("r4g_p%d", geo$point[k])
    for (nm in c("angle_deg", "b1", "b2", "est", "est_bc", "se_rob", "h", "n_left", "n_right", "truth", "covers"))
      put(paste0(pre, "_", nm), geo[[nm]][k])
    put_ci(paste0(pre, "_ci"), geo$ci[[k]])
    put(paste0(pre, "_n_within4_inside"), multi$sample$n_within4[k, "inside"])
    put(paste0(pre, "_n_within4_outside"), multi$sample$n_within4[k, "outside"])
  }
  put("r4g_share_covered", mean(geo$covers)); put("r4g_n_covered", sum(geo$covers))
  put_obj("r4g_pooled", multi$geo_pooled)
  put("r4g_naive", multi$geo_naive)
  put("r4g_min_n_within4", min(multi$sample$n_within4))

  nums
}
