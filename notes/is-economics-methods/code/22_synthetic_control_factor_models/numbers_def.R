# ------------------------------------------------------------------
# numbers_def.R -- the number registry of chapter 22.
#
# collect_numbers22(dir = ".") -> named list of finite numeric scalars, read from the
# res_*.rds files (02 to 06, the screen cache through res_stag.rds$screen), truth22.rds and
# the block metadata in aster_S*.rds, plus one constant of the legacy S1 block and one of S2
# (their effect plateaus) read from the 01_dgp.R lines that use them (aster_S1.rds keeps the
# legacy meta, which does not store them). Nothing is estimated here, no random numbers
# are drawn and no design value is typed in: the list is a pure function of those files. A
# field that is not a finite number stops the collection with its key. 08_summary.R writes
# it to numbers.json; tests/test_numbers_fresh.R recomputes it and compares. Task 11 may add
# keys, never rename them.
#
# Key scheme: <block>_<object>_<field>, block in s1, s2, s3, s4; truth_<block>_* is the
# simulation truth, dgp_<block>_* a design parameter or sample fact.
#   s1_core_*            the inherited S1 analysis (02_core.R, tidysynth, legacy conventions)
#   s1_*, s2_*           03_s1_modern.R and 04_hull.R; *_hand_* = sc_fns.R by hand
#   s3_<estimator>_*     05_staggered.R on the stored draw (seed seeds22[["S3"]])
#   s3_screen_*          05b_stag_screen.R, the 24-seed screen (Ruling S6); s3_screen_stored_*
#                        describes the stored draw (truth22$seeds[["S3"]]) within the screen:
#                        *_rank_est_<k> its estimate's rank among the draws (1 = lowest),
#                        *_rank_abs_err_<k> (1 = closest to the truth), *_rank_share_twfe_bias_
#                        removed_<k> (1 = largest share); s3_screen_n_est_above_truth_<k> counts
#                        draws whose estimate exceeds the truth
#   s3_scpi_*            05c_stag_scpi.R, scpi on the stored draw (Ruling S12); *_e00 ... *_e11
#                        one value per event month
#   s4_*                 06_design.R
#   *_est, *_se          point estimate and standard error
#   *_ci_lo, *_ci_hi     95% interval (for scpi: the prediction interval, coverage at least 90%)
#   *_err                estimate minus the truth;  *_covers  1 if the interval covers the truth
#   *_w_<unit>           a weight, unit names lower-cased without the hyphen (donor10)
#   *_n_w_above_0p01     the number of weights above 0.01 ("p" for the decimal point)
#   *_m25 ... *_m36      one value per post month
# Integers (sample sizes, seeds, counts) are registered too; the number checker ignores them.
# ------------------------------------------------------------------

collect_numbers22 <- function(dir = ".") {
  rd <- function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) stop("collect_numbers22: ", f, " is missing (run 01 to 06 first)", call. = FALSE)
    readRDS(p)
  }
  tr <- rd("truth22.rds")
  core <- rd("res_sc_core.rds"); s1 <- rd("res_s1_ext.rds"); s2 <- rd("res_hull.rds")
  s3 <- rd("res_stag.rds"); s4 <- rd("res_design.rds")
  A1 <- rd("aster_S1.rds"); A2 <- rd("aster_S2.rds"); A3 <- rd("aster_S3.rds"); A4 <- rd("aster_S4.rds")

  nums <- list()
  put <- function(name, value) {
    if (!is.character(name) || length(name) != 1L || !nzchar(name)) stop("collect_numbers22: bad key", call. = FALSE)
    if (name %in% names(nums)) stop("collect_numbers22: duplicate key ", name, call. = FALSE)
    if (is.logical(value)) value <- as.numeric(value)
    value <- unname(value)
    if (!(is.numeric(value) && length(value) == 1L && is.finite(value)))
      stop("collect_numbers22: ", name, " is not a finite numeric scalar", call. = FALSE)
    nums[[name]] <<- as.numeric(value)
    invisible(NULL)
  }
  put_ci <- function(name, ci) { put(paste0(name, "_ci_lo"), ci[[1]]); put(paste0(name, "_ci_hi"), ci[[2]]) }
  # est, se (when present), ci, err and covers of an estimator object list(est, ci, se)
  put_est <- function(name, o, truth = NULL) {
    put(paste0(name, "_est"), o$est)
    if (!is.null(o$se)) put(paste0(name, "_se"), o$se)
    if (!is.null(o$ci)) put_ci(name, o$ci)
    if (!is.null(truth)) {
      put(paste0(name, "_err"), o$est - truth)
      if (!is.null(o$ci)) put(paste0(name, "_covers"), o$ci[[1]] <= truth && truth <= o$ci[[2]])
    }
  }
  uid <- function(u) tolower(gsub("-", "", u))                    # Donor-10 -> donor10
  put_w <- function(prefix, w) {
    if (is.null(names(w)) || anyNA(names(w)) || any(!nzchar(names(w))) || anyDuplicated(names(w)))
      stop("collect_numbers22: the weights for ", prefix, " are not named by unit", call. = FALSE)
    for (u in names(w)) put(paste0(prefix, "_w_", uid(u)), w[[u]])
  }
  put_months <- function(prefix, v, months) {
    stopifnot(length(v) == length(months))
    for (k in seq_along(months)) put(sprintf("%s_m%02d", prefix, months[k]), v[[k]])
  }
  # a numeric literal from the 01_dgp.R line matching `pattern` (constants the legacy meta
  # does not store); the line must exist exactly once
  dgp_lines <- readLines(file.path(dir, "01_dgp.R"))
  lit <- function(pattern, regex) {
    hit <- grep(pattern, dgp_lines, value = TRUE)
    if (length(hit) != 1L) stop("collect_numbers22: expected one 01_dgp.R line matching ", pattern, call. = FALSE)
    as.numeric(sub(regex, "\\1", hit))
  }

  ## ==================================================================
  ## S1: design and truth
  ## ==================================================================
  m1 <- A1$meta
  put("dgp_s1_seed", tr$seeds[["S1"]]); put("dgp_s1_n_donors", m1$n_donors); put("dgp_s1_n_months", m1$Tt)
  put("dgp_s1_first_reform_month", m1$T0); put("dgp_s1_n_pre", m1$T0 - 1); put("dgp_s1_n_post", m1$Tt - m1$T0 + 1)
  put("dgp_s1_sigma_e", m1$sigma_e)
  for (k in seq_along(m1$anchors)) {
    put(sprintf("dgp_s1_anchor_%d", k), m1$anchors[k]); put(sprintf("dgp_s1_w_true_%d", k), m1$w_true[k])
  }
  put("dgp_s1_tau_plateau", lit("^tau_path\\[post\\] <- [0-9.]+ \\* \\(1 - exp", "^tau_path\\[post\\] <- ([0-9.]+) .*$"))
  put("truth_s1_att_post", tr$S1$att_post); put("truth_s1_tau_last", tr$S1$tau_path[12])
  put_months("truth_s1_tau", tr$S1$tau_path, 25:36)

  ## ---- S1 inherited analysis (02_core.R) ----------------------------
  put("s1_core_naive_did", core$naive_did); put("s1_core_did_bias", core$did_bias); put("s1_core_did_ratio", core$did_ratio)
  put("s1_core_tr_pre_slope", core$tr_pre_slope); put("s1_core_do_pre_slope", core$do_pre_slope)
  put("s1_core_tr_minus_do_slope", core$tr_pre_slope - core$do_pre_slope)
  put("s1_core_att_post", core$att_post); put("s1_core_err", core$att_post - tr$S1$att_post)
  put("s1_core_gap_last", core$gap_last); put("s1_core_pre_rmspe", core$pre_rmspe)
  put("s1_core_mspe_ratio", core$mspe_ratio); put("s1_core_p_perm", core$p_perm); put("s1_core_rank_perm", core$rank_perm)
  put("s1_core_n_perm", core$n_perm); put("s1_core_intime_gap", core$intime_gap); put("s1_core_t_fake", core$T_fake)
  put("s1_core_n_fig_placebos_kept", core$n_fig_placebos_kept)
  put_w("s1_core", core$w)
  put("s1_core_n_w_above_0p01", sum(core$w > 0.01)); put("s1_core_w_anchor_sum", sum(core$w[sprintf("Donor-%02d", m1$anchors)]))

  ## ---- S1 modern (03_s1_modern.R) -----------------------------------
  put("s1_naive_did", s1$naive_did); put("s1_did_uniform", s1$did_uniform); put("s1_twfe", s1$twfe)
  put("s1_did_identity_gap", s1$did_identity_gap); put("s1_twfe_identity_gap", s1$twfe_identity_gap)
  put_est("s1_did_placebo", s1$did_placebo, tr$S1$att_post)
  put_est("s1_sdid", s1$sdid, tr$S1$att_post)
  put("s1_sdid_n_omega_nonzero", s1$sdid$n_omega_nonzero); put("s1_sdid_n_lambda_nonzero", s1$sdid$n_lambda_nonzero)
  put("s1_sdid_placebo_reps", s1$sdid$placebo_reps)
  # the pre-month carrying SDID's largest time weight (on S1 a single month carries all of it)
  put("s1_sdid_lambda_month", which.max(s1$sdid$lambda)); put("s1_sdid_lambda_max", max(s1$sdid$lambda))
  put_est("s1_sc_synthdid", s1$sc_synthdid, tr$S1$att_post)
  put("s1_asc_est", s1$asc$est); put("s1_asc_err", s1$asc$est - tr$S1$att_post); put("s1_asc_p_conformal", s1$asc$p_conformal)
  put("s1_asc_n_neg", s1$asc$n_neg); put("s1_asc_min_w", s1$asc$min_w); put("s1_asc_ridge_lambda", s1$asc$ridge_lambda)
  put("s1_sc_hand_est", s1$sc_hand$est); put("s1_sc_hand_err", s1$sc_hand$est - tr$S1$att_post)
  put("s1_sc_hand_pre_rmspe", s1$sc_hand$pre_rmspe); put("s1_sc_hand_gap_last", s1$sc_hand$gap_last)
  put("s1_sc_hand_n_nonzero", s1$sc_hand$n_nonzero); put_w("s1_sc_hand", s1$sc_hand$w)
  put("s1_sc_pre_rmspe", s1$sc_pre_rmspe)
  put("s1_conformal_p", s1$conformal$p); put("s1_conformal_n_perm", s1$conformal$n_perm); put("s1_conformal_q", s1$conformal$q)
  put("s1_conformal_p_moving_block", s1$conformal$p_moving_block)
  put("s1_conformal_p_moving_block_q2", s1$conformal$p_moving_block_q2); put("s1_conformal_n_shifts", s1$conformal$n_shifts)
  put_months("s1_scpi_est", s1$scpi$est, 25:36); put_months("s1_scpi_ci_lo", s1$scpi$lo, 25:36)
  put_months("s1_scpi_ci_hi", s1$scpi$hi, 25:36)
  put("s1_scpi_n_covered", sum(s1$scpi$covered)); put("s1_scpi_n_covered_in_sample", sum(s1$scpi$covered_in))
  put("s1_scpi_mean_length", mean(s1$scpi$hi - s1$scpi$lo)); put("s1_scpi_sims", s1$scpi$sims)
  put("s1_scpi_u_alpha", s1$scpi$u_alpha); put("s1_scpi_e_alpha", s1$scpi$e_alpha)
  put("s1_perm_hand_ratio", s1$perm_hand$ratio); put("s1_perm_hand_p", s1$perm_hand$p)
  put("s1_perm_hand_rank", s1$perm_hand$rank); put("s1_perm_hand_n_units", s1$perm_hand$n_units)

  ## ==================================================================
  ## S2: design, truth, 04_hull.R
  ## ==================================================================
  m2 <- A2$meta
  put("dgp_s2_seed", tr$seeds[["S2"]]); put("dgp_s2_n_donors", m2$n_donors); put("dgp_s2_n_pre", m2$n_pre)
  put("dgp_s2_level_treated", m2$params["level", 1]); put("dgp_s2_loading1_treated", m2$params["loading1", 1])
  put("dgp_s2_loading2_treated", m2$params["loading2", 1]); put("dgp_s2_level_max_donor", max(m2$params["level", -1]))
  put("dgp_s2_level_gap", m2$level_gap); put("dgp_s2_hull_l1_gap", m2$hull_l1_gap)
  put("dgp_s2_outside_hull", m2$treated_outside_hull); put("dgp_s2_sigma_e", m2$sigma_e)
  put("dgp_s2_naive_did", m2$did_est)
  put("dgp_s2_tau_plateau", lit("^  tau_path\\[post\\] <- [0-9.]+ \\* \\(1 - exp", "^  tau_path\\[post\\] <- ([0-9.]+) .*$"))
  put("truth_s2_att_post", tr$S2$att_post); put("truth_s2_tau_last", tr$S2$tau_path[12])
  put_months("truth_s2_tau", tr$S2$tau_path, 25:36)
  a2 <- tr$S2$att_post
  put("s2_sc_est", s2$sc$est); put("s2_sc_err", s2$sc$est - a2); put("s2_sc_pre_rmspe", s2$sc$pre_rmspe)
  put("s2_sc_n_nonzero", s2$sc$n_nonzero); put_w("s2_sc", s2$sc$w[s2$sc$w > 1e-6])
  put("s2_ratio_pre_rmspe_s2_to_s1", s2$sc$pre_rmspe / s1$sc_pre_rmspe)
  put("s2_ratio_pre_rmspe_s2_to_s1core", s2$sc$pre_rmspe / core$pre_rmspe)
  put("s2_asc_est", s2$asc$est); put("s2_asc_err", s2$asc$est - a2); put("s2_asc_p_conformal", s2$asc$p_conformal)
  put("s2_asc_n_neg", s2$asc$n_neg); put("s2_asc_min_w", s2$asc$min_w); put("s2_asc_max_w", s2$asc$max_w)
  put("s2_asc_sum_w", s2$asc$sum_w); put("s2_asc_ridge_lambda", s2$asc$ridge_lambda); put("s2_asc_pre_rmspe", s2$asc$pre_rmspe)
  put_w("s2_asc", s2$asc$w)
  put_est("s2_sdid", s2$sdid, a2)
  put("s2_sdid_n_omega_nonzero", s2$sdid$n_omega_nonzero); put("s2_sdid_level_shift", s2$sdid$level_shift)
  put_w("s2_sdid_omega", s2$sdid$omega)                    # keys s2_sdid_omega_w_donorNN
  put("s2_did_est", s2$did$est); put("s2_did_err", s2$did$est - a2)
  put("s2_pen_pen", s2$pen$pen); put("s2_pen_est", s2$pen$est); put("s2_pen_err", s2$pen$est - a2)
  put("s2_pen_pre_rmspe", s2$pen$pre_rmspe); put("s2_pen_n_nonzero", s2$pen$n_nonzero)
  put("s2_pen_max_diff_across_starts", s2$pen$max_diff_across_starts); put("s2_pen_n_orderings", s2$pen$n_orderings)
  put_w("s2_pen", s2$pen$w[s2$pen$w > 1e-6])
  put("s2_pen_grid_min", min(s2$pen$grid)); put("s2_pen_grid_max", max(s2$pen$grid)); put("s2_pen_grid_n", length(s2$pen$grid))
  put("s2_pen_loo_loss_min", min(s2$pen$loo_loss))
  put("s2_pen_holdout_pen", s2$pen$holdout$pen); put("s2_pen_holdout_est", s2$pen$holdout$est)
  put("s2_pen_holdout_n_nonzero", s2$pen$holdout$n_nonzero); put("s2_pen_holdout_pre_rmspe", s2$pen$holdout$pre_rmspe)
  for (k in seq_len(nrow(s2$pen$path))) {
    key <- sprintf("s2_penpath_%02d", k)                    # grid point k = 10^(-4 + (k - 1) / 4)
    put(paste0(key, "_pen"), s2$pen$path$pen[k]); put(paste0(key, "_n_nonzero"), s2$pen$path$n_nonzero[k])
    put(paste0(key, "_est"), s2$pen$path$est[k])
  }

  ## ==================================================================
  ## S3: design, truth, 05_staggered.R, the screen
  ## ==================================================================
  m3 <- A3$meta; d3 <- s3$design; a3 <- tr$S3$att
  put("dgp_s3_seed", tr$seeds[["S3"]]); put("dgp_s3_n_units", d3$n_units); put("dgp_s3_n_months", d3$n_months)
  put("dgp_s3_n_adopters", d3$n_adopters); put("dgp_s3_n_never", d3$n_never); put("dgp_s3_n_contract", d3$n_contract)
  put("dgp_s3_n_wave", d3$n_wave); for (k in seq_along(d3$waves)) put(sprintf("dgp_s3_wave_%d", k), d3$waves[k])
  put("dgp_s3_n_treated_cells", d3$n_treated_cells); put("dgp_s3_tau", m3$tau); put("dgp_s3_sigma_e", m3$sigma_e)
  put("dgp_s3_ar_rho", m3$ar_rho); put("dgp_s3_ar_sd", m3$ar_sd); put("dgp_s3_prio_sd", m3$prio_sd)
  put("dgp_s3_loading1_adopters", d3$loading1_adopters); put("dgp_s3_loading1_never", d3$loading1_never)
  put("dgp_s3_loading1_gap", d3$loading1_adopters - d3$loading1_never)
  lc <- d3$loading1_by_cohort
  for (g in names(lc)) put(paste0("dgp_s3_loading1_cohort_", if (g == "Inf") "never" else g), lc[[g]])
  put("dgp_s3_cor_loading1_adopt_rank", d3$cor_loading1_adopt_rank)
  put("dgp_s3_cor_loading1_cohort_adopters", d3$cor_loading1_cohort_adopters)
  put("dgp_s3_cor_loading1_treated", d3$cor_loading1_treated)
  put("truth_s3_att", a3); put("truth_s3_n_att_gt_cells", nrow(tr$S3$att_gt))
  for (k in c("twfe", "cs", "cs_nyt", "ife", "mc", "mc_default", "multisynth", "sdid_cohort")) put_est(paste0("s3_", k), s3[[k]], a3)
  put("s3_ife_r_cv", s3$ife$r_cv); put("s3_mc_default_lambda_cv", s3$mc_default$lambda_cv)
  put("s3_mc_default_lambda_norm", s3$mc_default$lambda_norm); put("s3_mc_lambda", s3$mc$lambda)
  put("s3_mc_tol", s3$mc$tol); put("s3_mc_max_iter", s3$mc$max_iter); put("s3_mc_niter", s3$mc$niter)
  put("s3_mc_between_truth_twfe", s3$mc$est > min(a3, s3$twfe$est) && s3$mc$est < max(a3, s3$twfe$est))
  put("s3_mc_share_of_twfe_bias_removed", 1 - (s3$mc$est - a3) / (s3$twfe$est - a3))
  put("s3_mc_default_share_of_twfe_bias_removed", 1 - (s3$mc_default$est - a3) / (s3$twfe$est - a3))
  for (g in names(s3$sdid_cohort$by_cohort)) put(paste0("s3_sdid_cohort_est_", g), s3$sdid_cohort$by_cohort[[g]])
  # multisynth fits unit by unit (time_cohort = FALSE, the default): its summary gives each
  # treated metro's average post-period effect (Level = the metro's unit index); the cohort
  # value is the simple mean over the cohort's metros
  ms_bc <- s3$multisynth$by_cohort
  ms_unit <- as.integer(as.character(ms_bc$Level))
  if (anyNA(ms_unit) || !all(is.finite(m3$cohort[ms_unit])))
    stop("collect_numbers22: multisynth levels are not treated unit indices", call. = FALSE)
  for (k in seq_along(ms_unit)) put(sprintf("s3_multisynth_est_unit%02d", ms_unit[k]), ms_bc$Estimate[k])
  for (g in sort(unique(m3$cohort[ms_unit])))
    put(paste0("s3_multisynth_est_cohort", g), mean(ms_bc$Estimate[m3$cohort[ms_unit] == g]))
  put("s3_multisynth_nu", s3$multisynth$nu)
  # multisynth's default horizon (n_leads = NULL in augsynth 0.2.0) is the post window of the
  # last cohort, n_months - last wave + 1 months; its Average row (s3_multisynth_est) weights
  # each treated metro equally over its first n_leads post months, s3_multisynth_n_cells of
  # the dgp_s3_n_treated_cells treated metro-months (Task 14 final review)
  ms_leads <- d3$n_months - max(d3$waves) + 1
  put("s3_multisynth_n_leads", ms_leads); put("s3_multisynth_n_cells", d3$n_adopters * ms_leads)
  put("s3_sdid_cohort_B", s3$sdid_cohort$B); put("s3_fect_nboots", s3$settings$fect_nboots)
  put("s3_cs_biters", s3$settings$cs_biters)
  scr <- s3$screen
  if (is.null(scr)) stop("collect_numbers22: res_stag.rds has no screen (run 05b_stag_screen.R)", call. = FALSE)
  put("s3_screen_n_seeds", scr$n_seeds); put("s3_screen_seed_min", min(scr$seeds)); put("s3_screen_seed_max", max(scr$seeds))
  put("s3_screen_fect_nboots", scr$settings$fect_nboots)
  for (k in names(scr$covered)) {
    put(paste0("s3_screen_covered_", k), scr$covered[[k]])
    put(paste0("s3_screen_share_covered_", k), scr$covered[[k]] / scr$n_seeds)
    put(paste0("s3_screen_task9_covered_", k), scr$task9_record$covered[[k]])
    est <- scr$per_seed[[paste0(k, "_est")]]
    put(paste0("s3_screen_mean_est_", k), mean(est)); put(paste0("s3_screen_mean_err_", k), mean(est) - scr$truth)
    put(paste0("s3_screen_min_est_", k), min(est)); put(paste0("s3_screen_max_est_", k), max(est))
  }
  put("s3_screen_mc_between_mc_default", scr$mc_between$mc_default); put("s3_screen_mc_between_mc", scr$mc_between$mc)
  put("s3_screen_task9_mc_between_mc_default", scr$task9_record$mc_between$mc_default)
  put("s3_screen_task9_mc_between_mc", scr$task9_record$mc_between$mc)
  put("s3_screen_ife_n_r2", sum(scr$per_seed$ife_r_cv == 2))
  is_st <- scr$per_seed$seed == tr$seeds[["S3"]]
  if (sum(is_st) != 1L) stop("collect_numbers22: the screen has no row for the stored seed", call. = FALSE)
  put("s3_screen_stored_mc_covers", scr$per_seed$mc_cover[is_st])
  put("s3_screen_stored_mc_default_covers", scr$per_seed$mc_default_cover[is_st])
  # where the stored draw sits among the screened draws (Task 10 review, Important 1)
  twfe_e <- scr$per_seed$twfe_est
  for (k in names(scr$covered)) {
    e <- scr$per_seed[[paste0(k, "_est")]]
    put(paste0("s3_screen_stored_rank_est_", k), rank(e)[is_st])
    put(paste0("s3_screen_stored_rank_abs_err_", k), rank(abs(e - scr$truth))[is_st])
    put(paste0("s3_screen_n_est_above_truth_", k), sum(e > scr$truth))
  }
  for (k in c("mc", "mc_default")) {
    share <- 1 - (scr$per_seed[[paste0(k, "_est")]] - scr$truth) / (twfe_e - scr$truth)
    put(paste0("s3_screen_mean_share_twfe_bias_removed_", k), mean(share))
    put(paste0("s3_screen_min_share_twfe_bias_removed_", k), min(share))
    put(paste0("s3_screen_max_share_twfe_bias_removed_", k), max(share))
    put(paste0("s3_screen_stored_share_twfe_bias_removed_", k), share[is_st])
    put(paste0("s3_screen_stored_rank_share_twfe_bias_removed_", k), rank(-share)[is_st])
  }

  ## ==================================================================
  ## S4: design, truth, 06_design.R
  ## ==================================================================
  m4 <- A4$meta; sy <- s4$synthetic
  put("dgp_s4_seed", tr$seeds[["S4"]]); put("dgp_s4_n_metros", m4$N); put("dgp_s4_n_history", m4$T_pre)
  put("dgp_s4_n_pilot_months", m4$T_pilot); put("dgp_s4_m_pilot", m4$m_pilot); put("dgp_s4_reps", m4$R)
  put("dgp_s4_rho_g", m4$rho_g); put("dgp_s4_sd_g", m4$sd_g); put("dgp_s4_sigma_e", m4$sigma_e)
  put("truth_s4_tau", tr$S4$tau)
  put("s4_rmse_synthetic_design", s4$rmse_synthetic_design); put("s4_rmse_random_design", s4$rmse_random_design)
  put("s4_rmse_random_did", s4$rmse_random_did)
  put("s4_bias_synthetic", s4$bias_synthetic); put("s4_sd_synthetic", s4$sd_synthetic)
  put("s4_bias_random", s4$bias_random); put("s4_sd_random", s4$sd_random)
  put("s4_bias_random_did", s4$bias_random_did); put("s4_sd_random_did", s4$sd_random_did)
  put("s4_ratio_random_to_synthetic", s4$ratio_random_to_synthetic)
  put("s4_ratio_random_did_to_synthetic", s4$ratio_random_did_to_synthetic)
  put("s4_est_rep1", s4$est_rep1); put("s4_reps", s4$reps)
  put("s4_syn_n_treated", sy$n_treated); put("s4_syn_n_controls", sy$n_controls); put("s4_syn_obj", sy$obj)
  put_w("s4_syn", sy$w)
  put("s4_syn_rmse_fit_treated", sy$rmse_fit_treated); put("s4_syn_rmse_fit_control", sy$rmse_fit_control)
  put("s4_syn_rmse_fit_gap", sy$rmse_fit_gap)
  put("s4_syn_max_v", max(sy$v)); put("s4_syn_n_starts", sy$n_starts); put("s4_syn_n_starts_at_best", sy$n_starts_at_best)

  ## ==================================================================
  ## S3: scpi on the staggered rollout, 05c_stag_scpi.R (Ruling S12)
  ## ==================================================================
  # scdataMulti(effect = "time"): one effect and one prediction interval per event month
  # e = 0..11 (suffix _e00 .. _e11), no interval for their mean; s3_scpi_est is the mean of
  # the 12 effects (the same 180 treated metro-months as s3_multisynth_est)
  sc3 <- rd("res_stag_scpi.rds"); be <- sc3$by_event
  put("s3_scpi_est", sc3$est); put("s3_scpi_err", sc3$est - a3)
  put("s3_scpi_n_months", sc3$n_months); put("s3_scpi_n_cells", sc3$n_cells)
  put("s3_scpi_n_treated", sc3$n_treated); put("s3_scpi_n_donors", sc3$n_donors)
  for (i in seq_len(nrow(be))) {
    sfx <- sprintf("_e%02d", be$e[i])
    put(paste0("s3_scpi_est", sfx), be$est[i]); put(paste0("s3_scpi_ci_lo", sfx), be$lo[i])
    put(paste0("s3_scpi_ci_hi", sfx), be$hi[i]); put(paste0("s3_scpi_covers", sfx), be$covered[i])
  }
  put("s3_scpi_min_est", min(be$est)); put("s3_scpi_max_est", max(be$est))
  put("s3_scpi_ci_lo_min", sc3$lo_min); put("s3_scpi_ci_hi_max", sc3$hi_max); put("s3_scpi_mean_length", sc3$mean_length)
  put("s3_scpi_n_covered", sc3$n_covered); put("s3_scpi_n_covered_in_sample", sc3$n_covered_in)
  put("s3_scpi_n_below_truth", sc3$n_below_truth); put("s3_scpi_n_above_truth", sc3$n_above_truth)
  put("s3_scpi_pre_rmspe_mean", mean(sc3$pre_rmspe)); put("s3_scpi_pre_rmspe_max", max(sc3$pre_rmspe))
  put("s3_scpi_pre_gap_mean", sc3$pre_gap_mean); put("s3_scpi_q_pool12", sc3$q_pool12)   # Theorem 7 q^pool over the last 12 pre months
  put("s3_scpi_sims", sc3$settings$sims); put("s3_scpi_seed", sc3$settings$seed)
  put("s3_scpi_u_alpha", sc3$settings$u_alpha); put("s3_scpi_e_alpha", sc3$settings$e_alpha)

  nums
}
