# ------------------------------------------------------------------
# 06_summary.R -- Collect every quotable number of chapter 20 as a flat
# {name: number} object (numbers.json, full precision), and write sessionInfo.txt.
#
# Run this script LAST (after 01, 02, 03, 04): it reads truth20.rds, the four
# tessera_*.rds data sets and res_fuzzy.rds, res_dose.rds, res_switch.rds.
#
# Name scheme: <script prefix>_<object>_<field>. Prefixes: dgp (01: parameters,
# truth, realised sample facts), fuzzy (02), dose (03), switch (04), summary
# (this script: biases, gaps, coverage flags). Conventions:
#   *_est, *_se       point estimate and standard error
#   *_lb, *_ub        interval produced by the package or by the bootstrap
#                     (percentile 95% for fuzzy_*; the package's own CI for
#                     did_had, did_multiplegt_stat and did_multiplegt_dyn)
#   *_lb95, *_ub95    est -/+ qnorm(.975) * se (contdid) or the t-based interval
#                     of fixest (TWFE regressions)
#   *_truth           the true value (truth20.rds)
#   l<k>, e<k>        horizon l = k (DID_l); placebo k is pl<k>
#   d0p5, d1p5        dose 0.5, 1.5 (decimal point written p)
#   A, Arising        block A and variant A'; B; C
# A p-value below 1e-4 is not stored as a number: the key <name>_p_lt_0p0001 = 1
# is stored instead (also for a package p-value of exactly 0, which is a
# floating-point underflow of 1 - pchisq(), i.e. p < 1e-16).
# Flags (covers_*) are 1 when the interval contains the truth and 0 otherwise.
# Non-finite values are omitted.
# ------------------------------------------------------------------
library(data.table)
library(jsonlite)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
need <- c("truth20.rds", "tessera_A.rds", "tessera_A_rising.rds", "tessera_B.rds", "tessera_B_had.rds", "tessera_C.rds",
          "res_fuzzy.rds", "res_dose.rds", "res_switch.rds")
for (f in need) if (!file.exists(file.path(code_dir, f))) stop(f, " is missing: run 01-04 before 06_summary.R")
tr <- readRDS(file.path(code_dir, "truth20.rds"))
A  <- readRDS(file.path(code_dir, "tessera_A.rds")); AR <- readRDS(file.path(code_dir, "tessera_A_rising.rds"))
B  <- readRDS(file.path(code_dir, "tessera_B.rds")); H <- readRDS(file.path(code_dir, "tessera_B_had.rds"))
C  <- readRDS(file.path(code_dir, "tessera_C.rds"))
f  <- readRDS(file.path(code_dir, "res_fuzzy.rds")); dz <- readRDS(file.path(code_dir, "res_dose.rds"))
sw <- readRDS(file.path(code_dir, "res_switch.rds"))

nums <- list()
put <- function(name, value) {
  stopifnot(length(name) == 1, !name %in% names(nums))
  if (length(value) == 1 && is.numeric(value) && is.finite(value)) nums[[name]] <<- as.numeric(value)
  invisible(NULL)
}
put_p <- function(name, p) {                                  # name: key stem without the trailing _p
  if (length(p) != 1 || !is.finite(p)) return(invisible(NULL))
  if (p < 1e-4) put(paste0(name, "_p_lt_0p0001"), 1) else put(paste0(name, "_p"), p)
}
put_r <- function(prefix, r, ci = c("lb", "ub"), fields = c("est", "se")) {   # list(est, se, lo, hi, ...)
  for (k in fields) if (!is.null(r[[k]])) put(paste0(prefix, "_", k), r[[k]])
  if (!is.null(r$lo)) { put(paste0(prefix, "_", ci[1]), r$lo); put(paste0(prefix, "_", ci[2]), r$hi) }
}
pid <- function(x) gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE))     # 1.5 -> "1p5"
cov_flag <- function(r, v) as.numeric(r$lo <= v && v <= r$hi)
nearest <- function(grid, x) { i <- which.min(abs(grid - x)); stopifnot(abs(grid[i] - x) < 1e-9); i }

## ====================================================================
## dgp (01_dgp.R): parameters, truth, realised sample facts
## ====================================================================
for (nm in names(tr$params$seed)) put(paste0("dgp_seed_", nm), tr$params$seed[[nm]])
for (blk in c("A", "B")) for (nm in names(tr$params[[blk]])) put(sprintf("dgp_%s_%s", blk, nm), tr$params[[blk]][[nm]])
PC <- tr$params$C
for (nm in setdiff(names(PC), c("trans", "f", "omega"))) put(paste0("dgp_C_", nm), PC[[nm]])
for (i in 1:3) for (j in 1:3) put(sprintf("dgp_C_trans_%d_%d", i - 1, j - 1), PC$trans[i, j])    # from tier i-1 to tier j-1
for (k in seq_along(PC$f)) put(sprintf("dgp_C_f_%d", k - 1), PC$f[k])
for (k in seq_along(PC$omega)) put(sprintf("dgp_C_omega_%d", k - 1), PC$omega[k])

## block A
put("dgp_A_n_obs", nrow(A)); put("dgp_A_n_couriers", uniqueN(A$unit)); put("dgp_A_n_cities", uniqueN(A$city))
put("dgp_A_n_cities_treated", uniqueN(A[g == 1, city])); put("dgp_A_n_cities_control", uniqueN(A[g == 0, city]))
put("dgp_A_n_couriers_treated", uniqueN(A[g == 1, unit])); put("dgp_A_n_couriers_control", uniqueN(A[g == 0, unit]))
sh <- tr$A$shares
for (i in seq_len(nrow(sh))) {
  k <- sprintf("dgp_A_share_%s%d", sh$type_true[i], sh$g[i]); put(k, sh$share[i]); put(paste0(k, "_n"), sh$N[i])
  if ("share_pop" %in% names(sh)) put(paste0(k, "_pop"), sh$share_pop[i])
}
shr <- tr$A_rising$shares
for (i in seq_len(nrow(shr))) { k <- sprintf("dgp_Arising_share_%s%d", shr$type_true[i], shr$g[i]); put(k, shr$share[i]); put(paste0(k, "_n"), shr$N[i]) }
put("dgp_A_late", tr$A$late); put("dgp_A_late_pop", tr$A$late_pop); stopifnot(tr$A$lambda == tr$params$A$lambda)
put("dgp_A_w_did_pop", tr$A$w_did_pop); put("dgp_A_bias_did_pop", tr$A$bias_did_pop)
for (i in seq_len(nrow(tr$A$att_always))) put(sprintf("dgp_A_att_always_g%d_t%d", tr$A$att_always$g[i], tr$A$att_always$t[i]), tr$A$att_always$att[i])
for (nm in names(tr$A$att_always_pop)) put(paste0("dgp_A_att_always_pop_", nm), tr$A$att_always_pop[[nm]])
for (nm in setdiff(names(tr$A$pop_check), "n_city_per_group")) put(paste0("dgp_A_popcheck_", nm), tr$A$pop_check[[nm]])
put("dgp_A_popcheck_n_city_per_group", tr$A$pop_check$n_city_per_group)
## variant A'
put("dgp_Arising_late", tr$A_rising$late); put("dgp_Arising_late_pop", tr$A_rising$late_pop)
put("dgp_Arising_late_control", tr$A_rising$late_control)
for (nm in setdiff(names(tr$A_rising$pop_check), "n_city_per_group")) put(paste0("dgp_Arising_popcheck_", nm), tr$A_rising$pop_check[[nm]])
put("dgp_Arising_tc_bias_switchers_vs_always_Y1", tr$A_rising$tc_bias_parts[["switchers_vs_always_Y1"]])
put("dgp_Arising_tc_bias_never_vs_switchers_Y0", tr$A_rising$tc_bias_parts[["never_vs_switchers_Y0"]])
put("dgp_Arising_tc_bias_total", sum(tr$A_rising$tc_bias_parts))
put("dgp_Arising_popcheck_tc_bias", tr$A_rising$pop_check$w_tc - tr$A_rising$pop_check$late)
put("dgp_Arising_popcheck_cic_bias", tr$A_rising$pop_check$w_cic - tr$A_rising$pop_check$late)
put("dgp_Arising_popcheck_did_bias", tr$A_rising$pop_check$w_did - tr$A_rising$pop_check$late)
for (gg in 0:1) for (tt in 0:1) put(sprintf("dgp_Arising_rate_g%d_t%d", gg, tt), AR[AR$g == gg & AR$t == tt, mean(d)])
## block B
put("dgp_B_n_cities", uniqueN(B$unit)); put("dgp_B_n_obs", nrow(B))
put("dgp_B_n_treated", B[post == 1, sum(dose > 0)]); put("dgp_B_n_stayers", B[post == 1, sum(dose == 0)])
put("dgp_B_dose_mean_treated", B[post == 1 & dose > 0, mean(dose)]); put("dgp_B_dose_min_treated", B[post == 1 & dose > 0, min(dose)])
put("dgp_B_dose_max_treated", B[post == 1 & dose > 0, max(dose)])
put("dgp_B_had_dose_min", H[post == 1, min(dose)]); put("dgp_B_had_share_below_025", H[post == 1, mean(dose < 0.25)])
put("dgp_B_had_n_below_025", H[post == 1, sum(dose < 0.25)])
for (nm in setdiff(names(tr$B), c("att_curve", "definition", "share_zero"))) put(paste0("dgp_B_", nm), tr$B[[nm]])
put("dgp_B_share_zero_realised", tr$B$share_zero)                  # the parameter dgp_B_share_zero is the target share
cur <- tr$B$att_curve
for (dd in c(0.5, 1, 1.5, 2, 2.5, 3)) {
  i <- nearest(cur$dose, dd)
  for (nm in setdiff(names(cur), "dose")) put(sprintf("dgp_B_curve_d%s_%s", pid(dd), nm), cur[[nm]][i])
}
## block C
put("dgp_C_n_cities", uniqueN(C$unit)); put("dgp_C_n_obs", nrow(C)); put("dgp_C_n_cells_treated", C[, sum(dose > 0)])
put("dgp_C_share_cells_treated", C[, mean(dose > 0)])
for (k in 1:2) { put(sprintf("dgp_C_n_cells_tier%d", k), C[, sum(dose == k)]); put(sprintf("dgp_C_share_cells_tier%d", k), C[, mean(dose == k)]) }
put("dgp_C_n_never_switchers", C[, uniqueN(unit[g == 0])]); put("dgp_C_n_switchers", C[, uniqueN(unit[g > 0])])
put("dgp_C_first_switch_min", C[g > 0, min(g)]); put("dgp_C_first_switch_max", C[g > 0, max(g)])
put("dgp_C_first_tier2_share", C[g > 0 & t == g, mean(dose == 2)])
nsw <- C[, .(n = sum(diff(dose) != 0)), by = unit]
put("dgp_C_n_switch_events_mean_switchers", nsw[unit %in% C[g > 0, unique(unit)], mean(n)])
for (nm in c("avg_total_effect", "avg_effect", "avg_effect_unweighted", "share_switch_2plus")) put(paste0("dgp_C_", nm), tr$C[[nm]])
for (l in 1:8) {
  put(sprintf("dgp_C_delta_l%d", l), tr$C$delta[l]); put(sprintf("dgp_C_delta_norm_l%d", l), tr$C$delta_norm[l])
  put(sprintf("dgp_C_dose_cum_l%d", l), tr$C$dose_cum[l]); put(sprintf("dgp_C_dose_inc_l%d", l), tr$C$dose_inc[l])
  put(sprintf("dgp_C_n_switchers_l%d", l), tr$C$n_switchers[l])
}

## ====================================================================
## fuzzy (02_fuzzy.R)
## ====================================================================
put("fuzzy_boot_B", f$meta$B); put("fuzzy_boot_seed", f$meta$seed)
for (k in c("did", "tc", "cic")) {
  put_r(paste0("fuzzy_", k), f[[k]])                       # block A: est, se, lb, ub (percentile)
  put_r(paste0("fuzzy_Arising_", k), f$rising[[k]])        # variant A'
}
for (v in c("A", "rising")) {
  pr <- f$parts[[v]]; pre <- if (v == "A") "fuzzy_A_" else "fuzzy_Arising_"
  for (nm in names(pr)) put(paste0(pre, nm), pr[[nm]])
}
for (i in seq_len(nrow(f$cells))) {
  z <- f$cells[i, ]; pre <- sprintf("fuzzy_cell_%s_g%d_t%d", if (z$variant == "A") "A" else "Arising", z$g, z$t)
  put(paste0(pre, "_n"), z$n); put(paste0(pre, "_take_up"), z$take_up); put(paste0(pre, "_mean_y"), z$mean_y)
}
put("fuzzy_first_stage_A", f$cells[f$cells$variant == "A" & f$cells$g == 1 & f$cells$t == 1, "take_up"] -
                           f$cells[f$cells$variant == "A" & f$cells$g == 1 & f$cells$t == 0, "take_up"])

## ====================================================================
## dose (03_dose.R)
## ====================================================================
cn <- dz$cont
put_r("dose_cont_overall", cn$overall, ci = c("lb95", "ub95"))
put_r("dose_cont_overall_boot", cn$overall_boot)
put_r("dose_cont_acrt_overall", cn$acrt_overall, ci = c("lb95", "ub95"))
put_r("dose_cont_acrt_overall_boot", cn$acrt_overall_boot)
for (blk in c("att", "acrt_vs_slope", "acrt_vs_true_acrt")) for (nm in names(cn$fit[[blk]])) put(sprintf("dose_cont_fit_%s_%s", blk, nm), cn$fit[[blk]][[nm]])
for (nm in c("se_ratio_curve_att", "se_ratio_curve_acrt", "se_ratio_overall_acrt", "se_ratio_overall_att", "n_over_n_treated", "B_city"))
  put(paste0("dose_cont_check_", nm), cn$check[[nm]])
for (nm in names(cn$check$replication_error)) put(paste0("dose_cont_check_replication_error_", nm), cn$check$replication_error[[nm]])
for (nm in c("biters", "seed", "n_dvals", "n", "n_treated", "n_untreated", "degree", "num_knots")) put(paste0("dose_cont_spec_", nm), cn$spec[[nm]])
for (i in seq_len(nrow(cn$robust))) {
  z <- cn$robust[i, ]; pre <- sprintf("dose_cont_robust_deg%d_k%d", z$degree, z$num_knots)
  for (nm in c("att", "att_se", "acrt", "acrt_se", "curve_rmse", "curve_max_abs_err")) put(paste0(pre, "_", nm), z[[nm]])
}
cv <- cn$curve; ac <- cn$acrt_curve
put("dose_cont_dvals_min", min(cv$dose)); put("dose_cont_dvals_max", max(cv$dose))
for (dd in c(0.5, 1, 1.5, 2, 2.5)) {
  i <- nearest(cv$dose, dd); pre <- sprintf("dose_cont_curve_d%s", pid(dd))
  put_r(pre, list(est = cv$est[i], se = cv$se[i], lo = cv$lo[i], hi = cv$hi[i]), ci = c("lb95", "ub95"))
  put(paste0(pre, "_truth"), cv$truth[i]); put(paste0(pre, "_truth_strong"), cv$truth_strong[i])
  pre <- sprintf("dose_cont_acrt_curve_d%s", pid(dd))
  put_r(pre, list(est = ac$est[i], se = ac$se[i], lo = ac$lo[i], hi = ac$hi[i]), ci = c("lb95", "ub95"))
  put(paste0(pre, "_truth_slope"), ac$truth_slope[i]); put(paste0(pre, "_truth_acrt"), ac$truth_acrt[i])
}
for (end in c("min", "max")) {                                   # the two ends of the fitted curve (smallest and largest treated dose)
  i <- if (end == "min") which.min(cv$dose) else which.max(cv$dose); pre <- paste0("dose_cont_curve_", end)
  put(paste0(pre, "_dose"), cv$dose[i]); put_r(pre, list(est = cv$est[i], se = cv$se[i], lo = cv$lo[i], hi = cv$hi[i]), ci = c("lb95", "ub95"))
  put(paste0(pre, "_truth"), cv$truth[i])
}
put_r("dose_binary", dz$binary, ci = c("lb95", "ub95")); put_p("dose_binary", dz$binary$p); put("dose_binary_t", dz$binary$est / dz$binary$se)
put_r("dose_twfe", dz$twfe, ci = c("lb95", "ub95")); put_p("dose_twfe", dz$twfe$p); put("dose_twfe_t", dz$twfe$est / dz$twfe$se)
put("dose_twfe_df", dz$twfe$df); put("dose_twfe_oracle_est", dz$twfe$oracle_est)
for (nm in setdiff(names(dz$twfe$weights), character())) put(paste0("dose_twfe_weights_", nm), dz$twfe$weights[[nm]])
h <- dz$had
put_r("dose_had", h, ci = c("lb", "ub"))
for (nm in c("n", "bw", "n_bw", "qug_stat", "qug_p", "min_dose", "share_below_025")) put(paste0("dose_had_", nm), h[[nm]])
for (nm in names(h$yatchew)) put(paste0("dose_had_yatchew_", nm), h$yatchew[[nm]])
put_r("dose_had_twfe", h$twfe, ci = c("lb95", "ub95")); put_p("dose_had_twfe", h$twfe$p)
for (e in c("aoss", "waoss")) {
  put_r(paste0("dose_stat_", e), dz$stat[[e]], ci = c("lb", "ub"))
  put(paste0("dose_stat_", e, "_switchers"), dz$stat[[e]]$switchers); put(paste0("dose_stat_", e, "_stayers"), dz$stat[[e]]$stayers)
}
put_r("dose_stat_diff", dz$stat$diff, ci = c("lb", "ub")); put("dose_stat_diff_t", dz$stat$diff$t); put_p("dose_stat_diff", dz$stat$diff$p)

## ====================================================================
## switch (04_switching.R)
## ====================================================================
tw <- sw$twfe
put_r("switch_twfe", tw, ci = c("lb95", "ub95")); put_p("switch_twfe", tw$p); put("switch_twfe_t", tw$est / tw$se)
put("switch_twfe_df", tw$df); put("switch_twfe_oracle_est", tw$oracle_est); put("switch_twfe_n_obs", tw$n_obs)
W <- sw$weights
for (nm in setdiff(names(W), c("by_event", "by_dose"))) put(paste0("switch_weights_", nm), W[[nm]])
for (i in seq_len(nrow(W$by_dose))) { z <- W$by_dose[i, ]; for (nm in setdiff(names(z), "dose")) put(sprintf("switch_weights_tier%d_%s", z$dose, nm), z[[nm]]) }
for (i in seq_len(nrow(W$by_event))) { z <- W$by_event[i, ]; for (nm in c("n_cells", "n_negative", "sum_weight", "sum_negative_weight")) put(sprintf("switch_weights_e%d_%s", z$event_time, nm), z[[nm]]) }
put("switch_weights_share_negative_mass", -W$sum_negative / (W$sum_positive - W$sum_negative))          # |negative| / (positive + |negative|)
D <- sw$dyn
for (l in 1:8) {
  pre <- sprintf("switch_dyn_l%d", l); put_r(pre, D[[paste0("l", l)]], ci = c("lb", "ub")); put(paste0(pre, "_n"), D[[paste0("l", l)]]$n)
  put(paste0(pre, "_switchers"), D[[paste0("l", l)]]$switchers); put(paste0(pre, "_truth"), tr$C$delta[l])
  pre <- sprintf("switch_dynnorm_l%d", l); put_r(pre, D$norm[[paste0("l", l)]], ci = c("lb", "ub")); put(paste0(pre, "_truth"), tr$C$delta_norm[l])
}
for (l in 1:4) {
  for (v in c("", "norm")) {
    z <- if (v == "") D[[paste0("pl", l)]] else D$norm[[paste0("pl", l)]]
    pre <- sprintf("switch_dyn%s_pl%d", v, l); put_r(pre, z, ci = c("lb", "ub")); put(paste0(pre, "_n"), z$n)
  }
}
put_r("switch_dyn_avg_total", D$avg_total, ci = c("lb", "ub")); put("switch_dyn_avg_total_n", D$avg_total$n)
put("switch_dyn_avg_total_switchers", D$avg_total$switchers); put("switch_dyn_avg_total_truth", tr$C$avg_total_effect)
put("switch_dyn_delta_D_avg_total", D$delta_D_avg_total)
put_p("switch_dyn_jointplacebo", D$p_jointplacebo); put_p("switch_dyn_jointeffects", D$p_jointeffects)
put_p("switch_dynnorm_jointplacebo", D$norm$p_jointplacebo); put_p("switch_dynnorm_jointeffects", D$norm$p_jointeffects)
put("switch_dyn_n_placebos_excluding_zero", sum(D$placebos$lo > 0 | D$placebos$hi < 0))

## ====================================================================
## summary (this script): biases, gaps and coverage flags
## ====================================================================
for (k in c("did", "tc", "cic")) {
  put(sprintf("summary_fuzzy_%s_bias", k), f[[k]]$est - tr$A$late)
  put(sprintf("summary_fuzzy_%s_bias_in_se", k), (f[[k]]$est - tr$A$late) / f[[k]]$se)
  put(sprintf("summary_fuzzy_%s_covers_late", k), cov_flag(f[[k]], tr$A$late))
  put(sprintf("summary_fuzzy_Arising_%s_bias", k), f$rising[[k]]$est - tr$A_rising$late)
  put(sprintf("summary_fuzzy_Arising_%s_bias_in_se", k), (f$rising[[k]]$est - tr$A_rising$late) / f$rising[[k]]$se)
  put(sprintf("summary_fuzzy_Arising_%s_covers_late", k), cov_flag(f$rising[[k]], tr$A_rising$late))
}
put("summary_fuzzy_did_bias_closed_form", with(tr$params$A, (pA1 - pA0) * lambda / pS1))      # (pA1 - pA0) lambda / pS1
put("summary_fuzzy_tc_minus_cic", f$tc$est - f$cic$est)
put("summary_dose_cont_overall_minus_truth", cn$overall$est - tr$B$att_overall)
put("summary_dose_cont_overall_covers_truth", cov_flag(cn$overall, tr$B$att_overall))
put("summary_dose_cont_acrt_overall_minus_true_acrt", cn$acrt_overall$est - tr$B$acrt_overall)
put("summary_dose_cont_acrt_overall_minus_pt_slope", cn$acrt_overall$est - tr$B$acrt_overall_pt_pop)
put("summary_dose_cont_acrt_overall_covers_true_acrt", cov_flag(cn$acrt_overall, tr$B$acrt_overall))
put("summary_dose_cont_acrt_overall_covers_pt_slope", cov_flag(cn$acrt_overall, tr$B$acrt_overall_pt_pop))
put("summary_dose_cont_acrt_overall_boot_covers_true_acrt", cov_flag(cn$acrt_overall_boot, tr$B$acrt_overall))
put("summary_dose_cont_acrt_overall_boot_covers_pt_slope", cov_flag(cn$acrt_overall_boot, tr$B$acrt_overall_pt_pop))
put("summary_dose_twfe_minus_was", dz$twfe$est - tr$B$was); put("summary_dose_twfe_minus_as", dz$twfe$est - tr$B$as)
put("summary_dose_twfe_minus_true_acrt", dz$twfe$est - tr$B$acrt_overall)
put("summary_dose_twfe_minus_pt_slope", dz$twfe$est - tr$B$acrt_overall_pt_pop)
put("summary_dose_twfe_covers_was", cov_flag(dz$twfe, tr$B$was)); put("summary_dose_twfe_covers_as", cov_flag(dz$twfe, tr$B$as))
put("summary_dose_twfe_covers_true_acrt", cov_flag(dz$twfe, tr$B$acrt_overall))
put("summary_dose_had_minus_truth", h$est - tr$B$had_was); put("summary_dose_had_covers_truth", cov_flag(h, tr$B$had_was))
put("summary_dose_stat_aoss_minus_truth", dz$stat$aoss$est - tr$B$as); put("summary_dose_stat_aoss_covers_truth", cov_flag(dz$stat$aoss, tr$B$as))
put("summary_dose_stat_waoss_minus_truth", dz$stat$waoss$est - tr$B$was); put("summary_dose_stat_waoss_covers_truth", cov_flag(dz$stat$waoss, tr$B$was))
put("summary_dose_curve_share_covered", cn$fit$att$share_covered)
put("summary_switch_twfe_minus_avg_effect", tw$est - tr$C$avg_effect); put("summary_switch_twfe_bias_in_se", (tw$est - tr$C$avg_effect) / tw$se)
put("summary_switch_twfe_share_of_avg_effect", tw$est / tr$C$avg_effect)
put("summary_switch_twfe_minus_avg_effect_unweighted", tw$est - tr$C$avg_effect_unweighted)
put("summary_switch_twfe_covers_avg_effect", cov_flag(tw, tr$C$avg_effect))
put("summary_switch_twfe_minus_dyn_avg_total", tw$est - D$avg_total$est)
put("summary_switch_dyn_avg_total_minus_truth", D$avg_total$est - tr$C$avg_total_effect)
put("summary_switch_dyn_avg_total_covers_truth", cov_flag(D$avg_total, tr$C$avg_total_effect))
for (l in 1:8) {
  put(sprintf("summary_switch_dyn_l%d_minus_truth", l), D[[paste0("l", l)]]$est - tr$C$delta[l])
  put(sprintf("summary_switch_dyn_l%d_covers_truth", l), cov_flag(D[[paste0("l", l)]], tr$C$delta[l]))
  put(sprintf("summary_switch_dynnorm_l%d_minus_truth", l), D$norm[[paste0("l", l)]]$est - tr$C$delta_norm[l])
  put(sprintf("summary_switch_dynnorm_l%d_covers_truth", l), cov_flag(D$norm[[paste0("l", l)]], tr$C$delta_norm[l]))
}
put("summary_switch_dyn_n_covered", sum(vapply(1:8, function(l) cov_flag(D[[paste0("l", l)]], tr$C$delta[l]), 1)))

## ---- self-consistency with the res_*.rds files and write -------------------------------
stopifnot(length(nums) > 0, all(vapply(nums, function(z) is.numeric(z) && length(z) == 1 && is.finite(z), TRUE)))
stopifnot(nums$fuzzy_tc_est == f$tc$est, nums$fuzzy_did_lb == f$did$lo, nums$dose_cont_overall_est == dz$cont$overall$est,
          nums$switch_twfe_est == sw$twfe$est, nums$switch_dyn_l1_est == sw$dyn$l1$est, nums$switch_weights_share_negative == sw$weights$share_negative,
          nums$dgp_C_delta_l1 == tr$C$delta[1], nums$dgp_A_late == tr$A$late, nums$dgp_B_att_overall == tr$B$att_overall)
out_json <- file.path(code_dir, "numbers.json")
write_json(nums, out_json, auto_unbox = TRUE, digits = NA, pretty = TRUE)
chk <- fromJSON(out_json)                                          # round trip at the stored precision
stopifnot(identical(sort(names(chk)), sort(names(nums))),
          max(abs(unlist(chk)[names(nums)] - unlist(nums)) / pmax(1, abs(unlist(nums)))) < 1e-13)
cat(sprintf("\nSaved numbers.json (%d numbers)\n", length(nums)))

## ---- sessionInfo.txt: lists every package that 02-05 attach -------------------------------------
suppressPackageStartupMessages({
  library(fixest); library(splines2); library(contdid); library(DIDHAD); library(DIDmultiplegtSTAT); library(TwoWayFEWeights)
  library(polars); library(DIDmultiplegtDYN); library(ggplot2); library(svglite)
  # Rfuzzydid is used by tests/test_fuzzy_fns.R only; attach it here when installed
  if (requireNamespace("Rfuzzydid", quietly = TRUE)) library(Rfuzzydid)
})
writeLines(capture.output(sessionInfo()), file.path(code_dir, "sessionInfo.txt"))
cat("Saved sessionInfo.txt\n")
