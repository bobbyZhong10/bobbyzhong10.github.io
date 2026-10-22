# ------------------------------------------------------------------
# 06_summary.R -- Collect all estimates against the true simple ATT,
# and write numbers.json: every quotable number of chapter 19 as a flat
# {name: number} object at full precision.
#
# Run this script LAST (after 01-04 and 07-09): it reads every res_*.rds.
#
# Name scheme: <script prefix>_<object>_<field>. Prefixes: dgp (01),
# twfe (02), bacon (03), modern (04), dcdh (07), honest (08), family (09),
# summary (this script). Event-time suffixes: e6 = event time 6,
# em6 = event time -6; dcdh DID_l is dcdh_dyn_l<l>_*, placebo l is
# dcdh_pl_l<l>_* (DID_l <-> e = l-1, placebo l <-> e = -1-l).
# Intervals without a package CI are normal-approximation 95% intervals
# (est +/- qnorm(.975) * se); the static TWFE interval is fixest's
# t-based one (59 df). Non-finite values are omitted.
# ------------------------------------------------------------------
library(data.table)
library(fixest)
library(jsonlite)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
truth  <- readRDS(file.path(code_dir, "truth.rds"))
r_twfe <- readRDS(file.path(code_dir, "res_twfe.rds"))
r_bac  <- readRDS(file.path(code_dir, "res_bacon.rds"))
r_mod  <- readRDS(file.path(code_dir, "res_modern.rds"))
for (f in c("res_dcdh.rds", "res_honest.rds", "res_family.rds", "dgp_params.rds", "kestrel_panel.rds"))
  if (!file.exists(file.path(code_dir, f)))
    stop(f, " is missing: run 01-04 and 07-09 before 06_summary.R")
r_dcdh <- readRDS(file.path(code_dir, "res_dcdh.rds"))
r_hon  <- readRDS(file.path(code_dir, "res_honest.rds"))
r_fam  <- readRDS(file.path(code_dir, "res_family.rds"))
dgp_par <- readRDS(file.path(code_dir, "dgp_params.rds"))
dt     <- readRDS(file.path(code_dir, "kestrel_panel.rds"))

tab <- data.table(
  estimator = c("True simple ATT (DGP)",
                "Static TWFE (feols)",
                "CS simple ATT (att_gt, not-yet-treated)",
                "SA average ATT (sunab)",
                "BJS imputation (did_imputation)"),
  estimate = c(truth$true_att_simple,
               r_twfe$static$est,
               r_mod$cs_simple$est,
               r_mod$sa_att$est,
               r_mod$imp$est),
  se = c(NA_real_,
         r_twfe$static$se,
         r_mod$cs_simple$se,
         r_mod$sa_att$se,
         r_mod$imp$se)
)
tab[, `:=`(estimate = round(estimate, 4), se = round(se, 4))]
tab[, share_of_truth := round(estimate / truth$true_att_simple, 3)]

cat("== Module 11 summary: estimators vs truth ==\n")
print(tab)

cat("\n== Goodman-Bacon weights by comparison type ==\n")
print(r_bac$summary[, .(type, weight = round(weight, 4),
                        avg_estimate = round(avg_estimate, 4))])

cat("\n== CS group ATTs vs truth by cohort ==\n")
cmp <- merge(r_mod$cs_group, truth$true_att_by_g, by = "g")
print(cmp[, .(g, cs_att = round(att.x, 4), cs_se = round(se, 4),
              true_att = round(att.y, 4))])

fwrite(tab, file.path(code_dir, "summary_table.csv"))
cat("\nSaved summary_table.csv\n")

## ====================================================================
## numbers.json
## ====================================================================
nums <- list()
put <- function(name, value) {
  stopifnot(length(name) == 1, !name %in% names(nums))
  if (length(value) == 1 && is.numeric(value) && is.finite(value)) nums[[name]] <<- as.numeric(value)
  invisible(NULL)
}
eid <- function(e) ifelse(e < 0, paste0("m", -e), as.character(e))           # -6 -> "m6"
gridid <- function(x) gsub("\\.", "p", format(x, scientific = FALSE, drop0trailing = TRUE))
put_ci <- function(prefix, est, se) {
  put(paste0(prefix, "_est"), est); put(paste0(prefix, "_se"), se)
  put(paste0(prefix, "_lb95"), est - qnorm(0.975) * se)
  put(paste0(prefix, "_ub95"), est + qnorm(0.975) * se)
}
put_es <- function(prefix, es, truth_theta = NULL) {                          # es: data.frame(e, estimate, se)
  for (i in seq_len(nrow(es))) {
    k <- paste0(prefix, "_e", eid(es$e[i]))
    put(paste0(k, "_est"), es$estimate[i]); put(paste0(k, "_se"), es$se[i])
  }
}
true_att <- truth$true_att_simple

## ---- dgp (01_dgp.R) -------------------------------------------------
put("dgp_seed", dgp_par$seed); put("dgp_t_months", dgp_par$T_months)
put("dgp_n_cities", dgp_par$n_cities); put("dgp_n_never", dgp_par$n_never)
put("dgp_n_obs", nrow(dt)); put("dgp_n_treated_obs", sum(dt$D == 1))
put("dgp_share_treated_obs", mean(dt$D == 1))
for (k in seq_along(dgp_par$cohort_g)) {
  g <- dgp_par$cohort_g[k]
  put(sprintf("dgp_g%d_month", g), g); put(sprintf("dgp_g%d_size", g), dgp_par$cohort_size[k])
  put(sprintf("dgp_g%d_boost", g), dgp_par$cohort_boost[k])
  put(sprintf("dgp_g%d_level_shift", g), dgp_par$level_shift[k])
}
for (nm in c("tau0", "tau5", "plateau", "slope_late", "mult_sd", "sd_cityfe",
             "trend_common", "season_amp", "ar_rho", "ar_sd"))
  put(paste0("dgp_", nm), dgp_par[[nm]])
put("dgp_true_att_simple", true_att)
for (i in seq_len(nrow(truth$true_att_by_g))) {
  g <- truth$true_att_by_g$g[i]
  put(sprintf("dgp_true_att_g%d", g), truth$true_att_by_g$att[i])
  put(sprintf("dgp_true_att_g%d_cells", g), truth$true_att_by_g$n_cells[i])
}
for (i in seq_len(nrow(truth$true_theta_e))) {
  k <- paste0("dgp_true_theta_e", truth$true_theta_e$event_time[i])
  put(k, truth$true_theta_e$theta[i]); put(paste0(k, "_cells"), truth$true_theta_e$n_cells[i])
}

## ---- twfe (02_twfe.R) -----------------------------------------------
m_static <- feols(lgmv ~ D | city_id + month, data = dt, cluster = ~city_id)
stopifnot(abs(coef(m_static)[["D"]] - r_twfe$static$est) < 1e-12,
          abs(se(m_static)[["D"]]   - r_twfe$static$se)  < 1e-12)
ci_static <- confint(m_static)
put("twfe_static_est", r_twfe$static$est); put("twfe_static_se", r_twfe$static$se)
put("twfe_static_t", r_twfe$static$est / r_twfe$static$se)
put("twfe_static_p", pvalue(m_static)[["D"]])
put("twfe_static_lb95", ci_static[["2.5 %"]]); put("twfe_static_ub95", ci_static[["97.5 %"]])
put("twfe_static_df", degrees_freedom(m_static, "t"))
put("twfe_static_share_of_truth", r_twfe$static$est / true_att)
put_es("twfe_naive", r_twfe$naive_es)

## ---- bacon (03_bacon.R) ---------------------------------------------
bc <- r_bac$components
lab <- function(x) ifelse(x == 99999, "never", paste0("g", x))
for (i in seq_len(nrow(bc))) {
  k <- sprintf("bacon_%s_vs_%s", lab(bc$treated[i]), lab(bc$untreated[i]))
  put(paste0(k, "_weight"), bc$weight[i]); put(paste0(k, "_est"), bc$estimate[i])
  put(paste0(k, "_contribution"), bc$weight[i] * bc$estimate[i])
}
tkey <- c("Treated vs Untreated" = "treated_vs_untreated",
          "Earlier vs Later Treated" = "earlier_vs_later",
          "Later vs Earlier Treated" = "later_vs_earlier")
for (i in seq_len(nrow(r_bac$summary))) {
  k <- paste0("bacon_type_", tkey[[r_bac$summary$type[i]]])
  put(paste0(k, "_weight"), r_bac$summary$weight[i])
  put(paste0(k, "_avg_est"), r_bac$summary$avg_estimate[i])
  put(paste0(k, "_contribution"), r_bac$summary$contribution[i])
}
put("bacon_total_weight", sum(bc$weight)); put("bacon_total_est", sum(bc$weight * bc$estimate))

## ---- modern (04_modern.R) -------------------------------------------
put_ci("modern_sa_att", r_mod$sa_att$est, r_mod$sa_att$se)
put_ci("modern_cs_simple", r_mod$cs_simple$est, r_mod$cs_simple$se)
put_ci("modern_imp", r_mod$imp$est, r_mod$imp$se)
for (i in seq_len(nrow(r_mod$cs_group))) {
  g <- r_mod$cs_group$g[i]
  put_ci(sprintf("modern_cs_group_g%d", g), r_mod$cs_group$att[i], r_mod$cs_group$se[i])
  put(sprintf("modern_cs_group_g%d_minus_truth", g),
      r_mod$cs_group$att[i] - truth$true_att_by_g[g == r_mod$cs_group$g[i], att])
}
put_ci("modern_cs_group_overall", r_mod$cs_group_overall$est, r_mod$cs_group_overall$se)
put_es("modern_sa_es", r_mod$sa_es)
put_es("modern_cs_es", r_mod$cs_es)
put_es("modern_imp_es", r_mod$imp_es)

## ---- dcdh (07_dcdh.R) -----------------------------------------------
W <- r_dcdh$twfe_weights; S <- r_dcdh$twfe_summary
put("dcdh_w_beta", S$beta); put("dcdh_w_n_cells", S$n_cells)
put("dcdh_w_n_positive", S$n_positive); put("dcdh_w_n_negative", S$n_negative)
put("dcdh_w_share_negative", S$share_negative)
put("dcdh_w_sum_positive", S$sum_positive); put("dcdh_w_sum_negative", S$sum_negative)
put("dcdh_w_min_sd_beta_delta0", S$min_sd_beta_delta0)
put("dcdh_w_min_sd_opposite_sign", S$min_sd_opposite_sign)
put("dcdh_w_min", min(W$weight)); put("dcdh_w_max", max(W$weight))
put("dcdh_w_mean_negative", mean(W$weight[W$weight < 0]))
for (i in seq_len(nrow(S$by_cohort))) {
  k <- sprintf("dcdh_w_g%d", S$by_cohort$cohort[i])
  put(paste0(k, "_n_cells"), S$by_cohort$n_cells[i]); put(paste0(k, "_n_negative"), S$by_cohort$n_negative[i])
  put(paste0(k, "_sum_weight"), S$by_cohort$sum_weight[i])
  put(paste0(k, "_sum_negative"), S$by_cohort$sum_negative_weight[i])
}
for (i in seq_len(nrow(S$by_event))) {
  k <- sprintf("dcdh_w_e%d", S$by_event$event_time[i])
  put(paste0(k, "_n_negative"), S$by_event$n_negative[i]); put(paste0(k, "_sum_weight"), S$by_event$sum_weight[i])
}
dy <- r_dcdh$dyn
cs_by_e <- setNames(r_mod$cs_es$estimate, r_mod$cs_es$e)
for (i in seq_len(nrow(dy$effects))) {
  r <- dy$effects[i, ]; k <- sprintf("dcdh_dyn_l%d", r$l)
  put(paste0(k, "_est"), r$estimate); put(paste0(k, "_se"), r$se)
  put(paste0(k, "_lb"), r$lb); put(paste0(k, "_ub"), r$ub)
  put(paste0(k, "_n"), r$n); put(paste0(k, "_switchers"), r$switchers)
  put(paste0(k, "_truth"), truth$true_theta_e[event_time == r$e, theta])    # DID_l <-> e = l-1
  put(paste0(k, "_cs_est"), cs_by_e[[as.character(r$e)]])
}
for (i in seq_len(nrow(dy$placebos))) {
  r <- dy$placebos[i, ]; k <- sprintf("dcdh_pl_l%d", r$l)
  put(paste0(k, "_est"), r$estimate); put(paste0(k, "_se"), r$se)
  put(paste0(k, "_lb"), r$lb); put(paste0(k, "_ub"), r$ub)
  put(paste0(k, "_n"), r$n); put(paste0(k, "_switchers"), r$switchers)
}
put("dcdh_dyn_max_abs_diff_cs",                                                # max |DID_l - CS(e = l-1)|
    max(abs(dy$effects$estimate - cs_by_e[as.character(dy$effects$e)])))
put("dcdh_avg_total_est", dy$avg_total$estimate); put("dcdh_avg_total_se", dy$avg_total$se)
put("dcdh_avg_total_lb", dy$avg_total$lb); put("dcdh_avg_total_ub", dy$avg_total$ub)
put("dcdh_avg_total_n", dy$avg_total$n); put("dcdh_avg_total_switchers", dy$avg_total$switchers)
put("dcdh_delta_D_avg_total", dy$delta_D_avg_total)
put("dcdh_p_jointeffects", dy$p_jointeffects); put("dcdh_p_jointplacebo", dy$p_jointplacebo)

## ---- honest (08_honestdid.R) ----------------------------------------
put("honest_alpha", r_hon$meta$alpha)
put("honest_rm_theta_grid_step", r_hon$meta$rm_grid$step)
put("honest_rm_theta_grid_points", r_hon$meta$rm_grid$gridPoints[["e6"]])
put("honest_npre", r_hon$meta$npre); put("honest_npost", r_hon$meta$npost)
for (j in seq_along(r_hon$meta$Mvec))    put(sprintf("honest_sd_grid_%d", j), r_hon$meta$Mvec[j])
for (j in seq_along(r_hon$meta$Mbarvec)) put(sprintf("honest_rm_grid_%d", j), r_hon$meta$Mbarvec[j])
put_es("honest_es", r_hon$es)
## descriptive anchors for the scale of M and Mbar: first and second differences of the
## estimated pre-period coefficients (e = -6..-2) joined to the reference period e = -1 (= 0)
pre_path <- c(r_hon$es$estimate[r_hon$es$e < 0], 0)
put("honest_pre_est_max_abs_diff1", max(abs(diff(pre_path))))
put("honest_pre_est_max_abs_diff2", max(abs(diff(pre_path, differences = 2))))
for (nm in c("e6", "e12")) {
  tab <- r_hon[[nm]]
  put(sprintf("honest_%s_point", nm), tab$point[1])
  put(sprintf("honest_%s_rm_grid_limit", nm), tab$grid_limit[tab$method == "RM"][1])
  put(sprintf("honest_%s_orig_lb", nm), r_hon$orig[[nm]]$lb)
  put(sprintf("honest_%s_orig_ub", nm), r_hon$orig[[nm]]$ub)
  for (i in seq_len(nrow(tab))) {
    k <- sprintf("honest_%s_%s_M%s", nm, tolower(tab$method[i]), gridid(tab$M[i]))
    put(paste0(k, "_lb"), tab$lb[i]); put(paste0(k, "_ub"), tab$ub[i])
  }
  ## breakdown = smallest scanned M (Mbar) whose set contains 0; stored at the scan
  ## resolution (SD 1e-5, RM 0.005). _grid: first display-grid value containing 0.
  ## _fine: SD re-scan at 1e-6. _at_breakdown_*: the set at the breakdown (display settings).
  for (me in c("SD", "RM")) {
    k <- sprintf("honest_%s_%s", nm, tolower(me))
    put(paste0(k, "_breakdown"), r_hon$breakdown[[nm]][[me]])
    put(paste0(k, "_breakdown_resolution"), r_hon$breakdown_resolution[[nm]][[me]])
    put(paste0(k, "_breakdown_grid"), r_hon$breakdown_grid[[nm]][[me]])
    put(paste0(k, "_breakdown_fine"), r_hon$breakdown_fine[[nm]][[me]])
    put(paste0(k, "_breakdown_reversals"), r_hon$breakdown_reversals[[nm]][[me]])
    bs <- r_hon$breakdown_set[r_hon$breakdown_set$target == nm & r_hon$breakdown_set$method == me, ]
    put(paste0(k, "_breakdown_prev"), bs$prev_M)
    put(paste0(k, "_at_breakdown_lb"), bs$lb); put(paste0(k, "_at_breakdown_ub"), bs$ub)
  }
  jt <- r_hon$meta$sd_jitter[[nm]]
  put(sprintf("honest_%s_sd_bound_noise", nm), max(jt$max_lb_resid, jt$max_ub_resid))
}

## ---- family (09_family.R) -------------------------------------------
for (nm in c("did2s", "etwfe", "etwfe_never")) {
  r <- r_fam[[nm]]; k <- paste0("family_", nm)
  put(paste0(k, "_est"), r$est); put(paste0(k, "_se"), r$se)
  put(paste0(k, "_lb95"), r$lb); put(paste0(k, "_ub95"), r$ub)
}
put("family_did2s_minus_bjs", r_fam$compare$did2s_minus_bjs)
put("family_etwfe_minus_bjs", r_fam$compare$etwfe_minus_bjs)
put("family_etwfe_never_minus_sa", r_fam$compare$etwfe_never_minus_sa)

## ---- summary (this script): shares of truth and biases --------------
est_all <- c(twfe = r_twfe$static$est, cs = r_mod$cs_simple$est, sa = r_mod$sa_att$est,
             imp = r_mod$imp$est, did2s = r_fam$did2s$est, etwfe = r_fam$etwfe$est,
             etwfe_never = r_fam$etwfe_never$est)
for (nm in names(est_all)) {
  put(sprintf("summary_%s_share_of_truth", nm), est_all[[nm]] / true_att)
  put(sprintf("summary_%s_bias", nm), est_all[[nm]] - true_att)
}

stopifnot(length(nums) > 0, all(vapply(nums, function(z) is.numeric(z) && length(z) == 1 && is.finite(z), TRUE)))
write_json(nums, file.path(code_dir, "numbers.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)
cat(sprintf("\nSaved numbers.json (%d numbers)\n", length(nums)))
