# ============================================================================
# 07_summary.R  --  consolidate every quotable number for the chapter
# ============================================================================

d    <- readRDS("meridian_data.rds")
ols  <- readRDS("res_ols.rds")
iv   <- readRDS("res_iv.rds")
mte  <- readRDS("res_mte.rds")
wk   <- readRDS("res_weak.rds")
par  <- d$params; truth <- d$truth

line <- function() cat(strrep("-", 66), "\n")

cat("\n")
cat("================ MERIDIAN IV / LATE / MTE : SUMMARY ================\n")
cat(sprintf("seed = %d , n = %d\n", par$seed, par$n))
line()
cat("TRUE ESTIMANDS (oracle, from generated potential outcomes)\n")
cat(sprintf("  TRUE ATE            = %.3f\n", truth$ATE))
cat(sprintf("  TRUE ATT            = %.3f\n", truth$ATT))
cat(sprintf("  TRUE complier-LATE  = %.3f\n", truth$LATE))
cat(sprintf("  MTE(u) = %.3f - %.3f * u\n", par$t0, par$t1))
line()
cat("NAIVE ESTIMATORS (endogeneity bias)\n")
cat(sprintf("  OLS  b(D)           = %.3f  (SE %.3f)\n", ols$b_ols, ols$se_ols))
cat(sprintf("  Adjusted b(D|W)     = %.3f  (SE %.3f)\n", ols$b_adj, ols$se_adj))
cat(sprintf("  selection bias      = %.3f  (= E[Y0|D=1]-E[Y0|D=0])\n", ols$sel_bias))
line()
cat("INSTRUMENTAL VARIABLES (main design, strong instrument)\n")
cat(sprintf("  first stage         = %.3f  (SE %.3f)\n", iv$fs_coef, iv$fs_se))
cat(sprintf("  first-stage F (IID) = %.1f\n", iv$F_iid))
cat(sprintf("  effective F (robust)= %.1f\n", iv$eff_F))
cat(sprintf("  ITT (reduced form)  = %.3f  (SE %.3f)\n", iv$itt, iv$itt_se))
cat(sprintf("  Wald = ITT/1st      = %.3f  (SE %.3f)\n", iv$wald, iv$wald_se))
cat(sprintf("  2SLS (feols IV)     = %.3f  (SE %.3f)\n", iv$b_2sls, iv$se_2sls))
cat(sprintf("  --> 2SLS ~ TRUE complier-LATE (%.3f)\n", truth$LATE))
cat(sprintf("  complier share      = %.3f\n", iv$complier_share))
cat(sprintf("  E[W|complier] Abadie= %.3f  (oracle %.3f ; all firms %.3f)\n",
            iv$w_complier, iv$w_complier_oracle, iv$w_all))
line()
cat("MTE (continuous instrument)\n")
cat(sprintf("  local-IV MTE(u)     = %.3f + (%.3f) u   [slope SE %.3f]\n",
            mte$mte_a_int, mte$mte_a_slope, mte$mte_a_slope_se))
cat(sprintf("  control-fn t0,t1    = %.3f , %.3f\n", mte$cf_t0, mte$cf_t1))
cat(sprintf("  control-fn ATE      = %.3f  (bootstrap SE %.3f)\n",
            mte$cf_ate, mte$cf_ate_se))
cat("  ATE/ATT/LATE as weighted integrals of the MTE:\n")
cat(sprintf("    ATE  (uniform)          = %.3f  [true %.3f]\n", mte$ate_int, truth$ATE))
cat(sprintf("    ATT  (P(P>u)/E[P])      = %.3f  [oracle ATT(cont) %.3f]\n",
            mte$att_int, mte$att_oracle_cont))
cat(sprintf("    LATE (uniform on window)= %.3f  [2SLS %.3f, true %.3f]\n",
            mte$late_int, iv$b_2sls, truth$LATE))
line()
cat("WEAK-INSTRUMENT VARIANT\n")
cat(sprintf("  weak first stage    = %.3f  (complier share %.3f)\n",
            wk$fs_w_coef, par$complier_share_weak))
cat(sprintf("  first-stage F (IID) = %.2f   effective F = %.2f\n", wk$F_w_iid, wk$F_w_hc))
cat(sprintf("  weak 2SLS           = %.3f  (SE %.3f)\n", wk$b_w, wk$se_w))
cat(sprintf("  OLS (weak data)     = %.3f\n", wk$b_ols_w))
cat(sprintf("  conventional 95%% CI = [%.3f, %.3f]\n", wk$ci_w[1], wk$ci_w[2]))
cat(sprintf("  Anderson-Rubin CI   = [%.3f, %.3f]\n", wk$ar_ci[1], wk$ar_ci[2]))
cat(sprintf("  MC(1000): median 2SLS = %.3f , mean = %.3f\n", wk$mc_median, wk$mc_mean))
cat(sprintf("  MC conventional-CI coverage of LATE = %.3f  (nominal 0.95)\n", wk$mc_cov))
cat(sprintf("  MC median effective F = %.2f\n", wk$mc_medF))
cat("===================================================================\n")
