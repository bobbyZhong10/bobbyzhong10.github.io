# ============================================================================
# 02_ols_naive.R  --  the endogeneity problem, in numbers
# ----------------------------------------------------------------------------
# Regress observed engagement Y on adoption D, first raw, then adjusting for the
# observed covariate W (firm size).  Because adoption is driven by unobserved
# sophistication (low U_D), and sophistication ALSO lifts baseline Y0 and the
# gain, OLS overstates the causal effect.  Covariate adjustment on W removes
# only the part of sophistication that W proxies; residual selection remains.
# Benchmark against the TRUE ATT (what a "clean" adopter-vs-nonadopter contrast
# would target if adoption were random).
# ============================================================================

d <- readRDS("meridian_data.rds")
dat <- d$dat; truth <- d$truth

## naive OLS: Y ~ D
m_ols  <- lm(Y ~ D, data = dat)

## covariate-adjusted OLS: Y ~ D + W
m_adj  <- lm(Y ~ D + W_size, data = dat)

b_ols <- coef(m_ols)["D"]
b_adj <- coef(m_adj)["D"]
se_ols <- summary(m_ols)$coefficients["D", "Std. Error"]
se_adj <- summary(m_adj)$coefficients["D", "Std. Error"]

## the exact selection-bias decomposition (oracle):
##   E[Y|D=1]-E[Y|D=0] = ATT + { E[Y0|D=1] - E[Y0|D=0] }
sel_bias <- mean(dat$Y0[dat$D == 1]) - mean(dat$Y0[dat$D == 0])

cat("=== OLS / covariate-adjusted OLS ===\n")
cat(sprintf("Naive OLS  b(D)          = %.4f  (SE %.4f)\n", b_ols, se_ols))
cat(sprintf("Adjusted   b(D | W)      = %.4f  (SE %.4f)\n", b_adj, se_adj))
cat(sprintf("TRUE ATT                 = %.4f\n", truth$ATT))
cat(sprintf("TRUE ATE                 = %.4f\n", truth$ATE))
cat(sprintf("TRUE complier-LATE       = %.4f\n", truth$LATE))
cat("\n--- selection-bias check (oracle) ---\n")
cat(sprintf("E[Y0|D=1]-E[Y0|D=0]      = %.4f  (selection bias)\n", sel_bias))
cat(sprintf("ATT + selection bias     = %.4f  (should match naive OLS)\n",
            truth$ATT + sel_bias))
cat(sprintf("Naive OLS overstates ATT by %.4f; adjusted still overstates by %.4f\n",
            b_ols - truth$ATT, b_adj - truth$ATT))

saveRDS(list(b_ols = b_ols, se_ols = se_ols, b_adj = b_adj, se_adj = se_adj,
             sel_bias = sel_bias),
        "res_ols.rds")
