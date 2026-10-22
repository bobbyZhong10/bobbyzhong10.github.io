# ============================================================================
# 02_lpm_naive.R  --  the "obviously wrong number": a linear probability model
# on a binary outcome. OLS of D on x gives a single constant slope and, worse,
# predicted probabilities that fall outside [0,1] for a nontrivial share of
# firms. This motivates a proper likelihood (logit MLE).
# ============================================================================

dat <- readRDS("helix_data.rds")

lpm <- lm(D ~ x, data = dat)
b   <- coef(lpm)
phat <- fitted(lpm)

cat(sprintf("LPM intercept     = %.4f\n", b[1]))
cat(sprintf("LPM slope (x)     = %.4f (SE %.4f)\n",
            b[2], summary(lpm)$coefficients[2, 2]))
cat(sprintf("min fitted phat   = %.4f\n", min(phat)))
cat(sprintf("max fitted phat   = %.4f\n", max(phat)))
cat(sprintf("share phat<0 or >1= %.4f\n", mean(phat < 0 | phat > 1)))
cat(sprintf("  share phat > 1   = %.4f\n", mean(phat > 1)))
cat(sprintf("  share phat < 0   = %.4f\n", mean(phat < 0)))
