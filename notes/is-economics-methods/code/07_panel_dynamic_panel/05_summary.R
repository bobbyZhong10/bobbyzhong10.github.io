# ============================================================================
# 05_summary.R  --  master reconciliation of every number chapter 7 cites.
# ============================================================================
st <- readRDS("res_static.rds"); dy <- readRDS("res_dynamic.rds")
L <- function(...) cat(sprintf(...), "\n")
cat("================ CHAPTER 7 NUMBER RECONCILIATION ================\n\n")

cat("--- static panel (beta, truth 0.50) ---\n")
L("pooled OLS = %.4f (SE %.4f)", st$pool, st$se_pool)
L("random effects = %.4f", st$re)
L("fixed effects = %.4f (SE %.4f; cluster-by-merchant %.4f)", st$fe, st$se_fe, st$se_fe_cl)
L("Hausman chi2 = %.1f  p = %.3g", st$haus_stat, st$haus_p)

cat("\n--- dynamic panel (rho, truth 0.70) ---\n")
L("pooled OLS = %.4f (up)", dy$pool)
L("within/FE  = %.4f (Nickell, down)", dy$fe)
L("Arellano-Bond (diff GMM)  = %.4f", dy$ab)
L("Blundell-Bond (sys GMM)   = %.4f", dy$bb)
L("AB tests: Sargan p = %.3f | AR(1) p = %.4f | AR(2) p = %.3f", dy$sargan_p, dy$ar1_p, dy$ar2_p)

cat("\n--- Nickell bias vs T (FE / AB / -(1+rho)/(T-1) approx; truth 0.70) ---\n")
print(round(dy$nick, 4))

cat("\n--- rho near 1 (truth 0.90) ---\n")
L("difference GMM = %.4f (SE %.4f)", dy$ab9, dy$se_ab9)
L("system GMM     = %.4f (SE %.4f)", dy$bb9, dy$se_bb9)
cat("\n================================================================\n")
