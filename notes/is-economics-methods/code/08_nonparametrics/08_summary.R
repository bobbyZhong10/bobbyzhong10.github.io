# ============================================================================
# 08_summary.R  --  master reconciliation of every number chapter 8 cites.
# ============================================================================
S <- readRDS("solstice.rds"); truth <- S$truth
rl <- readRDS("res_linear.rds"); rk <- readRDS("res_kernel.rds")
rc <- readRDS("res_cv.rds"); rr <- readRDS("res_robinson.rds"); rcu <- readRDS("res_curse.rds")
L <- function(...) cat(sprintf(...), "\n")
cat("================ CHAPTER 8 NUMBER RECONCILIATION ================\n\n")

cat("--- the demand curve (truth: local slope/elasticity varies) ---\n")
L("true local slope g'(x) at price 2 / 5.5 / 8 = %.3f / %.3f / %.3f",
  truth$slope[1], truth$slope[2], truth$slope[3])
L("true point elasticity x*g'(x) at 2 / 5.5 / 8 = %.3f / %.3f / %.3f",
  truth$elas[1], truth$elas[2], truth$elas[3])
L("constant-elasticity fit (single number) = %.3f ; linear slope = %.3f",
  rl$elas_hat, rl$slope_lin)
L("linear R^2 = %.3f ; RESET p = %.2e (reject linear)", rl$r2, rl$reset_p)

cat("\n--- Nadaraya-Watson vs local linear ---\n")
L("single-fit MSE vs truth (h=%.1f): NW = %.5f  local-linear = %.5f", rk$h, rk$mse_nw, rk$mse_ll)
L("local-linear slope at 2/5.5/8 = %.3f / %.3f / %.3f",
  rk$ll_slopes[1], rk$ll_slopes[2], rk$ll_slopes[3])
cat("MC bias (NW row / LL row) at left / interior / right:\n"); print(round(rk$bias_tbl, 4))

cat("\n--- bias-variance and cross-validation ---\n")
cat("bias^2 / variance / MSE vs bandwidth h:\n"); print(round(rc$bv, 5))
L("MSE-optimal h on grid = %.2f ; least-squares CV bandwidth = %.3f", rc$h_opt, rc$cv_bw)
L("smoothing-spline GCV df = %.1f ; MSE oversmooth(df3)=%.5f GCV=%.5f overfit(df30)=%.5f",
  rc$gcv_df, rc$mse_over, rc$mse_gcv, rc$mse_under)

cat("\n--- partially linear (Robinson), truth beta = 0.80 ---\n")
L("y~z = %.4f ; y~z+x = %.4f ; y~z+poly3 = %.4f ; Robinson = %.4f (SE %.4f)",
  rr$b_naive1, rr$b_naive2, rr$b_naive3, rr$b_rob, rr$se_rob)
L("MC: linear-in-x = %.4f ; Robinson = %.4f ; Robinson SD = %.4f",
  rr$mc_lin, rr$mc_rob, rr$mc_sd)

cat("\n--- curse of dimensionality ---\n")
L("edge length e_d(0.1) for d=1/2/5/10 = %.3f / %.3f / %.3f / %.3f",
  rcu$edge[1], rcu$edge[2], rcu$edge[3], rcu$edge[4])
L("NW MSE at x=0 vs d=1/2/3/5 (n=2000) = %.4f / %.4f / %.4f / %.4f",
  rcu$mse_dim[1], rcu$mse_dim[2], rcu$mse_dim[3], rcu$mse_dim[4])
cat("\n================================================================\n")
