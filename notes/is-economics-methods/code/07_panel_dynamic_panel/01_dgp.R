# ============================================================================
# 01_dgp.R  --  Nimbus marketplace merchant panel (two facets, one seed)
# ----------------------------------------------------------------------------
# "Nimbus" is a fictional marketplace. We observe merchants i over months t.
#
#  (A) STATIC facet -- FE vs RE and the Hausman test.
#      y_it = beta * x_it + alpha_i + eps_it,   with Corr(alpha_i, x_it) > 0
#      (better merchants, high alpha_i, also run more promotions x_it). Pooled
#      OLS and random effects are biased for beta because they treat alpha_i as
#      uncorrelated with x; fixed effects (within) is consistent; Hausman
#      rejects RE. This teaches FE, strict exogeneity, and the Hausman test.
#
#  (B) DYNAMIC facet -- Nickell bias and Arellano-Bond / Blundell-Bond.
#      y_it = rho * y_i,t-1 + alpha_i + eps_it   (persistent GMV).
#      Pooled OLS of y on its lag is biased UP (alpha_i sits in both y and its
#      lag); the within/FE estimator is biased DOWN by the Nickell O(1/T) bias
#      (demeaning correlates the transformed lag with the transformed error);
#      the truth is bracketed between them, and difference/system GMM recover it.
#
# All TRUE parameters known. Seed fixed. Scripts use relative paths; run with
# the working directory set to this folder.
# ============================================================================

set.seed(737)

## ---- (A) static panel: FE vs RE --------------------------------------------
gen_static <- function(N = 600, T = 6, beta = 0.5, ca = 0.7) {
  alpha <- rnorm(N, 0, 1)                              # merchant quality (unobs.)
  id <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  a  <- alpha[id]
  x  <- ca * a + rnorm(N*T, 0, 1)                      # promotions, corr with alpha
  y  <- beta * x + a + rnorm(N*T, 0, 1)
  data.frame(id, t = tt, y, x)
}
static <- gen_static()

## ---- (B) dynamic panel: AR(1) with merchant effects ------------------------
# stationary start; burn-in to erase initial-condition transient
gen_dyn <- function(N, T, rho, burn = 50, sigma_a = 1, sigma_e = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  alpha <- rnorm(N, 0, sigma_a)
  Ttot  <- T + burn
  y <- matrix(0, N, Ttot)
  y[,1] <- alpha/(1-rho) + rnorm(N, 0, sigma_e/sqrt(1-rho^2))  # stationary draw
  for (s in 2:Ttot) y[,s] <- rho*y[,s-1] + alpha + rnorm(N, 0, sigma_e)
  y <- y[, (burn+1):Ttot]                              # keep last T
  data.frame(id = rep(1:N, each = T), t = rep(1:T, times = N),
             y = as.vector(t(y)))
}
dyn_main <- gen_dyn(N = 500, T = 8, rho = 0.70, seed = 808)    # main
dyn_persist <- gen_dyn(N = 500, T = 8, rho = 0.90, seed = 7372) # rho near 1

saveRDS(list(static = static, dyn_main = dyn_main, dyn_persist = dyn_persist,
             gen_dyn = gen_dyn,
             truth = list(beta = 0.5, ca = 0.7, rho = 0.70, rho_persist = 0.90)),
        "nimbus.rds")

cat(sprintf("static: N=600 T=6, beta=0.5, Corr(alpha,x) built via ca=0.7\n"))
cat(sprintf("  implied pooled-OLS beta plim = beta + Cov(x,a)/Var(x) = %.3f\n",
            0.5 + 0.7/(0.7^2*1 + 1)))
cat(sprintf("dynamic: N=500 T=8, rho=0.70 (main), rho=0.90 (persistent)\n"))
cat(sprintf("  Nickell leading bias -(1+rho)/(T-1) at rho=.7,T=8 = %.3f\n",
            -(1+0.7)/(8-1)))
