# ------------------------------------------------------------------
# 04_twostep.R -- block M3: two-step M-estimation with a generated regressor.
#
# Step 1: OLS of one week of log visits v on z = (1, z1, z2): gamma_hat, traffic index t_hat = z'gamma_hat.
# Step 2: NLS of monthly sales y on exp(theta0 + theta1 t_hat) (gauss_newton).
# Variances of theta_hat (twostep_fit in m_fns.R):
#   naive      treats t_hat as data: A2^-1 B2 A2^-1 / n (Wooldridge 2010 eq. 12.47 with the second-step scores);
#   corrected  A2^-1 Var(s2 + F2 r1) A2^-1 / n, Wooldridge 2010 eq. (12.41) with D_o of eq. (12.40), estimated by
#              eqs. (12.61) for F2 and (12.62) for D (section 12.5.2); r1 is the OLS influence function (eq. 12.38).
#              Checked in tests/test_fns.R against the stacked (joint) M-estimator sandwich.
# Because theta1 != 0, F_o = E[grad_gamma s] != 0 (eq. 12.36): condition (12.37) fails and the first-step
# error adds variance (the scores of the two steps are uncorrelated here, so D_o - B_o is p.s.d.,
# Wooldridge 2010 p. 413).
# Also reported: NLS on the noisy v itself (errors in variables: attenuated) and on the true index t
# (infeasible benchmark).
# Monte Carlo: R_mc = 2000 replications at the design value theta1 = 0.8 (seed M3_mc), and R_grid = 1000 at
# each other theta1 on the grid (seed M3_grid + j for grid point j); theta1 = 0 is the case where (12.37)
# holds and the two variances agree. Coverage of the nominal 95% interval est -/+ qnorm(0.975) se.
#
# Output: res_twostep.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("m_fns.R")

tr <- readRDS("truth05.rds")
P3 <- tr$params$M3
M3 <- readRDS("larkspur_M3.rds")
z <- qnorm(0.975)

## ---- the sample ------------------------------------------------------------------------
f <- twostep_fit(M3)
stopifnot(f$converged)
est <- f$theta[2]
se_naive <- sqrt(f$V_naive[2, 2])
se_corrected <- sqrt(f$V_corrected[2, 2])
nls_on <- function(xvar) {
  X <- cbind(1, xvar)
  g <- gauss_newton(M3$y, X, exp_mean, exp_jac, c(0, 0.5))
  stopifnot(g$converged)
  p <- nls_parts(M3$y, X, g$theta)
  se <- sqrt(m_vcov(p$S, p$A)[2, 2])
  list(est = g$theta[2], se = se, ci = ci95(g$theta[2], se), theta = g$theta)
}
direct <- nls_on(M3$v)
oracle <- nls_on(M3$t_true)
r2 <- summary(lm(v ~ z1 + z2, data = M3))$r.squared

## ---- Monte Carlo --------------------------------------------------------------------------
one_rep <- function(par) {
  r <- twostep_fit(gen_M3(par))
  c(est = r$theta[2], se_naive = sqrt(r$V_naive[2, 2]), se_corrected = sqrt(r$V_corrected[2, 2]), converged = r$converged)
}
summarise_mc <- function(sim, theta1) {
  ok <- sim[, "converged"] == 1
  s <- sim[ok, , drop = FALSE]
  list(R = nrow(sim), n_converged = sum(ok), theta1 = theta1,
       cover_naive = mean(abs(s[, "est"] - theta1) <= z * s[, "se_naive"]),
       cover_corrected = mean(abs(s[, "est"] - theta1) <= z * s[, "se_corrected"]),
       mean_est = mean(s[, "est"]), sd_est = sd(s[, "est"]),
       mean_se_naive = mean(s[, "se_naive"]), mean_se_corrected = mean(s[, "se_corrected"]))
}
run_mc <- function(par, R, seed) {
  set.seed(seed)
  t(vapply(seq_len(R), function(r) one_rep(par), numeric(4)))
}
t0 <- proc.time()[["elapsed"]]
sim_main <- run_mc(P3, P3$R_mc, tr$params$seed$M3_mc)
mc <- summarise_mc(sim_main, P3$theta[2])
grid <- do.call(rbind, lapply(seq_along(P3$grid_theta1), function(j) {
  th1 <- P3$grid_theta1[j]
  if (isTRUE(all.equal(th1, P3$theta[2]))) {
    s <- mc                                                     # the design point: the main Monte Carlo
  } else {
    par <- P3
    par$theta[2] <- th1
    s <- summarise_mc(run_mc(par, P3$R_grid, tr$params$seed$M3_grid + j), th1)
  }
  as.data.frame(s)
}))
mc_seconds <- proc.time()[["elapsed"]] - t0
stopifnot(mc$n_converged == mc$R, all(grid$n_converged == grid$R))

ts <- list(est = est, theta = f$theta, gamma = f$gamma, se_naive = se_naive, se_corrected = se_corrected,
           ci_naive = ci95(est, se_naive), ci_corrected = ci95(est, se_corrected), se_ratio = se_corrected / se_naive,
           F2 = f$F2, r2_first_stage = r2, iter = f$iter,
           covers_naive = as.numeric(abs(est - P3$theta[2]) <= z * se_naive),
           covers_corrected = as.numeric(abs(est - P3$theta[2]) <= z * se_corrected),
           direct = direct, oracle = oracle, mc = mc, grid = grid)   # run time is printed, not saved (byte-reproducible .rds)
saveRDS(ts, "res_twostep.rds")

cat(sprintf("M3: theta1_hat = %.4f; naive SE %.4f, corrected SE %.4f (ratio %.3f); direct NLS on v %.4f; on true t %.4f\n",
            est, se_naive, se_corrected, ts$se_ratio, direct$est, oracle$est))
cat(sprintf("M3 Monte Carlo (%d reps): coverage naive %.4f, corrected %.4f; sd of estimates %.4f; mean SE %.4f / %.4f (%.1f s incl. grid)\n",
            mc$R, mc$cover_naive, mc$cover_corrected, mc$sd_est, mc$mean_se_naive, mc$mean_se_corrected, mc_seconds))
print(grid[, c("theta1", "R", "cover_naive", "cover_corrected", "sd_est", "mean_se_naive", "mean_se_corrected")], row.names = FALSE)
cat("Saved res_twostep.rds\n")
