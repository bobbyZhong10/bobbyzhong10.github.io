# ------------------------------------------------------------------
# 02_md.R -- blocks G1 and G1h: classical minimum distance on the 16 market price tests.
#
# Reduced form: in market m, pi_hat_m = OLS slope of log weekly lessons booked on log price within the A/B test,
# sigma_hat_m its HC0 standard error (market_slopes). Structural restriction h(theta) = theta iota_r (r = 16,
# k = 1). V_pi = diag(sigma_hat_m^2) (the tests are independent), so the optimal MD weight solve(V_pi) is
# diagonal and optimal MD is the inverse-variance weighted mean (fixed-effect meta-analysis); its J statistic
# (Wooldridge 2010 eq. 14.69) is Cochran's Q, chi^2_15 under homogeneity.
#   G1   homogeneous elasticity: MD estimate, 95% interval, J; equal-weight MD (W = I_r) for comparison;
#        Monte Carlo (R_mc = 2000, seed G1_mc): size of the J test, coverage of the MD interval.
#   G1h  heterogeneous elasticities: J rejects; the MD estimate is near theta_star (precision-weighted mean of
#        theta_m), not the simple or household-weighted mean.
#   MD and GMM: two-step efficient GMM on the pooled household data with moments (market dummies, market
#        dummies x log price) x residual, parameters (a_1, ..., a_16, theta): 32 moments, 17 parameters, 15
#        overidentifying restrictions (the same count as MD). Asymptotically equivalent to optimal MD
#        (Wooldridge 2010 section 14.6.2); md_vs_gmm is the difference in the estimates of theta. With the GMM
#        weight built from the unrestricted (per-market) fits, GMM and MD are the same estimator (the
#        market-by-market criterion profiles out a_m and leaves (pi_hat_m - theta)^2 / sigma_hat_m^2):
#        md_vs_gmm_unrestricted, zero to rounding error.
#
# Output: res_md.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("gmm_fns.R")

tr <- readRDS("truth07.rds")
P1 <- tr$params$G1
G1 <- readRDS("tamarind_G1.rds")
G1h <- readRDS("tamarind_G1h.rds")
z <- qnorm(0.975)
r <- length(P1$markets)
h_hom <- function(th) rep(th, r)                     # h(theta) = theta iota_r

md_block <- function(d) {
  ms <- market_slopes(d)
  V <- diag(ms$se^2)                                  # V_pi: variance of pi_hat (independent tests)
  f <- md_fit(ms$pi_hat, V, h_hom, theta0 = 0)
  w <- 1 / ms$se^2
  se <- sqrt(f$vcov[1, 1])
  eq <- md_fit(ms$pi_hat, V, h_hom, theta0 = 0, W = diag(r))          # equal weights: the simple mean, sandwich SE
  se_eq <- sqrt(eq$vcov[1, 1])
  ms$lo <- ms$pi_hat - z * ms$se
  ms$hi <- ms$pi_hat + z * ms$se
  ms$weight <- w / sum(w)
  list(slopes = ms, theta = f$theta, se = se, ci = c(lo = f$theta - z * se, hi = f$theta + z * se),
       J = f$J, J_p = f$p, df = f$df, ivw_mean = sum(w * ms$pi_hat) / sum(w), converged = f$converged,
       equal = list(theta = eq$theta, se = se_eq, ci = c(lo = eq$theta - z * se_eq, hi = eq$theta + z * se_eq)))
}

## ---- G1: homogeneous elasticity -------------------------------------------------------------
md1 <- md_block(G1)
stopifnot(md1$converged)
md1$covers_truth <- as.numeric(md1$ci[["lo"]] <= tr$G1$theta && tr$G1$theta <= md1$ci[["hi"]])
md1$n_market_ci_cover <- sum(md1$slopes$lo <= tr$G1$theta & tr$G1$theta <= md1$slopes$hi)

## ---- G1: Monte Carlo of the 16 tests under homogeneity ---------------------------------------------
t0 <- proc.time()[["elapsed"]]
set.seed(tr$params$seed$G1_mc)
sim <- t(vapply(seq_len(P1$R_mc), function(i) {
  ms <- market_slopes(gen_G1(P1))
  f <- md_fit(ms$pi_hat, diag(ms$se^2), h_hom, theta0 = 0)
  c(theta = f$theta, se = sqrt(f$vcov[1, 1]), J = f$J, p = f$p, converged = f$converged)
}, numeric(5)))
mc_seconds <- proc.time()[["elapsed"]] - t0
stopifnot(all(sim[, "converged"] == 1))
md1$mc_R <- nrow(sim)
md1$mc_reject_J <- mean(sim[, "p"] < 0.05)
md1$mc_reject_J_01 <- mean(sim[, "p"] < 0.01)
md1$mc_reject_J_10 <- mean(sim[, "p"] < 0.10)
md1$mc_mean_J <- mean(sim[, "J"])
md1$mc_cover <- mean(abs(sim[, "theta"] - tr$G1$theta) <= z * sim[, "se"])
md1$mc_mean_theta <- mean(sim[, "theta"])
md1$mc_sd_theta <- sd(sim[, "theta"])
md1$mc_mean_se <- mean(sim[, "se"])

## ---- G1h: heterogeneous elasticities ---------------------------------------------------------------
md1h <- md_block(G1h)
stopifnot(md1h$converged)
th_m <- tr$G1h$theta_m
md1h$slopes$theta_m <- th_m
md1h$n_market_ci_cover <- sum(md1h$slopes$lo <= th_m & th_m <= md1h$slopes$hi)
cv <- function(ci, v) as.numeric(ci[["lo"]] <= v && v <= ci[["hi"]])
md1h$covers_star <- cv(md1h$ci, tr$G1h$theta_star)
md1h$covers_mean <- cv(md1h$ci, tr$G1h$theta_mean)
md1h$covers_hh <- cv(md1h$ci, tr$G1h$theta_hh)
md1h$equal_covers_mean <- cv(md1h$equal$ci, tr$G1h$theta_mean)

## ---- MD and efficient GMM on the pooled household data ---------------------------------------------
pdat <- pool_data(G1)
n <- length(pdat$y)
W1 <- solve(crossprod(pdat$Z) / n)                   # first step: pooled OLS with market intercepts and a common slope
g2 <- gmm_fit(g_pool, rep(0, r + 1), pdat, "two_step", W1 = W1)
stopifnot(g2$converged, g2$df == r - 1)
gmm2 <- list(theta = g2$theta[r + 1], se = sqrt(g2$vcov[r + 1, r + 1]), J = g2$J, J_p = g2$p, df = g2$df,
             theta1 = g2$theta1[r + 1], n = n, n_moments = ncol(pdat$Z), n_par = r + 1)
# weight from the unrestricted per-market fits (intercept and slope of each market's own regression)
ms <- md1$slopes
a_hat <- as.numeric(tapply(G1$y, G1$market, mean)) - ms$pi_hat * as.numeric(tapply(G1$log_price, G1$market, mean))
e_u <- pdat$y - a_hat[pdat$mk] - ms$pi_hat[pdat$mk] * pdat$lp
S_u <- moment_var(pdat$Z * e_u)
gu <- gmm_min(g_pool, g2$theta, pdat, solve(S_u))
stopifnot(gu$converged)
gmm_u <- list(theta = gu$theta[r + 1], se = sqrt(solve(crossprod(gu$D, solve(S_u, gu$D)))[r + 1, r + 1] / n),
              J = n * gu$obj)

res <- list(G1 = md1, G1h = md1h, gmm = list(two_step = gmm2, unrestricted = gmm_u),
            md_vs_gmm = md1$theta - gmm2$theta, md_vs_gmm_unrestricted = md1$theta - gmm_u$theta,
            J_md_vs_gmm_unrestricted = md1$J - gmm_u$J, se_md_vs_gmm_unrestricted = md1$se - gmm_u$se)
saveRDS(res, "res_md.rds")                           # run time is printed, not saved (byte-reproducible .rds)

cat(sprintf("G1: MD theta %.4f (SE %.4f, 95%% CI %.4f to %.4f; truth %.1f covered: %d); IVW mean %.4f; J %.3f (df %d, p %.4f)\n",
            md1$theta, md1$se, md1$ci[1], md1$ci[2], tr$G1$theta, md1$covers_truth, md1$ivw_mean, md1$J, md1$df, md1$J_p))
cat(sprintf("G1: equal-weight MD %.4f (SE %.4f); per-market intervals covering the truth: %d of %d\n",
            md1$equal$theta, md1$equal$se, md1$n_market_ci_cover, r))
cat(sprintf("G1 Monte Carlo (%d reps, %.1f s): J rejects at 5%% in %.4f; MD interval covers in %.4f; mean J %.3f\n",
            md1$mc_R, mc_seconds, md1$mc_reject_J, md1$mc_cover, md1$mc_mean_J))
cat(sprintf("G1h: MD theta %.4f (SE %.4f); theta_star %.4f, mean %.4f; J %.2f (p %.2e); covers theta_star %d, mean %d\n",
            md1h$theta, md1h$se, tr$G1h$theta_star, tr$G1h$theta_mean, md1h$J, md1h$J_p, md1h$covers_star, md1h$covers_mean))
cat(sprintf("Pooled GMM: first step %.4f, two-step %.4f (SE %.4f, J %.3f); MD - GMM %.2e; unrestricted-weight GMM %.10f (MD - GMM %.2e, J %.6f vs %.6f)\n",
            gmm2$theta1, gmm2$theta, gmm2$se, gmm2$J, res$md_vs_gmm, gmm_u$theta, res$md_vs_gmm_unrestricted, gmm_u$J, md1$J))
cat("Saved res_md.rds\n")
