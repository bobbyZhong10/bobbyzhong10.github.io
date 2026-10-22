# ------------------------------------------------------------------
# 04_weak.R -- block G3: weak moments (penny-sized price perturbations), Wald against the S confidence set.
#
# Moments g_conc (gmm_fns.R): g_i(theta) = z_i (y~_i - theta p~_i), the K = 4 perturbation moments with the
# intercept partialled out (y~, p~ demeaned; each z_j has mean exactly 0): r = 4, k = 1.
#   wald             two-step efficient GMM (first step 2SLS), conventional SE: est -/+ qnorm(0.975) SE
#   wald_windmeijer  the same estimate with the Windmeijer (2005) SE
#   wald_2sls        2SLS with the robust sandwich SE
#   wald_cue         CUE with the conventional SE (over the replications in which CUE converged; mc$n_cue_converged)
#   s_set            the Stock and Wright (2000) set {theta : S_n(theta) <= chi^2_{4, 0.95}}; it contains theta_0
#                    exactly when S_n(theta_0) <= chi^2_{4, 0.95}, so its Monte Carlo coverage is computed from
#                    S_n(theta_0) itself (no grid). The set is unbounded when the first-stage S statistic
#                    (the limit of S_n(theta) as |theta| -> infinity, moments z p~) is below the critical value.
# Sample: estimates, first-stage F (HC0 robust Wald / K), the S curve on a grid (figure 3) and the shape of
# the S set. Monte Carlo: R_mc = 2000 (seed G3_mc).
#
# Output: res_weak.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("gmm_fns.R")

tr <- readRDS("truth07.rds")
P3 <- tr$params$G3
G3 <- readRDS("tamarind_G3.rds")
z <- qnorm(0.975)
th0 <- tr$G3$theta
K <- length(P3$step)
crit <- qchisq(0.95, K)

## first-stage S statistic: moments z_i p~_i (the S statistic of theta as |theta| -> infinity)
s_inf <- function(dl) {
  G <- dl$Zs * dl$pc
  nrow(G) * drop(crossprod(colMeans(G), solve(crossprod(G) / nrow(G), colMeans(G))))
}
first_stage_F <- function(dl) {
  ZZi <- solve(crossprod(dl$Zs))
  b <- ZZi %*% crossprod(dl$Zs, dl$pc)
  e <- drop(dl$pc - dl$Zs %*% b)
  V <- ZZi %*% crossprod(dl$Zs * e) %*% ZZi
  drop(t(b) %*% solve(V, b)) / ncol(dl$Zs)
}
fit_weak <- function(dl) {
  n <- length(dl$y)
  W1 <- solve(crossprod(dl$Zs) / n)
  f2 <- gmm_fit(g_conc, 0, dl, "two_step", W1 = W1)
  wv <- windmeijer_vcov(g_conc, f2$theta1, f2$theta, dl, W1)
  # CUE under weak identification: its criterion can keep falling as |theta| grows (towards the first-stage S
  # statistic); such a run (optim not converged, or an error) is recorded as not converged and left out of
  # the CUE coverage only
  fc <- tryCatch(gmm_fit(g_conc, f2$theta1, dl, "cue", W1 = W1), error = function(e) NULL)
  cue_ok <- !is.null(fc) && fc$converged
  c(est_two_step = f2$theta, se_two_step = sqrt(f2$vcov[1, 1]), se_windmeijer = sqrt(wv[1, 1]),
    est_2sls = f2$theta1, se_2sls = sqrt(attr(wv, "V1")[1, 1]),
    est_cue = if (cue_ok) fc$theta else NA_real_, se_cue = if (cue_ok) sqrt(fc$vcov[1, 1]) else NA_real_,
    S0 = s_stat(g_conc, th0, dl), S_inf = s_inf(dl), F = first_stage_F(dl), J = f2$J,
    converged = as.numeric(f2$converged), cue_converged = as.numeric(cue_ok))
}

## ---- the sample -------------------------------------------------------------------------------
dl <- iv_data(G3)
v <- fit_weak(dl)
stopifnot(v[["converged"]] == 1, v[["cue_converged"]] == 1)
wald_ci <- function(est, se) c(lo = est - z * se, hi = est + z * se)
grid <- seq(-12, 8, by = 0.01)
ss <- s_set(g_conc, grid, dl)
runs <- rle(ss$in_set)
set_lo <- if (ss$in_set[1]) -Inf else min(ss$theta[ss$in_set])
set_hi <- if (ss$in_set[nrow(ss)]) Inf else max(ss$theta[ss$in_set])
# the bounds read off the grid ends agree with the S-at-infinity flag: both infinite when S_inf <= crit (the set
# contains the tails), both finite otherwise (as 03_gmm_mc.R asserts for G2)
stopifnot(is.infinite(set_lo) == (v[["S_inf"]] <= crit), is.infinite(set_hi) == (v[["S_inf"]] <= crit))
sample <- list(
  two_step = list(est = v[["est_two_step"]], se = v[["se_two_step"]], ci = wald_ci(v[["est_two_step"]], v[["se_two_step"]])),
  windmeijer = list(est = v[["est_two_step"]], se = v[["se_windmeijer"]], ci = wald_ci(v[["est_two_step"]], v[["se_windmeijer"]])),
  tsls = list(est = v[["est_2sls"]], se = v[["se_2sls"]], ci = wald_ci(v[["est_2sls"]], v[["se_2sls"]])),
  cue = list(est = v[["est_cue"]], se = v[["se_cue"]], ci = wald_ci(v[["est_cue"]], v[["se_cue"]])),
  ols = list(est = unname(coef(lm(y ~ log_price, data = G3))[2])),
  S0 = v[["S0"]], S0_p = pchisq(v[["S0"]], K, lower.tail = FALSE), S_inf = v[["S_inf"]],
  unbounded = as.numeric(v[["S_inf"]] <= crit), first_stage_F = v[["F"]], J = v[["J"]], crit = crit, df = K,
  s_curve = ss, s_set = c(lo = set_lo, hi = set_hi), s_set_pieces = sum(runs$values),
  s_min = min(ss$S), s_argmin = ss$theta[which.min(ss$S)], grid = range(grid))
sample$wald_covers <- as.numeric(abs(sample$two_step$est - th0) <= z * sample$two_step$se)
sample$s_covers <- as.numeric(sample$S0 <= crit)

## ---- Monte Carlo ------------------------------------------------------------------------------------
t0 <- proc.time()[["elapsed"]]
set.seed(tr$params$seed$G3_mc)
sim <- t(vapply(seq_len(P3$R_mc), function(i) fit_weak(iv_data(gen_G3(P3))), numeric(length(v))))
mc_seconds <- proc.time()[["elapsed"]] - t0
ok <- sim[, "converged"] == 1
s <- sim[ok, , drop = FALSE]
cov_w <- function(est, se) mean(abs(s[, est] - th0) <= z * s[, se], na.rm = TRUE)
cover <- list(wald = cov_w("est_two_step", "se_two_step"), wald_windmeijer = cov_w("est_two_step", "se_windmeijer"),
              wald_2sls = cov_w("est_2sls", "se_2sls"), wald_cue = cov_w("est_cue", "se_cue"),
              s_set = mean(s[, "S0"] <= crit))
mc <- list(R = nrow(sim), n_converged = sum(ok), n_cue_converged = sum(s[, "cue_converged"] == 1),
           share_unbounded = mean(s[, "S_inf"] <= crit),
           median_F = median(s[, "F"]), mean_F = mean(s[, "F"]), share_F_below_10 = mean(s[, "F"] < 10),
           median_two_step = median(s[, "est_two_step"]), median_2sls = median(s[, "est_2sls"]), median_cue = median(s[, "est_cue"], na.rm = TRUE),
           mean_se_two_step = mean(s[, "se_two_step"]), iqr_two_step = unname(diff(quantile(s[, "est_two_step"], c(0.25, 0.75)))),
           reject_S_at_truth = mean(s[, "S0"] > crit), reject_J = mean(s[, "J"] > qchisq(0.95, K - 1)))
stopifnot(mc$n_converged == mc$R)

res <- list(sample = sample, cover = cover, mc = mc)
saveRDS(res, "res_weak.rds")                         # run time is printed, not saved (byte-reproducible .rds)

cat(sprintf("G3 sample (n = %d): first-stage F %.2f; two-step %.4f (SE %.4f, CI %.3f to %.3f); 2SLS %.4f; CUE %.4f; OLS %.4f\n",
            length(dl$y), sample$first_stage_F, sample$two_step$est, sample$two_step$se, sample$two_step$ci[1], sample$two_step$ci[2],
            sample$tsls$est, sample$cue$est, sample$ols$est))
cat(sprintf("G3 sample: S(theta_0) = %.3f (p %.3f); S set on [%g, %g]: %g to %g in %d piece(s); S at infinity %.3f (crit %.3f, unbounded %d)\n",
            sample$S0, sample$S0_p, min(grid), max(grid), set_lo, set_hi, sample$s_set_pieces, sample$S_inf, crit, sample$unbounded))
cat(sprintf("G3 Monte Carlo (%d reps, %.1f s): coverage %s; unbounded S sets %.4f; median F %.2f\n", mc$R, mc_seconds,
            paste(sprintf("%s %.4f", names(cover), unlist(cover)), collapse = ", "), mc$share_unbounded, mc$median_F))
cat("Saved res_weak.rds\n")
