# ------------------------------------------------------------------
# 03_gmm_mc.R -- block G2: overidentified linear demand moments, experimental price shocks as instruments.
#
# Moments g_i(theta) = z_i (y_i - beta0 - theta log p_i), z_i = (1, z_i1, ..., z_i10): r = 11, k = 2, 9
# overidentifying restrictions (linear in the parameters: the case Windmeijer 2005 proves; ruling 16).
# Estimators of theta (the elasticity), each with its 95% Wald interval est -/+ qnorm(0.975) SE:
#   one_step        2SLS, W^(1) = (Z'Z / n)^-1, robust sandwich SE (gmm_sandwich)
#   identity        W = I_r with the experiments in their recorded units (USD, pp), robust sandwich SE
#   two_step_naive  two-step efficient GMM, weight S(theta^(1))^-1, conventional SE (D'WD)^-1 / n
#   windmeijer      the same estimate, Windmeijer (2005) corrected SE (windmeijer_vcov)
#   iterated        iterated GMM (Hansen, Heaton and Yaron 1996), conventional SE
#   cue             continuously updated GMM (Hansen, Heaton and Yaron 1996), conventional SE
# Also for the sample: OLS (biased by the discretionary promotions), the two-step J test, the AGS (2017)
# sensitivity Lambda-hat = -(D'WD)^-1 D'W of theta-hat to each moment, and the S statistic over a grid of
# theta (concentrated moments g_conc, intercept partialled out, 10 degrees of freedom; for figure 3).
# Monte Carlo: R_mc = 2000 replications at n = 200 (seed G2_mc); R_grid = 1000 at each other n on grid_n
# (seed G2_grid + j for grid point j).
#
# Output: res_gmm_mc.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("gmm_fns.R")

tr <- readRDS("truth07.rds")
P2 <- tr$params$G2
G2 <- readRDS("tamarind_G2.rds")
z <- qnorm(0.975)
th0 <- tr$G2$theta
est_names <- c("one_step", "identity", "two_step_naive", "windmeijer", "iterated", "cue")

## All six estimators on one data set (list from iv_data). Returns the fits and a named vector of estimates and SEs of theta.
fit_all <- function(dl) {
  n <- length(dl$y)
  r <- ncol(dl$Z)
  W1 <- solve(crossprod(dl$Z) / n)
  f2 <- gmm_fit(g_iv, c(0, 0), dl, "two_step", W1 = W1)
  wv <- windmeijer_vcov(g_iv, f2$theta1, f2$theta, dl, W1)
  fi <- gmm_fit(g_iv, f2$theta1, dl, "iterated", W1 = W1)
  fc <- gmm_fit(g_iv, f2$theta1, dl, "cue", W1 = W1)
  fid <- gmm_min(g_iv, f2$theta1, dl, diag(r))
  V_id <- gmm_sandwich(fid$D, diag(r), moment_var(g_iv(fid$theta, dl)), n)
  V1 <- attr(wv, "V1")
  v <- c(est_one_step = f2$theta1[2], se_one_step = sqrt(V1[2, 2]),
         est_identity = fid$theta[2], se_identity = sqrt(V_id[2, 2]),
         est_two_step_naive = f2$theta[2], se_two_step_naive = sqrt(f2$vcov[2, 2]),
         est_windmeijer = f2$theta[2], se_windmeijer = sqrt(wv[2, 2]),
         est_iterated = fi$theta[2], se_iterated = sqrt(fi$vcov[2, 2]),
         est_cue = fc$theta[2], se_cue = sqrt(fc$vcov[2, 2]),
         J_two_step = f2$J, J_iterated = fi$J, J_cue = fc$J, iter_iterated = fi$iter,
         converged = as.numeric(f2$converged && fi$converged && fc$converged && fid$converged))
  list(v = v, f2 = f2, wv = wv, fi = fi, fc = fc, fid = fid, V_id = V_id, W1 = W1)
}

## ---- the sample -------------------------------------------------------------------------------
dl <- iv_data(G2)
n <- length(dl$y)
fa <- fit_all(dl)
stopifnot(fa$v[["converged"]] == 1)
sample_est <- lapply(setNames(est_names, est_names), function(k) {
  est <- fa$v[[paste0("est_", k)]]
  se <- fa$v[[paste0("se_", k)]]
  list(est = est, se = se, ci = c(lo = est - z * se, hi = est + z * se), covers = as.numeric(abs(est - th0) <= z * se))
})
ols <- lm(y ~ log_price, data = G2)
se_ols <- sqrt(sandwich::vcovHC(ols, type = "HC0")[2, 2])
f2 <- fa$f2
Lambda <- -solve(crossprod(f2$D, f2$W %*% f2$D), crossprod(f2$D, f2$W))          # k x r sensitivity (AGS 2017)
fs <- lm(dl$X[, 2] ~ dl$Zs)
fs_b <- coef(fs)[-1]
fs_V <- sandwich::vcovHC(fs, type = "HC0")[-1, -1]
first_stage_F <- drop(t(fs_b) %*% solve(fs_V, fs_b)) / ncol(dl$Zs)               # robust (HC0) Wald / K
grid_s <- seq(-8, 2, by = 0.005)
ss <- s_set(g_conc, grid_s, dl)
stopifnot(!ss$in_set[1], !ss$in_set[nrow(ss)], sum(diff(ss$in_set) != 0) == 2)  # a bounded interval on the grid
sample <- list(est = sample_est, beta0_two_step = f2$theta[1], theta1 = f2$theta1, J = f2$J, J_p = f2$p, J_df = f2$df,
               J_iterated = fa$fi$J, J_cue = fa$fc$J, iter_iterated = fa$fi$iter,
               ols = list(est = unname(coef(ols)[2]), se = se_ols), first_stage_F = first_stage_F,
               first_stage_r2 = summary(fs)$r.squared,
               lambda_theta = Lambda[2, ], lambda_theta_std = Lambda[2, ] * sqrt(diag(f2$S)),
               windmeijer_D = attr(fa$wv, "dtheta2_dtheta1"),
               se_ratio_windmeijer = fa$v[["se_windmeijer"]] / fa$v[["se_two_step_naive"]],
               s_curve = ss, s_set = c(lo = min(ss$theta[ss$in_set]), hi = max(ss$theta[ss$in_set])),
               s_min = min(ss$S), s_argmin = ss$theta[which.min(ss$S)], s_df = ncol(dl$Zs))

## ---- Monte Carlo ------------------------------------------------------------------------------------
one_rep <- function(nn) fit_all(iv_data(gen_G2(P2, n = nn)))$v
run_mc <- function(nn, R, seed) {
  set.seed(seed)
  t(vapply(seq_len(R), function(i) one_rep(nn), numeric(length(fa$v))))
}
summarise_mc <- function(sim, nn) {
  ok <- sim[, "converged"] == 1
  s <- sim[ok, , drop = FALSE]
  out <- list(n = nn, R = nrow(sim), n_converged = sum(ok))
  for (k in est_names) {
    e <- s[, paste0("est_", k)]
    se <- s[, paste0("se_", k)]
    out[[paste0("cover_", k)]] <- mean(abs(e - th0) <= z * se)
    out[[paste0("mean_", k)]] <- mean(e)
    out[[paste0("median_", k)]] <- median(e)
    out[[paste0("sd_", k)]] <- sd(e)
    out[[paste0("mean_se_", k)]] <- mean(se)
  }
  df <- length(P2$step) - 1
  out$reject_J_two_step <- mean(s[, "J_two_step"] > qchisq(0.95, df))
  out$reject_J_iterated <- mean(s[, "J_iterated"] > qchisq(0.95, df))
  out$reject_J_cue <- mean(s[, "J_cue"] > qchisq(0.95, df))
  out$mean_iter_iterated <- mean(s[, "iter_iterated"])
  out
}
t0 <- proc.time()[["elapsed"]]
mc <- summarise_mc(run_mc(P2$n, P2$R_mc, tr$params$seed$G2_mc), P2$n)
grid <- do.call(rbind, lapply(seq_along(P2$grid_n), function(j) {
  nn <- P2$grid_n[j]
  s <- if (nn == P2$n) mc else summarise_mc(run_mc(nn, P2$R_grid, tr$params$seed$G2_grid + j), nn)
  as.data.frame(s)
}))
mc_seconds <- proc.time()[["elapsed"]] - t0
stopifnot(mc$n_converged == mc$R, all(grid$n_converged == grid$R))

res <- list(sample = sample, mc = mc, grid = grid,
            cover = lapply(setNames(est_names, est_names), function(k) mc[[paste0("cover_", k)]]),
            se_eff = sample$est$two_step_naive$se, se_identity = sample$est$identity$se, se_one_step = sample$est$one_step$se)
saveRDS(res, "res_gmm_mc.rds")                       # run time is printed, not saved (byte-reproducible .rds)

cat(sprintf("G2 sample (n = %d): OLS %.4f; 2SLS %.4f (SE %.4f); identity %.4f (SE %.4f); two-step %.4f (SE %.4f, Windmeijer %.4f); iterated %.4f; CUE %.4f; J %.3f (p %.3f)\n",
            n, sample$ols$est, sample_est$one_step$est, sample_est$one_step$se, sample_est$identity$est, sample_est$identity$se,
            sample_est$two_step_naive$est, sample_est$two_step_naive$se, sample_est$windmeijer$se,
            sample_est$iterated$est, sample_est$cue$est, sample$J, sample$J_p))
cat(sprintf("G2 sample: intervals covering theta = %.1f: %s; S set %.3f to %.3f; first-stage F %.1f\n", th0,
            paste(sprintf("%s %d", est_names, vapply(sample_est, `[[`, 1, "covers")), collapse = ", "),
            sample$s_set[1], sample$s_set[2], first_stage_F))
cat(sprintf("G2 Monte Carlo (%d reps at n = %d, %.1f s incl. grid): coverage %s\n", mc$R, P2$n, mc_seconds,
            paste(sprintf("%s %.4f", est_names, unlist(res$cover)), collapse = ", ")))
cat(sprintf("G2 Monte Carlo: sd two-step %.4f, mean SE naive %.4f, Windmeijer %.4f; sd 2SLS %.4f, identity %.4f; J rejects %.4f\n",
            mc$sd_two_step_naive, mc$mean_se_two_step_naive, mc$mean_se_windmeijer, mc$sd_one_step, mc$sd_identity, mc$reject_J_two_step))
print(grid[, c("n", "R", paste0("cover_", est_names))], row.names = FALSE)
cat("Saved res_gmm_mc.rds\n")
