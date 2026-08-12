# ------------------------------------------------------------------
# 03_cuped.R -- CUPED variance reduction, Module 14.
# Y_cv = Y - theta * (X - mean(X)), theta = Cov(Y,X)/Var(X).
# Randomization keeps X balanced, so the ATE is unchanged in
# expectation; variance falls by the factor (1 - rho^2).
#   - same n=15000 run, plain vs CUPED-adjusted (SE, t, p)
#   - empirical variance-reduction factor vs 1-rho^2
#   - empirical power at n=15000, plain vs CUPED
#   - n needed for 80% power under CUPED
# Writes res_cuped.rds
# ------------------------------------------------------------------
set.seed(1414)
p <- readRDS("lumen_params.rds")
gen <- p$gen; delta <- p$delta; sigma <- p$sd_Y0
za <- qnorm(0.975); zb <- qnorm(0.80)

cuped_estimate <- function(Y, X, Z) {
  theta <- cov(Y, X) / var(X)                 # pooled theta
  Ycv   <- Y - theta * (X - mean(X))
  est   <- mean(Ycv[Z == 1]) - mean(Ycv[Z == 0])
  n1 <- sum(Z == 1); n0 <- sum(Z == 0)
  se  <- sqrt(var(Ycv[Z == 1]) / n1 + var(Ycv[Z == 0]) / n0)
  c(est = est, se = se, t = est / se, theta = theta)
}

## ---- same opening run, plain vs CUPED ----------------------------
n_run <- 15000L
d <- gen(2 * n_run, seed = 24)
Z <- rep(0:1, each = n_run)
Y <- ifelse(Z == 1, d$Y1, d$Y0); X <- d$X

plain <- {
  est <- mean(Y[Z == 1]) - mean(Y[Z == 0])
  se  <- sqrt(var(Y[Z == 1]) / n_run + var(Y[Z == 0]) / n_run)
  c(est = est, se = se, t = est / se, p = 2 * pnorm(-abs(est / se)))
}
cv <- cuped_estimate(Y, X, Z)
cv_p <- 2 * pnorm(-abs(cv["t"]))

## ---- empirical variance-reduction factor -------------------------
R <- 3000
vr <- numeric(R); est_plain <- est_cv <- numeric(R)
set.seed(11)
for (r in 1:R) {
  dd <- gen(2 * n_run)
  zz <- rep(0:1, each = n_run)
  yy <- ifelse(zz == 1, dd$Y1, dd$Y0); xx <- dd$X
  est_plain[r] <- mean(yy[zz == 1]) - mean(yy[zz == 0])
  th <- cov(yy, xx) / var(xx)
  yc <- yy - th * (xx - mean(xx))
  est_cv[r] <- mean(yc[zz == 1]) - mean(yc[zz == 0])
}
vr_emp <- var(est_cv) / var(est_plain)       # should be ~ 1 - rho^2
rho2   <- p$rho_corr^2

## ---- empirical power, plain vs CUPED, at n=15000 -----------------
sim_power2 <- function(n, del, R = 3000, seed = 3) {
  set.seed(seed); rp <- rc <- logical(R)
  for (r in 1:R) {
    dd <- gen(2 * n); zz <- rep(0:1, each = n)
    yy <- ifelse(zz == 1, dd$Y0 + del, dd$Y0); xx <- dd$X
    ep <- mean(yy[zz == 1]) - mean(yy[zz == 0])
    sp <- sqrt(var(yy[zz == 1]) / n + var(yy[zz == 0]) / n)
    rp[r] <- abs(ep / sp) > za
    th <- cov(yy, xx) / var(xx); yc <- yy - th * (xx - mean(xx))
    ec <- mean(yc[zz == 1]) - mean(yc[zz == 0])
    sc <- sqrt(var(yc[zz == 1]) / n + var(yc[zz == 0]) / n)
    rc[r] <- abs(ec / sc) > za
  }
  c(plain = mean(rp), cuped = mean(rc))
}
pw <- sim_power2(n_run, delta)

## ---- n for 80% power under CUPED ---------------------------------
sigma_cv <- sigma * sqrt(1 - rho2)
n_cuped  <- ceiling(2 * sigma_cv^2 * (za + zb)^2 / delta^2)

res_cuped <- list(
  plain = plain, cuped = c(cv, p = unname(cv_p)),
  vr_emp = vr_emp, one_minus_rho2 = 1 - rho2, rho = p$rho_corr,
  bias_plain = mean(est_plain) - delta, bias_cv = mean(est_cv) - delta,
  power = pw, n_cuped = n_cuped, n_plain = 56516
)
saveRDS(res_cuped, "res_cuped.rds")

cat(sprintf("Opening run, plain : est %+.4f SE %.4f t %.2f p %.3f\n",
            plain["est"], plain["se"], plain["t"], plain["p"]))
cat(sprintf("Opening run, CUPED : est %+.4f SE %.4f t %.2f p %.3f (theta %.3f)\n",
            cv["est"], cv["se"], cv["t"], cv_p, cv["theta"]))
cat(sprintf("\nEmpirical variance-reduction factor = %.3f  (1-rho^2 = %.3f)\n",
            vr_emp, 1 - rho2))
cat(sprintf("Unbiasedness: bias(plain) = %+.5f, bias(CUPED) = %+.5f\n",
            mean(est_plain) - delta, mean(est_cv) - delta))
cat(sprintf("Empirical power at n=15000: plain %.3f, CUPED %.3f\n",
            pw["plain"], pw["cuped"]))
cat(sprintf("n for 80%% power: plain 56516/arm, CUPED %d/arm (%.2fx fewer)\n",
            n_cuped, 56516 / n_cuped))
