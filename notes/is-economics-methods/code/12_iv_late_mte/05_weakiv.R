# ============================================================================
# 05_weakiv.R  --  the weak-instrument pathology
# ----------------------------------------------------------------------------
# Same firms, same outcomes, but the encouragement barely moves adoption
# (a tiny a1).  The first stage is weak: F is low, 2SLS is biased BACK toward
# OLS, its conventional t-based confidence interval is distorted, and only a
# weak-instrument-robust interval (Anderson-Rubin) is trustworthy.
# ============================================================================

suppressMessages({library(fixest); library(ivmodel)})

d <- readRDS("meridian_data.rds")
dat <- d$dat; truth <- d$truth; par <- d$params

## ---- weak first stage ------------------------------------------------------
fs_w <- feols(D_weak ~ Z, data = dat)
fs_w_coef <- coef(fs_w)["Z"]
F_w_iid <- fitstat(fs_w, "wald", vcov = "iid")$wald$stat
F_w_hc  <- fitstat(fs_w, "wald", vcov = "hetero")$wald$stat   # effective F

## ---- weak 2SLS -------------------------------------------------------------
iv_w <- feols(Y_weak ~ 1 | D_weak ~ Z, data = dat)
b_w  <- coef(iv_w)["fit_D_weak"]
se_w <- se(iv_w)["fit_D_weak"]
ci_w <- confint(iv_w)["fit_D_weak", ]          # conventional (distorted) CI

## for reference: OLS on the weak-design data, and the strong-design 2SLS
b_ols_w <- coef(lm(Y_weak ~ D_weak, data = dat))["D_weak"]

## ---- Anderson-Rubin weak-IV-robust inference -------------------------------
ivm <- ivmodel(Y = dat$Y_weak, D = dat$D_weak, Z = matrix(dat$Z, ncol = 1))
ar <- ivm$AR
ar_ci <- ar$ci                                  # robust confidence set
ar_stat <- ar$Fstat
ar_p <- ar$p.value

cat("=== Weak-instrument variant ===\n")
cat(sprintf("Weak first stage E[D|Z=1]-E[D|Z=0] = %.4f  (complier share %.4f)\n",
            fs_w_coef, par$complier_share_weak))
cat(sprintf("First-stage F (IID)                = %.2f\n", F_w_iid))
cat(sprintf("Effective F (robust)               = %.2f   <-- well below 10\n", F_w_hc))
cat("\n--- point estimates ---\n")
cat(sprintf("Weak 2SLS                          = %.4f  (SE %.4f)\n", b_w, se_w))
cat(sprintf("OLS (weak-design data)             = %.4f\n", b_ols_w))
cat(sprintf("TRUE complier-LATE (strong design) = %.4f\n", truth$LATE))
cat(sprintf("--> weak 2SLS is pulled toward OLS, away from the LATE target\n"))
cat("\n--- confidence intervals ---\n")
cat(sprintf("Conventional 2SLS 95%% CI           = [%.3f, %.3f]  (width %.3f)\n",
            ci_w[1], ci_w[2], ci_w[2] - ci_w[1]))
cat(sprintf("Anderson-Rubin 95%% CI (robust)     = [%.3f, %.3f]  (width %.3f)\n",
            ar_ci[1], ar_ci[2], ar_ci[2] - ar_ci[1]))
cat(sprintf("AR test stat = %.3f , p = %.3f\n", ar_stat, ar_p))

## ---- Monte Carlo: the SAMPLING DISTRIBUTION of weak 2SLS -------------------
# One realized dataset is noisy; the pathology is a property of the sampling
# distribution.  Redraw the weak design B times (same parameters, same n) and
# record 2SLS, its conventional-CI coverage of the true LATE, and the F.
mc_weak <- function(B = 1000, n = par$n) {
  a0 <- par$a0; a1w <- par$a1_weak
  b2sls <- numeric(B); cov_conv <- logical(B); Fs <- numeric(B)
  P0 <- pnorm(a0); P1 <- pnorm(a0 + a1w)
  for (b in seq_len(B)) {
    U <- runif(n); lat <- qnorm(U)
    W <- par$rho_wu * lat + sqrt(1 - par$rho_wu^2) * rnorm(n)
    Y0 <- par$b_y0 + par$b_w * W - par$g0 * U + rnorm(n, 0, par$sd_e0)
    tau <- par$t0 - par$t1 * U + rnorm(n, 0, par$sd_tau)
    Z <- rbinom(n, 1, 0.5); D <- as.integer(U < pnorm(a0 + a1w * Z))
    Y <- Y0 + D * tau
    LATE_b <- mean(tau[U >= P0 & U < P1])
    m <- feols(Y ~ 1 | D ~ Z, data = data.frame(Y, D, Z), notes = FALSE)
    bb <- coef(m)["fit_D"]; ss <- se(m)["fit_D"]
    b2sls[b] <- bb
    cov_conv[b] <- (LATE_b >= bb - 1.96 * ss) & (LATE_b <= bb + 1.96 * ss)
    Fs[b] <- fitstat(feols(D ~ Z, data = data.frame(D, Z), notes = FALSE),
                     "wald", vcov = "hetero")$wald$stat
  }
  list(b2sls = b2sls, cov = mean(cov_conv), medF = median(Fs))
}
set.seed(2024)
mc <- mc_weak(B = 1000)
mc_median <- median(mc$b2sls)
mc_mean   <- mean(mc$b2sls)
mc_cov    <- mc$cov

cat("\n=== Weak-IV Monte Carlo (1000 redraws) ===\n")
cat(sprintf("median 2SLS  = %.4f\n", mc_median))
cat(sprintf("mean   2SLS  = %.4f  (OLS ~ %.1f , true LATE %.3f: pulled toward OLS)\n",
            mc_mean, b_ols_w, truth$LATE))
cat(sprintf("conventional 95%% CI coverage of true LATE = %.3f  (nominal 0.95)\n", mc_cov))
cat(sprintf("median first-stage effective F = %.2f\n", mc$medF))

saveRDS(list(fs_w_coef = fs_w_coef, F_w_iid = F_w_iid, F_w_hc = F_w_hc,
             b_w = b_w, se_w = se_w, ci_w = ci_w, b_ols_w = b_ols_w,
             ar_ci = ar_ci, ar_stat = ar_stat, ar_p = ar_p,
             mc_median = mc_median, mc_mean = mc_mean, mc_cov = mc_cov,
             mc_medF = mc$medF),
        "res_weak.rds")
