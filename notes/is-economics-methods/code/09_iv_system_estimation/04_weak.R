# ------------------------------------------------------------------
# 04_weak.R -- block S3: the weak fee pilot (one endogenous regressor, one instrument).
#
# Buyer equation buyers = const + beta n_sellers + ad_budget + coupon + log_pop + eps_B, instrument fee,
# with the fee nudged by only about USD 0.35 (concentration parameter mu2 = 1 at n = 1,500).
#   The sample: 2SLS with the HC0 Wald interval (tsls) and the homoskedastic one, first-stage F_fs and the
#   effective F_eff (Montiel Olea and Pflueger 2013; with one instrument it is the robust Wald F), the
#   Anderson-Rubin set (ar_ci, checked against ivmodel::AR.test), and from ivmodel the CLR set (equal to AR
#   with one instrument), LIML (equal to 2SLS when just identified) and Fuller (b = 1).
#   Monte Carlo: R_mc = 2000 replications at mu2 = 1 (seed S3_mc) and R_grid = 1000 at each other mu2 on the
#   grid (seed S3_grid + j for grid point j): coverage of the nominal 95% Wald (HC0 and homoskedastic) and AR
#   sets, the distribution of F_fs, the median of the 2SLS estimate (its mean does not exist when just
#   identified) and the shape of the AR set (bounded, two rays, the whole line).
# Each replication partials (1, ad_budget, coupon, log_pop) out once (Frisch-Waugh-Lovell) and uses the
# scalar formulas; on the sample they are checked against tsls(), first_stage_F() and ar_ci().
#
# Output: res_weak.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("iv_fns.R")
suppressPackageStartupMessages(library(ivmodel))

tr <- readRDS("truth09.rds")
P3 <- tr$params$S3
S3 <- readRDS("bramble_S3.rds")
b <- tr$S3$beta[["n_sellers"]]
z975 <- qnorm(0.975)

## ---- one replication (scalar FWL formulas) ------------------------------------------------
weak_stats <- function(d) {
  W <- cbind(1, d$ad_budget, d$coupon, d$log_pop)
  qw <- qr(W)
  rs <- function(v) v - qr.fitted(qw, v)
  yt <- rs(d$buyers); xt <- rs(d$n_sellers); zt <- rs(d$fee)
  n <- length(yt); df2 <- n - ncol(W) - 1
  zz <- sum(zt^2); zx <- sum(zt * xt); zy <- sum(zt * yt)
  est <- zy / zx
  e <- yt - est * xt
  se_r <- sqrt(sum(zt^2 * e^2)) / abs(zx)                             # HC0 (eq. 5.34) for the n_sellers coefficient
  se_h <- sqrt(sum(e^2) / (n - ncol(W) - 1) * zz) / abs(zx)          # eq. 5.27 (n - k, k = 5)
  pi <- zx / zz; v <- xt - pi * zt
  Ffs <- pi^2 * zz / (sum(v^2) / df2)
  a <- qf(0.95, 1, df2) / df2                                        # AR set: A b^2 + B b + C <= 0
  yy <- sum(yt^2); xy <- sum(xt * yt); xx <- sum(xt^2)
  A <- zx^2 / zz * (1 + a) - a * xx; B <- -2 * (zy * zx / zz * (1 + a) - a * xy); C <- zy^2 / zz * (1 + a) - a * yy
  disc <- B^2 - 4 * A * C
  r <- if (disc >= 0) sort(c((-B - sqrt(disc)) / (2 * A), (-B + sqrt(disc)) / (2 * A))) else c(NA, NA)
  type <- if (A > 0) 1 else if (disc > 0) 2 else 3                   # 1 bounded, 2 two rays, 3 whole line
  ar_cov <- switch(type, r[1] <= b && b <= r[2], b <= r[1] || b >= r[2], TRUE)
  c(est = est, se_r = se_r, se_h = se_h, F = Ffs, wald_r = abs(est - b) <= z975 * se_r, wald_h = abs(est - b) <= z975 * se_h,
    ar = ar_cov, ar_type = type, ar_lo = r[1], ar_hi = r[2])
}

## =================================================================
## The sample
## =================================================================
W3 <- cbind(const = 1, ad_budget = S3$ad_budget, coupon = S3$coupon, log_pop = S3$log_pop)
X3 <- cbind(const = 1, n_sellers = S3$n_sellers, W3[, -1]); Z3 <- cbind(W3, fee = S3$fee)
t3 <- tsls(S3$buyers, X3, Z3); t3h <- tsls(S3$buyers, X3, Z3, vcov = "homo")
fs3 <- first_stage_F(S3$n_sellers, cbind(S3$fee), W3)
Feff3 <- eff_F(S3$n_sellers, cbind(S3$fee), W3)
ar3 <- ar_ci(S3$buyers, S3$n_sellers, cbind(S3$fee), W3)
ws <- weak_stats(S3)
# (relative tolerance 1e-6: with an instrument this weak X_hat is nearly collinear and rounding grows)
stopifnot(abs(ws[["est"]] / t3$b[["n_sellers"]] - 1) < 1e-8, abs(ws[["se_r"]] / sqrt(t3$V["n_sellers", "n_sellers"]) - 1) < 1e-6,
          abs(ws[["se_h"]] / sqrt(t3h$V["n_sellers", "n_sellers"]) - 1) < 1e-6, abs(ws[["F"]] / fs3$F - 1) < 1e-8,
          ws[["ar_type"]] == match(ar3$type, c("bounded", "two_rays", "real_line")),
          ar3$type == "real_line" || (abs(ws[["ar_lo"]] - ar3$lo) < 1e-8 && abs(ws[["ar_hi"]] - ar3$hi) < 1e-8))
ivm <- ivmodel(Y = S3$buyers, D = S3$n_sellers, Z = cbind(S3$fee), X = W3[, -1])
ar_pkg <- AR.test(ivm); clr_pkg <- CLR(ivm); liml <- LIML(ivm); full <- Fuller(ivm, b = 1)
# the hand-written AR set equals ivmodel's (one row per piece: [lo, hi], or (-Inf, lo] and [hi, Inf))
pk <- ar_pkg$ci
stopifnot(switch(ar3$type,
                 bounded = nrow(pk) == 1 && abs(pk[1, 1] - ar3$lo) < 1e-6 && abs(pk[1, 2] - ar3$hi) < 1e-6,
                 two_rays = nrow(pk) == 2 && abs(pk[1, 2] - ar3$lo) < 1e-6 && abs(pk[2, 1] - ar3$hi) < 1e-6,
                 real_line = grepl("Whole Real Line", ar_pkg$ci.info), FALSE),
          abs(liml$point.est - t3$b[["n_sellers"]]) < 1e-6)                 # LIML = 2SLS when just identified
ols3 <- ols_fit(S3$buyers, X3)
sample <- list(tsls_est = t3$b[["n_sellers"]], tsls_se = sqrt(t3$V["n_sellers", "n_sellers"]),
               tsls_se_homo = sqrt(t3h$V["n_sellers", "n_sellers"]),
               wald_ci = ci95(t3$b[["n_sellers"]], sqrt(t3$V["n_sellers", "n_sellers"])),
               wald_covers = covers_num(ci95(t3$b[["n_sellers"]], sqrt(t3$V["n_sellers", "n_sellers"])), b),
               F_fs = fs3$F, F_eff = Feff3, pi_fee = fs3$pi, ar = ar3, ar_covers = as.numeric(ar_covers(ar3, b)),
               ar_pkg_ci = ar_pkg$ci, ar_pkg_info = ar_pkg$ci.info, ar_at_truth = ar_test(S3$buyers, S3$n_sellers, cbind(S3$fee), W3, b),
               ar_at_zero = ar_test(S3$buyers, S3$n_sellers, cbind(S3$fee), W3, 0),
               clr_ci = clr_pkg$ci, clr_info = clr_pkg$ci.info,
               liml_est = unname(liml$point.est), fuller_est = unname(full$point.est), fuller_se = unname(full$std.err),
               ols_est = ols3$b[["n_sellers"]], ols_se = sqrt(ols3$V["n_sellers", "n_sellers"]))

## =================================================================
## Monte Carlo
## =================================================================
run_mc <- function(par, R, seed) {
  set.seed(seed)
  as.data.frame(t(vapply(seq_len(R), function(r) weak_stats(gen_S3(par)), numeric(10))))
}
summ <- function(sim, mu2, plim_ols) {
  list(mu2 = mu2, R = nrow(sim), mean_F = mean(sim$F), median_F = median(sim$F), share_F_below_10 = mean(sim$F < 10),
       cover_wald = mean(sim$wald_r), cover_wald_homo = mean(sim$wald_h), cover_ar = mean(sim$ar),
       share_ar_bounded = mean(sim$ar_type == 1), share_ar_two_rays = mean(sim$ar_type == 2),
       share_ar_real_line = mean(sim$ar_type == 3),
       median_est = median(sim$est), q25_est = unname(quantile(sim$est, 0.25)), q75_est = unname(quantile(sim$est, 0.75)),
       median_bias = median(sim$est) - b, plim_ols = plim_ols,
       median_se_r = median(sim$se_r), share_wald_reject_zero = mean(abs(sim$est) > z975 * sim$se_r))
}
t0 <- proc.time()[["elapsed"]]
g3 <- tr$S3$grid
sim_main <- run_mc(P3, P3$R_mc, tr$params$seed$S3_mc)
main <- summ(sim_main, P3$mu2, tr$S3$plim_ols)
grid <- do.call(rbind, lapply(seq_len(nrow(g3)), function(j) {
  s <- if (isTRUE(all.equal(g3$mu2[j], P3$mu2))) main
       else summ(run_mc(modifyList(P3, list(fee_sd = g3$fee_sd[j])), P3$R_grid, tr$params$seed$S3_grid + j), g3$mu2[j], g3$plim_ols[j])
  as.data.frame(s)
}))
mc_seconds <- proc.time()[["elapsed"]] - t0

cover <- list(wald = main$cover_wald, wald_homo = main$cover_wald_homo, ar = main$cover_ar, R = main$R, mu2 = main$mu2)
res <- list(cover = cover, main = main, grid = grid, sample = sample,
            draws = sim_main[, c("est", "se_r", "F", "wald_r", "ar", "ar_type")])
saveRDS(res, "res_weak.rds")

cat(sprintf("S3 sample: 2SLS %.4f (HC0 se %.4f) Wald [%.3f, %.3f]; F_fs %.3f, F_eff %.3f; AR set %s [%.3f, %.3f]; OLS %.4f\n",
            sample$tsls_est, sample$tsls_se, sample$wald_ci[1], sample$wald_ci[2], sample$F_fs, sample$F_eff, ar3$type,
            ar3$lo, ar3$hi, sample$ols_est))
cat(sprintf("S3 Monte Carlo (mu2 = %g, R = %d): Wald coverage %.4f (homo %.4f), AR %.4f; median 2SLS %.4f; median F %.3f\n",
            main$mu2, main$R, main$cover_wald, main$cover_wald_homo, main$cover_ar, main$median_est, main$median_F))
print(grid[, c("mu2", "R", "mean_F", "cover_wald", "cover_wald_homo", "cover_ar", "median_est", "share_ar_bounded")], row.names = FALSE)
cat(sprintf("Monte Carlo time %.1f s\nSaved res_weak.rds\n", mc_seconds))
