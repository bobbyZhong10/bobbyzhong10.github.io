# ------------------------------------------------------------------
# 05_overid.R -- block S4: over-identification tests in the buyer equation.
#
# buyers = const + beta n_sellers + ad_budget + coupon + log_pop + eps_B with three excluded seller-side
# instruments (fee, payout_days, referral_bonus): r - k = 2 over-identifying restrictions.
#   S4  (invalid): the referral bonus also brings buyers directly (theta = 0.05), so E[referral_bonus eps] != 0.
#   S4v (valid):   theta = 0.
# Tests: the robust Hansen J (hansen_j: two-step efficient GMM, Wooldridge 2010 eq. 8.58) and the
# homoskedastic Sargan N R^2_u (sargan, eq. 6.31; equals ivreg's diagnostics, checked). Errors are
# heteroskedastic in the referral bonus, so only J has the right size.
#   Samples: 2SLS (HC0) with all three instruments and with the two valid ones only (J with one degree of
#   freedom), J, Sargan, first-stage F_fs and F_eff, ivreg's diagnostics.
#   Monte Carlo: R_mc = 2000 valid samples (seed S4_mc): rejection rates of J and Sargan at 5% (size), and
#   R_power = 1000 invalid samples (seed S4_power): rejection rates at 5% and 1% (power), the 2SLS bias.
#
# Output: res_overid.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("iv_fns.R")
suppressPackageStartupMessages(library(ivreg))

tr <- readRDS("truth09.rds")
S4 <- readRDS("bramble_S4.rds"); S4v <- readRDS("bramble_S4v.rds")
b <- tr$S4$beta[["n_sellers"]]

analyse <- function(d) {
  m <- s4_mats(d)
  Z <- cbind(m$W, m$Zx); Z2 <- cbind(m$W, m$Zx[, c("fee", "payout_days")])
  t3 <- tsls(m$y, m$X, Z); t2 <- tsls(m$y, m$X, Z2)
  J3 <- hansen_j(m$y, m$X, Z); J2 <- hansen_j(m$y, m$X, Z2)
  sg <- sargan(m$y, m$X, Z)
  se3 <- sqrt(t3$V["n_sellers", "n_sellers"]); se2 <- sqrt(t2$V["n_sellers", "n_sellers"])
  list(tsls_est = t3$b[["n_sellers"]], tsls_se = se3, tsls_ci = ci95(t3$b[["n_sellers"]], se3),
       tsls_covers = covers_num(ci95(t3$b[["n_sellers"]], se3), b),
       tsls_valid2_est = t2$b[["n_sellers"]], tsls_valid2_se = se2, tsls_valid2_ci = ci95(t2$b[["n_sellers"]], se2),
       tsls_valid2_covers = covers_num(ci95(t2$b[["n_sellers"]], se2), b),
       J = J3$J, J_p = J3$p, J_df = J3$df, gmm_est = J3$b[["n_sellers"]], sargan = sg$stat, sargan_p = sg$p,
       J_valid2 = J2$J, J_valid2_p = J2$p, J_valid2_df = J2$df,
       F_fs = first_stage_F(m$x, m$Zx, m$W)$F, F_eff = eff_F(m$x, m$Zx, m$W))
}
invalid <- analyse(S4)
valid <- analyse(S4v)
iv4 <- ivreg(buyers ~ n_sellers + ad_budget + coupon + log_pop | ad_budget + coupon + log_pop + fee + payout_days + referral_bonus,
             data = S4)
dg4 <- summary(iv4, diagnostics = TRUE)$diagnostics
stopifnot(abs(dg4["Sargan", "statistic"] - invalid$sargan) < 1e-8, abs(coef(iv4)[["n_sellers"]] - invalid$tsls_est) < 1e-10)
invalid$diagnostics <- dg4

## ---- Monte Carlo ------------------------------------------------------------------------------
one_rep <- function(par) {
  m <- s4_mats(gen_S4(par))
  Z <- cbind(m$W, m$Zx)
  J <- hansen_j(m$y, m$X, Z); sg <- sargan(m$y, m$X, Z)
  t3 <- tsls(m$y, m$X, Z)
  c(J = J$J, sargan = sg$stat, est = t3$b[["n_sellers"]])
}
run_mc <- function(par, R, seed) {
  set.seed(seed)
  as.data.frame(t(vapply(seq_len(R), function(r) one_rep(par), numeric(3))))
}
t0 <- proc.time()[["elapsed"]]
P4 <- tr$params$S4; P4v <- tr$params$S4v
simv <- run_mc(P4v, P4$R_mc, tr$params$seed$S4_mc)
simi <- run_mc(P4, P4$R_power, tr$params$seed$S4_power)
mc_seconds <- proc.time()[["elapsed"]] - t0
c05 <- qchisq(0.95, 2); c01 <- qchisq(0.99, 2)
valid$R <- nrow(simv); valid$mc_reject_J <- mean(simv$J > c05); valid$mc_reject_sargan <- mean(simv$sargan > c05)
valid$mc_reject_J_01 <- mean(simv$J > c01); valid$mc_reject_sargan_01 <- mean(simv$sargan > c01)
valid$mc_mean_est <- mean(simv$est); valid$mc_sd_est <- sd(simv$est)
invalid$R_power <- nrow(simi); invalid$mc_reject_J <- mean(simi$J > c05); invalid$mc_reject_J_01 <- mean(simi$J > c01)
invalid$mc_reject_sargan <- mean(simi$sargan > c05); invalid$mc_mean_est <- mean(simi$est); invalid$mc_sd_est <- sd(simi$est)
invalid$mc_bias <- invalid$mc_mean_est - b

res <- list(valid = valid, invalid = invalid, crit = list(chi2_2_05 = c05, chi2_2_01 = c01),
            draws = list(valid = simv, invalid = simi))
saveRDS(res, "res_overid.rds")

for (k in c("valid", "invalid")) with(res[[k]], cat(sprintf(
  "S4 %-7s: 2SLS %.4f [%.4f, %.4f]; J %.3f (p %.4g), Sargan %.3f (p %.4g); two valid instruments: 2SLS %.4f, J %.3f (p %.3g); F_fs %.1f, F_eff %.1f\n",
  k, tsls_est, tsls_ci[1], tsls_ci[2], J, J_p, sargan, sargan_p, tsls_valid2_est, J_valid2, J_valid2_p, F_fs, F_eff)))
cat(sprintf("Monte Carlo valid (R = %d): reject at 5%%: J %.4f, Sargan %.4f; at 1%%: J %.4f, Sargan %.4f; mean 2SLS %.4f\n",
            valid$R, valid$mc_reject_J, valid$mc_reject_sargan, valid$mc_reject_J_01, valid$mc_reject_sargan_01, valid$mc_mean_est))
cat(sprintf("Monte Carlo invalid (R = %d): reject at 5%%: J %.4f (1%%: %.4f), Sargan %.4f; mean 2SLS %.4f (bias %.4f)\n",
            invalid$R_power, invalid$mc_reject_J, invalid$mc_reject_J_01, invalid$mc_reject_sargan, invalid$mc_mean_est, invalid$mc_bias))
cat(sprintf("Monte Carlo time %.1f s\nSaved res_overid.rds\n", mc_seconds))
