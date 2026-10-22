# ------------------------------------------------------------------
# 05_twostep.R -- block C4: two-step MLE with a generated usage index (Wooldridge 2010 section 13.10.1).
#
# Step 1 (CMLE): Poisson regression of first-week check-ins v on z = (1, hh_size, quiz_std): gamma_hat,
#   usage-intensity index u_hat = z'gamma_hat.
# Step 2 (MLE): Weibull PH hazard with x = (1, fee_usd, u_hat); the coefficient of interest is beta_u
#   (truth06$C4$beta = -0.80).
# Standard errors of beta_u_hat (twostep_weibull in ml_fns.R):
#   naive      A2^-1 / n, the second-step Hessian treating u_hat as data (what software prints);
#   naive_sand the second-step sandwich, still treating u_hat as data;
#   corrected  A2^-1 D A2^-1 / n, D = n^-1 sum (s_i + Gamma_hat psi_i)(s_i + Gamma_hat psi_i)',
#              Wooldridge 2010 eqs. (12.39) to (12.41), estimated by (12.61) and (12.62); psi_i is the
#              Poisson CMLE influence function (eq. 13.65). Checked in tests/test_fns.R against the
#              stacked (joint) M-estimator sandwich.
#   oracle     the infeasible second step on the true index u = z'gamma (gamma known).
# Here beta_u != 0 makes Gamma_0 = E[grad_gamma s] != 0, and the first-step and second-step scores are
# uncorrelated (given z, v and t* are independent), so D_0 = B_0 + Gamma_0 Avar(gamma_hat) Gamma_0' and the
# corrected SE is the larger one. Section 13.10.2's result (a CMLE first step can make the two-step
# estimator MORE efficient than plugging in gamma_0) needs the second-step objective itself to depend on
# the first-step data v (eq. 13.67, as in inverse probability weighting); here it does not.
# Monte Carlo: R_mc = 1000 samples of n = 5000 (seed C4_mc); coverage of est -/+ qnorm(0.975) se of the
# nominal 95% intervals for beta_u.
#
# Output: res_twostep.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("ml_fns.R")

tr <- readRDS("truth06.rds")
P4 <- tr$params$C4
C4 <- readRDS("sorrel_C4.rds")
z <- qnorm(0.975)
bu <- tr$C4$beta
j <- 3L                                                           # beta_u in (beta_0, beta_fee, beta_u, log kappa)

## ---- the sample ---------------------------------------------------------------------------------
f <- twostep_weibull(C4)
stopifnot(f$converged)
est <- f$est
se_naive <- sqrt(f$V_naive[j, j])
se_naive_sand <- sqrt(f$V_naive_sand[j, j])
se_corrected <- sqrt(f$V_corrected[j, j])
fo <- fit_weibull(C4$t_months, C4$delta, cbind(1, C4$fee_usd, C4$u_true), engine = "newton")
stopifnot(fo$converged)
covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])
naive <- list(est = est, se = se_naive, ci = ci95(est, se_naive))
corrected <- list(est = est, se = se_corrected, ci = ci95(est, se_corrected))
oracle <- list(est = fo$par[j], se = sqrt(fo$V_hess[j, j]), ci = ci95(fo$par[j], sqrt(fo$V_hess[j, j])), par = fo$par)
g1 <- glm(usage_wk1 ~ hh_size + quiz_std, family = poisson(), data = C4, control = glm.control(epsilon = 1e-12))
stopifnot(max(abs(unname(coef(g1)) - f$gamma)) < 1e-8)

## ---- Monte Carlo ----------------------------------------------------------------------------------
one_rep <- function() {
  d <- gen_C4(P4)
  r <- twostep_weibull(d, start = c(P4$beta, log(P4$kappa)))
  o <- fit_weibull(d$t_months, d$delta, cbind(1, d$fee_usd, d$u_true), start = c(P4$beta, log(P4$kappa)), engine = "newton")
  c(est = r$est, se_naive = sqrt(r$V_naive[j, j]), se_naive_sand = sqrt(r$V_naive_sand[j, j]),
    se_corrected = sqrt(r$V_corrected[j, j]), oracle = o$par[j], se_oracle = sqrt(o$V_hess[j, j]),
    converged = r$converged && o$converged)
}
t0 <- proc.time()[["elapsed"]]
set.seed(tr$params$seed$C4_mc)
sim <- t(vapply(seq_len(P4$R_mc), function(r) one_rep(), numeric(7)))
mc_seconds <- proc.time()[["elapsed"]] - t0
stopifnot(all(sim[, "converged"] == 1))
cov <- function(e, s) mean(abs(sim[, e] - bu) <= z * sim[, s])
mc <- list(R = nrow(sim), n_converged = sum(sim[, "converged"] == 1), n = P4$n,
           cover_naive = cov("est", "se_naive"), cover_naive_sand = cov("est", "se_naive_sand"),
           cover_corrected = cov("est", "se_corrected"), cover_oracle = cov("oracle", "se_oracle"),
           mean_est = mean(sim[, "est"]), sd_est = sd(sim[, "est"]),
           mean_se_naive = mean(sim[, "se_naive"]), mean_se_naive_sand = mean(sim[, "se_naive_sand"]),
           mean_se_corrected = mean(sim[, "se_corrected"]),
           mean_oracle = mean(sim[, "oracle"]), sd_oracle = sd(sim[, "oracle"]), mean_se_oracle = mean(sim[, "se_oracle"]),
           share_corrected_ge_naive = mean(sim[, "se_corrected"] >= sim[, "se_naive"]))

res <- list(n = nrow(C4), share_cens = 1 - mean(C4$delta), mean_usage = mean(C4$usage_wk1), share_zero_usage = mean(C4$usage_wk1 == 0),
            gamma = f$gamma, gamma_se = sqrt(diag(crossprod(f$R1) / nrow(C4)^2)), par = f$par, est = est,
            se_naive = se_naive, se_naive_sand = se_naive_sand, se_corrected = se_corrected, se_ratio = se_corrected / se_naive,
            naive = c(naive, list(covers = covers_num(naive$ci, bu))),
            corrected = c(corrected, list(covers = covers_num(corrected$ci, bu))),
            oracle = c(oracle, list(covers = covers_num(oracle$ci, bu))),
            F2 = f$F2, mc = mc)   # run time is printed, not saved (byte-reproducible .rds)
saveRDS(res, "res_twostep.rds")

cat(sprintf("C4: gamma_hat %s; beta_u_hat %.4f (truth %.2f); SE naive %.4f, naive sandwich %.4f, corrected %.4f (ratio %.3f)\n",
            paste(sprintf("%.4f", f$gamma), collapse = " "), est, bu, se_naive, se_naive_sand, se_corrected, res$se_ratio))
cat(sprintf("  corrected CI [%.4f, %.4f]; naive CI [%.4f, %.4f]; oracle (true index) %.4f (SE %.4f)\n",
            corrected$ci[1], corrected$ci[2], naive$ci[1], naive$ci[2], oracle$est, oracle$se))
cat(sprintf("C4 Monte Carlo (%d reps, %.1f s): coverage naive %.4f, naive sandwich %.4f, corrected %.4f, oracle %.4f; sd %.4f, mean SE %.4f / %.4f\n",
            mc$R, mc_seconds, mc$cover_naive, mc$cover_naive_sand, mc$cover_corrected, mc$cover_oracle, mc$sd_est,
            mc$mean_se_naive, mc$mean_se_corrected))
cat("Saved res_twostep.rds\n")
