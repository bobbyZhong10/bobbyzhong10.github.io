# ------------------------------------------------------------------
# 02_censoring.R -- block C1: Weibull tenure under independent right censoring.
#
# x = (1, fee_usd, referred); the coefficient of interest is beta_fee in the PH hazard
# h(t | x) = exp(x'beta) kappa t^(kappa - 1) (truth06$C1$beta = 0.02).
#   ols_observed    the naive regression: OLS of log t_i on x with every censored spell treated as an
#                   exit, mapped to the PH scale through the Weibull log-linear form (ols_weibull):
#                   log t* = -x'beta / kappa + eps / kappa, so beta_fee = -kappa_hat b_fee with
#                   kappa_hat = pi / sqrt(6 s2). Censoring cuts the long spells, so the slope and the
#                   residual variance are both wrong.
#   ols_uncensored  the same regression on the uncensored spells only (selects short spells).
#   ols_latent      the same regression on t_i* (infeasible benchmark: consistent on complete data).
#   mle             Weibull MLE with the censored-data likelihood (Wooldridge 2010 eq. 22.24), maxLik
#                   "NR" with analytic derivatives; 95% interval from the Hessian (model-based) SE.
#   survreg         survival::survreg(dist = "weibull") on the same data: an AFT fit. Conversion to PH
#                   (lit_software 13.5, question 11): kappa = 1 / fit$scale, beta_PH = -coef / fit$scale
#                   (intercept included: baseline log scale = -intercept / scale). Checked here to 1e-5.
#   se_hess, se_opg the fee coefficient's SE from A_hat^-1 / n and from B_hat^-1 / n, both by info_hat
#                   (numDeriv): the model is right, so A_0 = B_0 and the two agree.
#   trinity         W, LR, LM for H_0: kappa = 1 (exponential tenure: no duration dependence), i.e.
#                   R(theta) = log kappa = 0, nu = 1; the restricted MLE is the exponential model.
#                   The LR is also checked against survreg's exponential and Weibull log likelihoods and
#                   against lmtest::lrtest on the two survreg fits.
#   aft             the same comparison on the log-tenure (AFT) scale, -beta_fee / kappa, where OLS
#                   lives: OLS slope on observed and latent tenure against the survreg coefficient.
#
# Output: res_cens.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("ml_fns.R")
suppressPackageStartupMessages({ library(survival); library(lmtest) })

tr <- readRDS("truth06.rds")
C1 <- readRDS("sorrel_C1.rds")
X <- cbind(1, C1$fee_usd, C1$referred)
t <- C1$t_months
d <- C1$delta
n <- nrow(X)
j <- 2L                                                          # position of beta_fee in (beta_0, beta_fee, beta_ref, log kappa)

## ---- naive regressions (ols_weibull) -------------------------------------------------------
ols_block <- function(o) {
  est <- o$beta[j]
  se <- o$se_beta[j]
  list(est = est, se = se, ci = ci95(est, se), kappa = o$kappa, se_kappa = o$se_kappa, beta = o$beta, se_beta = o$se_beta,
       slope = o$b[j], slope_se = o$se_b[j], slope_ci = ci95(o$b[j], o$se_b[j]), s2 = o$s2)
}
ols_observed <- ols_block(ols_weibull(log(t), X))
u <- d == 1
ols_uncensored <- c(ols_block(ols_weibull(log(t[u]), X[u, ])), list(n = sum(u)))
ols_latent <- ols_block(ols_weibull(log(C1$t_star_true), X))

## ---- MLE ---------------------------------------------------------------------------------------
fm <- fit_weibull(t, d, X)                                     # maxLik, method "NR"
stopifnot(fm$converged)
p <- fm$par
kap <- exp(p[4])
ih <- info_hat(p, nll_weibull_cens, t = t, delta = d, X = X)    # numDeriv: A_hat and B_hat
V_h <- solve(ih$A) / n
V_o <- solve(ih$B) / n
V_s <- solve(ih$A) %*% ih$B %*% solve(ih$A) / n
se_hess <- sqrt(V_h[j, j])
se_opg <- sqrt(V_o[j, j])
se_sand <- sqrt(V_s[j, j])
stopifnot(abs(se_hess / sqrt(fm$V_hess[j, j]) - 1) < 1e-5, abs(se_opg / sqrt(fm$V_opg[j, j]) - 1) < 1e-5)   # numDeriv = analytic
mle <- list(est = p[j], se = se_hess, ci = ci95(p[j], se_hess), par = p,
            se_hess_vec = sqrt(diag(V_h)), se_opg_vec = sqrt(diag(V_o)), se_sand_vec = sqrt(diag(V_s)),
            ci_opg = ci95(p[j], se_opg), ci_sand = ci95(p[j], se_sand),
            kappa = kap, se_kappa = kap * sqrt(V_h[4, 4]), ci_kappa = ci95(kap, kap * sqrt(V_h[4, 4])),
            ref = list(est = p[3], se = sqrt(V_h[3, 3]), ci = ci95(p[3], sqrt(V_h[3, 3]))),
            nll = fm$nll, iter = fm$iter, grad_max = fm$grad_max,
            median_fee60 = (log(2) / exp(p[1] + 60 * p[2]))^(1 / kap))

## ---- survreg (AFT) and the conversion ------------------------------------------------------------
sv <- survreg(Surv(t_months, delta) ~ fee_usd + referred, data = C1, dist = "weibull",
              control = survreg.control(rel.tolerance = 1e-12, maxiter = 100))
ph <- survreg_to_ph(sv)
stopifnot(max(abs(ph$par - p)) < 1e-5, abs(sv$loglik[2] + fm$nll) < 1e-6)
svr <- list(coef = unname(coef(sv)), scale = sv$scale, se_coef = sqrt(diag(vcov(sv)))[1:3], ph_par = ph$par,
            ph_se = sqrt(diag(ph$V)), max_abs_diff = max(abs(ph$par - p)), loglik = sv$loglik[2])

## ---- trinity: H_0 kappa = 1 ------------------------------------------------------------------------
fr <- fit_weibull(t, d, X, log_kappa = 0)
stopifnot(fr$converged)
tri <- trinity(p, fr$par, nll_weibull_cens, function(q) q[4], t = t, delta = d, X = X)
se_exp <- survreg(Surv(t_months, delta) ~ fee_usd + referred, data = C1, dist = "exponential",
                  control = survreg.control(rel.tolerance = 1e-12, maxiter = 100))
LR_sv <- 2 * (sv$loglik[2] - se_exp$loglik[2])
lrt <- lmtest::lrtest(se_exp, sv)
stopifnot(abs(LR_sv - tri$LR) < 1e-5, abs(lrt$Chisq[2] - tri$LR) < 1e-5, lrt$Df[2] == tri$df)
trinity_res <- c(tri, list(p_W = pchisq(tri$W, tri$df, lower.tail = FALSE), p_LR = pchisq(tri$LR, tri$df, lower.tail = FALSE),
                           p_LM = pchisq(tri$LM, tri$df, lower.tail = FALSE), crit_5pct = qchisq(0.95, tri$df),
                           par_r = fr$par, LR_survreg = LR_sv))

## ---- the AFT (log-tenure) scale ----------------------------------------------------------------------
aft <- list(truth = tr$C1$beta_aft,
            ols_observed = ols_observed$slope, ols_observed_se = ols_observed$slope_se, ols_observed_ci = ols_observed$slope_ci,
            ols_latent = ols_latent$slope, ols_latent_se = ols_latent$slope_se, ols_latent_ci = ols_latent$slope_ci,
            mle = svr$coef[j], mle_se = svr$se_coef[j], mle_ci = ci95(svr$coef[j], svr$se_coef[j]))

res <- list(n = n, n_exit = sum(d), share_cens = 1 - mean(d), mean_t_observed = mean(t),
            ols_observed = ols_observed, ols_uncensored = ols_uncensored, ols_latent = ols_latent, mle = mle,
            se_hess = se_hess, se_opg = se_opg, se_sand = se_sand, A = ih$A, B = ih$B,
            survreg = svr, trinity = trinity_res, aft = aft)
saveRDS(res, "res_cens.rds")

cat(sprintf("C1: censored %.3f. beta_fee: truth %.4f\n", res$share_cens, tr$C1$beta))
cat(sprintf("  OLS on observed log tenure  %.5f (SE %.5f) CI [%.5f, %.5f]; kappa_hat %.3f\n",
            ols_observed$est, ols_observed$se, ols_observed$ci[1], ols_observed$ci[2], ols_observed$kappa))
cat(sprintf("  OLS on uncensored spells    %.5f (SE %.5f) CI [%.5f, %.5f]; kappa_hat %.3f\n",
            ols_uncensored$est, ols_uncensored$se, ols_uncensored$ci[1], ols_uncensored$ci[2], ols_uncensored$kappa))
cat(sprintf("  OLS on latent tenure t*     %.5f (SE %.5f) CI [%.5f, %.5f]; kappa_hat %.3f\n",
            ols_latent$est, ols_latent$se, ols_latent$ci[1], ols_latent$ci[2], ols_latent$kappa))
cat(sprintf("  Weibull MLE                 %.5f (SE %.5f) CI [%.5f, %.5f]; kappa_hat %.4f (SE %.4f)\n",
            mle$est, mle$se, mle$ci[1], mle$ci[2], kap, mle$se_kappa))
cat(sprintf("  SE Hessian %.6f, OPG %.6f, sandwich %.6f (Hessian / OPG = %.4f); survreg max |diff| %.2e\n",
            se_hess, se_opg, se_sand, se_hess / se_opg, svr$max_abs_diff))
cat(sprintf("  H0 kappa = 1: W %.3f, LR %.3f, LM %.3f (df %d); LR from survreg %.3f\n",
            tri$W, tri$LR, tri$LM, tri$df, LR_sv))
cat(sprintf("  AFT scale (truth %.5f): OLS observed %.5f, OLS latent %.5f, survreg %.5f\n",
            aft$truth, aft$ols_observed, aft$ols_latent, aft$mle))
cat("Saved res_cens.rds\n")
