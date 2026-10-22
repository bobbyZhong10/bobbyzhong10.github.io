# ------------------------------------------------------------------
# 03_qmle.R -- block C2: a Weibull fitted to a log-logistic mixture (QMLE).
#
# The Weibull PH model with x = (1, fee_usd, referred) is misspecified: the truth is a mixture of
# deal-seekers and regulars with log-logistic tenure. The QMLE theta_hat converges to the pseudo-true
# value theta_* (truth06$C2: the Weibull MLE on 2,000,000 draws; quadrature check in 01_dgp.R), and
# sqrt(n)(theta_hat - theta_*) -> N(0, A_0^-1 B_0 A_0^-1) with A_0 != B_0 (Wooldridge 2010 section 13.11,
# sandwich estimate eq. 13.71; White 1982). The quantity of interest is the pseudo-true shape kappa_*
# (truth06$C2$pseudo): a single Weibull summarises two segments by a hazard that falls with tenure.
# Three standard errors for kappa_hat (delta method from log kappa: se(kappa) = kappa se(log kappa)):
#   model     (A_hat^-1 / n)  the usual MLE variance, right only if A_0 = B_0;
#   opg       (B_hat^-1 / n)  the outer product of the scores (BHHH), also right only if A_0 = B_0;
#   sandwich  A_hat^-1 B_hat A_hat^-1 / n.
# A_hat and B_hat by info_hat (numDeriv), checked against the analytic pieces.
# Also: the same three SEs for the fee coefficient beta_fee,*; the Vuong (1989) test of the Weibull
# against a single log-logistic AFT model (survreg dist = "loglogistic"; both have four parameters,
# both are misspecified), Wooldridge 2010 section 13.11.2.
# Monte Carlo: R_mc = 2000 samples of n = 8000 (seed C2_mc), QMLE by Newton from theta_*; coverage of
# kappa_hat -/+ qnorm(0.975) se for each SE, against kappa_*; the same for beta_fee.
#
# Output: res_qmle.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("ml_fns.R")
suppressPackageStartupMessages(library(survival))

tr <- readRDS("truth06.rds")
P2 <- tr$params$C2
C2 <- readRDS("sorrel_C2.rds")
X <- cbind(1, C2$fee_usd, C2$referred)
t <- C2$t_months
d <- C2$delta
n <- nrow(X)
z <- qnorm(0.975)
k_star <- tr$C2$pseudo
b_star <- tr$C2$pseudo_fee

## ---- the sample --------------------------------------------------------------------------------
fm <- fit_weibull(t, d, X)
stopifnot(fm$converged)
p <- fm$par
kap <- exp(p[4])
ih <- info_hat(p, nll_weibull_cens, t = t, delta = d, X = X)
stopifnot(max(abs(ih$A - fm$A)) / max(abs(fm$A)) < 1e-5, max(abs(ih$B - fm$B)) / max(abs(fm$B)) < 1e-5)
Ai <- solve(ih$A)
V_h <- Ai / n
V_o <- solve(ih$B) / n
V_s <- Ai %*% ih$B %*% Ai / n
se3 <- function(i, scale = 1) scale * sqrt(c(model = V_h[i, i], opg = V_o[i, i], sandwich = V_s[i, i]))
sk <- se3(4, kap)
sb <- se3(2)
# seed C2 6102 (fixed before results, ruling 25): this interval misses kappa_* (kappa_hat 2.14 sandwich SEs below);
# Ruling 29 tests the Monte Carlo below instead (mean of kappa_hat within 0.003 of kappa_*, sandwich coverage >= 0.93)
qmle <- list(est = kap, se = sk[["sandwich"]], ci = ci95(kap, sk[["sandwich"]]), ci_model = ci95(kap, sk[["model"]]),
             ci_opg = ci95(kap, sk[["opg"]]), par = p, log_kappa = p[4],
             se_model_vec = sqrt(diag(V_h)), se_opg_vec = sqrt(diag(V_o)), se_sand_vec = sqrt(diag(V_s)),
             nll = fm$nll, iter = fm$iter)
fee <- list(est = p[2], se_model = sb[["model"]], se_opg = sb[["opg"]], se_sandwich = sb[["sandwich"]],
            ci = ci95(p[2], sb[["sandwich"]]), ci_model = ci95(p[2], sb[["model"]]))
covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])

## ---- Vuong: Weibull against a single log-logistic ------------------------------------------------
sl <- survreg(Surv(t_months, delta) ~ fee_usd + referred, data = C2, dist = "loglogistic",
              control = survreg.control(rel.tolerance = 1e-12, maxiter = 100))
ll_ll <- loglogistic_ll_obs(sl, t, d, X)
ll_wb <- -nll_weibull_cens(p, t, d, X, by_obs = TRUE)
stopifnot(abs(sum(ll_ll) - sl$loglik[2]) < 1e-6)
vu <- vuong(ll_wb, ll_ll)
vuong_res <- c(vu, list(loglik_weibull = sum(ll_wb), loglik_loglogistic = sum(ll_ll),
                        loglogistic_coef = unname(coef(sl)), loglogistic_scale = sl$scale))

## ---- Monte Carlo -----------------------------------------------------------------------------------
one_rep <- function() {
  dd <- gen_C2(P2)
  Xr <- cbind(1, dd$fee_usd, dd$referred)
  f <- fit_weibull(dd$t_months, dd$delta, Xr, start = tr$C2$pseudo_vec, engine = "newton")
  k <- exp(f$par[4])
  c(kappa = k, se_model = k * sqrt(f$V_hess[4, 4]), se_opg = k * sqrt(f$V_opg[4, 4]), se_sand = k * sqrt(f$V_sand[4, 4]),
    fee = f$par[2], fee_se_model = sqrt(f$V_hess[2, 2]), fee_se_opg = sqrt(f$V_opg[2, 2]), fee_se_sand = sqrt(f$V_sand[2, 2]),
    converged = f$converged)
}
t0 <- proc.time()[["elapsed"]]
set.seed(tr$params$seed$C2_mc)
sim <- t(vapply(seq_len(P2$R_mc), function(r) one_rep(), numeric(9)))
mc_seconds <- proc.time()[["elapsed"]] - t0
ok <- sim[, "converged"] == 1
stopifnot(all(ok))
cov <- function(est, se, v) mean(abs(sim[, est] - v) <= z * sim[, se])
mc <- list(R = nrow(sim), n_converged = sum(ok), n = P2$n,
           cover_model = cov("kappa", "se_model", k_star), cover_opg = cov("kappa", "se_opg", k_star),
           cover_sandwich = cov("kappa", "se_sand", k_star),
           mean_est = mean(sim[, "kappa"]), sd_est = sd(sim[, "kappa"]),
           mean_se_model = mean(sim[, "se_model"]), mean_se_opg = mean(sim[, "se_opg"]), mean_se_sandwich = mean(sim[, "se_sand"]),
           share_ratio_gt_0p3 = mean(abs(sim[, "se_model"] / sim[, "se_opg"] - 1) > 0.30),
           mean_ratio_model_opg = mean(sim[, "se_model"] / sim[, "se_opg"]),
           fee = list(cover_model = cov("fee", "fee_se_model", b_star), cover_opg = cov("fee", "fee_se_opg", b_star),
                      cover_sandwich = cov("fee", "fee_se_sand", b_star),
                      mean_est = mean(sim[, "fee"]), sd_est = sd(sim[, "fee"]),
                      mean_se_model = mean(sim[, "fee_se_model"]), mean_se_opg = mean(sim[, "fee_se_opg"]),
                      mean_se_sandwich = mean(sim[, "fee_se_sand"])))

res <- list(n = n, n_exit = sum(d), share_cens = 1 - mean(d), share_deal = mean(C2$deal_true),
            qmle = qmle, se_hess = sk[["model"]], se_opg = sk[["opg"]], se_sand = sk[["sandwich"]],
            ratio_hess_opg = sk[["model"]] / sk[["opg"]], ratio_sand_hess = sk[["sandwich"]] / sk[["model"]],
            covers_model = covers_num(qmle$ci_model, k_star), covers_opg = covers_num(qmle$ci_opg, k_star),
            covers_sandwich = covers_num(qmle$ci, k_star),
            fee = c(fee, list(covers_sandwich = covers_num(fee$ci, b_star), covers_model = covers_num(fee$ci_model, b_star))),
            A = ih$A, B = ih$B, vuong = vuong_res, mc = mc)   # run time is printed, not saved (byte-reproducible .rds)
saveRDS(res, "res_qmle.rds")

cat(sprintf("C2: censored %.3f, deal-seekers %.3f. QMLE kappa_hat %.4f (pseudo-true %.4f); beta_fee %.5f (pseudo-true %.5f)\n",
            res$share_cens, res$share_deal, kap, k_star, p[2], b_star))
cat(sprintf("  SE of kappa_hat: model %.5f, OPG %.5f, sandwich %.5f (model / OPG = %.3f); sandwich CI [%.4f, %.4f]\n",
            sk[["model"]], sk[["opg"]], sk[["sandwich"]], res$ratio_hess_opg, qmle$ci[1], qmle$ci[2]))
cat(sprintf("  SE of beta_fee: model %.6f, OPG %.6f, sandwich %.6f\n", sb[["model"]], sb[["opg"]], sb[["sandwich"]]))
cat(sprintf("  Vuong (Weibull vs log-logistic): V_n = %.3f, p = %.4f; log likelihoods %.2f / %.2f\n",
            vu$stat, vu$p_two_sided, sum(ll_wb), sum(ll_ll)))
cat(sprintf("C2 Monte Carlo (%d reps, %.1f s): coverage for kappa_*: model %.4f, OPG %.4f, sandwich %.4f; sd %.5f, mean SE %.5f / %.5f / %.5f\n",
            mc$R, mc_seconds, mc$cover_model, mc$cover_opg, mc$cover_sandwich, mc$sd_est, mc$mean_se_model, mc$mean_se_opg, mc$mean_se_sandwich))
cat(sprintf("  coverage for beta_fee,*: model %.4f, OPG %.4f, sandwich %.4f\n", mc$fee$cover_model, mc$fee$cover_opg, mc$fee$cover_sandwich))
cat("Saved res_qmle.rds\n")
