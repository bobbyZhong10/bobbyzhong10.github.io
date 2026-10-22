# ------------------------------------------------------------------
# 04_grouped.R -- block C3: grouped durations, pooled cloglog = grouped-duration MLE.
#
# Data: sorrel_C3.rds, C1's subscribers as monthly person-period rows for months 1 to 24 of tenure.
# Two estimators of (eta_1, ..., eta_24, beta_fee, beta_ref):
#   grouped_mle     Newton-Raphson on the interval form of the likelihood, Wooldridge 2010 eq. (22.45),
#                   with alpha_m(x) = exp(-exp(eta_m + x'beta)), eq. (22.47) with
#                   eta_m = log(lambda_m (a_m - a_{m-1})) (grouped_mle in ml_fns.R; it never forms the rows).
#   cloglog_pooled  glm(y ~ 0 + factor(m) + fee_usd + referred, binomial("cloglog")) on the person-period
#                   rows: a pooled binary-response (partial) MLE, Wooldridge 2010 section 13.8.
# Why they are the same estimator (the derivation behind the identity; Wooldridge never says "cloglog"):
# take Wooldridge's convention: m_i is the interval in which spell i ends, by an observed exit
# (delta_i = 1) or by censoring (delta_i = 0; the spell is then known to survive intervals 1, ..., m_i - 1).
# Spell i has person-period rows m = 1, ..., m_i - 1 + delta_i, with y_im = delta_i 1{m = m_i}. By (22.45)
# and (22.47), with alpha_m(x) = exp(-exp(eta_m + x'beta)), the log likelihood of spell i is
#   sum_{m < m_i} log alpha_m(x_i) + delta_i log[1 - alpha_m_i(x_i)]
#   = sum_{m = 1}^{m_i - 1 + delta_i} { y_im log[1 - alpha_m(x_i)] + (1 - y_im) log alpha_m(x_i) } = sum_m ell_im,
# for both censored spells (rows 1, ..., m_i - 1, all y_im = 0) and uncensored ones (rows 1, ..., m_i, y_im = 1
# in the last). ell_im (the registry's ℓ_im) is the Bernoulli log likelihood of row (i, m), success probability
# 1 - exp(-exp(eta_m + x_i'beta)), the inverse complementary log-log link (Prentice and Gloeckler 1978).
# (In ml_fns.R, grouped_data's mi is the number of rows, m_i - 1 + delta_i: the exit interval for an
# observed exit, the number of intervals survived completely for a censored spell.)
# Same objective, same maximiser: the two vectors agree to the optimisers' tolerance.
# Because each row's density is conditional on survival to the start of its interval, the pooled
# likelihood is the full likelihood (the scores of a spell's rows are uncorrelated: the "dynamically
# complete" case of Wooldridge 2010 section 13.8.3), so the usual and the subscriber-clustered SEs agree
# up to sampling error. Both are reported; no ordering is tested (controller decision).
# The continuous Weibull MLE of C1 estimates the same beta: with a PH hazard, grouping with one free
# eta_m per interval is exact, eta_m = beta_0 + log(m^kappa - (m - 1)^kappa) (truth06$C3$eta).
#
# Output: res_grouped.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("ml_fns.R")
suppressPackageStartupMessages(library(sandwich))

tr <- readRDS("truth06.rds")
P3 <- tr$params$C3
C1 <- readRDS("sorrel_C1.rds")
C3 <- readRDS("sorrel_C3.rds")
c1 <- readRDS("res_cens.rds")
M <- P3$M

## ---- grouped-duration MLE (interval form) --------------------------------------------------------
gd <- grouped_data(C1$t_months, C1$delta, P3$cuts)
Xg <- cbind(fee_usd = C1$fee_usd, referred = C1$referred)
gm <- grouped_mle(gd$mi, gd$yi, Xg, M)
stopifnot(gm$converged, max(abs(gm$grad)) < 1e-6)
nm <- c(sprintf("eta_%02d", seq_len(M)), "fee_usd", "referred")
grouped_mle_par <- setNames(gm$par, nm)

## ---- pooled cloglog on the person-period rows --------------------------------------------------------
gl <- glm(y ~ 0 + factor(m) + fee_usd + referred, family = binomial(link = "cloglog"), data = C3,
          control = glm.control(epsilon = 1e-14, maxit = 500))
stopifnot(gl$converged)
cloglog_par <- setNames(unname(coef(gl)), nm)
max_abs_diff <- max(abs(cloglog_par - grouped_mle_par))
# the same objective: minus the cloglog log likelihood at the glm estimates = the grouped nll there
nll_gap <- grouped_nll(unname(cloglog_par), gd$mi, gd$yi, Xg) + as.numeric(logLik(gl))

## ---- standard errors --------------------------------------------------------------------------------
V_glm <- unname(vcov(gl))                                         # expected information (IRLS weights)
V_cl <- unname(vcovCL(gl, cluster = ~id, type = "HC0", cadjust = FALSE))  # clustered by subscriber
V_obs <- gm$V                                                     # observed information of (22.45)
k <- M + 1:2
se <- function(V) sqrt(diag(V))[k]
fee <- list(est = grouped_mle_par[["fee_usd"]], se = se(V_obs)[1], ci = ci95(grouped_mle_par[["fee_usd"]], se(V_obs)[1]),
            se_glm = se(V_glm)[1], se_cluster = se(V_cl)[1])
ref <- list(est = grouped_mle_par[["referred"]], se = se(V_obs)[2], ci = ci95(grouped_mle_par[["referred"]], se(V_obs)[2]),
            se_glm = se(V_glm)[2], se_cluster = se(V_cl)[2])
covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])

res <- list(cloglog_pooled = cloglog_par, grouped_mle = grouped_mle_par, max_abs_diff = max_abs_diff, nll_gap = nll_gap,
            iter_grouped = gm$iter, iter_glm = gl$iter, n_rows = nrow(C3), n_exits = sum(C3$y), n_subscribers = nrow(C1),
            exits_by_month = as.numeric(tapply(C3$y, C3$m, sum)), at_risk_by_month = as.numeric(tapply(C3$y, C3$m, length)),
            fee = c(fee, list(covers_truth = covers_num(fee$ci, tr$C3$beta[1]))),
            ref = c(ref, list(covers_truth = covers_num(ref$ci, tr$C3$beta[2]))),
            eta = grouped_mle_par[seq_len(M)], eta_se = sqrt(diag(V_obs))[seq_len(M)], eta_true = tr$C3$eta,
            ratio_cluster_glm_fee = se(V_cl)[1] / se(V_glm)[1], ratio_cluster_glm_ref = se(V_cl)[2] / se(V_glm)[2],
            fee_minus_continuous = grouped_mle_par[["fee_usd"]] - c1$mle$est)
saveRDS(res, "res_grouped.rds")

cat(sprintf("C3: %d person-period rows, %d exits; grouped MLE (%d Newton steps) vs pooled cloglog (%d IRLS steps): max |diff| %.2e\n",
            res$n_rows, res$n_exits, gm$iter, gl$iter, max_abs_diff))
cat(sprintf("  beta_fee %.5f (SE observed %.5f, glm %.5f, clustered %.5f), CI [%.5f, %.5f]; continuous Weibull MLE %.5f; truth %.4f\n",
            fee$est, fee$se, fee$se_glm, fee$se_cluster, fee$ci[1], fee$ci[2], c1$mle$est, tr$C3$beta[1]))
cat(sprintf("  beta_ref %.4f (SE %.4f, clustered %.4f); truth %.2f\n", ref$est, ref$se, ref$se_cluster, tr$C3$beta[2]))
cat("Saved res_grouped.rds\n")
