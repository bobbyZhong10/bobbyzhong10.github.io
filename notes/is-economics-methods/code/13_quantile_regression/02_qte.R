# ------------------------------------------------------------------
# 02_qte.R -- blocks Q1 and Q1r: mean effect, marginal QTE, conditional QR (chapter 13)
#   ols      OLS of revenue on D (difference in means), HC2 standard error (Neyman's variance)
#   ols_x    OLS of revenue on D, X and D X (regression adjustment), HC2
#   qte      marginal QTE at taus (difference of arm quantiles), bootstrap that resamples each arm
#            (999 draws, seed 13101); qte_fine: the same on the figure grid 0.05, ..., 0.95
#   qte_log  marginal QTE of log revenue (equivariance: the log of the arm quantiles)
#   cqr_log  QR of log revenue on D, X, D X: the coefficient on D is the conditional QTE at X = 0
#            (correctly specified: Quant_tau(log Y | D, X) is linear in D, X, D X), nid standard errors
#   Q1r      the same estimators on the re-paired world (bootstrap seed 13102)
# Output: res_qte.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("qr_fns.R")
suppressPackageStartupMessages(library(quantreg))
tr <- readRDS("truth13.rds"); P <- tr$params; taus <- P$taus; taus_fine <- P$taus_fine
Q1 <- readRDS("marigold_Q1.rds"); Q1r <- readRDS("marigold_Q1r.rds")
z975 <- qnorm(0.975)

ols_hc2 <- function(fit, term) {
  V <- sandwich::vcovHC(fit, type = "HC2")
  est <- coef(fit)[[term]]; se <- sqrt(V[term, term])
  list(est = est, se = se, ci = c(est - z975 * se, est + z975 * se))
}

estimate_block <- function(d, boot_seed) {
  out <- list()
  out$ols <- ols_hc2(lm(y ~ D, data = d), "D")
  out$ols_x <- ols_hc2(lm(y ~ D * x, data = d), "D")
  set.seed(boot_seed)
  fine <- qte_marginal(d$y, d$D, taus_fine, B = 999)
  out$qte_fine <- fine
  out$qte <- fine[fine$tau %in% taus, ]
  rownames(out$qte) <- NULL
  out$qte_log <- data.frame(tau = taus, est = emp_quantile(log(d$y[d$D == 1]), taus) - emp_quantile(log(d$y[d$D == 0]), taus))
  out$cqr_log <- do.call(rbind, lapply(taus, function(t) {
    fit <- .muffle_nonunique(rq(log(y) ~ D * x, tau = t, data = d))
    cf <- summary(fit, se = "nid")$coefficients
    data.frame(tau = t, est = cf["D", 1], se = cf["D", 2], lo = cf["D", 1] - z975 * cf["D", 2], hi = cf["D", 1] + z975 * cf["D", 2])
  }))
  out$sample <- list(n = nrow(d), n_treat = sum(d$D), mean_y1 = mean(d$y[d$D == 1]), mean_y0 = mean(d$y[d$D == 0]),
                     median_y0 = emp_quantile(d$y[d$D == 0], 0.5),
                     share_lose_true = mean(d$y1_true < d$y0_true))       # simulation truth, not observable
  out
}

q <- estimate_block(Q1, P$seed$Q1_boot)
q$Q1r <- estimate_block(Q1r, P$seed$Q1r_boot)
# Q1 and Q1r share the control arm, so their QTE estimates differ only through the treated arm's draws
q$Q1r$max_abs_qte_diff_vs_Q1 <- max(abs(q$Q1r$qte$est - q$qte$est))

saveRDS(q, "res_qte.rds")
cat("Saved res_qte.rds\n")
cat(sprintf("OLS %.3f [%.3f, %.3f]\n", q$ols$est, q$ols$ci[1], q$ols$ci[2]))
print(cbind(q$qte, truth = tr$Q1$qte))
print(cbind(q$cqr_log, truth = tr$Q1$cqte_log, marginal_log = q$qte_log$est, truth_marg = tr$Q1$qte_log))
print(cbind(q$Q1r$qte, truth = tr$Q1r$qte))
