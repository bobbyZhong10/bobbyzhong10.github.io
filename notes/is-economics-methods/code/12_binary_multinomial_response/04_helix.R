# ------------------------------------------------------------------
# 04_helix.R -- block B1 (inherited Helix adoption data): res_helix.rds, the numbers of the legacy
# scripts 02_lpm_naive.R (printed only) and 03_mle_logit.R (res_mle.rds) under the field names of the
# chapter card, plus the v3 additions on the same data (no random numbers):
#   lpm      OLS of D on x (Y_i on X_i in the text): gamma_0, gamma_1, SEs, fitted-value range and shares
#   logit    hand-rolled Newton logit (legacy 03), its iterations and distance to glm, Fisher and sandwich
#            SEs, the AME with a standard error (delta method plus the sampling variance of X), the PEA
#   probit   probit on the same data: coefficients, APE, PEA; logit / probit coefficient ratio
#   consistency   the legacy Monte Carlo table (03_mle_logit.R, set.seed(778))
# The last part prints the legacy reconciliation table (the logic of the legacy 04_summary.R, unchanged).
# Reads helix_data.rds, oracle.rds, res_mle.rds, truth12.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("bin_fns.R")
suppressPackageStartupMessages(library(sandwich))
dat <- readRDS("helix_data.rds"); oracle <- readRDS("oracle.rds"); mle <- readRDS("res_mle.rds")
tr <- readRDS("truth12.rds")

## ---- LPM (legacy 02_lpm_naive.R, the same lm call) ----------------------------------------------------
lpm <- lm(D ~ x, data = dat)
b <- coef(lpm); phat <- fitted(lpm)
lpm_res <- list(intercept = unname(b[1]), slope = unname(b[2]),
                se = summary(lpm)$coefficients[2, 2], se_intercept = summary(lpm)$coefficients[1, 2],
                se_hc1 = sqrt(vcovHC(lpm, type = "HC1")[2, 2]), se_hc0 = sqrt(vcovHC(lpm, type = "HC0")[2, 2]),
                min_phat = min(phat), max_phat = max(phat), share_out = mean(phat < 0 | phat > 1),
                share_above = mean(phat > 1), share_below = mean(phat < 0),
                n_above = sum(phat > 1), n_below = sum(phat < 0))
lpm_res$ci_hc1 <- ci95(lpm_res$slope, lpm_res$se_hc1)

## ---- logit MLE (legacy 03_mle_logit.R; the Newton iterations are re-run without random numbers) ------
X <- cbind(1, dat$x); y <- dat$D
newton <- function(X, y, tol = 1e-10, maxit = 100) {         # the legacy newton_logit, with the step size kept
  b <- rep(0, ncol(X)); steps <- numeric(0)
  for (it in 1:maxit) {
    F <- plogis(drop(X %*% b)); score <- crossprod(X, y - F); H <- -crossprod(X, X * (F * (1 - F)))
    step <- solve(H, score); b <- b - step; steps <- c(steps, max(abs(step)))
    if (max(abs(step)) < tol) break
  }
  list(coef = drop(b), iter = it, steps = steps)
}
nt <- newton(X, y)
g <- glm(D ~ x, data = dat, family = binomial())
stopifnot(max(abs(nt$coef - drop(mle$coef))) < 1e-12)          # = the legacy res_mle.rds
bl <- drop(mle$coef)
ame_inf <- avg_effect_inference(function(t) dlogis(drop(X %*% t)) * t[2], bl, solve(crossprod(X, X * (dlogis(drop(X %*% bl))))))
logit_res <- list(beta0 = bl[1], beta1 = bl[2], se0_fisher = mle$se_fisher[1], se1_fisher = mle$se_fisher[2],
                  se0_sand = mle$se_sand[1], se1_sand = mle$se_sand[2], ame = mle$ame,
                  iter = nt$iter, glm_maxdiff = max(abs(nt$coef - coef(g))),
                  se1_glm = sqrt(vcov(g)[2, 2]), se1_hc0 = sqrt(vcovHC(g, type = "HC0")[2, 2]),
                  sand_over_fisher = mle$se_sand[2] / mle$se_fisher[2],
                  ame_se = ame_inf$se, ame_se_delta = ame_inf$se_delta, ame_ci = ame_inf$ci,
                  pea = bl[2] * dlogis(bl[1] + bl[2] * mean(dat$x)), mean_x = mean(dat$x),
                  ci1 = ci95(bl[2], mle$se_fisher[2]), loglik = as.numeric(logLik(g)))
stopifnot(abs(ame_inf$est - mle$ame) < 1e-12)

## ---- probit on the same data: coefficients differ by a scale factor, APEs agree --------------------
gp <- glm(D ~ x, data = dat, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
bp <- coef(gp)
ape_p <- avg_effect_inference(function(t) dnorm(drop(X %*% t)) * t[2], bp, vcov(gp))
probit_res <- list(beta0 = unname(bp[1]), beta1 = unname(bp[2]), se1 = sqrt(vcov(gp)[2, 2]),
                   ape = ape_p$est, ape_se = ape_p$se, ape_ci = ape_p$ci,
                   pea = unname(bp[2] * dnorm(bp[1] + bp[2] * mean(dat$x))),
                   ratio_logit_probit = unname(bl[2] / bp[2]), loglik = as.numeric(logLik(gp)))

h <- list(lpm = lpm_res, logit = logit_res, probit = probit_res, consistency = mle$consistency,
          sample = list(n = nrow(dat), adoption_rate = mean(dat$D), mean_x = mean(dat$x), sd_x = sd(dat$x)))
saveRDS(h, "res_helix.rds")
cat(sprintf("LPM slope %.4f (SE %.4f, HC1 %.4f), fitted range [%.4f, %.4f], outside [0, 1] %.4f\n",
            lpm_res$slope, lpm_res$se, lpm_res$se_hc1, lpm_res$min_phat, lpm_res$max_phat, lpm_res$share_out))
cat(sprintf("logit beta1 %.4f, AME %.4f (SE %.4f), PEA %.4f; probit beta1 %.4f, APE %.4f; ratio %.3f\n",
            logit_res$beta1, logit_res$ame, logit_res$ame_se, logit_res$pea, probit_res$beta1, probit_res$ape,
            probit_res$ratio_logit_probit))
cat("Saved res_helix.rds\n\n")

## ---- legacy 04_summary.R (reconciliation print), logic unchanged ----------------------------------
line <- function(...) cat(sprintf(...), "\n")
cat("================ CHAPTER 12 NUMBER RECONCILIATION ================\n\n")

cat("--- adoption / logit MLE ---\n")
line("logit MLE  b0 = %.4f  b1 = %.4f   (truth -0.40 / 1.10)", mle$coef[1], mle$coef[2])
line("AME (logit) = %.4f   (oracle %.4f)", mle$ame, oracle$ame_logit)
line("b1 Fisher SE = %.4f   sandwich SE = %.4f  (agree: correct spec)",
     mle$se_fisher[2], mle$se_sand[2])
cat("consistency of b1_hat (mean / sd) across n = 100..102400:\n")
print(round(mle$consistency[c("mean","sd"),], 4))
cat("\n================================================================\n")
