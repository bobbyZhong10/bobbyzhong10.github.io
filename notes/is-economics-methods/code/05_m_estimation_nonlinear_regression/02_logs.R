# ------------------------------------------------------------------
# 02_logs.R -- block M1: log(1 + y) against Poisson PML on the Spotlight experiment.
#
# Every regression has the same right-hand side: D, category dummies, tenure_std.
#   log(1 + a y): OLS of log(1 + a * y_usd) on the regressors for a unit multiplier a
#                 (a = 0.001: thousands of USD, 1: USD, 100: cents) and on a grid of a;
#                 the coefficient on D moves with a (Chen and Roth 2024, Proposition 1).
#   asinh(a y):   the same with the inverse hyperbolic sine (USD and cents).
#   PPML:         glm(family = quasipoisson) on USD and on cents; the coefficient on D estimates
#                 log(E[y(1)] / E[y(0)]) whatever the unit (Gourieroux, Monfort and Trognon 1984).
#   Margins:      linear probability model for 1{y > 0}; OLS of log y among merchants with y > 0.
#                 The seller-only comparison is a comparison of selected groups and is not a causal
#                 effect in general. In THIS DGP it identifies delta = 0.10 exactly: the selling
#                 threshold U is independent of the intensive-margin shock eps (and of tenure), and
#                 Spotlight scales p0[c] by the same k in every category, so sellers under Spotlight
#                 and under control have the same distribution of (category, tenure, eps). There is
#                 no selection bias here to point at (numbers: logs_logpos_*, logs_logpos_covers_delta).
#
# Variance convention for the whole chapter: HC0, i.e. the sandwich A^-1 B A^-1 / n with no
# small-sample factor (sandwich::sandwich(), equal to vcovHC(type = "HC0") for lm and glm, and to
# m_vcov() on the hand-coded scores). fixest::fepois(vcov = "hetero") applies a small-sample factor;
# it is run once below as a cross-check of the point estimate, and its standard error is recorded
# only to show the size of that factor.
#
# Output: res_logs.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("m_fns.R")
suppressPackageStartupMessages({ library(sandwich); library(fixest) })

tr <- readRDS("truth05.rds")
M1 <- readRDS("larkspur_M1.rds")
n <- nrow(M1)

ols_D <- function(resp, data = M1) {
  d <- data
  d$resp <- resp
  f <- lm(resp ~ D + category + tenure_std, data = d)
  est <- unname(coef(f)[["D"]])
  se <- sqrt(sandwich::sandwich(f)["D", "D"])
  list(est = est, se = se, ci = ci95(est, se))
}
ppml_D <- function(resp, rhs = "D + category + tenure_std") {
  d <- M1
  d$resp <- resp
  f <- glm(as.formula(paste("resp ~", rhs)), family = quasipoisson(), data = d,
           control = glm.control(epsilon = 1e-12, maxit = 100))
  stopifnot(f$converged)
  est <- unname(coef(f)[["D"]])
  V <- sandwich::sandwich(f)
  se <- sqrt(V["D", "D"])
  list(est = est, se = se, ci = ci95(est, se), fit = f)
}

## ---- log(1 + a y) and asinh ---------------------------------------------------------
l_kusd <- ols_D(log1p(M1$usd / 1000))
l_usd <- ols_D(log1p(M1$usd))
l_cents <- ols_D(log1p(M1$cents))
a_usd <- ols_D(asinh(M1$usd))
a_cents <- ols_D(asinh(M1$cents))
sw <- tr$M1$pop$sweep
sweep <- do.call(rbind, lapply(seq_len(nrow(sw)), function(j) {
  r <- ols_D(log1p(sw$a[j] * M1$usd))
  data.frame(a = sw$a[j], est = r$est, se = r$se, lo = r$ci[1], hi = r$ci[2], ate_pop = sw$ate_log1p[j])
}))
stopifnot(abs(sweep$est[abs(sweep$a - 1) < 1e-12] - l_usd$est) < 1e-12,
          abs(sweep$est[abs(sweep$a - 100) < 1e-9] - l_cents$est) < 1e-9)

## ---- PPML -------------------------------------------------------------------------------
p_usd <- ppml_D(M1$usd)
p_cents <- ppml_D(M1$cents)
p_raw <- ppml_D(M1$usd, rhs = "D")
stopifnot(abs(p_raw$est - log(mean(M1$usd[M1$D == 1]) / mean(M1$usd[M1$D == 0]))) < 1e-10)   # saturated: log of the ratio of means
# hand-coded scores and Hessian (q = mu - y log mu) give the same sandwich as sandwich::sandwich(glm)
Xd <- model.matrix(p_usd$fit)
pp <- poisson_parts(M1$usd, Xd, coef(p_usd$fit))
V_hand <- m_vcov(pp$S, pp$A)
check_mvcov <- max(abs(V_hand - unname(sandwich::sandwich(p_usd$fit)))) / max(abs(V_hand))
stopifnot(check_mvcov < 1e-6)     # glm keeps the IRLS weights of its last iteration: agreement to about 1e-7
# the Poisson model-based variance A^-1 / n (Var(y | x) = E[y | x]), for comparison only
se_model_usd <- sqrt(solve(pp$A)[2, 2] / n)
# fixest cross-check: category fixed effects instead of dummies; same point estimate
fe <- fixest::fepois(usd ~ D + tenure_std | category, data = M1, vcov = "hetero", glm.tol = 1e-10, fixef.tol = 1e-10)
fe_est <- unname(coef(fe)[["D"]])
fe_se <- unname(fixest::se(fe)[["D"]])
stopifnot(abs(fe_est - p_usd$est) < 1e-6)

## ---- margins ----------------------------------------------------------------------------
ext <- ols_D(as.numeric(M1$usd > 0))
# sellers only: E[log y | y > 0, D = 1, c, tenure] - E[log y | y > 0, D = 0, c, tenure] = delta in this DGP
# (U independent of eps, common k), so lp$est estimates delta = 0.10 without selection bias
pos <- M1[M1$usd > 0, ]
lp <- ols_D(log(pos$usd), data = pos)

lg <- list(
  log1p = list(usd = l_usd$est, cents = l_cents$est, kusd = l_kusd$est,
               se_usd = l_usd$se, se_cents = l_cents$se, se_kusd = l_kusd$se,
               ci_usd = l_usd$ci, ci_cents = l_cents$ci, ci_kusd = l_kusd$ci),
  asinh = list(usd = a_usd$est, cents = a_cents$est, se_usd = a_usd$se, se_cents = a_cents$se),
  ppml = list(usd = list(est = p_usd$est, se = p_usd$se, ci = p_usd$ci, se_model = se_model_usd,
                         pct = exp(p_usd$est) - 1, pct_ci = exp(p_usd$ci) - 1),
              cents = list(est = p_cents$est, se = p_cents$se, ci = p_cents$ci),
              raw = list(est = p_raw$est, se = p_raw$se, ci = p_raw$ci),
              fepois = list(est = fe_est, se_hetero = fe_se, se_ratio = fe_se / p_usd$se),
              check_mvcov = check_mvcov),
  ext = ext,
  logpos = list(est = lp$est, se = lp$se, ci = lp$ci, n = nrow(pos)),
  sweep = sweep,
  sample = list(n = n, n_treated = sum(M1$D), share_zero0 = mean(M1$usd[M1$D == 0] == 0),
                share_zero1 = mean(M1$usd[M1$D == 1] == 0), mean_usd0 = mean(M1$usd[M1$D == 0]),
                mean_usd1 = mean(M1$usd[M1$D == 1]), mean_pos_usd = mean(pos$usd), median_pos_usd = median(pos$usd),
                p99_usd = unname(quantile(M1$usd, 0.99)), max_usd = max(M1$usd))
)
saveRDS(lg, "res_logs.rds")

cat(sprintf("log(1 + y) coefficient on D: %.4f (USD thousands) %.4f (USD) %.4f (cents)\n", l_kusd$est, l_usd$est, l_cents$est))
cat(sprintf("asinh coefficient on D: %.4f (USD) %.4f (cents)\n", a_usd$est, a_cents$est))
cat(sprintf("PPML: %.6f (USD) %.6f (cents), HC0 SE %.4f, 95%% CI [%.4f, %.4f]; truth %.4f; fepois %.6f (SE %.4f)\n",
            p_usd$est, p_cents$est, p_usd$se, p_usd$ci[1], p_usd$ci[2], tr$M1$log_ratio, fe_est, fe_se))
cat(sprintf("Extensive margin %.4f (SE %.4f); log y among sellers %.4f (SE %.4f; delta = %.2f)\n", ext$est, ext$se, lp$est, lp$se,
            tr$params$M1$delta))
cat("Saved res_logs.rds\n")
