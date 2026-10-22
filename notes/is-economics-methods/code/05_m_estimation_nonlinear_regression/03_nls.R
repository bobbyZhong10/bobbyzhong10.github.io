# ------------------------------------------------------------------
# 03_nls.R -- blocks M2 and M2m: NLS with an exponential mean.
#
# M2 (mean correctly specified, Var(y | x) proportional to the squared mean):
#   Gauss-Newton (gauss_newton, Wooldridge 2010 section 12.7.3) from theta = 0 against stats::nls;
#   three variance estimators:
#     model    sigma^2 (sum grad m' grad m)^-1, valid only under homoskedasticity (vcov(nls));
#     robust   the heteroskedasticity-robust NLS sandwich, eq. 12.52 = m_vcov(S, A) with
#              A = mean grad m' grad m (eq. 12.44); equal to sandwich::sandwich(nls) and to
#              sandwich::vcovCL(type = "HC0", cadjust = FALSE) (Task 2e: vcovHC() errors on nls);
#     full     the same with the full Hessian (eq. 12.42, eq. 12.48); same limit as robust here.
# M2m (working mean exp(b0 + b1 x), true mean exp(a0 + a1 x + a2 x^2)):
#   NLS converges to the pseudo-true value theta_* (White 1981). Inference on theta_* needs the
#   fully robust sandwich with the full Hessian (eq. 12.48): with a wrong mean, E[H | x] is no
#   longer grad m' grad m, so the semirobust eq. 12.52 is not valid (Wooldridge 2010 p. 416).
#   Poisson PML on the same data converges to a different pseudo-true value.
# Variance convention: HC0 (see 02_logs.R).
#
# Output: res_nls.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("m_fns.R")
source("../../tools/R/chapter_kit.R")   # covers()
suppressPackageStartupMessages({ library(sandwich); library(numDeriv) })

tr <- readRDS("truth05.rds")
M2 <- readRDS("larkspur_M2.rds")
M2m <- readRDS("larkspur_M2m.rds")
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
se_of <- function(V) sqrt(diag(V))

## =================================================================
## M2
## =================================================================
n2 <- nrow(M2)
X2 <- cbind(1, M2$x1, M2$free_ship)
gn <- gauss_newton(M2$y, X2, exp_mean, exp_jac, c(0, 0, 0))
stopifnot(gn$converged)
fit <- nls(y ~ exp(b0 + b1 * x1 + b2 * free_ship), data = M2, start = list(b0 = 0, b1 = 0, b2 = 0),
           control = nls.control(tol = 1e-8, maxiter = 200))
nls_coef <- coef(fit)
p_out <- nls_parts(M2$y, X2, gn$theta)
p_full <- nls_parts(M2$y, X2, gn$theta, hessian = "full")
V_robust <- m_vcov(p_out$S, p_out$A)
V_full <- m_vcov(p_full$S, p_full$A)
V_model <- unname(vcov(fit))
V_pkg <- unname(sandwich::sandwich(fit))
V_cl <- unname(sandwich::vcovCL(fit, cluster = seq_len(n2), type = "HC0", cadjust = FALSE))
stopifnot(rel(V_robust, V_pkg) < 1e-5, rel(V_robust, V_cl) < 1e-5)
Q2 <- function(b) mean((M2$y - exp_mean(X2, b))^2) / 2
stopifnot(rel(p_full$A, numDeriv::hessian(Q2, gn$theta)) < 1e-6)
se_r <- se_of(V_robust)
M2res <- list(
  gn = gn, nls_coef = nls_coef, nls_iter = fit$convInfo$finIter,
  max_abs_diff = max(abs(gn$theta - nls_coef)),
  theta = gn$theta, se_robust = se_r, se_model = se_of(V_model), se_full = se_of(V_full), se_pkg = se_of(V_pkg),
  ci = lapply(1:3, function(j) ci95(gn$theta[j], se_r[j])),
  sigma2 = sum(p_out$u^2) / (n2 - 3),
  # the intercept interval misses theta0 with seed 5102: an expected 5% event, seed not tuned (see 01_dgp.R seeds)
  covers = vapply(1:3, function(j) as.numeric(covers(ci95(gn$theta[j], se_r[j]), tr$M2$theta[j])), numeric(1))
)

## =================================================================
## M2m
## =================================================================
nm <- nrow(M2m)
Xm <- cbind(1, M2m$x)
gm <- gauss_newton(M2m$y, Xm, exp_mean, exp_jac, c(0, 0.5))
stopifnot(gm$converged)
fit_m <- nls(y ~ exp(b0 + b1 * x), data = M2m, start = list(b0 = 0, b1 = 0.5), control = nls.control(tol = 1e-8, maxiter = 200))
stopifnot(max(abs(gm$theta - coef(fit_m))) < 1e-6)
pm_full <- nls_parts(M2m$y, Xm, gm$theta, hessian = "full")
pm_out <- nls_parts(M2m$y, Xm, gm$theta)
Vm_full <- m_vcov(pm_full$S, pm_full$A)
Vm_semi <- m_vcov(pm_out$S, pm_out$A)
stopifnot(rel(Vm_semi, unname(sandwich::sandwich(fit_m))) < 1e-5)
Qm <- function(b) mean((M2m$y - exp_mean(Xm, b))^2) / 2
stopifnot(rel(pm_full$A, numDeriv::hessian(Qm, gm$theta)) < 1e-6)
b1 <- gm$theta[2]
se_full <- sqrt(Vm_full[2, 2])
se_semi <- sqrt(Vm_semi[2, 2])
# Poisson PML on the same data (y is continuous and positive; quasipoisson only uses the mean)
gpm <- glm(y ~ x, family = quasipoisson(), data = M2m, control = glm.control(epsilon = 1e-12, maxit = 100))
pp_est <- unname(coef(gpm)[2])
pp_se <- sqrt(sandwich::sandwich(gpm)[2, 2])
misspec <- list(theta = gm$theta, est = b1, se = se_full, ci = ci95(b1, se_full),
                se_semi = se_semi, ci_semi = ci95(b1, se_semi), se_model = sqrt(vcov(fit_m)[2, 2]),
                iter = gm$iter, covers = as.numeric(covers(ci95(b1, se_full), tr$M2m$pseudo)),
                covers_a1 = as.numeric(covers(ci95(b1, se_full), tr$params$M2m$a1)),
                ppml = list(theta = unname(coef(gpm)), est = pp_est, se = pp_se, ci = ci95(pp_est, pp_se),
                            covers = as.numeric(covers(ci95(pp_est, pp_se), tr$M2m$pseudo_ppml)),
                            covers_nls_pseudo = as.numeric(covers(ci95(pp_est, pp_se), tr$M2m$pseudo))))

nl <- c(M2res, list(misspec = misspec))
saveRDS(nl, "res_nls.rds")

cat(sprintf("M2: Gauss-Newton (%d iterations) %s; nls %s; max difference %.2e\n", gn$iter,
            paste(sprintf("%.5f", gn$theta), collapse = " "), paste(sprintf("%.5f", nls_coef), collapse = " "), M2res$max_abs_diff))
cat(sprintf("M2: SE model %s | robust %s\n", paste(sprintf("%.4f", M2res$se_model), collapse = " "),
            paste(sprintf("%.4f", se_r), collapse = " ")))
cat(sprintf("M2m: b1 = %.4f, fully robust SE %.4f (semirobust %.4f), CI [%.4f, %.4f]; pseudo-true %.4f; PPML %.4f (pseudo-true %.4f)\n",
            b1, se_full, se_semi, misspec$ci[1], misspec$ci[2], tr$M2m$pseudo, pp_est, tr$M2m$pseudo_ppml))
cat("Saved res_nls.rds\n")
