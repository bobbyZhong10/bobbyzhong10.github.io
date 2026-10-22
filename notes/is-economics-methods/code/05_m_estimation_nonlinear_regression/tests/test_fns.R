# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of m_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# The first block is the chapter card's assertions, verbatim; the rest checks the
# helpers that 02 to 05 use against packages (stats::nls, glm, sandwich) and
# against numerical derivatives (numDeriv).
# ------------------------------------------------------------------
source("m_fns.R")

## ---- chapter card (verbatim) ---------------------------------------------------
x <- seq(0, 1, length.out = 200); y <- exp(0.5 + 1.2 * x)                       # noise-free: GN converges exactly
gn <- gauss_newton(y, cbind(1, x), function(X, b) exp(drop(X %*% b)), function(X, b) exp(drop(X %*% b)) * X, c(0, 0))
stopifnot(gn$converged, max(abs(gn$theta - c(0.5, 1.2))) < 1e-8)
set.seed(1); S <- matrix(rnorm(400), 200); A <- diag(2)
stopifnot(max(abs(two_step_vcov(S, A, matrix(0, 2, 1), matrix(rnorm(200), 200)) - m_vcov(S, A))) < 1e-12)  # F2 = 0: the first step has no effect

## ---- gauss_newton: return value, non-convergence flag, same minimiser as stats::nls ----
stopifnot(identical(names(gn), c("theta", "iter", "converged")), gn$iter >= 1, gn$iter < 200,
          is.numeric(gn$theta), length(gn$theta) == 2)
gn1 <- gauss_newton(y, cbind(1, x), exp_mean, exp_jac, c(0, 0), maxit = 1)
stopifnot(!gn1$converged, gn1$iter == 1)
set.seed(2); n <- 500
X <- cbind(1, rnorm(n), rbinom(n, 1, 0.4)); yy <- exp(drop(X %*% c(0.2, 0.5, -0.3))) * rexp(n)
fit <- nls(yy ~ exp(b0 + b1 * x1 + b2 * x2), data = data.frame(yy, x1 = X[, 2], x2 = X[, 3]),
           start = list(b0 = 0, b1 = 0, b2 = 0), control = nls.control(tol = 1e-8, maxiter = 200))
g2 <- gauss_newton(yy, X, exp_mean, exp_jac, c(0, 0, 0))
stopifnot(g2$converged, max(abs(g2$theta - coef(fit))) < 1e-6)

## ---- m_vcov with NLS scores: the heteroskedasticity-robust NLS variance (Wooldridge 2010 eq. 12.52) ----
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
p2 <- nls_parts(yy, X, g2$theta)
stopifnot(all(dim(p2$S) == c(n, 3)), all(dim(p2$A) == c(3, 3)))
V2 <- m_vcov(p2$S, p2$A)
stopifnot(rel(V2, unname(sandwich::sandwich(fit))) < 1e-5,
          rel(V2, unname(sandwich::vcovCL(fit, cluster = seq_len(n), type = "HC0", cadjust = FALSE))) < 1e-5)
# full Hessian (eq. 12.42) = numerical Hessian of the sample objective (1/2n) sum (y - m)^2
Q <- function(b) mean((yy - exp(drop(X %*% b)))^2) / 2
stopifnot(rel(nls_parts(yy, X, g2$theta, hessian = "full")$A, numDeriv::hessian(Q, g2$theta)) < 1e-6)
# mean score at the minimiser is zero
stopifnot(max(abs(colMeans(p2$S))) < 1e-8)

## ---- m_vcov with Poisson scores = sandwich::sandwich(glm); model-based = vcov(glm) ----
set.seed(4); xp <- rnorm(800); Xp <- cbind(1, xp); cnt <- rnbinom(800, size = 2, mu = exp(1 + 0.4 * xp))
gp <- glm(cnt ~ xp, family = poisson(), control = glm.control(epsilon = 1e-12, maxit = 100))
pp <- poisson_parts(cnt, Xp, coef(gp))
stopifnot(rel(m_vcov(pp$S, pp$A), unname(sandwich::sandwich(gp))) < 1e-8,
          rel(solve(pp$A) / 800, unname(vcov(gp))) < 1e-6)

## ---- two_step_vcov = theta block of the stacked (joint) M-estimator sandwich ----
## First step OLS of v on z, second step exponential NLS of y on the generated regressor z'gamma_hat.
## The stacked scores psi_i = (s1_i, s2_i) have a block-triangular mean Jacobian G; the theta block of
## G^-1 Omega G^-1' / n is A2^-1 Var(s2 + F2 r1) A2^-1 / n with A2 = G22, F2 = G21, r1 = -G11^-1 s1
## (Wooldridge 2010 eqs. 12.36 to 12.41).
set.seed(3)
d3 <- gen_M3(list(gamma = c(0, 0.3, 0.2), theta = c(log(1.5), 0.8), sigma_e = 0.6, sigma_eta = 0.35), n = 800)
f3 <- twostep_fit(d3)
Z <- cbind(1, d3$z1, d3$z2); n3 <- nrow(Z)
psi <- function(phi) {
  g <- phi[1:3]; th <- phi[4:5]; t <- drop(Z %*% g); e <- d3$v - t
  m <- exp(th[1] + th[2] * t); u <- d3$y - m
  cbind(-Z * e, -(m * u) * cbind(1, t))
}
phi_hat <- c(f3$gamma, f3$theta); P <- psi(phi_hat)
stopifnot(f3$converged, max(abs(colMeans(P))) < 1e-8)                 # both steps solve their first-order conditions
G <- numDeriv::jacobian(function(p) colMeans(psi(p)), phi_hat)
stopifnot(max(abs(G[1:3, 4:5])) < 1e-8)                                # the first step does not depend on theta
Om <- crossprod(P) / n3; Gi <- solve(G); V_joint <- Gi %*% Om %*% t(Gi) / n3
R1n <- -P[, 1:3] %*% t(solve(G[1:3, 1:3]))
V_two <- two_step_vcov(P[, 4:5], G[4:5, 4:5], G[4:5, 1:3, drop = FALSE], R1n)
stopifnot(rel(V_two, V_joint[4:5, 4:5]) < 1e-9)
# the analytic pieces returned by twostep_fit: F2 (eq. 12.61), r1 (OLS influence function), s2
stopifnot(all(dim(f3$F2) == c(2, 3)), all(dim(f3$R1) == c(n3, 3)),
          rel(f3$F2, G[4:5, 1:3]) < 1e-6, rel(f3$R1, R1n) < 1e-8, max(abs(f3$S2 - P[, 4:5])) < 1e-12)
stopifnot(rel(f3$V_naive, m_vcov(f3$S2, f3$A2)) < 1e-14,
          rel(f3$V_corrected, two_step_vcov(f3$S2, f3$A2, f3$F2, f3$R1)) < 1e-14)

## ---- pseudo_true_exp: a correctly specified exponential mean returns its own coefficients ----
for (crit in c("nls", "ppml")) {
  pt <- pseudo_true_exp(function(x) exp(0.3 + 0.5 * x), sd_x = 1, criterion = crit)
  stopifnot(max(abs(pt - c(0.3, 0.5))) < 1e-8)
}

## ---- ci95 ----
stopifnot(isTRUE(all.equal(ci95(1, 0.5), c(1 - qnorm(0.975) * 0.5, 1 + qnorm(0.975) * 0.5))))
