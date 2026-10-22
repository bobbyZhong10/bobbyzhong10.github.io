# ------------------------------------------------------------------
# m_fns.R -- functions for chapter 5 (M-estimation and nonlinear regression).
#
# Convention (notation_v3.md 1.2 and 1.4): an M-estimator MINIMISES the sample mean of
# q(w_i, theta); the score is s_i = grad_theta q_i' (k x 1), A_0 = E[grad^2_theta q] and
# B_0 = E[s s']; Avar sqrt(n)(theta_hat - theta_0) = A_0^-1 B_0 A_0^-1 (Wooldridge 2010,
# Theorem 12.3). For NLS q = (y - m(x, theta))^2 / 2, for Poisson PML q = mu - y log(mu).
#
# Chapter card interface:
#   gauss_newton(y, X, mfun, jac, theta0, tol = 1e-10, maxit = 200) -> list(theta, iter, converged)
#   m_vcov(S, A)                  -> A^-1 B A^-1' / n, B = S'S / n           (Wooldridge 2010 eq. 12.47)
#   two_step_vcov(S2, A2, F2, R1) -> A2^-1 D A2^-1' / n, D = n^-1 sum g_i g_i',
#                                    g_i = s2_i + F2 r1_i                        (eqs. 12.39 to 12.41, 12.61, 12.62)
# Helpers used by 01 to 05:
#   exp_mean(X, b), exp_jac(X, b)          m = exp(X b) and its n x k gradient
#   nls_parts(y, X, b, hessian)            NLS scores S (n x k) and mean Hessian A for an exponential mean
#   poisson_parts(y, X, b)                 Poisson PML scores and mean Hessian for an exponential mean
#   ci95(est, se)                          est -/+ qnorm(0.975) se
#   pseudo_true_exp(m0, sd_x, criterion)   NLS or PPML pseudo-true value theta_* of exp(b0 + b1 x)
#                                          when E[y | x] = m0(x), x ~ N(0, sd_x^2) (quadrature)
#   twostep_fit(d, theta_start)            block M3: OLS first step, NLS second step on the generated
#                                          regressor, naive and corrected variances
#   gen_M1, gen_M2, gen_M2m, gen_M3, gen_M4   the Larkspur data generators (parameters and seeds in 01_dgp.R)
# Only base R (stats) is used here.
# ------------------------------------------------------------------

## =================================================================
## Card functions
## =================================================================

#' Gauss-Newton for NLS (Wooldridge 2010 section 12.7.3) with step halving.
#' mfun(X, b) returns the n-vector m(x_i, b); jac(X, b) the n x k matrix grad_theta m(x_i, b).
#' Iteration j: the Gauss-Newton direction is the OLS coefficient of the residuals y - m on the
#' gradient matrix; the step is halved until the sum of squared residuals (SSR) does not rise by more
#' than rounding error (a relative 1e-12). A direction with max_j |step_j| < sqrt(machine eps) *
#' (1 + max_j |theta_j|) is taken in full. Near the minimum the SSR comparison is at rounding level,
#' and step halving only serves to bring a poor start into the region where Gauss-Newton contracts.
#' Converged when the full direction satisfies max_j |step_j| < tol * (1 + max_j |theta_j|). Stops
#' without convergence when the gradient matrix is rank deficient, when no halved step is acceptable,
#' or after maxit iterations.
gauss_newton <- function(y, X, mfun, jac, theta0, tol = 1e-10, maxit = 200) {
  theta <- as.numeric(theta0)
  ssr <- function(th) sum((y - mfun(X, th))^2)
  s_cur <- ssr(theta)
  if (!is.finite(s_cur)) stop("gauss_newton: the objective is not finite at theta0")
  small <- sqrt(.Machine$double.eps)
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(maxit)) {
    J <- jac(X, theta)
    u <- y - mfun(X, theta)
    step <- as.numeric(qr.coef(qr(J), u))
    if (anyNA(step)) break
    scale <- 1 + max(abs(theta))
    lambda <- 1
    repeat {
      cand <- theta + lambda * step
      s_new <- ssr(cand)
      if (is.finite(s_new) && (s_new <= s_cur * (1 + 1e-12) || max(abs(step)) < small * scale)) break
      lambda <- lambda / 2
      if (lambda < 2^-30) break
    }
    if (!is.finite(s_new) || lambda < 2^-30) break
    theta <- cand
    s_cur <- s_new
    if (max(abs(step)) < tol * (1 + max(abs(theta)))) {
      converged <- TRUE
      break
    }
  }
  list(theta = theta, iter = iter, converged = converged)
}

#' Sandwich variance of an M-estimator: S is the n x k matrix of scores s_i' at theta_hat, A the
#' k x k mean Hessian (eq. 12.42) or mean expected Hessian (eq. 12.44). B_hat = S'S / n (eq. 12.45);
#' returns A^-1 B A^-1' / n (eq. 12.47). The transpose keeps the formula right when A is a numerical,
#' not exactly symmetric, Jacobian; for a symmetric A it is the textbook A^-1 B A^-1 / n.
m_vcov <- function(S, A) {
  S <- as.matrix(S)
  A <- as.matrix(A)
  if (nrow(A) != ncol(A) || ncol(S) != nrow(A)) stop("m_vcov: S is n x k and A is k x k")
  Ai <- solve(A)
  Ai %*% (crossprod(S) / nrow(S)) %*% t(Ai) / nrow(S)
}

#' Variance of a two-step M-estimator (Wooldridge 2010 section 12.4.2 and 12.5.2).
#' S2: n x k second-step scores s_i(theta_hat, gamma_hat); A2: k x k mean Hessian of the second step;
#' F2: k x k_gamma, F2 = n^-1 sum grad_gamma s_i(theta_hat, gamma_hat) (eq. 12.61; Wooldridge's F_o of
#' eq. 12.36, the book's Gamma_0); R1: n x k_gamma, row i is the first-step influence function r_i'
#' with sqrt(n)(gamma_hat - gamma_*) = n^-1/2 sum r_i + o_p(1) (eq. 12.38; the book's psi_i).
#' g_i = s_i + F2 r_i (eq. 12.39), D_hat = n^-1 sum g_i g_i' (eq. 12.62, uncentred as in the book;
#' both s_i and r_i average to zero at the estimates), Avar = A^-1 D A^-1 (eq. 12.41).
#' F2 = 0 (condition 12.37) gives back m_vcov(S2, A2) exactly.
two_step_vcov <- function(S2, A2, F2, R1) {
  S2 <- as.matrix(S2)
  F2 <- as.matrix(F2)
  R1 <- as.matrix(R1)
  if (nrow(R1) != nrow(S2) || nrow(F2) != ncol(S2) || ncol(F2) != ncol(R1))
    stop("two_step_vcov: S2 is n x k, F2 is k x k_gamma and R1 is n x k_gamma")
  m_vcov(S2 + R1 %*% t(F2), A2)
}

## =================================================================
## Exponential-mean pieces
## =================================================================
exp_mean <- function(X, b) exp(drop(X %*% b))
exp_jac <- function(X, b) exp(drop(X %*% b)) * X

#' NLS with m = exp(X b): q_i = (y_i - m_i)^2 / 2, s_i = -grad m_i' u_i (eq. 12.50 up to sign),
#' A = n^-1 sum grad m_i' grad m_i (the expected Hessian, eq. 12.44: valid when the mean is right) or,
#' with hessian = "full", n^-1 sum [grad m_i' grad m_i - grad^2 m_i u_i] (eq. 12.42, the fully robust
#' choice when the mean may be misspecified, White 1981; grad^2 m_i = m_i x_i x_i').
nls_parts <- function(y, X, b, hessian = c("outer", "full")) {
  hessian <- match.arg(hessian)
  m <- exp_mean(X, b)
  J <- m * X
  u <- y - m
  A <- crossprod(J) / length(y)
  if (hessian == "full") A <- A - crossprod(X * (m * u), X) / length(y)
  list(S = -J * u, A = A, m = m, u = u)
}

#' Poisson PML with mu = exp(X b): q_i = mu_i - y_i log(mu_i), s_i = (mu_i - y_i) x_i,
#' A = n^-1 sum mu_i x_i x_i' (the Hessian does not involve y). The model-based (information-equality)
#' variance is A^-1 / n; the sandwich is m_vcov(S, A).
poisson_parts <- function(y, X, b) {
  mu <- exp_mean(X, b)
  list(S = (mu - y) * X, A = crossprod(X * mu, X) / length(y), mu = mu)
}

ci95 <- function(est, se) unname(c(est - qnorm(0.975) * se, est + qnorm(0.975) * se))

#' Pseudo-true value theta_* of the working mean exp(b0 + b1 x) when E[y | x] = m0(x) and
#' x ~ N(0, sd_x^2). criterion = "nls": argmin E[(m0(x) - exp(b0 + b1 x))^2] (White 1981);
#' criterion = "ppml": the root of E[(m0(x) - exp(b0 + b1 x)) (1, x)'] = 0 (Poisson PML).
#' The expectation over x is a Riemann sum on a uniform grid over +/- range sd_x with normal weights
#' (spectrally accurate for these smooth, fast-decaying integrands).
pseudo_true_exp <- function(m0, sd_x, criterion = c("nls", "ppml"), n_grid = 40001, range = 10) {
  criterion <- match.arg(criterion)
  xg <- seq(-range * sd_x, range * sd_x, length.out = n_grid)
  w <- dnorm(xg, 0, sd_x)
  w <- w / sum(w)
  X <- cbind(1, xg)
  y0 <- m0(xg)
  b <- as.numeric(qr.coef(qr(sqrt(w) * X), sqrt(w) * log(y0)))      # start: weighted OLS of log m0 on (1, x)
  if (criterion == "nls") {
    rw <- sqrt(w)
    gn <- gauss_newton(rw * y0, X, function(X, b) rw * exp_mean(X, b), function(X, b) rw * exp_jac(X, b), b, tol = 1e-12)
    if (!gn$converged) stop("pseudo_true_exp: Gauss-Newton did not converge")
    return(gn$theta)
  }
  for (it in 1:200) {                                                  # Newton on the weighted Poisson FOC
    mu <- exp_mean(X, b)
    step <- solve(crossprod(X * (w * mu), X), crossprod(X, w * (y0 - mu)))
    b <- b + as.numeric(step)
    if (max(abs(step)) < 1e-13 * (1 + max(abs(b)))) return(b)
  }
  stop("pseudo_true_exp: Newton did not converge")
}

## =================================================================
## Block M3: two-step estimator with a generated regressor
## =================================================================

#' d: data frame with z1, z2, v, y (gen_M3). Step 1: OLS of v on z = (1, z1, z2) gives gamma_hat and the
#' traffic index t_hat = z'gamma_hat; r_i = (Z'Z/n)^-1 z_i e_i is its influence function (eq. 12.38).
#' Step 2: NLS of y on exp(theta0 + theta1 t_hat) by gauss_newton; S2, A2 from nls_parts (A2: eq. 12.44,
#' the mean is correctly specified in block M3). F2 is eq. 12.61 computed analytically:
#'   s2_i = -m_i u_i x_i with x_i = (1, t_i)', t_i = z_i'gamma, u_i = y_i - m_i, so
#'   grad_gamma s2_i = -theta1 m_i u_i x_i z_i' - m_i u_i (0, z_i')' + theta1 m_i^2 x_i z_i'.
#' V_naive treats t_hat as data (eq. 12.47 with s2); V_corrected is eq. 12.41 via two_step_vcov.
twostep_fit <- function(d, theta_start = c(0, 0.5)) {
  Z <- cbind(1, d$z1, d$z2)
  n <- nrow(Z)
  ZZ <- crossprod(Z) / n
  gamma <- as.numeric(solve(ZZ, crossprod(Z, d$v) / n))
  e <- d$v - drop(Z %*% gamma)
  R1 <- (Z * e) %*% solve(ZZ)
  t_hat <- drop(Z %*% gamma)
  X <- cbind(1, t_hat)
  gn <- gauss_newton(d$y, X, exp_mean, exp_jac, theta_start)
  theta <- gn$theta
  p <- nls_parts(d$y, X, theta)
  mu <- p$m * p$u
  F2 <- (-theta[2] * crossprod(X * mu, Z) - rbind(0, crossprod(mu, Z)) + theta[2] * crossprod(X * p$m^2, Z)) / n
  dimnames(F2) <- NULL
  dimnames(R1) <- NULL
  list(gamma = gamma, theta = theta, t_hat = t_hat, S2 = unname(p$S), A2 = unname(p$A), F2 = F2, R1 = R1,
       V_naive = m_vcov(p$S, p$A), V_corrected = two_step_vcov(p$S, p$A, F2, R1),
       converged = gn$converged, iter = gn$iter)
}

## =================================================================
## Larkspur data generators (parameters and seeds: 01_dgp.R, stored in truth05.rds$params)
## =================================================================

#' M1: merchant cross-section, Spotlight experiment. Category c with probabilities pi; tenure_std ~ N(0, 1);
#' D assigned by complete randomisation (half the merchants). A merchant makes at least one sale in the
#' month under control when U < p0[c] and under Spotlight when U < k p0[c] (same U: monotone extensive
#' margin); log sales given any sale are log(median0[c]) + beta_tenure tenure + sigma eps, plus delta
#' under Spotlight (same eps). Sales are recorded in whole cents; usd = cents / 100.
#' Hence E[y(d) | c, tenure] = exp(alpha_c + tau d + beta_tenure tenure) with tau = log(k) + delta.
gen_M1 <- function(par, n = par$n) {
  category <- sample.int(length(par$pi), n, replace = TRUE, prob = par$pi)
  tenure_std <- rnorm(n)
  D <- sample(rep(0:1, length.out = n))
  U <- runif(n)
  eps <- rnorm(n)
  p0 <- par$p0[category]
  log_pos <- log(par$median0[category]) + par$beta_tenure * tenure_std + par$sigma * eps
  y0c <- round(100 * ifelse(U < p0, exp(log_pos), 0))
  y1c <- round(100 * ifelse(U < par$k * p0, exp(log_pos + par$delta), 0))
  cents <- ifelse(D == 1, y1c, y0c)
  data.frame(merchant = seq_len(n), category = factor(par$categories[category], levels = par$categories),
             tenure_std = tenure_std, D = D, cents = cents, usd = cents / 100,
             y0_cents_true = y0c, y1_cents_true = y1c)
}

#' M2: active merchants. listings ~ round(listings_median exp(sd_log_listings N(0, 1))) (at least 1),
#' x1 = log(listings / listings_median), free_ship ~ Bernoulli(p_free_ship);
#' y (monthly sales, USD thousands) = exp(theta0 + theta1 x1 + theta2 free_ship) eta with eta a
#' mean-one lognormal, sd of log eta = sigma_eta: E[y | x] is exponential, Var(y | x) is proportional to
#' the square of the mean.
gen_M2 <- function(par, n = par$n) {
  listings <- pmax(1, round(par$listings_median * exp(par$sd_log_listings * rnorm(n))))
  x1 <- log(listings / par$listings_median)
  free_ship <- rbinom(n, 1, par$p_free_ship)
  m <- exp(par$theta[1] + par$theta[2] * x1 + par$theta[3] * free_ship)
  y <- m * exp(par$sigma_eta * rnorm(n) - par$sigma_eta^2 / 2)
  data.frame(listings = listings, x1 = x1, free_ship = free_ship, y = y, m_true = m)
}

#' M2m: the misspecified-mean variant. x ~ N(0, sd_x^2) (log catalogue size relative to the median),
#' E[y | x] = exp(a0 + a1 x + a2 x^2) (a2 < 0: diminishing returns), y = E[y | x] eta as in M2.
gen_M2m <- function(par, n = par$n) {
  x <- rnorm(n, 0, par$sd_x)
  m <- exp(par$a0 + par$a1 * x + par$a2 * x^2)
  y <- m * exp(par$sigma_eta * rnorm(n) - par$sigma_eta^2 / 2)
  data.frame(x = x, y = y, m_true = m)
}

#' M3: z1 (search-placement score), z2 (log followers), both N(0, 1); the traffic index
#' t = gamma0 + gamma1 z1 + gamma2 z2; v = t + sigma_e e is one week of log storefront visits;
#' y (monthly sales, USD thousands) = exp(theta0 + theta1 t) eta, eta mean-one lognormal (sigma_eta).
gen_M3 <- function(par, n = par$n) {
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  t <- par$gamma[1] + par$gamma[2] * z1 + par$gamma[3] * z2
  v <- t + par$sigma_e * rnorm(n)
  y <- exp(par$theta[1] + par$theta[2] * t) * exp(par$sigma_eta * rnorm(n) - par$sigma_eta^2 / 2)
  data.frame(z1 = z1, z2 = z2, v = v, y = y, t_true = t)
}

#' M4: x ~ N(0, 1) (standardised seller rating), mu = exp(beta0 + beta1 x), y = monthly orders,
#' negative binomial with mean mu and size `size` (Var = mu + mu^2 / size); size = Inf gives Poisson.
gen_M4 <- function(par, n = par$n, size = par$size) {
  x <- rnorm(n)
  mu <- exp(par$beta[1] + par$beta[2] * x)
  y <- if (is.infinite(size)) rpois(n, mu) else rnbinom(n, size = size, mu = mu)
  data.frame(x = x, y = as.numeric(y), mu_true = mu)
}
