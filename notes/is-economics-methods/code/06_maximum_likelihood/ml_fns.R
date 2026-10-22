# ------------------------------------------------------------------
# ml_fns.R -- functions for chapter 6 (maximum likelihood).
#
# Convention (notation_v3.md 1.1, 1.2, 1.4, 1.5): an MLE MINIMISES the sample mean of
# q(w_i, theta) = -ln f(w_i; theta), so the log likelihood is l_i(theta) = -q_i. The score is
# s_i = grad_theta q_i' (k x 1; it is minus the log-likelihood score grad_theta l_i'), A_0 = E[grad^2 q],
# B_0 = E[s s']; under correct specification A_0 = B_0 = I(theta_0) (information matrix equality,
# Wooldridge 2010 eqs. 13.26 and 13.27), otherwise Avar = A_0^-1 B_0 A_0^-1 (QMLE, section 13.11).
#
# Durations (registry 1.5): t_i* latent duration, c_i censoring time, t_i = min(t_i*, c_i) (argument t),
# delta_i = 1{t_i* <= c_i} (argument delta, 1 = exit observed). Weibull in proportional-hazard (PH) form:
#     h(t | x) = exp(x'beta) kappa t^(kappa - 1),   S(t | x) = exp(-exp(x'beta) t^kappa),
# parameter vector par = theta = (beta', log kappa)' (the shape enters on the log scale so that it is
# unrestricted). Under independent right censoring (Wooldridge 2010 eq. 22.22) the log likelihood of
# (t_i, delta_i) given (x_i, c_i) is eq. (22.24) with the Weibull density (22.25):
#     l_i = delta_i [x_i'beta + log kappa + (kappa - 1) log t_i] - exp(x_i'beta) t_i^kappa.
# Grouped durations (Wooldridge 2010 section 22.4): intervals m = 1, ..., M with cut points
# a_0 = 0 < a_1 < ... < a_M (argument cuts); y_im = 1 if the exit falls in interval m; with a PH hazard
# the probability of surviving interval m is alpha_m(x) = exp(-exp(eta_m + x'beta)) (eq. 22.47 with
# eta_m = log of the integrated baseline hazard over the interval), so the exit probability in
# interval m given survival to its start is 1 - exp(-exp(eta_m + x'beta)): a complementary log-log
# model with interval dummies.
#
# Chapter card interface:
#   nll_weibull_cens(par, t, delta, X)        -> numeric: sum of q_i (by_obs = TRUE: the n-vector q_i)
#   info_hat(par, nll, ...)                   -> list(A, B): A_hat = n^-1 sum grad^2 q_i (numDeriv::hessian),
#                                                B_hat = n^-1 sum s_i s_i' (numDeriv::jacobian of q_i)
#   trinity(par_u, par_r, nll, restr, ...)    -> list(W, LR, LM, df)
#   person_period(t, delta, X, cuts)          -> data.frame(id, m, y, columns of X)
# Helpers used by 01 to 06: weibull_parts, fit_weibull, survreg_to_ph, ols_weibull, newton_min,
#   grouped_data, grouped_nll, grouped_grad, grouped_hess, grouped_mle, two_step_vcov, twostep_weibull,
#   c2_surv, c2_dens, c2_hazard, c2_cells, c2_power_moments, pop_nll_weibull, pop_nll_weibull_grad,
#   pseudo_true_weibull, gauss_legendre, loglogistic_ll_obs, vuong, ci95,
#   gen_C1, gen_C2, gen_C4 (the Sorrel data generators; parameters and seeds in 01_dgp.R).
# Packages: stats, numDeriv (info_hat, trinity), maxLik (fit_weibull's default engine; its hessian()
# is always called as maxLik::hessian because numDeriv::hessian masks it, lit_software 13.5).
# ------------------------------------------------------------------

ci95 <- function(est, se) unname(c(est - qnorm(0.975) * se, est + qnorm(0.975) * se))

## =================================================================
## Card functions
## =================================================================

#' Negative log likelihood of the right-censored Weibull PH model (minimisation convention).
#' par is theta = (beta', log kappa)'; t is t_i (observed duration, here months), delta is delta_i
#' (1 = exit observed, 0 = censored at t_i = c_i); X is the n x k matrix with rows x_i' (first column
#' the intercept when there is one). Returns sum_i q_i with q_i = -l_i (Wooldridge 2010 eq. 22.24 with
#' the density 22.25 in PH form), or the n-vector of q_i when by_obs = TRUE (needed for B_hat).
nll_weibull_cens <- function(par, t, delta, X, by_obs = FALSE) {
  X <- as.matrix(X)
  k <- ncol(X)
  if (length(par) != k + 1L) stop("nll_weibull_cens: par must be (beta', log kappa)', length ncol(X) + 1")
  rho <- par[k + 1L]
  eta <- drop(X %*% par[seq_len(k)])
  lt <- log(t)
  q <- -(delta * (eta + rho + (exp(rho) - 1) * lt) - exp(eta + exp(rho) * lt))
  if (by_obs) q else sum(q)
}

#' A_hat and B_hat of an MLE or QMLE at par (Wooldridge 2010 eqs. 13.31 to 13.33 for the CMLE,
#' 13.71 for the QMLE sandwich). nll(par, ..., by_obs = FALSE) must return sum_i q_i and
#' nll(par, ..., by_obs = TRUE) the n-vector q_i (q_i = -l_i). A_hat = n^-1 sum_i grad^2 q_i is the
#' numerical Hessian of the sample objective (numDeriv::hessian); B_hat = n^-1 sum_i s_i s_i' with the
#' scores s_i = grad q_i' from numDeriv::jacobian. Model-based variance A_hat^-1 / n, OPG (BHHH)
#' variance B_hat^-1 / n, sandwich A_hat^-1 B_hat A_hat^-1 / n.
#' numDeriv::hessian is called with eps = 0.01, its step for a coordinate that is exactly zero (log kappa = 0
#' under H_0: kappa = 1); with the default 1e-4 rounding costs about four digits there.
.hess_args <- list(eps = 1e-2, d = 0.1, zero.tol = sqrt(.Machine$double.eps / 7e-7), r = 4, v = 2)
info_hat <- function(par, nll, ...) {
  n <- length(nll(par, ..., by_obs = TRUE))
  A <- numDeriv::hessian(function(p) nll(p, ...), par, method.args = .hess_args) / n
  S <- numDeriv::jacobian(function(p) nll(p, ..., by_obs = TRUE), par)
  list(A = (A + t(A)) / 2, B = crossprod(S) / n)
}

#' Wald, likelihood-ratio and score (LM) statistics for H_0: R(theta_0) = 0 (Wooldridge 2010
#' section 12.6 for M-estimators, section 13.6 for MLE; LR = 2(l_u - l_r) is eq. 13.35).
#' par_u is the unrestricted estimate theta_hat; par_r the restricted estimate theta_hat_r (a full-length
#' vector satisfying the restriction); nll(par, ...) returns sum_i q_i = -l_n(theta); restr(par) returns
#' R(theta), a nu-vector. All three use the Hessian form of the information (A_hat):
#'   W  = R(theta_hat)' [grad R V_u grad R']^-1 R(theta_hat),  V_u = [sum_i grad^2 q_i(theta_hat)]^-1
#'   LR = 2 [nll(theta_hat_r) - nll(theta_hat)] = 2 (l_u - l_r)
#'   LM = g_r' [sum_i grad^2 q_i(theta_hat_r)]^-1 g_r,  g_r = sum_i s_i(theta_hat_r)
#' Under H_0 and correct specification all three are asymptotically chi^2_nu (df = nu); LR relies on the
#' information matrix equality, W and LM could be made robust by replacing A_hat^-1 with the sandwich.
trinity <- function(par_u, par_r, nll, restr, ...) {
  f <- function(p) nll(p, ...)
  r_u <- restr(par_u)
  Rj <- numDeriv::jacobian(restr, par_u)
  V_u <- solve(numDeriv::hessian(f, par_u, method.args = .hess_args))
  W <- drop(t(r_u) %*% solve(Rj %*% V_u %*% t(Rj), r_u))
  LR <- 2 * (f(par_r) - f(par_u))
  g_r <- numDeriv::grad(f, par_r)
  LM <- drop(t(g_r) %*% solve(numDeriv::hessian(f, par_r, method.args = .hess_args), g_r))
  list(W = W, LR = LR, LM = LM, df = length(r_u))
}

#' Person-period (person-interval) rows for grouped duration analysis (Wooldridge 2010 section 22.4).
#' t is t_i, delta is delta_i, X the time-invariant covariates (columns kept with their names), cuts the
#' interval end points a_0 = 0 < a_1 < ... < a_M. Subscriber i contributes one row per interval m that
#' he or she is observed to enter: y_im = 1 in the interval (a_{m-1}, a_m] that contains an observed exit,
#' 0 otherwise. A censored spell contributes the intervals it survives completely (a_m <= t_i); an exit
#' after a_M is censored at a_M. This assumes that censoring happens at cut points, as Wooldridge's
#' section 22.4 does (in the Sorrel data c_i is a whole number of months); a censoring time between
#' two cuts drops the partially observed interval.
person_period <- function(t, delta, X, cuts) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  g <- grouped_data(t, delta, cuts)
  id <- rep(seq_along(t), g$mi)
  m <- sequence(g$mi)
  out <- data.frame(id = id, m = m, y = as.numeric(g$yi[id] == 1 & m == g$mi[id]))
  out <- cbind(out, as.data.frame(X[id, , drop = FALSE]))
  row.names(out) <- NULL
  out
}

## =================================================================
## Weibull pieces, fitting, survreg conversion
## =================================================================

#' Analytic pieces of the Weibull PH objective at par = (beta', log kappa)': q (n-vector q_i), S (n x (k+1)
#' scores s_i' = grad q_i) and H (sum_i grad^2 q_i, (k+1) x (k+1)). With Lambda_i = exp(x_i'beta) t_i^kappa
#' (the cumulative hazard at t_i; the text writes it in words) and rho = log kappa:
#'   s_beta = -x_i (delta_i - Lambda_i),  s_rho = -[delta_i (1 + kappa log t_i) - Lambda_i kappa log t_i]
#'   H_bb = sum Lambda_i x_i x_i',  H_brho = sum Lambda_i kappa log t_i x_i,
#'   H_rhorho = sum [-delta_i kappa log t_i + Lambda_i (kappa log t_i + kappa^2 (log t_i)^2)].
weibull_parts <- function(par, t, delta, X) {
  X <- as.matrix(X)
  k <- ncol(X)
  kap <- exp(par[k + 1L])
  eta <- drop(X %*% par[seq_len(k)])
  lt <- log(t)
  Lam <- exp(eta + kap * lt)
  klt <- kap * lt
  q <- -(delta * (eta + par[k + 1L] + (kap - 1) * lt) - Lam)
  S <- cbind(-X * (delta - Lam), -(delta * (1 + klt) - Lam * klt))
  Hbr <- colSums(X * (Lam * klt))
  H <- rbind(cbind(crossprod(X * Lam, X), Hbr), c(Hbr, sum(-delta * klt + Lam * (klt + klt^2))))
  dimnames(H) <- NULL
  dimnames(S) <- NULL
  list(q = q, S = S, H = H)
}

#' Newton-Raphson minimisation with step halving (used for the Weibull, the grouped-duration and the
#' population objectives). fn, gr, hs: objective, gradient, Hessian. A step that is not a descent
#' direction (Hessian not positive definite far from the optimum) is replaced by a scaled gradient step.
#' Converged when the accepted step is below tol * (1 + max |par|).
newton_min <- function(fn, gr, hs, start, tol = 1e-11, maxit = 200) {
  p <- as.numeric(start)
  f0 <- fn(p)
  if (!is.finite(f0)) stop("newton_min: objective not finite at the start")
  conv <- FALSE
  it <- 0L
  for (it in seq_len(maxit)) {
    g <- gr(p)
    H <- hs(p)
    step <- tryCatch(solve(H, g), error = function(e) rep(NA_real_, length(g)))
    if (anyNA(step) || sum(step * g) <= 0) step <- g / max(1, max(abs(diag(H))))
    lam <- 1
    repeat {
      cand <- p - lam * step
      f1 <- fn(cand)
      if (is.finite(f1) && f1 <= f0 + 1e-13 * abs(f0)) break
      lam <- lam / 2
      if (lam < 2^-40) break
    }
    if (lam < 2^-40) break
    p <- cand
    f0 <- f1
    if (max(abs(lam * step)) < tol * (1 + max(abs(p)))) {
      conv <- TRUE
      break
    }
  }
  list(par = p, value = f0, grad = gr(p), iter = it, converged = conv)
}

#' Weibull PH MLE (or QMLE) under right censoring. Returns par = theta_hat = (beta', log kappa)', nll
#' (minimised sum of q_i), A = H / n and B = S'S / n from weibull_parts, and three variances:
#' V_hess = A^-1 / n (model-based, information matrix equality), V_opg = B^-1 / n (outer product of
#' the scores, BHHH), V_sand = A^-1 B A^-1 / n (sandwich, valid under misspecification).
#' engine = "maxLik": maxLik::maxLik(method = "NR") with the analytic per-observation gradient and the
#' analytic Hessian (fit is the maxLik object); it stops on the gradient only (tol and reltol switched off)
#' and uses qac = "marquardt": with the fee in dollars the Hessian has condition number near 1e5 and the
#' default step halving stalls at a gradient of about 1e-4. engine = "newton": newton_min (Monte Carlo).
#' log_kappa = a number fixes the shape (the restricted model; log_kappa = 0 is the exponential model)
#' and always uses newton_min over beta. start defaults to the exponential rate for the intercept.
fit_weibull <- function(t, delta, X, start = NULL, log_kappa = NULL, engine = c("maxLik", "newton")) {
  engine <- match.arg(engine)
  X <- as.matrix(X)
  k <- ncol(X)
  n <- length(t)
  if (is.null(start)) start <- c(log(sum(delta) / sum(t)), rep(0, k - 1L), 0)
  fit <- NULL
  if (!is.null(log_kappa)) {
    fb <- function(b) nll_weibull_cens(c(b, log_kappa), t, delta, X)
    gb <- function(b) colSums(weibull_parts(c(b, log_kappa), t, delta, X)$S)[seq_len(k)]
    hb <- function(b) weibull_parts(c(b, log_kappa), t, delta, X)$H[seq_len(k), seq_len(k), drop = FALSE]
    o <- newton_min(fb, gb, hb, start[seq_len(k)])
    par <- c(o$par, log_kappa)
    conv <- o$converged
    iter <- o$iter
    grad_max <- max(abs(o$grad))
    engine <- "newton"
  } else if (engine == "maxLik") {
    fit <- maxLik::maxLik(logLik = function(p) -weibull_parts(p, t, delta, X)$q,
                          grad = function(p) -weibull_parts(p, t, delta, X)$S,
                          hess = function(p) -weibull_parts(p, t, delta, X)$H,
                          start = start, method = "NR",
                          control = list(tol = -1, reltol = -1, gradtol = 1e-7, iterlim = 200, qac = "marquardt"))
    par <- unname(fit$estimate)
    grad_max <- max(abs(colSums(weibull_parts(par, t, delta, X)$S)))
    conv <- fit$code %in% c(1L, 2L, 8L) && grad_max < 1e-6
    iter <- fit$iterations
  } else {
    fn <- function(p) nll_weibull_cens(p, t, delta, X)
    gr <- function(p) colSums(weibull_parts(p, t, delta, X)$S)
    hs <- function(p) weibull_parts(p, t, delta, X)$H
    o <- newton_min(fn, gr, hs, start)
    if (!o$converged) {                                         # fall back: BFGS from the start, then Newton
      ob <- optim(start, fn, gr, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))
      o <- newton_min(fn, gr, hs, ob$par)
    }
    par <- o$par
    conv <- o$converged
    iter <- o$iter
    grad_max <- max(abs(o$grad))
  }
  wp <- weibull_parts(par, t, delta, X)
  A <- wp$H / n
  B <- crossprod(wp$S) / n
  Ai <- solve(A)
  list(par = par, nll = sum(wp$q), converged = conv, iter = iter, grad_max = grad_max, n = n,
       A = A, B = B, V_hess = Ai / n, V_opg = solve(B) / n, V_sand = Ai %*% B %*% Ai / n,
       engine = engine, fit = fit)
}

#' survreg(dist = "weibull") fits the AFT form log t = x'mu + sigma W, W standard minimum extreme value
#' (fit$scale = sigma). The PH form h(t | x) = exp(x'beta) kappa t^(kappa - 1) has
#'   kappa = 1 / sigma,  beta = -mu / sigma (the intercept included: baseline log scale = -mu_0 / sigma),
#' so AFT and PH coefficients have opposite signs (lit_software 13.5, question 11). Returns
#' par = (beta', log kappa)' and its variance V by the delta method from vcov(fit), whose order is
#' (mu', log sigma).
survreg_to_ph <- function(fit) {
  mu <- unname(coef(fit))
  k <- length(mu)
  s <- fit$scale
  par <- c(-mu / s, -log(s))
  J <- matrix(0, k + 1L, k + 1L)
  J[cbind(seq_len(k), seq_len(k))] <- -1 / s                  # d beta_j / d mu_j
  J[seq_len(k), k + 1L] <- mu / s                              # d beta_j / d log sigma
  J[k + 1L, k + 1L] <- -1                                      # d log kappa / d log sigma
  V <- J %*% unname(vcov(fit)) %*% t(J)
  list(par = par, V = V)
}

#' The log-linear regression implied by the Weibull PH model, fitted by OLS as if every spell were
#' complete. With E ~ Exp(1), exp(x'beta) t*^kappa = E, so
#'   log t* = -x'beta / kappa + eps / kappa,  eps = log E standard minimum-Gumbel: E[eps] = -Euler
#'   (0.5772...), Var(eps) = pi^2 / 6.
#' OLS of log t on x gives b -> -beta / kappa (slopes) and the mean squared residual s2 -> pi^2 / (6 kappa^2), so
#'   kappa_hat = pi / sqrt(6 s2),  beta_j = -kappa_hat b_j (slopes),  beta_0 = -kappa_hat b_0 - Euler.
#' On complete data this is consistent; on censored spells treated as exits it is the "naive regression".
#' logt: the log durations used; X: n x k with the intercept first. V_bs2 is the stacked M-estimator
#' sandwich of (b, s2) (moments x_i e_i and e_i^2 - s2), V_beta its delta-method transform.
ols_weibull <- function(logt, X) {
  X <- as.matrix(X)
  n <- nrow(X)
  k <- ncol(X)
  b <- as.numeric(qr.coef(qr(X), logt))
  e <- logt - drop(X %*% b)
  s2 <- mean(e^2)
  kap <- pi / sqrt(6 * s2)
  euler <- 0.5772156649015329
  beta <- c(-kap * b[1] - euler, -kap * b[-1])
  G <- cbind(X * e, e^2 - s2)
  J <- rbind(cbind(-crossprod(X) / n, 0), c(-2 * colMeans(X * e), -1))
  Ji <- solve(J)
  V_bs2 <- Ji %*% (crossprod(G) / n) %*% t(Ji) / n
  D <- matrix(0, k, k + 1L)
  D[cbind(seq_len(k), seq_len(k))] <- -kap
  D[, k + 1L] <- b * kap / (2 * s2)                            # d(-kappa b_j) / d s2, kappa = pi / sqrt(6 s2)
  V_beta <- D %*% V_bs2 %*% t(D)
  list(b = b, se_b = sqrt(diag(V_bs2))[seq_len(k)], s2 = s2, kappa = kap,
       se_kappa = sqrt(V_bs2[k + 1L, k + 1L]) * kap / (2 * s2),
       beta = beta, V_bs2 = V_bs2, V_beta = V_beta, se_beta = sqrt(diag(V_beta)))
}

## =================================================================
## Grouped durations (block C3)
## =================================================================

#' For each spell: mi, the last interval observed (the exit interval when the exit is observed by a_M,
#' otherwise the number of intervals survived completely, at most M), and yi = 1 when the spell ends with
#' an observed exit in interval mi. Wooldridge's section 22.4 writes the exit interval as m_i.
grouped_data <- function(t, delta, cuts) {
  M <- length(cuts) - 1L
  if (cuts[1] != 0 || any(diff(cuts) <= 0)) stop("grouped_data: cuts must start at 0 and increase")
  exit_in <- delta == 1 & t <= cuts[M + 1L]
  mi <- ifelse(exit_in, findInterval(t, cuts, left.open = TRUE), pmin(findInterval(t, cuts) - 1L, M))
  data.frame(mi = as.integer(mi), yi = as.numeric(exit_in))
}

.g1 <- function(b) { e <- exp(b); e * exp(-e) / (-expm1(-e)) }                         # d/db log(1 - exp(-e^b))
.g2 <- function(b) { e <- exp(b); e * exp(-e) * (-expm1(-e) - e) / expm1(-e)^2 }       # its derivative

#' Grouped-duration negative log likelihood, Wooldridge 2010 eq. (22.45) with the piecewise-constant PH
#' hazard (22.46) and (22.47): par = (eta_1, ..., eta_M, beta')', eta_m = log(lambda_m (a_m - a_{m-1}))
#' (for a continuous PH baseline: the log of the baseline hazard integrated over interval m), X without an
#' intercept (the eta_m play that role). For spell i with K_i = mi - yi intervals survived:
#'   q_i = exp(x_i'beta) sum_{m <= K_i} exp(eta_m) - yi log(1 - exp(-exp(eta_mi + x_i'beta))),
#' which equals -sum_m ell_im over the spell's person-period rows (see 04_grouped.R).
grouped_nll <- function(par, mi, yi, X, by_obs = FALSE) {
  X <- as.matrix(X)
  k <- ncol(X)
  M <- length(par) - k
  eta <- par[seq_len(M)]
  a <- drop(X %*% par[M + seq_len(k)])
  CE <- c(0, cumsum(exp(eta)))
  q <- exp(a) * CE[mi - yi + 1L]
  ex <- yi == 1
  q[ex] <- q[ex] - log(-expm1(-exp(eta[mi[ex]] + a[ex])))
  if (by_obs) q else sum(q)
}

.grouped_pieces <- function(par, mi, yi, X) {
  X <- as.matrix(X)
  k <- ncol(X)
  M <- length(par) - k
  eta <- par[seq_len(M)]
  a <- drop(X %*% par[M + seq_len(k)])
  E <- exp(eta)
  ea <- exp(a)
  K <- mi - yi
  Ind <- outer(K, seq_len(M), ">=") * ea                       # ea_i 1{m <= K_i}, column m
  Ex <- outer(mi, seq_len(M), "==") * (yi == 1)                # 1{exit in interval m}, column m
  bi <- eta[pmax(mi, 1L)] + a
  G1 <- ifelse(yi == 1, .g1(bi), 0)
  G2 <- ifelse(yi == 1, .g2(bi), 0)
  list(X = X, M = M, E = E, ea = ea, CE = c(0, cumsum(E))[K + 1L], Ind = Ind, Ex = Ex, G1 = G1, G2 = G2)
}

grouped_grad <- function(par, mi, yi, X) {
  p <- .grouped_pieces(par, mi, yi, X)
  c(p$E * colSums(p$Ind) - colSums(p$Ex * p$G1), colSums(p$X * (p$ea * p$CE - p$G1)))
}

grouped_hess <- function(par, mi, yi, X) {
  p <- .grouped_pieces(par, mi, yi, X)
  Hee <- diag(p$E * colSums(p$Ind) - colSums(p$Ex * p$G2), p$M)
  Heb <- p$E * crossprod(p$Ind, p$X) - crossprod(p$Ex * p$G2, p$X)
  Hbb <- crossprod(p$X * (p$ea * p$CE - p$G2), p$X)
  H <- rbind(cbind(Hee, Heb), cbind(t(Heb), Hbb))
  dimnames(H) <- NULL
  H
}

#' Grouped-duration MLE by Newton-Raphson on the interval-form likelihood (22.45), started from the
#' life-table complementary log-log rates and beta = 0. Returns par, grad, H (observed information, the
#' Hessian of the summed q_i), V = H^-1, converged, iter.
grouped_mle <- function(mi, yi, X, M, start = NULL) {
  X <- as.matrix(X)
  if (is.null(start)) {
    at_risk <- vapply(seq_len(M), function(h) sum(mi >= h), numeric(1))
    exits <- vapply(seq_len(M), function(h) sum(yi == 1 & mi == h), numeric(1))
    start <- c(log(-log1p(-(exits + 0.5) / (at_risk + 1))), rep(0, ncol(X)))
  }
  o <- newton_min(function(p) grouped_nll(p, mi, yi, X), function(p) grouped_grad(p, mi, yi, X),
                  function(p) grouped_hess(p, mi, yi, X), start, tol = 1e-12)
  H <- grouped_hess(o$par, mi, yi, X)
  list(par = o$par, nll = o$value, grad = o$grad, H = H, V = solve(H), converged = o$converged, iter = o$iter)
}

## =================================================================
## Two-step estimation (block C4)
## =================================================================

#' Variance of a two-step M-estimator (Wooldridge 2010 eqs. 12.39 to 12.41, estimated by 12.61 and 12.62).
#' S2 is the n x k matrix of second-step scores s_i(theta_hat, gamma_hat)'; A2 the k x k mean Hessian of
#' the second step (A_0); F2 = n^-1 sum_i grad_gamma s_i (k x k_gamma, Gamma_hat; Wooldridge's F_o);
#' R1 = n x k_gamma, row i the first-step influence function psi_i' (Wooldridge's r_i). Returns
#' A2^-1 D A2^-1' / n with D = n^-1 sum_i (s_i + F2 psi_i)(s_i + F2 psi_i)'. F2 = 0 gives the sandwich.
two_step_vcov <- function(S2, A2, F2, R1) {
  S2 <- as.matrix(S2)
  F2 <- as.matrix(F2)
  R1 <- as.matrix(R1)
  if (nrow(R1) != nrow(S2) || nrow(F2) != ncol(S2) || ncol(F2) != ncol(R1))
    stop("two_step_vcov: S2 is n x k, F2 is k x k_gamma and R1 is n x k_gamma")
  G <- S2 + R1 %*% t(F2)
  Ai <- solve(A2)
  Ai %*% (crossprod(G) / nrow(G)) %*% t(Ai) / nrow(G)
}

#' Block C4. d: data frame with usage_wk1 (v_i, first-week in-app check-ins), hh_size, quiz_std,
#' fee_usd, t_months, delta. Step 1 (CMLE): Poisson regression of v on z = (1, hh_size, quiz_std),
#' gamma_hat; usage-intensity index u_hat = z'gamma_hat; psi_i = J1^-1 z_i (v_i - mu_i), J1 = n^-1 sum
#' mu_i z_i z_i' (eq. 13.65: the CMLE influence function). Step 2 (MLE): Weibull PH with
#' x_i = (1, fee_usd, u_hat_i), theta = (beta_0, beta_fee, beta_u, log kappa).
#' Gamma_hat = F2 = n^-1 sum grad_gamma s_i, analytically: with Lambda_i the cumulative hazard (the text
#' writes it in words),
#'   grad_gamma s_beta,i = beta_u Lambda_i x_i z_i' - (delta_i - Lambda_i) e_3 z_i',
#'   grad_gamma s_rho,i  = beta_u kappa log t_i Lambda_i z_i'.
#' V_naive = A2^-1 / n treats u_hat as data (what software prints); V_naive_sand is its sandwich;
#' V_corrected = two_step_vcov(S2, A2, F2, R1) (eq. 12.41). The second-step objective does not involve
#' v_i, so the section 13.10.2 argument (eq. 13.67) does not apply here: F2 != 0 adds variance.
twostep_weibull <- function(d, start = NULL) {
  Z <- cbind(1, d$hh_size, d$quiz_std)
  n <- nrow(Z)
  g <- glm.fit(Z, d$usage_wk1, family = poisson(), control = list(epsilon = 1e-12, maxit = 100))
  gamma <- unname(g$coefficients)
  uh <- drop(Z %*% gamma)
  mu <- exp(uh)
  J1 <- crossprod(Z * mu, Z) / n
  R1 <- (Z * (d$usage_wk1 - mu)) %*% solve(J1)
  X <- cbind(1, d$fee_usd, uh)
  f <- fit_weibull(d$t_months, d$delta, X, start = start, engine = "newton")
  p <- f$par
  kap <- exp(p[4])
  lt <- log(d$t_months)
  Lam <- exp(drop(X %*% p[1:3]) + kap * lt)
  wp <- weibull_parts(p, d$t_months, d$delta, X)
  F2 <- rbind(p[3] * crossprod(X * Lam, Z) / n, p[3] * colSums(Z * (kap * lt * Lam)) / n)
  F2[3, ] <- F2[3, ] - colSums(Z * (d$delta - Lam)) / n
  dimnames(F2) <- NULL
  dimnames(R1) <- NULL
  A2 <- wp$H / n
  list(gamma = gamma, u_hat = uh, par = p, est = p[3], S2 = wp$S, A2 = A2, F2 = F2, R1 = R1,
       V_naive = solve(A2) / n, V_naive_sand = f$V_sand, V_corrected = two_step_vcov(wp$S, A2, F2, R1),
       converged = isTRUE(g$converged) && f$converged, iter = f$iter)
}

## =================================================================
## Block C2: log-logistic mixture, population objective, pseudo-true value
## =================================================================

# Mixture parameters for a subscriber with monthly fee `fee` and referral indicator `ref`: deal-seeker
# probability pi(ref) = plogis(g0 + g_ref ref); each segment is log-logistic with S(t) = 1 / (1 + (t / a)^p),
# median a, and the fee rescales time in both segments (AFT): a = median * exp(-theta_fee (fee - fee_ref)).
.c2_type <- function(fee, ref, P) {
  sh <- exp(-P$theta_fee * (fee - P$fee_ref))
  list(pi = plogis(P$g0 + P$g_ref * ref), a_d = P$med_deal * sh, a_r = P$med_reg * sh, p_d = P$p_deal, p_r = P$p_reg)
}
c2_surv <- function(t, fee, ref, P) {
  ty <- .c2_type(fee, ref, P)
  ty$pi * plogis(-ty$p_d * (log(t) - log(ty$a_d))) + (1 - ty$pi) * plogis(-ty$p_r * (log(t) - log(ty$a_r)))
}
c2_dens <- function(t, fee, ref, P) {
  ty <- .c2_type(fee, ref, P)
  (ty$pi * dlogis(log(t), log(ty$a_d), 1 / ty$p_d) + (1 - ty$pi) * dlogis(log(t), log(ty$a_r), 1 / ty$p_r)) / t
}
c2_hazard <- function(t, fee, ref, P) c2_dens(t, fee, ref, P) / c2_surv(t, fee, ref, P)

#' Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
gauss_legendre <- function(n) {
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = rev(e$values), w = rev(2 * e$vectors[1, ]^2))
}

#' The population of block C2 as a finite set of cells (fee, ref, c): every covariate value and every
#' censoring time c = 1, ..., months has a known probability w (fee shares p_fee, referral share p_ref,
#' signup month k with probability proportional to growth^k, c = months - k + 1). For each cell the
#' expectations that the Weibull objective needs are one-dimensional integrals over the true mixture
#' density f_0(t | x) on (0, c]: F = int f_0 (closed form), L = int log t f_0, and, for a given kappa,
#' P = int t^kappa f_0 and Plog = int t^kappa log t f_0 (c2_power_moments). They are computed in
#' u = log t (where each segment's density is logistic) by n_gl-point Gauss-Legendre on [-120, -60],
#' [-60, -30], [-30, -10], [-10, 0] (t <= 1) and on [log(m - 1), log m] for m = 2, ..., months.
c2_cells <- function(P, n_gl = 40) {
  gl <- gauss_legendre(n_gl)
  seg <- rbind(c(-120, -60), c(-60, -30), c(-30, -10), c(-10, 0), cbind(log(seq_len(P$months - 1)), log(2:P$months)))
  seg_month <- c(1, 1, 1, 1, 2:P$months)
  mid <- rowMeans(seg)
  half <- (seg[, 2] - seg[, 1]) / 2
  u <- as.vector(outer(gl$x, half) + rep(mid, each = n_gl))
  wu <- as.vector(outer(gl$w, half))
  mon <- rep(seg_month, each = n_gl)
  pc <- P$growth^(P$months - seq_len(P$months) + 1)            # P(c), c = 1, ..., months
  pc <- pc / sum(pc)
  types <- expand.grid(fee = P$fees, ref = 0:1)
  cells <- list()
  nodes <- list()
  for (j in seq_len(nrow(types))) {
    fee <- types$fee[j]
    ref <- types$ref[j]
    ty <- .c2_type(fee, ref, P)
    g <- ty$pi * dlogis(u, log(ty$a_d), 1 / ty$p_d) + (1 - ty$pi) * dlogis(u, log(ty$a_r), 1 / ty$p_r)
    wg <- wu * g
    w_type <- P$p_fee[match(fee, P$fees)] * (if (ref == 1) P$p_ref else 1 - P$p_ref)
    cs <- seq_len(P$months)
    cells[[j]] <- data.frame(fee = fee, ref = ref, c = cs, w = w_type * pc,
                             F = cumsum(rowsum(wg, mon)[, 1]), L = cumsum(rowsum(wg * u, mon)[, 1]),
                             S = c2_surv(cs, fee, ref, P))
    nodes[[j]] <- list(u = u, wg = wg, mon = mon)
  }
  list(cells = do.call(rbind, cells), nodes = nodes, months = P$months)
}

c2_power_moments <- function(cobj, kappa) {
  P <- unlist(lapply(cobj$nodes, function(nd) cumsum(rowsum(nd$wg * exp(kappa * nd$u), nd$mon)[, 1])), use.names = FALSE)
  Plog <- unlist(lapply(cobj$nodes, function(nd) cumsum(rowsum(nd$wg * exp(kappa * nd$u) * nd$u, nd$mon)[, 1])), use.names = FALSE)
  list(P = P, Plog = Plog)
}

#' Population objective E[q(w, theta)] of the Weibull fit to block C2 (the expected negative log likelihood
#' over the cells), par = (beta_0, beta_fee, beta_ref, log kappa). For a cell with eta = x'beta:
#'   E[q] = -(eta + log kappa) F - (kappa - 1) L + exp(eta) (P + c^kappa S).
#' Its minimiser is the pseudo-true value theta_* (White 1982; Wooldridge 2010 section 13.11): the Weibull
#' closest to the truth in expected Kullback-Leibler divergence, for this censoring design.
pop_nll_weibull <- function(par, cobj) {
  cl <- cobj$cells
  kap <- exp(par[4])
  eta <- par[1] + par[2] * cl$fee + par[3] * cl$ref
  pm <- c2_power_moments(cobj, kap)
  sum(cl$w * (-(eta + par[4]) * cl$F - (kap - 1) * cl$L + exp(eta) * (pm$P + cl$c^kap * cl$S)))
}
pop_nll_weibull_grad <- function(par, cobj) {
  cl <- cobj$cells
  kap <- exp(par[4])
  eta <- par[1] + par[2] * cl$fee + par[3] * cl$ref
  pm <- c2_power_moments(cobj, kap)
  r <- cl$w * (-cl$F + exp(eta) * (pm$P + cl$c^kap * cl$S))
  c(sum(r), sum(r * cl$fee), sum(r * cl$ref),
    sum(cl$w * (-cl$F - kap * cl$L + exp(eta) * kap * (pm$Plog + cl$c^kap * log(cl$c) * cl$S))))
}

#' Pseudo-true value of the Weibull fit to block C2 by quadrature (a check on the large-sample value that
#' 01_dgp.R stores as the truth): BFGS on pop_nll_weibull, then Newton with the numerical Jacobian of the
#' analytic gradient.
pseudo_true_weibull <- function(P, n_gl = 40, start = c(-3, 0.015, -0.2, log(0.8))) {
  cobj <- c2_cells(P, n_gl)
  o <- optim(start, pop_nll_weibull, pop_nll_weibull_grad, cobj = cobj, method = "BFGS",
             control = list(maxit = 5000, reltol = 1e-15))
  hs <- function(p) { H <- numDeriv::jacobian(pop_nll_weibull_grad, p, cobj = cobj); (H + t(H)) / 2 }
  nt <- newton_min(function(p) pop_nll_weibull(p, cobj), function(p) pop_nll_weibull_grad(p, cobj), hs, o$par, tol = 1e-12)
  list(par = nt$par, kappa = exp(nt$par[4]), value = nt$value, grad = nt$grad,
       converged = nt$converged && max(abs(nt$grad)) < 1e-9)
}

#' Per-observation log likelihood of a survreg(dist = "loglogistic") fit on the time scale:
#' z = (log t - x'mu) / sigma, uncensored log f(t) = log dlogis(z) - log sigma - log t, censored
#' log S(t) = log plogis(-z). Sums to fit$loglik[2].
loglogistic_ll_obs <- function(fit, t, delta, X) {
  z <- (log(t) - drop(as.matrix(X) %*% coef(fit))) / fit$scale
  ifelse(delta == 1, dlogis(z, log = TRUE) - log(fit$scale) - log(t), plogis(-z, log.p = TRUE))
}

#' Vuong (1989) test for two non-nested models with the same number of parameters (Wooldridge 2010
#' section 13.11.2): ll1, ll2 are the per-observation log likelihoods, each at its own (Q)MLE.
#' V_n = n^-1/2 sum_i Delta l_i / sigma_hat (sigma_hat with divisor n); asymptotically N(0, 1) when the two
#' models are equally close to the truth in KL divergence; large positive values favour model 1.
vuong <- function(ll1, ll2) {
  dl <- ll1 - ll2
  n <- length(dl)
  s <- sqrt(mean((dl - mean(dl))^2))
  stat <- sqrt(n) * mean(dl) / s
  list(stat = stat, p_two_sided = 2 * pnorm(-abs(stat)), mean = mean(dl), sd = s, sum = sum(dl), n = n)
}

## =================================================================
## Sorrel data generators (parameters and seeds: 01_dgp.R, stored in truth06.rds$params)
## =================================================================
## Common design: a subscription starts on the first day of the signup month k = 1, ..., months
## (probability proportional to growth^k: a growing platform), and the data are pulled at the end of
## month `months`, so the censoring time is c = months - k + 1 whole months, drawn independently of the
## latent duration (independent censoring, Wooldridge 2010 eq. 22.22). Monthly fee (USD) from `fees`
## with shares p_fee. The order of the random draws is part of each generator's definition.

#' C1: Weibull PH tenure. referred ~ Bernoulli(p_ref); x = (1, fee_usd, referred);
#' t* = (E / exp(x'beta))^(1 / kappa), E ~ Exp(1), so S(t | x) = exp(-exp(x'beta) t^kappa).
gen_C1 <- function(par, n = par$n) {
  signup <- sample.int(par$months, n, replace = TRUE, prob = par$growth^seq_len(par$months))
  fee <- sample(par$fees, n, replace = TRUE, prob = par$p_fee)
  referred <- rbinom(n, 1, par$p_ref)
  c_months <- par$months - signup + 1
  eta <- par$beta[1] + par$beta[2] * fee + par$beta[3] * referred
  t_star <- (rexp(n) / exp(eta))^(1 / par$kappa)
  data.frame(id = seq_len(n), signup_month = signup, fee_usd = fee, referred = referred, c_months = c_months,
             t_months = pmin(t_star, c_months), delta = as.numeric(t_star <= c_months), t_star_true = t_star)
}

#' C2: two latent segments. deal_true ~ Bernoulli(plogis(g0 + g_ref referred)) marks a deal-seeker
#' (joined on a first-box discount); t* is log-logistic with median med_deal or med_reg times
#' exp(-theta_fee (fee - fee_ref)) and shape p_deal or p_reg: t* = a (U / (1 - U))^(1 / p), U ~ U(0, 1).
gen_C2 <- function(par, n = par$n) {
  signup <- sample.int(par$months, n, replace = TRUE, prob = par$growth^seq_len(par$months))
  fee <- sample(par$fees, n, replace = TRUE, prob = par$p_fee)
  referred <- rbinom(n, 1, par$p_ref)
  c_months <- par$months - signup + 1
  deal <- rbinom(n, 1, plogis(par$g0 + par$g_ref * referred))
  a <- ifelse(deal == 1, par$med_deal, par$med_reg) * exp(-par$theta_fee * (fee - par$fee_ref))
  p <- ifelse(deal == 1, par$p_deal, par$p_reg)
  U <- runif(n)
  t_star <- a * (U / (1 - U))^(1 / p)
  data.frame(id = seq_len(n), signup_month = signup, fee_usd = fee, referred = referred, c_months = c_months,
             t_months = pmin(t_star, c_months), delta = as.numeric(t_star <= c_months), deal_true = deal,
             t_star_true = t_star)
}

#' C4: z = (1, hh_size, quiz_std): household size from the onboarding form and a standardised cooking-quiz
#' score; usage-intensity index u = z'gamma (log of expected first-week check-ins); v = usage_wk1 ~
#' Poisson(exp(u)); x = (1, fee_usd, u) in a Weibull PH hazard. Given z, v and t* are independent.
gen_C4 <- function(par, n = par$n) {
  signup <- sample.int(par$months, n, replace = TRUE, prob = par$growth^seq_len(par$months))
  fee <- sample(par$fees, n, replace = TRUE, prob = par$p_fee)
  hh <- sample(par$hh_sizes, n, replace = TRUE, prob = par$p_hh)
  quiz <- rnorm(n)
  u <- par$gamma[1] + par$gamma[2] * hh + par$gamma[3] * quiz
  v <- as.numeric(rpois(n, exp(u)))
  c_months <- par$months - signup + 1
  eta <- par$beta[1] + par$beta[2] * fee + par$beta[3] * u
  t_star <- (rexp(n) / exp(eta))^(1 / par$kappa)
  data.frame(id = seq_len(n), signup_month = signup, fee_usd = fee, hh_size = hh, quiz_std = quiz, usage_wk1 = v,
             c_months = c_months, t_months = pmin(t_star, c_months), delta = as.numeric(t_star <= c_months),
             u_true = u, t_star_true = t_star)
}
