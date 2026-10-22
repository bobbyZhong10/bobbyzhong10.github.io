# ------------------------------------------------------------------
# gmm_fns.R -- functions for chapter 7 (GMM and minimum distance).
#
# Notation (notation_v3.md 1.3 and 1.6): moment function g(w, theta) (r x 1), sample moment
# gbar_n(theta) = n^-1 sum_i g(w_i, theta); r moments, k = dim(theta), r - k overidentifying
# restrictions; weight matrix W; Jacobian D_0 = E[dg(w, theta_0) / dtheta'] (r x k; a matrix of
# derivatives, not a treatment); moment variance S_0 = E[g g']; S_n(theta) = n^-1 sum_i g_i(theta) g_i(theta)'
# (uncentred); first-step estimate theta^(1) with weight W^(1), two-step theta^(2); J = n gbar' S^-1 gbar.
# Minimum distance: reduced-form estimate pi_hat (r x 1), pi_0 = h(theta_0), Sigma_pi its variance;
# D_0 = -dh(theta_0)/dtheta' when MD is read as GMM with moment pi_hat - h(theta).
#
# A moment function `gfun(theta, data)` returns the n x r matrix whose row i is g(w_i, theta)'.
# `data` can be anything gfun understands (a data frame or a list of matrices).
#
# Chapter card interface:
#   md_fit(pi_hat, V_pi, h, theta0, W = solve(V_pi))         -> list(theta, vcov, J, df, ...)
#   gmm_fit(gfun, theta0, data, type = c("two_step", "iterated", "cue"))  -> list(theta, vcov, J, df, ...)
#   windmeijer_vcov(gfun, theta1, theta2, data, W1)           -> k x k matrix
#   s_set(gfun, grid, data, level = 0.95)                     -> data.frame(theta, S, p, in_set)
# Arguments and the registered symbols they stand for (notation_v3.md 1.3, 1.6):
#   pi_hat is $\hat\pi$ (r x 1); V_pi is $\hat\Sigma_\pi / n$, the estimated variance of $\hat\pi$ itself;
#   h is $h(\theta)$; W is $W$; md_fit's H is $\partial h / \partial\theta' = -D_0$;
#   gfun is $g(w, \theta)$, returned as the n x r matrix of $g(w_i, \theta)'$; colMeans(gfun) is $\bar g_n(\theta)$;
#   theta0 is a starting value (not the true $\theta_0$); W1 is $W^{(1)}$; theta1 is $\hat\theta^{(1)}$,
#   theta2 is $\hat\theta^{(2)}$; the returned D is $\hat D$ (sample $D_0$), S is $\hat S$ (sample $S_0$),
#   moment_var(G) is $\hat S_n(\theta)$; s_stat and the S column of s_set are $\mathcal{S}_n(\theta)$; J is $J$;
#   df is $r - k$; Lambda (03_gmm_mc.R) is $\hat\Lambda = -(\hat D' W \hat D)^{-1} \hat D' W$ (AGS 2017, the registry's $\Lambda_0$).
# Helpers: gbar_jac, moment_var, gmm_min, gmm_sandwich, s_stat; the moment functions g_iv, g_conc,
# g_pool of blocks G1 to G3 and their data builders iv_data, pool_data; market_slopes; the Tamarind
# generators gen_G1, gen_G2, gen_G3 (parameters and seeds in 01_dgp.R); pop_moments (closed-form D_0, S_0
# and asymptotic standard errors of the G2 and G3 designs).
# Packages: stats and numDeriv (jacobian, grad).
# ------------------------------------------------------------------

## =================================================================
## Building blocks
## =================================================================

#' Jacobian of the sample moment, the sample analogue of D_0: the r x k matrix d gbar_n(theta) / d theta'
#' (numDeriv::jacobian, Richardson extrapolation; exact up to rounding when g is linear in theta).
gbar_jac <- function(gfun, theta, data) {
  numDeriv::jacobian(function(t) colMeans(as.matrix(gfun(t, data))), as.numeric(theta))
}

#' Moment variance from the n x r matrix G of moments g_i(theta)': S_n = n^-1 sum_i g_i g_i' (uncentred, the
#' registry's S_n(theta) and Windmeijer's W_N(theta)); centered = TRUE subtracts gbar gbar' (Hansen 2014, eq. 13.4).
moment_var <- function(G, centered = FALSE) {
  G <- as.matrix(G)
  if (centered) G <- sweep(G, 2L, colMeans(G))
  crossprod(G) / nrow(G)
}

#' Minimise the GMM criterion gbar_n(theta)' W gbar_n(theta) for a fixed weight W (r x r) by Gauss-Newton with
#' step halving: step = -(D'WD)^-1 D'W gbar at the current theta. For moments linear in theta the first step
#' lands on the minimiser and the second confirms it. Converged when max_j |step_j| < tol (1 + max_j |theta_j|), or
#' when a later step no longer lowers the criterion beyond rounding error (relative 1e-12).
#' Returns list(theta, obj = the criterion at theta, D = the Jacobian at the last iterate, iter, converged).
gmm_min <- function(gfun, theta0, data, W, tol = 1e-10, maxit = 100) {
  theta <- as.numeric(theta0)
  Q <- function(th) { gb <- colMeans(as.matrix(gfun(th, data))); drop(crossprod(gb, W %*% gb)) }
  q_cur <- Q(theta)
  if (!is.finite(q_cur)) stop("gmm_min: the criterion is not finite at theta0")
  converged <- FALSE
  D <- NULL
  iter <- 0L
  for (iter in seq_len(maxit)) {
    D <- gbar_jac(gfun, theta, data)
    gb <- colMeans(as.matrix(gfun(theta, data)))
    step <- -as.numeric(solve(crossprod(D, W %*% D), crossprod(D, W %*% gb)))
    scale <- 1 + max(abs(theta))
    lambda <- 1
    repeat {
      cand <- theta + lambda * step
      q_new <- Q(cand)
      if (is.finite(q_new) && (q_new <= q_cur * (1 + 1e-12) || max(abs(step)) < sqrt(.Machine$double.eps) * scale)) break
      lambda <- lambda / 2
      if (lambda < 2^-30) break
    }
    if (!is.finite(q_new) || lambda < 2^-30) break
    # after the first iteration a step that does not lower the criterion by more than rounding error is noise from
    # the numerical Jacobian (under weak identification D is tiny and that noise is amplified): stop where we are
    if (iter > 1L && q_cur - q_new <= 1e-12 * q_cur) {
      converged <- TRUE
      break
    }
    theta <- cand
    q_cur <- q_new
    if (max(abs(step)) < tol * (1 + max(abs(theta)))) {
      converged <- TRUE
      break
    }
  }
  list(theta = theta, obj = q_cur, D = D, iter = iter, converged = converged)
}

#' Sandwich variance of a GMM estimator with a fixed weight W (Wooldridge 2010, Theorem 14.2 and eqs. 14.13 to
#' 14.15): (D'WD)^-1 D'W S W D (D'WD)^-1 / n, with D the r x k Jacobian (D_0-hat) and S the moment variance (S_0-hat).
gmm_sandwich <- function(D, W, S, n) {
  B <- solve(crossprod(D, W %*% D))
  B %*% crossprod(D, W %*% S %*% W %*% D) %*% B / n
}

## =================================================================
## Card functions
## =================================================================

#' Classical minimum distance (Wooldridge 2010, section 14.5, eqs. 14.63 to 14.69).
#' pi_hat: the r-vector of reduced-form estimates, the book's pi-hat. V_pi: r x r, the estimated variance matrix of
#'   pi_hat itself, i.e. the registry's Sigma_pi-hat / n (with markets of different sizes there is no common n, so
#'   the variance of pi_hat is passed directly). h: function(theta) -> r-vector, the map h(theta) from structural to
#'   reduced-form parameters. theta0: start. W: r x r weight; the default solve(V_pi) is the optimal MD weight.
#' Minimises (pi_hat - h(theta))' W (pi_hat - h(theta)) by Gauss-Newton with step halving; H = dh/dtheta' by
#' numDeriv::jacobian (H = -D_0 in the registry's GMM reading of MD).
#' vcov = (H'WH)^-1 H'W V_pi W H (H'WH)^-1, which is (H' V_pi^-1 H)^-1 (eq. 14.67) when W = solve(V_pi).
#' J = (pi_hat - h(theta_hat))' W (pi_hat - h(theta_hat)), the overidentification statistic of eq. (14.69), chi^2 with
#' df = r - k under H0 when W = solve(V_pi) (for another W the value is reported but is not chi^2); p its p-value.
#' Returns list(theta, vcov, J, df, p, H, resid = pi_hat - h(theta_hat), iter, converged).
md_fit <- function(pi_hat, V_pi, h, theta0, W = solve(V_pi), tol = 1e-12, maxit = 200) {
  pi_hat <- as.numeric(pi_hat)
  V_pi <- as.matrix(V_pi)
  W <- as.matrix(W)
  r <- length(pi_hat)
  if (!all(dim(V_pi) == r) || !all(dim(W) == r)) stop("md_fit: V_pi and W must be r x r with r = length(pi_hat)")
  theta <- as.numeric(theta0)
  k <- length(theta)
  crit <- function(th) { e <- pi_hat - h(th); drop(crossprod(e, W %*% e)) }
  q_cur <- crit(theta)
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(maxit)) {
    H <- numDeriv::jacobian(h, theta)
    e <- pi_hat - h(theta)
    step <- as.numeric(solve(crossprod(H, W %*% H), crossprod(H, W %*% e)))
    lambda <- 1
    repeat {
      cand <- theta + lambda * step
      q_new <- crit(cand)
      if (is.finite(q_new) && (q_new <= q_cur * (1 + 1e-12) || max(abs(step)) < sqrt(.Machine$double.eps) * (1 + max(abs(theta))))) break
      lambda <- lambda / 2
      if (lambda < 2^-30) break
    }
    if (!is.finite(q_new) || lambda < 2^-30) break
    theta <- cand
    q_cur <- q_new
    if (max(abs(step)) < tol * (1 + max(abs(theta)))) {
      converged <- TRUE
      break
    }
  }
  H <- numDeriv::jacobian(h, theta)
  B <- solve(crossprod(H, W %*% H))
  vcov <- B %*% crossprod(H, W %*% V_pi %*% W %*% H) %*% B
  resid <- pi_hat - h(theta)
  J <- drop(crossprod(resid, W %*% resid))
  df <- r - k
  list(theta = theta, vcov = vcov, J = J, df = df, p = if (df > 0) pchisq(J, df, lower.tail = FALSE) else NA_real_,
       H = H, resid = resid, iter = iter, converged = converged)
}

#' GMM: two-step, iterated, or continuously updated (Hansen 1982; Hansen, Heaton and Yaron 1996; Wooldridge 2010,
#' section 14.1; Hansen 2014, sections 13.4 and 13.5).
#' gfun(theta, data): n x r moments. theta0: start of the first step. W1: the first-step weight W^(1) (r x r;
#' NULL = I_r). centered: whether S is centred (FALSE: S_n(theta) = n^-1 sum g g', the registry's definition and
#' Windmeijer's W_N; TRUE: Hansen 2014 eq. 13.4).
#'   first step   theta^(1) = argmin gbar' W1 gbar
#'   two_step     theta^(2) = argmin gbar' S(theta^(1))^-1 gbar
#'   iterated     theta^(j+1) = argmin gbar' S(theta^(j))^-1 gbar from theta^(2), until max |change| < tol (1 + max |theta|)
#'   cue          argmin gbar(theta)' S(theta)^-1 gbar(theta) (optim BFGS with numDeriv gradients from theta^(2)); with
#'                the uncentred S this is the minimiser of the Stock-Wright statistic, and the minimiser is the same
#'                with the centred S (J_c = J_u / (1 - J_u / n))
#' Conventional variance vcov = (D'WD)^-1 / n with D = D_0-hat at the estimate and W the weight used in the last
#' minimisation: S(theta^(1))^-1 for two_step (Windmeijer's "naive" variance; the gmm package instead re-evaluates S at
#' theta^(2)), S at the converged estimate for iterated and cue. J = n gbar(theta_hat)' W gbar(theta_hat) with that W;
#' df = r - k; p its chi^2 p-value. Returns list(theta, vcov, J, df, p, theta1, W, S (moment variance at theta_hat),
#' D, iter, converged, type, n).
gmm_fit <- function(gfun, theta0, data, type = c("two_step", "iterated", "cue"), W1 = NULL, centered = FALSE,
                    tol = 1e-10, maxit = 500) {
  type <- match.arg(type)
  G0 <- as.matrix(gfun(as.numeric(theta0), data))
  n <- nrow(G0)
  r <- ncol(G0)
  k <- length(theta0)
  if (r < k) stop("gmm_fit: fewer moments than parameters")
  if (is.null(W1)) W1 <- diag(r)
  s1 <- gmm_min(gfun, theta0, data, W1, tol)
  theta1 <- s1$theta
  W <- solve(moment_var(gfun(theta1, data), centered))
  s2 <- gmm_min(gfun, theta1, data, W, tol)
  theta <- s2$theta
  D <- s2$D
  iter <- 2L
  converged <- s1$converged && s2$converged
  if (type == "iterated") {
    converged <- FALSE
    for (j in seq_len(maxit)) {
      W <- solve(moment_var(gfun(theta, data), centered))
      s <- gmm_min(gfun, theta, data, W, tol)
      change <- max(abs(s$theta - theta))
      theta <- s$theta
      D <- s$D
      iter <- iter + 1L
      if (s$converged && change < tol * (1 + max(abs(theta)))) {
        converged <- TRUE
        break
      }
    }
  }
  if (type == "cue") {
    crit <- function(th) {
      G <- as.matrix(gfun(th, data))
      gb <- colMeans(G)
      drop(crossprod(gb, solve(moment_var(G, centered), gb)))
    }
    op <- optim(theta, crit, gr = function(th) numDeriv::grad(crit, th), method = "BFGS",
                control = list(reltol = 1e-15, maxit = 5000))
    theta <- op$par
    converged <- op$convergence == 0
    iter <- iter + as.integer(op$counts[["function"]])
    W <- solve(moment_var(gfun(theta, data), centered))
    D <- gbar_jac(gfun, theta, data)
  }
  G <- as.matrix(gfun(theta, data))
  gb <- colMeans(G)
  J <- n * drop(crossprod(gb, W %*% gb))
  df <- r - k
  list(theta = theta, vcov = solve(crossprod(D, W %*% D)) / n, J = J, df = df,
       p = if (df > 0) pchisq(J, df, lower.tail = FALSE) else NA_real_,
       theta1 = theta1, W = W, S = moment_var(G, centered), D = D, iter = iter, converged = converged, type = type, n = n)
}

#' Windmeijer (2005) finite-sample corrected variance of the two-step GMM estimator, for moments LINEAR in theta
#' (the case the paper proves; Phase 3 ruling 16). No R GMM package implements it for general moments.
#' gfun: moments; theta1: first-step estimate theta^(1) (obtained with weight W1 = W^(1)); theta2: two-step estimate
#' theta^(2) (weight S(theta^(1))^-1); W1: the first-step weight matrix; centered as in gmm_fit.
#'   V2 = (D'S1^-1 D)^-1 / n            naive two-step variance, S1 = S(theta^(1))
#'   V1 = gmm_sandwich(D, W1, S1, n)    robust variance of the first step
#'   column j of Dw = (D'S1^-1 D)^-1 D'S1^-1 (dS(theta)/dtheta_j at theta^(1)) S1^-1 gbar(theta^(2)),
#'                    the derivative of theta^(2) with respect to the point at which the weight is evaluated
#'                    (dS/dtheta_j = n^-1 sum (dg_i/dtheta_j g_i' + g_i dg_i'/dtheta_j), per-observation derivatives
#'                    by numDeriv::jacobian)
#'   corrected      V2 + Dw V2 + V2 Dw' + Dw V1 Dw'
#' Returns the corrected k x k matrix with attributes "naive" (V2), "V1" and "dtheta2_dtheta1" (Dw).
#' Checked in tests/test_fns.R against plm::vcovHC on a two-step Arellano-Bond pgmm fit.
windmeijer_vcov <- function(gfun, theta1, theta2, data, W1, centered = FALSE) {
  theta1 <- as.numeric(theta1)
  theta2 <- as.numeric(theta2)
  G1 <- as.matrix(gfun(theta1, data))
  n <- nrow(G1)
  r <- ncol(G1)
  k <- length(theta1)
  S1 <- moment_var(G1, centered)
  S1i <- solve(S1)
  D2 <- gbar_jac(gfun, theta2, data)
  D1 <- gbar_jac(gfun, theta1, data)
  gb1 <- colMeans(G1)
  gb2 <- colMeans(as.matrix(gfun(theta2, data)))
  V2 <- solve(crossprod(D2, S1i %*% D2)) / n
  V1 <- gmm_sandwich(D1, as.matrix(W1), S1, n)
  dG <- numDeriv::jacobian(function(t) as.vector(as.matrix(gfun(t, data))), theta1)   # (n r) x k, column-major rows
  A <- solve(crossprod(D2, S1i %*% D2), crossprod(D2, S1i))                          # k x r
  Dw <- matrix(0, k, k)
  for (j in seq_len(k)) {
    Gj <- matrix(dG[, j], n, r)
    dS <- (crossprod(Gj, G1) + crossprod(G1, Gj)) / n
    if (centered) {
      gj <- colMeans(Gj)
      dS <- dS - (gj %o% gb1 + gb1 %o% gj)
    }
    Dw[, j] <- A %*% dS %*% S1i %*% gb2
  }
  Vc <- V2 + Dw %*% V2 + V2 %*% t(Dw) + Dw %*% V1 %*% t(Dw)
  Vc <- (Vc + t(Vc)) / 2
  attr(Vc, "naive") <- V2
  attr(Vc, "V1") <- V1
  attr(Vc, "dtheta2_dtheta1") <- Dw
  Vc
}

#' Stock and Wright (2000) S statistic at theta: S_n(theta) = n gbar(theta)' S_n(theta)^-1 gbar(theta), with the
#' uncentred moment variance S_n(theta) evaluated at the same theta (the registry's calligraphic S_n). Under
#' H0: theta = theta_0 it is asymptotically chi^2 with r degrees of freedom whatever the strength of identification.
s_stat <- function(gfun, theta, data) {
  G <- as.matrix(gfun(theta, data))
  gb <- colMeans(G)
  nrow(G) * drop(crossprod(gb, solve(crossprod(G) / nrow(G), gb)))
}

#' S confidence set by test inversion over a grid of scalar theta values (k = 1): the set
#' {theta : S_n(theta) <= chi^2_{r, level}}. Returns data.frame(theta, S, p, in_set) with p = P(chi^2_r > S).
s_set <- function(gfun, grid, data, level = 0.95) {
  grid <- as.numeric(grid)
  r <- ncol(as.matrix(gfun(grid[1], data)))
  S <- vapply(grid, function(th) s_stat(gfun, th, data), numeric(1))
  data.frame(theta = grid, S = S, p = pchisq(S, r, lower.tail = FALSE), in_set = S <= qchisq(level, r))
}

## =================================================================
## Moment functions and data of blocks G1 to G3
## =================================================================

#' Linear IV moments g_i(theta) = z_i (y_i - x_i'theta): d$y (n), d$X (n x k regressors), d$Z (n x r instruments).
g_iv <- function(theta, d) d$Z * drop(d$y - d$X %*% theta)

#' Concentrated scalar moments of the demand equation with the intercept partialled out:
#' g_i(theta) = z_i (y~_i - theta p~_i), with z_i the mean-zero experimental shocks (d$Zs, n x K; the sample mean of
#' each column is exactly 0 because z_j = step_j (T_j - share_j) is centred at the exact share: G3's experiments are 50/50,
#' G2's deep-discount test (z2) is 20/80), y~ and p~ the outcome and log price minus their sample means (d$yc, d$pc).
#' Equivalent, by Frisch-Waugh, to including the intercept and concentrating it out.
g_conc <- function(theta, d) d$Zs * (d$yc - theta * d$pc)

#' G2 and G3 data frame -> list for g_iv and g_conc. Columns y, log_price, z1 ... zK.
iv_data <- function(df) {
  Zs <- as.matrix(df[, grep("^z[0-9]+$", names(df)), drop = FALSE])
  dimnames(Zs) <- NULL
  list(y = df$y, X = cbind(1, df$log_price), Z = cbind(1, Zs), Zs = Zs,
       yc = df$y - mean(df$y), pc = df$log_price - mean(df$log_price))
}

#' Block G1 market slopes: per market m, the OLS slope pi_hat_m of y (log weekly lessons booked) on log price and its
#' HC0 standard error sigma_hat_m (equal to sandwich::vcovHC(lm, type = "HC0")). d: market (factor), log_price, y.
market_slopes <- function(d) {
  out <- lapply(split(d, d$market, drop = TRUE), function(s) {
    x <- s$log_price - mean(s$log_price)
    b <- sum(x * s$y) / sum(x^2)
    e <- (s$y - mean(s$y)) - b * x
    data.frame(n = nrow(s), pi_hat = b, se = sqrt(sum(x^2 * e^2)) / sum(x^2))
  })
  data.frame(market = factor(names(out), levels = levels(d$market)), do.call(rbind, out), row.names = NULL)
}

#' Block G1 pooled micro data for GMM: list(y, lp, mk = market index, r = number of markets,
#' Z = (market dummies, market dummies x log price), n x 2r).
pool_data <- function(d) {
  mk <- as.integer(d$market)
  r <- nlevels(d$market)
  Dm <- outer(mk, seq_len(r), "==") * 1
  list(y = d$y, lp = d$log_price, mk = mk, r = r, Z = cbind(Dm, Dm * d$log_price))
}

#' Block G1 pooled moments: g_i(b) = z_i (y_i - a_{m(i)} - theta log p_i), b = (a_1, ..., a_r, theta): two moments
#' per market (its intercept and its price variation), one common elasticity; r + 1 parameters, 2r moments.
g_pool <- function(b, d) d$Z * (d$y - b[d$mk] - b[d$r + 1L] * d$lp)

## =================================================================
## Tamarind data generators (parameters and seeds: 01_dgp.R, stored in truth07.rds$params)
## =================================================================

#' G1 and G1h: one price A/B test per metro market. In market m, n[m] households (an even number) are split at
#' random into two equal arms: A pays the base price p0[m] (USD per lesson credit), B pays p0[m] (1 + delta[m]).
#' y = log of the household's average weekly lessons booked over the 8-week test,
#'   y = log(q0[m]) + theta_m (log price - log p0[m]) + sigma[m] s e,  e ~ N(0, 1), s = 1 in arm A and het in arm B,
#' so q0[m] is the median weekly lessons booked at the base price and theta_m the price elasticity (G1: theta_m = theta for
#' every market; G1h: theta_m varies). Draw order per market: assignment, then e.
gen_G1 <- function(par, theta_m = rep(par$theta, length(par$markets))) {
  r <- length(par$markets)
  out <- lapply(seq_len(r), function(m) {
    nm <- par$n[m]
    B <- sample(rep(0:1, length.out = nm))
    dlp <- B * log(1 + par$delta[m])
    s <- ifelse(B == 1, par$het, 1)
    y <- log(par$q0[m]) + theta_m[m] * dlp + par$sigma[m] * s * rnorm(nm)
    data.frame(market = m, B = B, log_price = log(par$p0[m]) + dlp, y = y)
  })
  d <- do.call(rbind, out)
  d$market <- factor(par$markets[d$market], levels = par$markets)
  d
}

#' G2 and G3: n regions in one week. K concurrent pricing experiments; experiment j puts exactly n share[j]
#' regions in its upper arm (T_j = 1) and the rest in its lower arm, the two arms differing by step[j] in the
#' experiment's own unit (USD or percentage points). The instrument is the centred assignment
#' z_j = step[j] (T_j - share[j]), so each z_j has sample mean exactly 0 (a balanced test, share 0.5, gives +/- step / 2).
#' Log price paid (USD per lesson credit):
#'   log_price = log(p_base) + sum_j gamma[j] z_j + v,   v = sigma_v (rho e + sqrt(1 - rho^2) xi)
#' (gamma[j]: change in log price per unit; v: the regional manager's discretionary promotions, correlated with the demand
#' shock). Log weekly lesson bookings per 1,000 households:
#'   y = log(q_base) + theta (log_price - log(p_base)) + sigma_i e,   sigma_i = sigma_eps (1 + het T_{het_on}),
#' with e a unit-variance draw: normal, or when mix_p > 0 a scale mixture (with probability mix_p the draw is mix_scale
#' times larger), rescaled to unit variance. Draw order: the K assignments, then e, then xi, then the mixture indicators.
gen_G2 <- function(par, n = par$n) {
  K <- length(par$step)
  n_up <- round(n * par$share)
  if (any(abs(n_up - n * par$share) > 1e-9)) stop("gen_G2: n * share must be whole numbers")
  Tm <- vapply(seq_len(K), function(j) sample(rep(c(1, 0), c(n_up[j], n - n_up[j]))), numeric(n))
  Zs <- sweep(sweep(Tm, 2L, par$share), 2L, par$step, "*")
  e <- rnorm(n)
  xi <- rnorm(n)
  if (par$mix_p > 0) {
    big <- runif(n) < par$mix_p
    e <- e * ifelse(big, par$mix_scale, 1) / sqrt(1 - par$mix_p + par$mix_p * par$mix_scale^2)
  }
  v <- par$sigma_v * (par$rho * e + sqrt(1 - par$rho^2) * xi)
  dlp <- drop(Zs %*% par$gamma) + v
  sig <- par$sigma_eps * (1 + par$het * Tm[, par$het_on])
  d <- data.frame(y = log(par$q_base) + par$theta * dlp + sig * e, log_price = log(par$p_base) + dlp)
  for (j in seq_len(K)) d[[paste0("z", j)]] <- Zs[, j]
  d$v_true <- v
  d$eps_true <- sig * e
  d
}
gen_G3 <- function(par, n = par$n) gen_G2(par, n)

#' Population moments of the G2 / G3 design (closed form; z_j independent across experiments, e independent of z with
#' E[e^2] = 1, sigma_i depends on T_{het_on} only), for moments g = z (y - beta0 - theta log p), z = (1, z_1, ..., z_K):
#'   D0 = -E[z x'] with x = (1, log p);  S0 = E[z z' sigma_i^2]  (the registry's D_0 and S_0)
#' and the asymptotic standard errors (at sample size par$n) of GMM with W = E[zz']^-1 (2SLS), S0^-1 (efficient)
#' and I_r (identity), the AGS sensitivity Lambda_0 = -(D0'WD0)^-1 D0'W of theta (efficient W), the concentration
#' parameter mu^2 = n sum_j Var(gamma_j z_j) / sigma_v^2, the first-stage R^2 and the OLS probability limit.
pop_moments <- function(par) {
  K <- length(par$step)
  q <- par$share
  vz <- par$step^2 * q * (1 - q)                                   # Var(z_j)
  eff2 <- par$gamma^2 * vz                                         # variance of experiment j's price effect
  s2_lo <- par$sigma_eps^2
  s2_hi <- (par$sigma_eps * (1 + par$het))^2
  h <- par$het_on
  Es2 <- q[h] * s2_hi + (1 - q[h]) * s2_lo
  Ezz <- diag(c(1, vz))
  S0 <- Es2 * Ezz
  S0[1, h + 1] <- S0[h + 1, 1] <- par$step[h] * q[h] * (1 - q[h]) * (s2_hi - s2_lo)
  S0[h + 1, h + 1] <- par$step[h]^2 * (q[h] * (1 - q[h])^2 * s2_hi + (1 - q[h]) * q[h]^2 * s2_lo)
  D0 <- -cbind(c(1, rep(0, K)), c(log(par$p_base), par$gamma * vz))
  avar <- function(W) { B <- solve(t(D0) %*% W %*% D0); B %*% t(D0) %*% W %*% S0 %*% W %*% D0 %*% B }
  Esig <- par$sigma_eps * (1 + par$het * q[h])
  var_lp <- sum(eff2) + par$sigma_v^2
  We <- solve(S0)
  list(beta0 = log(par$q_base) - par$theta * log(par$p_base), eff_sd = sqrt(eff2), D0 = D0, S0 = S0, Ezz = Ezz,
       mu2 = par$n * sum(eff2) / par$sigma_v^2, r2_first_stage = sum(eff2) / var_lp,
       ols_plim = par$theta + par$sigma_v * par$rho * Esig / var_lp,
       se_2sls = sqrt(avar(solve(Ezz))[2, 2] / par$n), se_eff = sqrt(avar(We)[2, 2] / par$n),
       se_identity = sqrt(avar(diag(K + 1))[2, 2] / par$n),
       sensitivity = -(solve(t(D0) %*% We %*% D0) %*% t(D0) %*% We)[2, ])
}
