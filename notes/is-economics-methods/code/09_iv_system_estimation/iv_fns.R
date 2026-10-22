# ------------------------------------------------------------------
# iv_fns.R -- functions for chapter 9 (instrumental variables and system estimation).
#
# Notation (notation_v3.md 1.4 and 1.7): structural equation y_i = x_i' beta + eps_i with instruments
# z_i (r x 1), E[z_i eps_i] = 0 and rank E[z_i x_i'] = k; stacked over observations y, X (n x k), Z (n x r),
# P_Z = Z (Z'Z)^-1 Z', X_hat = P_Z X. The single-equation first stage is x_i' = z_i' Pi_fs + v_i'
# (Pi_fs is r x k). A system has M equations (index m); Omega = E[eps_i eps_i'] (M x M) is the error
# variance matrix, Omega_hat its estimate (the card's `Sigma`). A simultaneous system is
# y_i' Gamma + z_i' Delta + eps_i' = 0 with reduced form y_i' = z_i' Pi + v_i', Pi = -Delta Gamma^-1
# (r x M), v_i' = -eps_i' Gamma^-1; Pi_(m) is the identification submatrix of equation m (the card's
# `Pi_sub`). The Anderson-Rubin null value is beta_circ (H0: beta = beta_circ); F_fs and F_eff are the
# first-stage and Montiel Olea-Pflueger (2013) effective F. Wooldridge (2010) writes G equations, u for
# errors and Sigma for the structural error variance; the book writes M, eps and Omega.
#
# Chapter card interface:
#   tsls(y, X, Z, vcov = c("robust", "homo"))       -> list(b, V)          2SLS, Wooldridge eqs. 5.22, 5.27, 5.34
#   sur_fgls(ylist, Xlist, Sigma = NULL)            -> list(b, V, Sigma)   SUR by FGLS, eqs. 7.40, 7.54 (Theorem 7.4)
#   three_sls(ylist, Xlist, Z, Sigma = NULL, vcov = c("homo", "robust"))
#                                                   -> list(b, V, Sigma)   3SLS, eqs. 8.39, 8.40, 8.45 (Theorem 8.4)
#   order_condition(n_excluded_exog, n_included_endog) -> logical          Theorem 9.1 (eq. 9.11)
#   rank_condition(Pi_sub, tol = 1e-8)              -> list(rank, ok)      eq. 9.10 (full column rank of Pi_(m))
# (sur_fgls and three_sls take an optional Sigma = Omega_hat, used for Theorems 7.5 and section 9.3.2.)
# Helpers used by 01 to 05 and the tests:
#   ols_fit, partial_out, ci95, covers_num         OLS with HC0 or homoskedastic variance; residuals on W; est -/+ 1.96 se; 0/1 coverage
#   first_stage_F, eff_F                            F_fs (homoskedastic) and F_eff (Montiel Olea and Pflueger 2013, eq. 3)
#   ar_test, ar_ci, ar_covers                       Anderson-Rubin statistic AR(beta_circ) and its confidence set (one endogenous regressor)
#   hansen_j, sargan                                robust Hansen (1982) J (two-step GMM, eq. 8.58) and Sargan's N R^2_u (eq. 6.31)
#   cf_test                                         control-function form of the endogeneity test (eqs. 6.15, 6.22)
#   reduced_form, pop_avar_sys                      Pi = -Delta Gamma^-1; population Avar of 2SLS and 3SLS
#   s1_system, s4_mats                              design matrices of the Bramble blocks
#   gen_S1, gen_S1m, gen_S2, gen_S3, gen_S4, gen_S4v  the Bramble data generators (parameters and seeds in 01_dgp.R)
# Only base R (stats) is used here.
#
# Arguments and the registered symbols they stand for (notation_v3.md 1.7):
#   y            $\mathbf{y}$ (n x 1)                X       $\mathbf{X}$ (n x k), rows $\mathbf{x}_i'$
#   Z            $\mathbf{Z}$ (n x r), rows $\mathbf{z}_i'$; in three_sls the common instruments, $\mathbf{Z}_i = I_M \otimes \mathbf{z}_i'$
#   ylist[[m]]   column m of the n x M matrix with rows $\mathbf{y}_i'$     Xlist[[m]]  n x k_m, rows $\mathbf{x}_{im}'$
#   Sigma        $\hat\Omega$ (M x M)                b       $\hat\beta$ ($\hat\beta_{2SLS}$, $\hat\beta_{FGLS}$, $\hat\beta_{3SLS}$)
#   V            $\widehat{\mathrm{Avar}}(\hat\beta)$, the estimated variance of $\hat\beta$ (already divided by n)
#   Pi_sub       $\Pi_{(m)}$                         n_excluded_exog, n_included_endog  rows and columns of $\Pi_{(m)}$
#   Gamma, Delta $\Gamma$ (M x M), $\Delta$ (r x M)   reduced_form() returns $\Pi$ (r x M)
#   Pilist[[m]]  $\Pi_{\mathrm{fs}}$ of equation m   Q       $\mathbb{E}[\mathbf{z}_i\mathbf{z}_i']$     Omega  $\Omega$
#   x, Zx, W     the endogenous regressor (one column of $\mathbf{X}$), the excluded instruments and the included
#                exogenous columns of $\mathbf{Z}$; first_stage_F() returns $F_{\mathrm{fs}}$ and pi = the Zx block of $\Pi_{\mathrm{fs}}$,
#                eff_F() returns $F_{\mathrm{eff}}$
#   beta0        $\beta^{\circ}$ in $\mathrm{AR}(\beta^{\circ})$ (H0: $\beta = \beta^{\circ}$)
#   cf_test      v_hat = $\hat{\mathbf{v}}_i$, rho = $\rho$ in $y_i = \mathbf{x}_i'\beta + \hat{\mathbf{v}}_i'\rho + e_i$
# ------------------------------------------------------------------

## =================================================================
## Card functions
## =================================================================

#' Two-stage least squares. y: the n-vector y; X: the n x k regressor matrix X (exogenous columns
#' included); Z: the n x r instrument matrix Z (r >= k, the exogenous columns of X included).
#' b = (X_hat'X)^-1 X_hat'y with X_hat = P_Z X (Wooldridge 2010 eq. 5.22): beta_hat_2SLS.
#' V is Avar_hat(beta_hat_2SLS): "robust" = (X_hat'X_hat)^-1 (sum_i x_hat_i x_hat_i' eps_hat_i^2)
#' (X_hat'X_hat)^-1 (eq. 5.34, HC0); "homo" = sigma_hat^2 (X_hat'X_hat)^-1 with sigma_hat^2 = sum eps_hat_i^2 / (n - k)
#' (eqs. 5.26 and 5.27; the same as ivreg's vcov). eps_hat = y - X b uses X, not X_hat.
#' Computed by QR (b is the least-squares coefficient of y on X_hat, since X_hat'X = X_hat'X_hat), which keeps
#' rounding at the level of ivreg's when the regressors have large means.
tsls <- function(y, X, Z, vcov = c("robust", "homo")) {
  vcov <- match.arg(vcov)
  y <- as.numeric(y); X <- as.matrix(X); Z <- as.matrix(Z)
  n <- length(y); k <- ncol(X)
  if (nrow(X) != n || nrow(Z) != n) stop("tsls: y, X and Z need the same number of rows")
  if (ncol(Z) < k) stop("tsls: fewer instruments (", ncol(Z), ") than regressors (", k, "): the order condition fails")
  qz <- qr(Z)
  if (qz$rank < ncol(Z)) stop("tsls: Z is not of full column rank")
  Xh <- qr.fitted(qz, X)
  qx <- qr(Xh)
  if (qx$rank < k) stop("tsls: P_Z X is rank deficient (the rank condition fails in the sample)")
  Ai <- qr_inv_crossprod(qx)
  b <- drop(qr.coef(qx, y))
  e <- drop(y - X %*% b)
  V <- if (vcov == "robust") Ai %*% crossprod(Xh * e) %*% Ai else sum(e^2) / (n - k) * Ai
  names(b) <- colnames(X)
  dimnames(V) <- list(colnames(X), colnames(X))
  list(b = b, V = V)
}

#' Seemingly unrelated regressions by feasible GLS (Zellner 1962; Wooldridge 2010 section 7.5).
#' ylist: M outcome vectors (y_m, n each); Xlist: M regressor matrices X_m (n x k_m), so that the
#' system is y_i = X_i beta + eps_i with X_i block diagonal (k = sum k_m). Sigma = Omega_hat (M x M); by
#' default n^-1 sum_i eps_hat_i eps_hat_i' from equation-by-equation OLS residuals (eq. 7.40, no
#' degrees-of-freedom correction, systemfit's methodResidCov = "noDfCor").
#' b = [X'(Omega_hat^-1 (x) I_n) X]^-1 X'(Omega_hat^-1 (x) I_n) y, computed block by block;
#' V = [X'(Omega_hat^-1 (x) I_n) X]^-1, the usual FGLS variance under SGLS.3 (Theorem 7.4, eq. 7.54).
sur_fgls <- function(ylist, Xlist, Sigma = NULL) {
  M <- length(ylist)
  if (length(Xlist) != M || M < 1) stop("sur_fgls: ylist and Xlist need the same length")
  Xlist <- lapply(Xlist, as.matrix); ylist <- lapply(ylist, as.numeric)
  n <- length(ylist[[1]])
  if (is.null(Sigma)) {
    E <- vapply(seq_len(M), function(m) ylist[[m]] - qr.fitted(qr(Xlist[[m]]), ylist[[m]]), numeric(n))
    Sigma <- crossprod(E) / n
  }
  gls_blocks(ylist, Xlist, Xlist, Sigma)
}

#' Three-stage least squares (Zellner and Theil 1962; Wooldridge 2010 section 8.3.4) with common
#' instruments Z (n x r) for all M equations. Sigma = Omega_hat; by default n^-1 sum eps_hat_i eps_hat_i'
#' from equation-by-equation 2SLS residuals (eq. 8.44 with the first-step 2SLS estimates).
#' b = [X'(Omega_hat^-1 (x) P_Z) X]^-1 X'(Omega_hat^-1 (x) P_Z) y (eqs. 8.39 and 8.40 for common
#' instruments; systemfit's method3sls = "GLS"). V = [X'(Omega_hat^-1 (x) P_Z) X]^-1 (eq. 8.45, valid
#' under SIV.5, system homoskedasticity) or, with vcov = "robust", the GMM sandwich of the 3SLS
#' weight W = (Omega_hat (x) Z'Z / n)^-1: A^-1 (sum_i s_i s_i') A^-1 with A = X'(Omega_hat^-1 (x) P_Z) X and
#' block m of s_i equal to x_hat_im sum_l [Omega_hat^-1]_ml eps_hat_il (residuals at the 3SLS estimate);
#' it equals sandwich::vcovCL(systemfit_fit, cluster = market, type = "HC0", cadjust = FALSE, residFit = FALSE).
three_sls <- function(ylist, Xlist, Z, Sigma = NULL, vcov = c("homo", "robust")) {
  vcov <- match.arg(vcov)
  M <- length(ylist)
  if (length(Xlist) != M || M < 1) stop("three_sls: ylist and Xlist need the same length")
  Xlist <- lapply(Xlist, as.matrix); ylist <- lapply(ylist, as.numeric); Z <- as.matrix(Z)
  n <- nrow(Z)
  qz <- qr(Z)
  if (qz$rank < ncol(Z)) stop("three_sls: Z is not of full column rank")
  if (is.null(Sigma)) {
    E <- vapply(seq_len(M), function(m) {
      t2 <- tsls(ylist[[m]], Xlist[[m]], Z, vcov = "homo")
      ylist[[m]] - drop(Xlist[[m]] %*% t2$b)
    }, numeric(n))
    Sigma <- crossprod(E) / n
  }
  Xh <- lapply(Xlist, function(X) qr.fitted(qz, X))
  out <- gls_blocks(ylist, Xh, Xlist, Sigma)
  if (vcov == "robust") {
    Si <- solve(Sigma)
    idx <- block_index(Xlist)
    E <- vapply(seq_len(M), function(m) ylist[[m]] - drop(Xlist[[m]] %*% out$b[idx[[m]]]), numeric(n))
    W <- E %*% Si                                                    # row i: sum_l [Omega_hat^-1]_ml eps_hat_il
    S <- do.call(cbind, lapply(seq_len(M), function(m) Xh[[m]] * W[, m]))
    Ai <- out$V
    out$V <- Ai %*% crossprod(S) %*% Ai
    dimnames(out$V) <- dimnames(Ai)
  }
  out
}

#' Order condition for equation m under exclusion restrictions (Wooldridge 2010 Theorem 9.1, eq. 9.11):
#' at least as many excluded exogenous variables as included right-hand-side endogenous variables.
#' n_excluded_exog: number of exogenous variables excluded from equation m (rows of Pi_(m));
#' n_included_endog: number of endogenous variables on its right-hand side (columns of Pi_(m)).
order_condition <- function(n_excluded_exog, n_included_endog) {
  ok_count <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0 && x == round(x)
  if (!ok_count(n_excluded_exog) || !ok_count(n_included_endog))
    stop("order_condition: the arguments are non-negative whole numbers")
  n_excluded_exog >= n_included_endog
}

#' Rank condition for equation m under exclusion restrictions (Wooldridge 2010 eq. 9.10; registry 1.7):
#' Pi_sub = Pi_(m), the block of the reduced form Pi with the exogenous variables excluded from equation m
#' as rows and the right-hand-side endogenous variables of equation m as columns. The equation is
#' identified when Pi_(m) has full column rank. rank = number of singular values above
#' tol * max(1, largest singular value) (an absolute tolerance for matrices with entries of order one or less);
#' ok = (rank == ncol(Pi_sub)), which fails automatically when the order condition fails (fewer rows than columns).
rank_condition <- function(Pi_sub, tol = 1e-8) {
  Pi_sub <- as.matrix(Pi_sub)
  if (!is.numeric(Pi_sub) || anyNA(Pi_sub)) stop("rank_condition: Pi_sub must be a numeric matrix without NA")
  d <- svd(Pi_sub, nu = 0, nv = 0)$d
  rank <- sum(d > tol * max(1, d))
  list(rank = as.integer(rank), ok = rank == ncol(Pi_sub))
}

## =================================================================
## Shared algebra
## =================================================================

#' (A'A)^-1 from the QR decomposition of a full-column-rank A (undoing LINPACK's column pivoting).
qr_inv_crossprod <- function(qx) {
  k <- ncol(qx$qr)
  Ri <- chol2inv(qr.R(qx))
  out <- matrix(0, k, k)
  out[qx$pivot, qx$pivot] <- Ri
  out
}

## GLS for systems (SUR and 3SLS)
block_index <- function(Xlist) {
  k <- vapply(Xlist, ncol, integer(1))
  split(seq_len(sum(k)), rep(seq_along(k), k))
}

#' GLS normal equations for a system: A[m, l] = [Omega_hat^-1]_ml L_m' R_l, c_m = sum_l [Omega_hat^-1]_ml L_m' y_l,
#' b = A^-1 c, V = A^-1. SUR: L = R = X; 3SLS: L = X_hat (= P_Z X), R = X. Coefficient names are
#' <equation>:<column> when ylist is named, otherwise eq<m>:<column>.
gls_blocks <- function(ylist, Lefts, Rights, Sigma) {
  M <- length(ylist)
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(M, M))) stop("Sigma must be M x M")
  Si <- solve(Sigma)
  idx <- block_index(Rights)
  K <- sum(lengths(idx))
  A <- matrix(0, K, K); cc <- numeric(K)
  for (m in seq_len(M)) for (l in seq_len(M)) {
    A[idx[[m]], idx[[l]]] <- Si[m, l] * crossprod(Lefts[[m]], Rights[[l]])
    cc[idx[[m]]] <- cc[idx[[m]]] + Si[m, l] * drop(crossprod(Lefts[[m]], ylist[[l]]))
  }
  A <- (A + t(A)) / 2                                                 # symmetric up to rounding (X_hat'X = X_hat'X_hat)
  V <- solve(A)
  b <- drop(V %*% cc)
  eq <- if (is.null(names(ylist))) paste0("eq", seq_len(M)) else names(ylist)
  cn <- unlist(lapply(seq_len(M), function(m) {
    nm <- colnames(Rights[[m]]); if (is.null(nm)) nm <- paste0("x", seq_len(ncol(Rights[[m]])))
    paste0(eq[m], ":", nm)
  }))
  names(b) <- cn
  dimnames(V) <- list(cn, cn)
  list(b = b, V = V, Sigma = Sigma)
}

## =================================================================
## Single-equation helpers
## =================================================================

#' OLS of y on X: b and V (HC0 sandwich, or sigma_hat^2 (X'X)^-1 with n - k degrees of freedom).
ols_fit <- function(y, X, vcov = c("robust", "homo")) {
  vcov <- match.arg(vcov)
  X <- as.matrix(X); y <- as.numeric(y)
  qx <- qr(X)
  if (qx$rank < ncol(X)) stop("ols_fit: X is rank deficient")
  XXi <- qr_inv_crossprod(qx)
  b <- drop(qr.coef(qx, y))
  e <- drop(y - X %*% b)
  V <- if (vcov == "robust") XXi %*% crossprod(X * e) %*% XXi else sum(e^2) / (length(y) - ncol(X)) * XXi
  names(b) <- colnames(X)
  dimnames(V) <- list(colnames(X), colnames(X))
  list(b = b, V = V)
}

#' Residuals of each column of A on W (Frisch-Waugh-Lovell).
partial_out <- function(A, W) {
  A <- as.matrix(A)
  A - qr.fitted(qr(as.matrix(W)), A)
}

ci95 <- function(est, se) unname(c(est - qnorm(0.975) * se, est + qnorm(0.975) * se))
covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])          # 1 when c(lo, hi) contains v

#' First-stage F_fs for one endogenous regressor x: H0 that the excluded instruments Zx (n x r_x) have
#' zero coefficients in the regression of x on (W, Zx); W holds the included exogenous variables
#' (constant included). Homoskedastic F with (r_x, n - ncol(W) - r_x) degrees of freedom (the same as
#' ivreg's "Weak instruments" diagnostic). Also returns pi_hat (the excluded-instrument block of Pi_fs).
first_stage_F <- function(x, Zx, W) {
  xt <- partial_out(x, W); Zt <- partial_out(Zx, W)
  n <- length(xt); kz <- ncol(Zt); df2 <- n - ncol(as.matrix(W)) - kz
  pi <- qr.coef(qr(Zt), xt)
  v <- xt - Zt %*% pi
  Fs <- (sum((Zt %*% pi)^2) / kz) / (sum(v^2) / df2)
  list(F = Fs, df1 = kz, df2 = df2, p = pf(Fs, kz, df2, lower.tail = FALSE), pi = drop(pi))
}

#' Montiel Olea and Pflueger (2013) effective F, eq. (3): F_eff = n pi_hat' Q pi_hat / tr(W_hat Q) after
#' partialling W out of x and Zx, with Q = Zx~'Zx~ / n and W_hat the estimated asymptotic variance of
#' sqrt(n)(pi_hat - pi): HC0 ("robust") or sigma_hat^2 Q^-1 with n - ncol(W) - r_x degrees of freedom ("homo",
#' where F_eff equals first_stage_F). With one instrument the robust F_eff is the robust Wald statistic t^2.
#' Invariant to nonsingular linear transformations of the instruments.
eff_F <- function(x, Zx, W, vcov = c("robust", "homo")) {
  vcov <- match.arg(vcov)
  xt <- partial_out(x, W); Zt <- partial_out(Zx, W)
  n <- length(xt); kz <- ncol(Zt)
  ZZ <- crossprod(Zt); ZZi <- solve(ZZ)
  pi <- ZZi %*% crossprod(Zt, xt)
  v <- drop(xt - Zt %*% pi)
  Wh <- if (vcov == "robust") n * ZZi %*% crossprod(Zt * v) %*% ZZi
        else sum(v^2) / (n - ncol(as.matrix(W)) - kz) * n * ZZi
  Q <- ZZ / n
  drop(n * t(pi) %*% Q %*% pi) / sum(diag(Wh %*% Q))
}

#' Anderson-Rubin statistic AR(beta_circ) for H0: beta = beta_circ with one endogenous regressor x,
#' excluded instruments Zx (n x r_x) and included exogenous W (constant included): the F statistic of Zx
#' in the regression of y - beta_circ x on (W, Zx), with (r_x, n - ncol(W) - r_x) degrees of freedom
#' (homoskedastic; exact under normal errors; the S statistic of chapter 7 in the linear IV case).
#' Equals ivmodel::AR.test(..., beta0 = beta_circ)$Fstat.
ar_test <- function(y, x, Zx, W, beta0) {
  yt <- partial_out(y, W); xt <- partial_out(x, W); Zt <- partial_out(Zx, W)
  n <- length(yt); kz <- ncol(Zt); df2 <- n - ncol(as.matrix(W)) - kz
  u <- yt - beta0 * xt
  pu <- qr.fitted(qr(Zt), u)
  stat <- (sum(pu^2) / kz) / ((sum(u^2) - sum(pu^2)) / df2)
  list(stat = drop(stat), df1 = kz, df2 = df2, p = pf(drop(stat), kz, df2, lower.tail = FALSE))
}

#' Anderson-Rubin confidence set {beta_circ : AR(beta_circ) <= F_{r_x, df2, 1 - alpha}} in closed form.
#' With u = y~ - b x~ (W partialled out) the condition is u'P u - a u'(I - P) u <= 0, a = c r_x / df2, a quadratic
#' A b^2 + B b + C <= 0 in b. type: "bounded" [lo, hi] (A > 0), "two_rays" (-Inf, lo] U [hi, Inf) (A < 0 and two
#' real roots), "real_line" (A < 0, no real root), "empty" (A > 0, no real root; possible only when r_x > 1).
ar_ci <- function(y, x, Zx, W, alpha = 0.05) {
  yt <- partial_out(y, W); xt <- partial_out(x, W); Zt <- partial_out(Zx, W)
  n <- length(yt); kz <- ncol(Zt); df2 <- n - ncol(as.matrix(W)) - kz
  a <- qf(1 - alpha, kz, df2) * kz / df2
  qz <- qr(Zt)
  Y <- cbind(yt, xt); PY <- qr.fitted(qz, Y)
  Mp <- crossprod(Y, PY); Mt <- crossprod(Y)
  Mr <- Mt - Mp
  A <- Mp[2, 2] - a * Mr[2, 2]; B <- -2 * (Mp[1, 2] - a * Mr[1, 2]); C <- Mp[1, 1] - a * Mr[1, 1]
  disc <- B^2 - 4 * A * C
  if (disc >= 0) r <- sort(c((-B - sqrt(disc)) / (2 * A), (-B + sqrt(disc)) / (2 * A)))
  if (A > 0) {
    if (disc >= 0) list(type = "bounded", lo = r[1], hi = r[2]) else list(type = "empty", lo = NA_real_, hi = NA_real_)
  } else {
    if (disc > 0) list(type = "two_rays", lo = r[1], hi = r[2]) else list(type = "real_line", lo = -Inf, hi = Inf)
  }
}

#' Does an ar_ci() set contain b?
ar_covers <- function(ci, b) {
  switch(ci$type,
         bounded = ci$lo <= b & b <= ci$hi,
         two_rays = b <= ci$lo | b >= ci$hi,
         real_line = rep(TRUE, length(b)),
         empty = rep(FALSE, length(b)))
}

#' Hansen (1982) J statistic, heteroskedasticity robust (Wooldridge 2010 eq. 8.58): the first step is 2SLS,
#' S_hat = n^-1 sum z_i z_i' eps_hat_i^2 (uncentered, at the 2SLS residuals); the second step is efficient GMM
#' b = (X'Z S_hat^-1 Z'X)^-1 X'Z S_hat^-1 Z'y and J = n g_bar(b)' S_hat^-1 g_bar(b), g_bar(b) = n^-1 Z'(y - X b),
#' asymptotically chi-square with r - k degrees of freedom. Equals gmm::gmm(type = "twoStep", vcov = "MDS",
#' centeredVcov = FALSE)'s specTest. y, X, Z as in tsls().
hansen_j <- function(y, X, Z) {
  y <- as.numeric(y); X <- as.matrix(X); Z <- as.matrix(Z)
  n <- length(y); df <- ncol(Z) - ncol(X)
  t1 <- tsls(y, X, Z, vcov = "homo")
  e1 <- drop(y - X %*% t1$b)
  Si <- solve(crossprod(Z * e1) / n)
  ZX <- crossprod(Z, X) / n; Zy <- crossprod(Z, y) / n
  b <- drop(solve(t(ZX) %*% Si %*% ZX, t(ZX) %*% Si %*% Zy))
  g <- drop(Zy - ZX %*% b)
  J <- n * drop(t(g) %*% Si %*% g)
  names(b) <- colnames(X)
  list(J = J, df = df, p = if (df > 0) pchisq(J, df, lower.tail = FALSE) else NA_real_, b = b)
}

#' Sargan (1958) over-identification statistic in Wooldridge's (2010) eq. 6.31 form: n R^2_u from the
#' regression of the 2SLS residuals on all instruments Z; chi-square with r - k degrees of freedom under
#' homoskedasticity (2SLS.3). Equals ivreg's summary(diagnostics = TRUE) "Sargan". NA when just identified.
sargan <- function(y, X, Z) {
  y <- as.numeric(y); X <- as.matrix(X); Z <- as.matrix(Z)
  df <- ncol(Z) - ncol(X)
  if (df == 0) return(list(stat = NA_real_, df = 0, p = NA_real_))
  t1 <- tsls(y, X, Z, vcov = "homo")
  e <- drop(y - X %*% t1$b)
  f <- qr.fitted(qr(Z), e)
  stat <- length(y) * sum(f^2) / sum(e^2)
  list(stat = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
}

#' Control-function (regression) form of the endogeneity test (Wooldridge 2010 eqs. 6.15 and 6.22):
#' regress each endogenous column of X (names in `endog`) on Z to get the first-stage residuals v_hat,
#' then regress y on (X, v_hat): y_i = x_i' beta + v_hat_i' rho + e_i. H0: rho = 0 (exogeneity); HC0 t (one
#' endogenous regressor) or robust Wald (several). The coefficients on X equal the 2SLS estimates.
cf_test <- function(y, X, Z, endog) {
  X <- as.matrix(X)
  V <- partial_out(X[, endog, drop = FALSE], Z)
  colnames(V) <- paste0("vhat_", endog)
  f <- ols_fit(y, cbind(X, V))
  j <- colnames(V)
  rho <- f$b[j]; Vr <- f$V[j, j, drop = FALSE]
  W <- drop(t(rho) %*% solve(Vr) %*% rho)
  list(rho = unname(rho), se = unname(sqrt(diag(Vr))), t = if (length(j) == 1) unname(rho / sqrt(Vr[1, 1])) else NA_real_,
       wald = W, p = pchisq(W, length(j), lower.tail = FALSE), b = f$b[colnames(X)])
}

## =================================================================
## Simultaneous-system helpers
## =================================================================

#' Reduced form of y_i' Gamma + z_i' Delta + eps_i' = 0: Pi = -Delta Gamma^-1 (r x M; Wooldridge 2010 eq. 9.14).
#' Gamma: M x M, rows = the endogenous variables y, columns = the equations; Delta: r x M, rows = the exogenous
#' variables z. Pi has rows = z and columns = y (column j of Pi is the reduced form of y_j).
reduced_form <- function(Gamma, Delta) {
  Pi <- -Delta %*% solve(Gamma)
  dimnames(Pi) <- list(rownames(Delta), rownames(Gamma))
  Pi
}

#' Population asymptotic variances of sqrt(n)(beta_hat - beta) under system homoskedasticity (SIV.5) when every
#' equation uses the instruments z_i with Q = E[z_i z_i'] (r x r). Pilist[[m]] = Pi_fs of equation m
#' (r x k_m), the linear projection coefficients of x_im on z_i, so x_hat_im = Pi_fs' z_i; Omega = Omega (M x M).
#' 2SLS of equation m: Omega_mm (Pi_m' Q Pi_m)^-1 (Theorem 5.2, eq. 5.24); V_2sls_joint is the joint variance of
#' the equation-by-equation 2SLS estimators (block [m, l] = Omega_ml H_m^-1 Pi_m' Q Pi_l H_l^-1, H_m = Pi_m' Q Pi_m);
#' 3SLS: A^-1 with A[m, l] = [Omega^-1]_ml Pi_m' Q Pi_l (Theorem 8.4). V_2sls_joint - V_3sls is positive semidefinite.
pop_avar_sys <- function(Pilist, Q, Omega) {
  M <- length(Pilist)
  H <- lapply(seq_len(M), function(m) solve(t(Pilist[[m]]) %*% Q %*% Pilist[[m]]))
  V2 <- lapply(seq_len(M), function(m) Omega[m, m] * H[[m]])
  Oi <- solve(Omega)
  idx <- block_index(Pilist)
  K <- sum(lengths(idx)); A <- matrix(0, K, K); VJ <- matrix(0, K, K)
  for (m in seq_len(M)) for (l in seq_len(M)) {
    G <- t(Pilist[[m]]) %*% Q %*% Pilist[[l]]
    A[idx[[m]], idx[[l]]] <- Oi[m, l] * G
    VJ[idx[[m]], idx[[l]]] <- Omega[m, l] * H[[m]] %*% G %*% H[[l]]
  }
  list(V_2sls = V2, V_2sls_joint = (VJ + t(VJ)) / 2, V_3sls = solve((A + t(A)) / 2))
}

#' Design of block S1 (and S1m, S3): instruments z = (const, ad_budget, coupon, fee, log_pop); buyer equation
#' regressors (const, n_sellers, ad_budget, coupon, log_pop), seller equation (const, buyers, fee, log_pop).
s1_system <- function(d) {
  Z <- cbind(const = 1, ad_budget = d$ad_budget, coupon = d$coupon, fee = d$fee, log_pop = d$log_pop)
  XB <- cbind(const = 1, n_sellers = d$n_sellers, ad_budget = d$ad_budget, coupon = d$coupon, log_pop = d$log_pop)
  XS <- cbind(const = 1, buyers = d$buyers, fee = d$fee, log_pop = d$log_pop)
  list(y = list(buyer = d$buyers, seller = d$n_sellers), X = list(buyer = XB, seller = XS), Z = Z)
}

#' Design of block S4: the buyer equation with the three seller-side programs as excluded instruments.
s4_mats <- function(d) {
  W <- cbind(const = 1, ad_budget = d$ad_budget, coupon = d$coupon, log_pop = d$log_pop)
  list(y = d$buyers, X = cbind(W[, 1, drop = FALSE], n_sellers = d$n_sellers, W[, -1]),
       W = W, Zx = cbind(fee = d$fee, payout_days = d$payout_days, referral_bonus = d$referral_bonus),
       x = d$n_sellers)
}

## =================================================================
## Bramble data generators (parameters and seeds: 01_dgp.R, stored in truth09.rds$params)
## =================================================================

#' S1 (and S1m with kappa != 0, S3 with r = 0 and a small fee_sd): one month in n local markets.
#' Draw order: log_pop, ad_budget, u, w, then the two structural errors.
#'   log_pop ~ N(0, 1); ad_budget = ad0 + ad_sd N(0, 1) (USD thousands); coupon = cp0 + cp_sd u (USD);
#'   fee = fee0 + fee_sd (r u + sqrt(1 - r^2) w) (USD per month), so corr(coupon, fee) = r;
#'   (eps_B, eps_S) bivariate normal, sd sB and sS, correlation rho.
#'   Buyers (thousands):  buyers    = aB + bBS n_sellers + gad ad_budget + gcp coupon + lB log_pop + eps_B
#'   Sellers (hundreds):  n_sellers = aS + (bSB + kappa (coupon - cp0)) buyers + dfee fee + lS log_pop + eps_S
#' solved for the equilibrium market by market (1 - bBS (bSB + kappa (coupon - cp0)) > 0).
gen_S1 <- function(par, n = par$n, kappa = 0) {
  log_pop <- rnorm(n)
  ad_budget <- par$ad0 + par$ad_sd * rnorm(n)
  u <- rnorm(n); w <- rnorm(n)
  coupon <- par$cp0 + par$cp_sd * u
  fee <- par$fee0 + par$fee_sd * (par$r * u + sqrt(1 - par$r^2) * w)
  Om <- matrix(c(par$sB^2, par$rho * par$sB * par$sS, par$rho * par$sB * par$sS, par$sS^2), 2)
  e <- matrix(rnorm(2 * n), n) %*% chol(Om)
  cB <- par$aB + par$gad * ad_budget + par$gcp * coupon + par$lB * log_pop + e[, 1]
  slope <- par$bSB + kappa * (coupon - par$cp0)
  cS <- par$aS + par$dfee * fee + par$lS * log_pop + e[, 2]
  buyers <- (cB + par$bBS * cS) / (1 - par$bBS * slope)
  n_sellers <- cS + slope * buyers
  data.frame(market = seq_len(n), buyers = buyers, n_sellers = n_sellers, ad_budget = ad_budget, coupon = coupon,
             fee = fee, log_pop = log_pop, eps_buyer_true = e[, 1], eps_seller_true = e[, 2])
}
gen_S1m <- function(par, n = par$n) gen_S1(par, n, kappa = par$kappa)
gen_S3 <- function(par, n = par$n) gen_S1(par, n)

#' S2: the same kind of markets, two outcome equations with correlated errors.
#'   free_ship_min ~ N(ship0, ship_sd^2) (USD); premium_share = prem0 + prem_ship (free_ship_min - ship0) + prem_sd N(0, 1);
#'   loyalty = 1{loy_ship (free_ship_min - ship0) + N(0, 1) > loy_cut};
#'   basket_usd = b1[1] + b1[2] free_ship_min + b1[3] premium_share + eps_1;
#'   orders_pb  = b2[1] + b2[2] free_ship_min + b2[3] loyalty + eps_2; (eps_1, eps_2) normal, sd s1, s2, correlation rho.
gen_S2 <- function(par, n = par$n) {
  free_ship_min <- par$ship0 + par$ship_sd * rnorm(n)
  premium_share <- par$prem0 + par$prem_ship * (free_ship_min - par$ship0) + par$prem_sd * rnorm(n)
  loyalty <- as.integer(par$loy_ship * (free_ship_min - par$ship0) + rnorm(n) > par$loy_cut)
  Om <- matrix(c(par$s1^2, par$rho * par$s1 * par$s2, par$rho * par$s1 * par$s2, par$s2^2), 2)
  e <- matrix(rnorm(2 * n), n) %*% chol(Om)
  basket_usd <- par$b1[1] + par$b1[2] * free_ship_min + par$b1[3] * premium_share + e[, 1]
  orders_pb <- par$b2[1] + par$b2[2] * free_ship_min + par$b2[3] * loyalty + e[, 2]
  data.frame(market = seq_len(n), free_ship_min = free_ship_min, premium_share = premium_share, loyalty = loyalty,
             basket_usd = basket_usd, orders_pb = orders_pb, eps_basket_true = e[, 1], eps_orders_true = e[, 2])
}

#' S4 (theta > 0) and S4v (theta = 0): one quarter in which three seller-side programs are randomized
#' independently: fee ~ fee0 + fee_sd N(0, 1), payout_days ~ U(pay_lo, pay_hi), referral_bonus ~ U(ref_lo, ref_hi).
#' Draw order: log_pop, ad_budget, coupon, fee, payout_days, referral_bonus, then the errors.
#'   buyers    = aB + bBS n_sellers + gad ad_budget + gcp coupon + lB log_pop + theta (referral_bonus - ref_mid) + eps_B
#'   n_sellers = aS + bSB buyers + dfee fee + dpay payout_days + dref referral_bonus + lS log_pop + eps_S
#' Both errors are scaled by h_i = exp(eta zr_i) / sqrt(E[exp(2 eta zr)]), zr = standardized referral bonus
#' (uniform on +/- sqrt(3)), so E[h^2] = 1: heteroskedasticity tied to an instrument.
gen_S4 <- function(par, n = par$n) {
  log_pop <- rnorm(n)
  ad_budget <- par$ad0 + par$ad_sd * rnorm(n)
  coupon <- par$cp0 + par$cp_sd * rnorm(n)
  fee <- par$fee0 + par$fee_sd * rnorm(n)
  payout_days <- runif(n, par$pay_lo, par$pay_hi)
  referral_bonus <- runif(n, par$ref_lo, par$ref_hi)
  Om <- matrix(c(par$sB^2, par$rho * par$sB * par$sS, par$rho * par$sB * par$sS, par$sS^2), 2)
  e <- matrix(rnorm(2 * n), n) %*% chol(Om)
  zr <- (referral_bonus - par$ref_mid) / ((par$ref_hi - par$ref_lo) / sqrt(12))
  t2 <- 2 * par$eta * sqrt(3)
  h <- exp(par$eta * zr) / sqrt(sinh(t2) / t2)                       # E[exp(2 eta zr)] = sinh(2 eta sqrt 3) / (2 eta sqrt 3)
  eB <- e[, 1] * h; eS <- e[, 2] * h
  cB <- par$aB + par$gad * ad_budget + par$gcp * coupon + par$lB * log_pop + par$theta * (referral_bonus - par$ref_mid) + eB
  cS <- par$aS + par$dfee * fee + par$dpay * payout_days + par$dref * referral_bonus + par$lS * log_pop + eS
  buyers <- (cB + par$bBS * cS) / (1 - par$bBS * par$bSB)
  n_sellers <- cS + par$bSB * buyers
  data.frame(market = seq_len(n), buyers = buyers, n_sellers = n_sellers, ad_budget = ad_budget, coupon = coupon,
             fee = fee, payout_days = payout_days, referral_bonus = referral_bonus, log_pop = log_pop,
             eps_buyer_true = eB, eps_seller_true = eS, h_true = h)
}
gen_S4v <- gen_S4
