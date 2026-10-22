# ------------------------------------------------------------------
# sc_fns.R -- hand-written synthetic-control tools for chapter 22
#
#   sc_weights(X1, X0, V)                        -> numeric J-vector (simplex weights)
#   penalized_sc(X1, X0, pen, V = diag(length(X1)), tol = 1e-6)
#                                                -> list(w, n_nonzero)
#   conformal_p(y1, y0hat, T0, scheme = c("moving_block", "iid"), q = 1, n_perm = 999)
#                                                -> numeric p-value
#   mspe_ratio_perm(Y, treated, T0, fit_fun, placebo_donors = c("all_others", "exclude_treated"))
#                                                -> list(ratio, p, rank, all_ratios, n_units)
#   sc_fit(y1, Y0, T0, V = NULL)                 -> SC counterfactual path (a ready-made fit_fun)
#   in_hull_lp(x1, X0, tol = 1e-9)               -> TRUE / FALSE (attribute "l1_gap")
#
# Conventions (shared by every function):
#   * Donors are COLUMNS: X0 is k x J (k predictors, J donors), X1 has length k;
#     in time-series use X0 = t(Y_donors[, 1:T0]) is T0 x J. Y in mspe_ratio_perm
#     is N x T with units in rows (the layout of `Y` in aster_S*.rds).
#   * T0 is the chapter's T_0 (notation_v3.md 2.9(a)): the LAST untreated period,
#     i.e. the number of pre-treatment periods; the post periods are T0+1, ..., T.
#     The legacy S1 objects store meta$T0 = 25, the FIRST reform month, so on
#     aster_S1.rds call these functions with T0 = meta$T0 - 1 = 24.
#   * V is the chapter's positive semi-definite predictor weighting matrix Omega
#     (ADH and Synth's custom.v call it V); it may be given as a k x k matrix or
#     as the vector of its diagonal.
#
# Sources:
#   Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic control methods
#     for comparative case studies. JASA 105(490), 493-505. (weights on the
#     simplex; placebo MSPE-ratio permutation inference: the intervention is
#     reassigned to each control unit in turn, "shifting California to the
#     donor pool", so the treated unit stays in every placebo's donor pool;
#     mspe_ratio_perm's default)
#   Abadie, A. and L'Hour, J. (2021). A penalized synthetic control estimator for
#     disaggregated data. JASA 116(536), 1817-1834. Objective, for pen = the
#     paper's lambda (the chapter's zeta_pen, notation_v3.md 2.9(f)):
#       min_w ||X1 - X0 w||_V^2 + pen * sum_j w_j ||X1 - X0[, j]||_V^2,
#       w_j >= 0, sum_j w_j = 1.
#     pen = 0 is plain SC; pen -> infinity puts all weight on the nearest donor.
#     For pen > 0 the solution is unique (and sparse) under the paper's general
#     position condition.
#   Chernozhukov, V., Wuthrich, K. and Zhu, Y. (2021). An exact and robust
#     conformal inference method for counterfactual and synthetic controls.
#     JASA 116(536), 1849-1864. Test of H0: tau_t = 0 for every post period.
#     Under H0 the post outcomes are untreated outcomes, so the counterfactual
#     model is fitted on ALL T periods (pre and post) and the residuals
#     u_t = y1_t - y0hat_t are compared under permutations of time:
#       S(u) = ( (T - T0)^(-1/2) * sum_{t > T0} |u_t|^q )^(1/q),
#       p = (1 / |Pi|) * sum_{pi in Pi} 1{ S(u_pi) >= S(u) }.
#     "moving_block": Pi = the T cyclic shifts t -> t + j mod T (j = 0, ..., T-1;
#     p in {1/T, ..., 1}); exact under exchangeable errors (CWZ Theorem D.1) and
#     only approximately valid under stationary, weakly dependent errors (CWZ
#     Theorem 1; chapter 22, Theorem 8 (ii)).
#     "iid": Pi = all T! permutations (exact under exchangeable errors),
#     approximated here by the identity plus n_perm uniform random permutations,
#     which keeps the test exact (draws use the session RNG).
#     conformal_p takes y0hat as given: fitting it under H0 is the caller's job,
#     e.g. for SC on S1 (T = 36, T0 = 24):
#       w0 <- sc_weights(y1, t(Y[-1, ]), diag(36)); y0hat <- drop(t(Y[-1, ]) %*% w0)
#       conformal_p(y1, y0hat, T0 = 24)
#     A null tau_t = theta is tested by passing y1 - theta * (t > T0).
#
# Numerical method: the SC and penalized-SC problems are convex quadratic
# programs on the simplex. They are solved by quadprog::solve.QP (Goldfarb-
# Idnani) on the objective rescaled by the mean of the diagonal of X0'VX0 (the
# argmin does not change); when X0'VX0 is singular (more donors than
# predictors, or a singular V) a ridge eps * I is added, eps escalating over
# 1e-12, 1e-10, 1e-8, 1e-6 of that scale until the Cholesky step succeeds. With
# a near-singular X0'VX0 the Goldfarb-Idnani start (the unconstrained minimiser)
# is huge and that answer can be off by 1e-4 or miss the sparse optimum, so the
# same problem is also solved by an exact primal active-set method (vertex
# start, ratio tests, zero-curvature steps on singular faces, KKT check at the
# end), and the answer with the lower objective is returned (the active-set one
# on a tie). Checked against CVXR/Clarabel on S1, S2 and random singular and
# non-singular problems (Task 9 report). When the optimum is not unique (pen = 0
# with more donors than predictors and X1 inside the hull) any optimal w may be
# returned; for pen > 0 the optimum is unique in general position. Returned
# weights are clipped at 0 and renormalised (the clip only removes rounding).
# ------------------------------------------------------------------

.scf_stop <- function(...) stop(..., call. = FALSE)

# V as a k x k symmetric matrix
.scf_V <- function(V, k) {
  if (is.null(dim(V))) {
    if (length(V) != k) .scf_stop("V: a vector V must have length ", k, " (one weight per predictor)")
    V <- diag(as.numeric(V), nrow = k)
  }
  if (!all(dim(V) == c(k, k))) .scf_stop("V must be ", k, " x ", k)
  if (anyNA(V)) .scf_stop("V has missing values")
  (V + t(V)) / 2
}

.scf_check_X <- function(X1, X0) {
  X0 <- as.matrix(X0)
  X1 <- as.numeric(X1)
  if (nrow(X0) != length(X1)) .scf_stop("X0 must have one row per element of X1 (donors in columns)")
  if (anyNA(X0) || anyNA(X1)) .scf_stop("X1 and X0 must not contain missing values")
  list(X1 = X1, X0 = X0)
}

# Moore-Penrose pseudo-inverse (for singular faces in the active-set solver)
.scf_pinv <- function(A) {
  sv <- svd(A)
  keep <- sv$d > 1e-12 * max(sv$d)
  sv$v[, keep, drop = FALSE] %*% (t(sv$u[, keep, drop = FALSE]) / sv$d[keep])
}

# Exact primal active-set solver for min f(w) = (1/2) w'Dw - d'w on the simplex
# (D positive semi-definite). Starts at the best vertex. On the current face (the
# donors in S) it takes the Newton step p (sum(p) = 0) to the face minimum,
# shortened by a ratio test when a weight would turn negative (that donor leaves
# S). On a singular face it first moves along a zero-curvature descent direction
# (D_SS p = 0, sum(p) = 0, g'p < 0), which must hit a bound, and otherwise uses
# the pseudo-inverse step. At a face minimum it checks the KKT multipliers of the
# donors off the face, (Dw - d)_j - mu >= 0, and lets the most violated one in.
# Returns NULL if it has not settled after maxit steps.
.scf_active_set <- function(D, d, tol = 1e-11, maxit = 2000L) {
  J <- length(d)
  j0 <- which.min(0.5 * diag(D) - d)
  w <- numeric(J); w[j0] <- 1
  S <- j0
  for (it in seq_len(maxit)) {
    m <- length(S)
    DS <- D[S, S, drop = FALSE]
    gS <- drop(DS %*% w[S]) - d[S]                       # gradient on the face (w = 0 off S)
    M <- rbind(cbind(DS, -1), c(rep(1, m), 0))
    p <- NULL
    if (rcond(M) > 1e-12) {
      p <- solve(M, c(-gS, 0))[seq_len(m)]
    } else {
      sv <- svd(rbind(DS, rep(1, m)), nu = 0)
      null <- sv$v[, sv$d <= 1e-10 * max(sv$d), drop = FALSE]
      if (m > length(sv$d)) null <- cbind(null, sv$v[, -seq_along(sv$d), drop = FALSE])
      if (ncol(null) > 0) {
        q <- -drop(null %*% crossprod(null, gS))
        if (sum(q * q) > 1e-28 && sum(gS * q) < -1e-14 * max(1, sqrt(sum(gS * gS)))) {
          neg <- which(q < 0)
          r <- -w[S][neg] / q[neg]
          k <- neg[which.min(r)]
          w[S] <- w[S] + min(r) * q
          w[S[k]] <- 0
          S <- S[-k]
          next
        }
      }
      p <- drop(.scf_pinv(M) %*% c(-gS, 0))[seq_len(m)]
    }
    neg <- which(p < -1e-15)
    if (length(neg)) {
      r <- -w[S][neg] / p[neg]
      if (min(r) < 1) {                                  # blocked: step to the bound, drop that donor
        k <- neg[which.min(r)]
        w[S] <- w[S] + min(r) * p
        w[S[k]] <- 0
        S <- S[-k]
        next
      }
    }
    w[S] <- w[S] + p
    g <- drop(D %*% w) - d
    mu <- mean(g[S])
    out <- setdiff(seq_len(J), S)
    if (!length(out)) return(w)
    lam <- g[out] - mu
    if (min(lam) >= -tol * max(1, abs(mu))) return(w)
    S <- c(S, out[which.min(lam)])
  }
  NULL
}

# min (1/2) w'Dw - d'w  s.t.  sum(w) = 1, w >= 0: quadprog (escalating ridge when
# D is singular) and the exact active-set solver; the answer with the lower
# objective wins, the active-set one on a tie (it satisfies the KKT conditions
# exactly; quadprog's ridge answer can be off by 1e-4 when D is near-singular).
.scf_qp_simplex <- function(D, d) {
  J <- ncol(D)
  if (J == 1L) return(1)
  s <- mean(diag(D))
  if (!is.finite(s) || s <= 0) s <- 1
  D <- (D + t(D)) / (2 * s)
  d <- d / s
  fobj <- function(w) 0.5 * sum(w * drop(D %*% w)) - sum(d * w)
  A <- cbind(rep(1, J), diag(J))
  b <- c(1, rep(0, J))
  w_qp <- NULL
  for (eps in c(0, 1e-12, 1e-10, 1e-8, 1e-6)) {
    sol <- tryCatch(quadprog::solve.QP(D + diag(eps, J), d, A, b, meq = 1), error = function(e) NULL)
    if (!is.null(sol)) {
      w_qp <- pmax(sol$solution, 0)
      w_qp <- w_qp / sum(w_qp)
      break
    }
  }
  w_as <- .scf_active_set(D, d)
  if (!is.null(w_as)) w_as <- pmax(w_as, 0) / sum(pmax(w_as, 0))
  if (is.null(w_qp) && is.null(w_as)) .scf_stop("could not solve the simplex QP (check X0 and V for non-finite values)")
  if (is.null(w_as)) return(w_qp)
  if (is.null(w_qp)) return(w_as)
  if (fobj(w_as) <= fobj(w_qp) + 1e-12 * max(1, abs(fobj(w_qp)))) w_as else w_qp
}

#' Synthetic-control weights: argmin over the simplex of (X1 - X0 w)' V (X1 - X0 w).
#' V is the chapter's predictor weighting matrix Omega (notation 2.9(c)); V = diag(k)
#' on k pre-period outcomes is SC on outcomes alone. Returns a numeric vector of
#' length ncol(X0), named by colnames(X0) when present.
sc_weights <- function(X1, X0, V) {
  z <- .scf_check_X(X1, X0)
  V <- .scf_V(V, length(z$X1))
  D <- 2 * crossprod(z$X0, V %*% z$X0)
  d <- 2 * drop(crossprod(z$X0, V %*% z$X1))
  w <- .scf_qp_simplex(D, d)
  names(w) <- colnames(z$X0)
  w
}

#' Penalized synthetic control (Abadie and L'Hour 2021), weights on the simplex.
#' pen >= 0 multiplies sum_j w_j ||X1 - X0[, j]||_V^2. n_nonzero counts weights above tol.
penalized_sc <- function(X1, X0, pen, V = diag(length(X1)), tol = 1e-6) {
  z <- .scf_check_X(X1, X0)
  if (!is.numeric(pen) || length(pen) != 1L || is.na(pen) || pen < 0) .scf_stop("pen must be a number >= 0")
  V <- .scf_V(V, length(z$X1))
  diffs <- z$X1 - z$X0                                   # k x J: X1 - X0[, j]
  dist2 <- colSums(diffs * (V %*% diffs))                # ||X1 - X0[, j]||_V^2
  D <- 2 * crossprod(z$X0, V %*% z$X0)
  d <- 2 * drop(crossprod(z$X0, V %*% z$X1)) - pen * dist2
  w <- .scf_qp_simplex(D, d)
  names(w) <- colnames(z$X0)
  list(w = w, n_nonzero = sum(w > tol))
}

#' Conformal p-value of Chernozhukov, Wuthrich and Zhu (2021) for H0: no effect in
#' periods T0+1, ..., T, given residuals y1 - y0hat from a fit under H0 (see header).
conformal_p <- function(y1, y0hat, T0, scheme = c("moving_block", "iid"), q = 1, n_perm = 999) {
  scheme <- match.arg(scheme)
  if (length(y1) != length(y0hat)) .scf_stop("y1 and y0hat must have the same length")
  u <- as.numeric(y1) - as.numeric(y0hat)
  Tn <- length(u)
  if (anyNA(u)) .scf_stop("y1 and y0hat must not contain missing values")
  if (!isTRUE(length(T0) == 1L && is.numeric(T0) && T0 >= 1 && T0 < Tn && T0 == round(T0)))
    .scf_stop("T0 must be an integer in 1, ..., length(y1) - 1")
  a <- abs(u)^q
  post <- (T0 + 1):Tn
  S0 <- sum(a[post])                                     # monotone in S(u): same ranking
  if (scheme == "moving_block") {
    S <- vapply(0:(Tn - 1), function(j) sum(a[((post - 1 + j) %% Tn) + 1]), numeric(1))
  } else {
    # the statistic only depends on which residuals land in the post positions
    n_post <- length(post)
    S <- c(S0, vapply(seq_len(n_perm), function(b) sum(a[sample.int(Tn, n_post)]), numeric(1)))
  }
  mean(S >= S0 - 1e-12 * max(1, abs(S0)))
}

#' SC counterfactual from pre-period outcomes: y1 is the length-T outcome path of the
#' (placebo) treated unit, Y0 the T x J donor outcomes. Weights fit periods 1..T0
#' with V = I (or V); returns the length-T path Y0 %*% w with attribute "weights".
sc_fit <- function(y1, Y0, T0, V = NULL) {
  Y0 <- as.matrix(Y0)
  if (is.null(V)) V <- diag(T0)
  w <- sc_weights(y1[seq_len(T0)], Y0[seq_len(T0), , drop = FALSE], V)
  out <- drop(Y0 %*% w)
  attr(out, "weights") <- w
  out
}

#' Placebo (in-space) permutation inference with the post/pre MSPE ratio
#' (Abadie, Diamond and Hainmueller 2010; the chapter's R_i and p_perm, notation 2.9(c)).
#' Y: N x T outcomes, units in rows; treated: row index (or rowname) of the treated unit;
#' T0: number of pre-periods; fit_fun(y1, Y0, T0) -> length-T counterfactual for y1
#' from donor matrix Y0 (T x J, donors in columns), e.g. sc_fit.
#' placebo_donors: "all_others" (default; ADH 2010, who shift the treated unit into
#' the donor pool of every placebo; also tidysynth 0.2.1): when unit i plays the
#' treated role, its donors are all other units, the treated unit included.
#' "exclude_treated" (alternative): a placebo's donors are all units except i and
#' the treated unit, which keeps the treated unit's post-period outcomes (and so
#' the effect) out of the placebo fits. The treated unit's own fit is the same
#' under both.
#' Returns ratio = R_treated, p = mean(R_i >= R_treated) over all N units (minimum
#' 1/N), rank = 1 + #{i: R_i > R_treated} (1 = largest), all_ratios (named), n_units.
#' Stops (with a message, also on NA input) when T0 is not an integer in 1, ..., T - 1,
#' when fit_fun returns missing or non-finite values, or when some unit's pre-period
#' MSPE is zero (below 1e-24 times its mean squared pre-period outcome), since its
#' ratio would then be Inf or NaN.
mspe_ratio_perm <- function(Y, treated, T0, fit_fun,
                            placebo_donors = c("all_others", "exclude_treated")) {
  placebo_donors <- match.arg(placebo_donors)
  Y <- as.matrix(Y)
  N <- nrow(Y); Tn <- ncol(Y)
  if (is.character(treated)) treated <- match(treated, rownames(Y))
  if (!isTRUE(length(treated) == 1L && !is.na(treated) && treated >= 1 && treated <= N)) .scf_stop("treated must name one row of Y")
  if (!isTRUE(length(T0) == 1L && is.numeric(T0) && T0 >= 1 && T0 < Tn && T0 == round(T0)))
    .scf_stop("T0 must be an integer in 1, ..., ncol(Y) - 1")
  pre <- seq_len(T0); post <- (T0 + 1):Tn
  ratios <- vapply(seq_len(N), function(i) {
    donors <- setdiff(seq_len(N), i)
    if (placebo_donors == "exclude_treated" && i != treated) donors <- setdiff(donors, treated)
    yhat <- fit_fun(Y[i, ], t(Y[donors, , drop = FALSE]), T0)
    g <- Y[i, ] - as.numeric(yhat)
    if (length(g) != Tn || !all(is.finite(g)))
      .scf_stop("mspe_ratio_perm: fit_fun returned missing, non-finite or too few values for unit ", i,
                if (!is.null(rownames(Y))) paste0(" (", rownames(Y)[i], ")"), " (or Y has missing values)")
    pre_mspe <- mean(g[pre]^2)
    if (!isTRUE(pre_mspe > 1e-24 * max(mean(Y[i, pre]^2), .Machine$double.xmin)))
      .scf_stop("mspe_ratio_perm: the pre-period MSPE of unit ", i,
                if (!is.null(rownames(Y))) paste0(" (", rownames(Y)[i], ")"),
                " is zero (an exact pre-period fit), so its post/pre ratio is undefined")
    mean(g[post]^2) / pre_mspe
  }, numeric(1))
  names(ratios) <- if (is.null(rownames(Y))) as.character(seq_len(N)) else rownames(Y)
  r1 <- ratios[[treated]]
  list(ratio = r1, p = mean(ratios >= r1), rank = 1L + sum(ratios > r1),
       all_ratios = ratios, n_units = N)
}

#' Is x1 in the convex hull of the columns of X0? LP feasibility (lpSolve):
#'   min sum(s_plus + s_minus)  s.t.  X0 w + s_plus - s_minus = x1, sum(w) = 1, w, s >= 0.
#' x1 is in the hull iff the minimum (the L1 distance from x1 to the hull, i.e. to
#' {X0 w : w in the simplex}, returned as attribute "l1_gap") is 0, up to
#' tol * max(1, max|x1|, max|X0|).
in_hull_lp <- function(x1, X0, tol = 1e-9) {
  z <- .scf_check_X(x1, X0)
  k <- length(z$X1); J <- ncol(z$X0)
  A <- rbind(cbind(z$X0, diag(k), -diag(k)),
             c(rep(1, J), rep(0, 2 * k)))
  sol <- lpSolve::lp("min", c(rep(0, J), rep(1, 2 * k)), A, rep("=", k + 1), c(z$X1, 1))
  if (sol$status != 0) .scf_stop("in_hull_lp: lpSolve status ", sol$status)
  gap <- sol$objval
  out <- gap <= tol * max(1, abs(z$X1), abs(z$X0))
  attr(out, "l1_gap") <- gap
  out
}
