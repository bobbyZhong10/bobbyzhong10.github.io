# ------------------------------------------------------------------
# qr_fns.R -- hand-written quantile estimators and the Marigold data generators for chapter 13
#             (quantile regression)
#
# Chapter card interface:
#   check_loss(u, tau)                          -> numeric      rho_tau(u) = (tau - 1{u < 0}) u
#   gen_inverse(support, cdf, tau)              -> numeric      F^{-1}(tau) = inf{y : F(y) >= tau}
#   qte_marginal(y, d, taus, B, level)          -> data.frame(tau, est, lo, hi, se)
#   rif_quantile(y, tau, bw = NULL)             -> numeric      RIF(Y_i; tau), one value per observation
#   ivqr_grid(y, d, z, X, tau, grid, level, se) -> list(est, ci, ...)   Chernozhukov-Hansen inverse QR
#   dist_reg(y, X, thresholds, link = "logit")  -> matrix       gamma(y), one row per threshold
#   canay_qr(y, X, id, tau)                     -> list(beta, beta_mean, alpha_hat)
#   cic_qtt(y00, y01, y10, taus, y11)           -> data.frame(tau, qtt, q11, cf)
#   cic_lqte(df, y, d, group, time, tau, B, level) -> list(est, ci, se, num, first_stage, wald, ...)
# Added here (the card is silent; each is used by an estimation script or a test):
#   emp_quantile(y, tau)                        F_n^{-1}(tau), the generalized inverse of the empirical cdf
#   qq_map(yeval, y_d00, y_d01)                 the changes-in-changes map Q_d(y) in exact integer arithmetic
#   qdid(y00, y01, y10, y11, taus)              quantile difference-in-differences (Athey and Imbens 2006)
#   lqte_iv(y, d, z, taus)                      complier QTE from the Wald-type complier cdfs (no covariates)
#   dr_cdf(gamma, Xnew, link)                   distribution-regression cdf averaged over the rows of Xnew
#   fe_dummy_qr(y, X, id, tau)                  QR with one dummy per unit (sparse Frisch-Newton)
#   save_extra(name, value, file)               add one block's results to res_extra.rds
#   gen_Q1 ... gen_Q5f, mix_cdf_Q2, mix_pdf_Q2, h_Q5   data generators and closed-form truth helpers
#
# Notation (notation_v3.md 1.11, 2.1, 2.2, 3.2): tau is the quantile level only (never an effect);
# individual effects are tau_i; F_Y, F_Y^{-1}, Quant_tau(Y | X = x); beta(tau); rho_tau(u); the CIC
# cells Y_st have the group s (S_i) first and the period t (T_i) second, as in chapter 20; Q_d is the
# CIC map, Q_{d=0} for the untreated; the panel length is \mathcal{T} (the code writes Tn).
# Sample quantiles are generalized inverses of empirical cdfs, F_n^{-1}(tau) = emp_quantile(y, tau) =
# gen_inverse(sort(y), (1:n) / n, tau) exactly. This is quantile(type = 1) in exact arithmetic, but R's
# quantile(type = 1) can return the next order statistic when n * tau rounds just above an integer
# (n = 800, tau = 0.55; 145 of the 2500 levels k / n for n = 2500), so the two are not always identical.
# rif_quantile is the exception: it follows rifreg 1.1.0 (type 7) so that the two agree exactly.
# ------------------------------------------------------------------

## =================================================================
## Check function, generalized inverse, sample quantiles
## =================================================================

#' rho_tau(u) = (tau - 1{u < 0}) u. u is the argument of the check function (a residual
#' Y_i - X_i'beta in a QR objective); tau is the quantile level tau in (0, 1).
check_loss <- function(u, tau) (tau - (u < 0)) * u

#' F^{-1}(tau) = inf{y in support : F(y) >= tau}, the generalized inverse (left-continuous quantile
#' function). support: the points y (strictly increasing); cdf: F_Y(y) at those points; tau: levels in
#' (0, 1]. Properties used in the chapter (Wacker 2023, arXiv 2306.06989, Lemma 1 (g)(3) and (i)):
#' F^{-1}(F(y)) <= y and, F right-continuous, F(F^{-1}(tau)) >= tau. The definition needs no
#' monotone cdf, so the same function inverts estimated cdfs that are not monotone in finite samples
#' (cic_lqte, lqte_iv). Returns NA where F never reaches tau.
gen_inverse <- function(support, cdf, tau) {
  support <- as.numeric(support); cdf <- as.numeric(cdf); tau <- as.numeric(tau)
  if (!length(support) || length(support) != length(cdf)) stop("gen_inverse: support and cdf must have the same, positive length")
  if (anyNA(support) || anyNA(cdf)) stop("gen_inverse: missing values in support or cdf")
  if (is.unsorted(support, strictly = TRUE)) stop("gen_inverse: support must be strictly increasing")
  if (!length(tau) || anyNA(tau) || any(tau <= 0) || any(tau > 1)) stop("gen_inverse: tau must lie in (0, 1]")
  vapply(tau, function(t) {
    j <- which(cdf >= t)[1L]
    if (is.na(j)) NA_real_ else support[j]
  }, numeric(1))
}

#' F_n^{-1}(tau): the generalized inverse of the empirical cdf of y (cdf k / n at the k-th order statistic),
#' equal to gen_inverse(sort(y), (1:n) / n, tau) exactly. The index is the smallest j with j / n >= tau,
#' computed with the same floating-point comparison. It agrees with quantile(y, tau, type = 1) except where
#' n * tau rounds just above an integer, e.g. n = 800 at tau = 0.55 (n * tau = 440.00000000000006), where
#' quantile(type = 1) returns the 441st order statistic and this function the 440th, as the definition
#' inf{y : F_n(y) >= tau} requires.
emp_quantile <- function(y, tau) {
  s <- sort(as.numeric(y)); n <- length(s)
  if (!n) stop("emp_quantile: y is empty")
  if (!length(tau) || anyNA(tau) || any(tau <= 0) || any(tau > 1)) stop("emp_quantile: tau must lie in (0, 1]")
  j <- ceiling(n * tau)
  j <- j - ((j - 1) / n >= tau)
  j <- j + (j / n < tau)
  s[pmin(pmax(j, 1), n)]
}

## =================================================================
## Quantile treatment effects
## =================================================================

#' Marginal QTE: QTE(tau) = F^{-1}_{Y | D = 1}(tau) - F^{-1}_{Y | D = 0}(tau) with sample quantiles
#' (emp_quantile: the generalized inverse of each arm's empirical cdf).
#' Under random assignment of D this estimates F^{-1}_{Y(1)}(tau) - F^{-1}_{Y(0)}(tau) (Firpo 2007), a
#' difference of two marginal distributions, not a quantile of the individual effects tau_i.
#' y: Y_i; d: D_i in {0, 1}; taus: levels tau. Standard errors: nonparametric bootstrap that resamples
#' each arm separately (B draws; set.seed before the call); interval est -/+ z se.
qte_marginal <- function(y, d, taus, B = 999, level = 0.95) {
  if (!all(d %in% c(0, 1))) stop("qte_marginal: d must take values in {0, 1}")
  y1 <- y[d == 1]; y0 <- y[d == 0]
  est <- emp_quantile(y1, taus) - emp_quantile(y0, taus)
  se <- rep(NA_real_, length(taus))
  if (B > 0) {
    bs <- vapply(seq_len(B), function(b) {
      emp_quantile(y1[sample.int(length(y1), replace = TRUE)], taus) -
        emp_quantile(y0[sample.int(length(y0), replace = TRUE)], taus)
    }, numeric(length(taus)))
    se <- apply(matrix(bs, nrow = length(taus)), 1, sd)
  }
  z <- qnorm(1 - (1 - level) / 2)
  data.frame(tau = taus, est = est, lo = est - z * se, hi = est + z * se, se = se)
}

#' Complier QTE without covariates (Abadie, Angrist and Imbens 2002; Frolich and Melly 2013 with X empty):
#' with a binary instrument Z_i (random), monotone D_i(z), the complier cdfs are
#'   F_{Y(1) | C}(y) = [E(1{Y <= y} D | Z = 1) - E(1{Y <= y} D | Z = 0)] / [E(D | Z = 1) - E(D | Z = 0)],
#'   F_{Y(0) | C}(y) = [E(1{Y <= y}(1 - D) | Z = 0) - E(1{Y <= y}(1 - D) | Z = 1)] / [E(D | Z = 1) - E(D | Z = 0)],
#' and LQTE(tau) = F^{-1}_{Y(1) | C}(tau) - F^{-1}_{Y(0) | C}(tau) (generalized inverse over the sample points;
#' the estimated cdfs need not be monotone). y: Y_i; d: D_i; z: Z_i; taus: levels. Point estimates only.
lqte_iv <- function(y, d, z, taus) {
  if (!all(d %in% c(0, 1)) || !all(z %in% c(0, 1))) stop("lqte_iv: d and z must take values in {0, 1}")
  n1 <- sum(z == 1); n0 <- sum(z == 0)
  fs <- mean(d[z == 1]) - mean(d[z == 0])
  if (abs(fs) < 1e-8) stop("lqte_iv: first stage is (near) zero")
  supp <- sort(unique(y))
  cnt <- function(sel) if (any(sel)) findInterval(supp, sort(y[sel])) else numeric(length(supp))
  F1 <- (cnt(z == 1 & d == 1) / n1 - cnt(z == 0 & d == 1) / n0) / fs
  F0 <- (cnt(z == 0 & d == 0) / n0 - cnt(z == 1 & d == 0) / n1) / fs
  q1 <- gen_inverse(supp, F1, taus); q0 <- gen_inverse(supp, F0, taus)
  data.frame(tau = taus, lqte = q1 - q0, q1c = q1, q0c = q0, first_stage = fs)
}

## =================================================================
## Unconditional quantiles: RIF and distribution regression
## =================================================================

#' Recentered influence function of the tau-quantile (Firpo, Fortin and Lemieux 2009, eq. 15):
#'   RIF(Y; tau) = F_Y^{-1}(tau) + (tau - 1{Y <= F_Y^{-1}(tau)}) / f_Y(F_Y^{-1}(tau)).
#' y: Y_i; tau: level; bw: the bandwidth h_n of the Gaussian kernel K for f_Y (NULL: bw.nrd0, Silverman's
#' rule; any value stats::density accepts, e.g. "SJ"). Sample quantile type 7 and density() evaluated at
#' the quantile, as in rifreg 1.1.0 (rifreg:::get_rif_quantile), so OLS of this RIF on X reproduces
#' rifreg() exactly. Returns the n RIF values with attributes quantile, density and bw.
rif_quantile <- function(y, tau, bw = NULL) {
  if (length(tau) != 1L || tau <= 0 || tau >= 1) stop("rif_quantile: tau must be one level in (0, 1)")
  q <- quantile(y, tau, type = 7, names = FALSE)
  dens <- if (is.null(bw)) stats::density(y, from = q, to = q) else stats::density(y, bw = bw, from = q, to = q)
  f <- dens$y[1]
  structure(q + (tau - as.numeric(y <= q)) / f, quantile = q, density = f, bw = dens$bw)
}

#' Distribution regression (Foresi and Peracchi 1995; Chernozhukov, Fernandez-Val and Melly 2013):
#' F_{Y | X}(y | x) = Lambda(x'gamma(y)), one binary-response GLM of 1{Y <= y} on (1, X) per threshold y.
#' y: Y_i; X: covariate matrix without intercept (named columns); thresholds: the values y; link: the
#' binomial link of Lambda ("logit", "probit", ...). glm.fit with glm's default control, so each row equals
#' coef(glm(...)). glm.fit's warning "fitted probabilities numerically 0 or 1 occurred" (expected at
#' thresholds in the tails, where some merchants' fitted probabilities are within 1e-8 of 0 or 1) is
#' muffled and counted in the attribute n_extreme; non-convergence is always reported.
#' Returns a matrix: one row per threshold, columns "(Intercept)" and colnames(X).
dist_reg <- function(y, X, thresholds, link = "logit") {
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  M <- cbind("(Intercept)" = 1, X)
  fam <- stats::binomial(link)
  n_extreme <- 0L
  out <- t(vapply(thresholds, function(cval) {
    g <- withCallingHandlers(stats::glm.fit(M, as.numeric(y <= cval), family = fam), warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w), fixed = TRUE)) {
        n_extreme <<- n_extreme + 1L
        invokeRestart("muffleWarning")
      }
    })
    if (!g$converged) warning("dist_reg: glm.fit did not converge at threshold ", format(cval))
    g$coefficients
  }, numeric(ncol(M))))
  dimnames(out) <- list(format(thresholds, digits = 8), colnames(M))
  attr(out, "n_extreme") <- n_extreme
  out
}

#' Distribution-regression cdf at each threshold, averaged over the rows of Xnew:
#' F(y_j) = n^{-1} sum_i Lambda((1, Xnew_i)'gamma(y_j)). gamma: the matrix from dist_reg; Xnew: covariate
#' matrix (same columns as in the fit, no intercept); link: the same link. A counterfactual covariate
#' distribution enters through Xnew (or through a weighted average of two such calls).
dr_cdf <- function(gamma, Xnew, link = "logit") {
  M <- cbind(1, as.matrix(Xnew))
  inv <- stats::binomial(link)$linkinv
  colMeans(inv(M %*% t(gamma)))
}

## =================================================================
## IV quantile regression (Chernozhukov and Hansen 2005, 2006, 2008)
## =================================================================

## internal: evaluate expr, muffling only quantreg's "Solution may be nonunique" warning (the simplex
## reports it when n tau is an integer; any minimiser then is a solution and other warnings pass through)
.muffle_nonunique <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("nonunique", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}

#' Inverse quantile regression on a grid. Structural model Y_i = D_i beta_D(U_i) + X_i'beta_X(U_i), U_i | Z_i
#' uniform, rank similarity; then P(Y <= D beta_D(tau) + X'beta_X(tau) | Z, X) = tau. For each candidate b in
#' `grid` (values of beta_D(tau)) the QR of Y - D b on (1, X, Z) at level tau gives a coefficient on Z;
#' beta_D(tau) is estimated by the grid point where that coefficient's t statistic is closest to zero, and
#' the confidence set collects the grid points where |t| <= z_{1 - (1 - level)/2} (the dual inference of
#' Chernozhukov and Hansen 2008, robust to a weak instrument).
#' y: Y_i; d: D_i (endogenous); z: Z_i (instrument); X: exogenous covariates without intercept (matrix or
#' NULL); tau: level; grid: candidate values of beta_D(tau); se: the summary.rq option for the Z coefficient
#' ("nid": Hendricks-Koenker sandwich with the Hall-Sheather bandwidth h_n).
#' Returns est, ci = c(lo, hi) (range of the accepted points), connected (accepted points contiguous),
#' at_boundary (an end of the grid accepted: widen the grid), and per grid point coef_z, se_z, tstat.
ivqr_grid <- function(y, d, z, X, tau, grid, level = 0.95, se = "nid") {
  M <- cbind("(Intercept)" = 1, X, z = z)
  jz <- ncol(M)
  crit <- qnorm(1 - (1 - level) / 2)
  res <- vapply(grid, function(b) {
    yb <- y - d * b
    fit <- .muffle_nonunique(quantreg::rq(yb ~ M - 1, tau = tau))
    cf <- .muffle_nonunique(summary(fit, se = se))$coefficients
    cf[jz, 1:2]
  }, numeric(2))
  coef_z <- res[1, ]; se_z <- res[2, ]; tstat <- coef_z / se_z
  acc <- which(abs(tstat) <= crit)
  if (!length(acc)) warning("ivqr_grid: no grid point accepted; the confidence set is empty on this grid")
  list(est = grid[which.min(abs(tstat))],
       ci = if (length(acc)) range(grid[acc]) else c(NA_real_, NA_real_),
       connected = length(acc) > 0 && all(diff(acc) == 1),
       at_boundary = length(acc) > 0 && (1L %in% acc || length(grid) %in% acc),
       grid = grid, coef_z = coef_z, se_z = se_z, tstat = tstat, tau = tau, level = level)
}

## =================================================================
## Panel QR with fixed effects
## =================================================================

#' Canay (2011) two-step estimator for Quant_tau(Y_it | X_it, alpha_i) = alpha_i + X_it'beta(tau) (the unit
#' effect a pure location shift). Step 1: within (FE) estimator beta_mean of the conditional-mean model and
#' alpha_hat_i = Tn^{-1} sum_t (Y_it - X_it'beta_mean). Step 2: QR of Y_it - alpha_hat_i on (1, X_it).
#' y: Y_it; X: X_it (matrix, no intercept, named columns); id: unit i; tau: level.
#' Canay's consistency result needs n / Tn^s -> 0 for some s > 1 (Theorem 4.1) and his rate condition is
#' disputed (Besstremyannaya and Golovan 2019): with Tn in {5, 20} the chapter's Monte Carlo is a finite-
#' sample illustration only, not a check of an asymptotic theorem.
#' Returns beta (named: "(Intercept)", colnames(X)), beta_mean, alpha_hat (one value per observation).
canay_qr <- function(y, X, id, tau) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  Xd <- X - apply(X, 2, function(col) ave(col, id))
  yd <- y - ave(y, id)
  bmean <- qr.coef(qr(Xd), yd)
  ah <- ave(drop(y - X %*% bmean), id)
  fit <- .muffle_nonunique(quantreg::rq.fit(cbind(1, X), y - ah, tau = tau, method = "br"))
  beta <- fit$coefficients
  names(beta) <- c("(Intercept)", colnames(X))
  list(beta = beta, beta_mean = setNames(as.numeric(bmean), colnames(X)), alpha_hat = ah)
}

#' QR with one dummy per unit (Koenker 2004 without the penalty): minimise sum rho_tau(Y_it - X_it'b - a_i)
#' over b and every a_i. y: Y_it; X: X_it (matrix, no intercept); id: unit i; tau: level. The design
#' [X, unit dummies] is built as a SparseM matrix.csr and solved by quantreg's sparse Frisch-Newton
#' interior-point method (rq.fit.sfn). With Tn fixed the incidental-parameter problem makes the slope
#' inconsistent (Wooldridge 2010, section 12.10.3). Returns beta (named slopes) and the unit effects.
fe_dummy_qr <- function(y, X, id, tau) {
  X <- as.matrix(X); n <- nrow(X); k <- ncol(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(k))
  f <- factor(id); nu <- nlevels(f)
  A <- methods::new(methods::getClass("matrix.csr", where = asNamespace("SparseM")), ra = as.numeric(t(cbind(X, 1))),
           ja = as.integer(t(cbind(matrix(seq_len(k), n, k, byrow = TRUE), k + as.integer(f)))),
           ia = as.integer(seq(1L, by = k + 1L, length.out = n + 1L)), dimension = as.integer(c(n, k + nu)))
  fit <- quantreg::rq.fit.sfn(A, y, tau = tau)
  if (!is.null(fit$ierr) && fit$ierr != 0) warning("fe_dummy_qr: rq.fit.sfn returned ierr = ", fit$ierr)
  list(beta = setNames(fit$coef[seq_len(k)], colnames(X)), alpha = fit$coef[k + seq_len(nu)])
}

## =================================================================
## Changes-in-changes (Athey and Imbens 2006; de Chaisemartin and D'Haultfoeuille 2018)
## =================================================================

#' Q_d(y) = F^{-1}_{Y_d01}(F_{Y_d00}(y)) with empirical cdfs and the exact generalized inverse.
#' yeval: points y; y_d00, y_d01: the control group's outcomes (with D = d) in periods 0 and 1.
#' With k0 = #{y_d00 <= y} = n0 F_{Y_d00}(y), the generalized inverse at k0 / n0 is the j-th order statistic
#' of y_d01 with j the smallest integer such that j n0 >= n1 k0; this is computed in integer arithmetic
#' (exact in double precision below 2^53), so equal or proportional cell sizes cannot produce the
#' floating-point off-by-one of quantile(type = 1) at probs = k0 / n0. y below min(y_d00) (k0 = 0) maps
#' to min(y_d01), the convention of chapter 20's Wald-CIC.
qq_map <- function(yeval, y_d00, y_d01) {
  s00 <- sort(y_d00); s01 <- sort(y_d01)
  n0 <- length(s00); n1 <- length(s01)
  if (!n0 || !n1) stop("qq_map: empty control cell")
  k0 <- findInterval(yeval, s00)
  target <- as.double(n1) * k0
  j <- ceiling(target / n0)
  j <- j + (j * n0 < target) - ((j - 1) * n0 >= target)
  s01[pmax(j, 1)]
}

#' CIC quantile treatment effect on the treated (Athey and Imbens 2006, Theorem 3.1 and the quantile effect
#' after eq. 3.16): QTT(tau) = F^{-1}_{Y11}(tau) - F^{-1}_{Y01}(F_{Y00}(F^{-1}_{Y10}(tau))).
#' y00, y01: control group (S_i = 0) in periods T_i = 0, 1; y10: treated group (S_i = 1) in period 0;
#' taus: levels; y11: treated group in period 1. (The card's signature has no y11; it is needed for the
#' effect and is added as the fifth argument.) The counterfactual quantile cf = Q_{d=0}(F^{-1}_{Y10}(tau))
#' equals the tau-quantile of Q_{d=0}(Y10) because Q_{d=0} is non-decreasing.
cic_qtt <- function(y00, y01, y10, taus, y11) {
  if (missing(y11)) stop("cic_qtt: y11 (treated group, period 1) is needed for the QTT")
  cf <- qq_map(emp_quantile(y10, taus), y00, y01)
  q11 <- emp_quantile(y11, taus)
  data.frame(tau = taus, qtt = q11 - cf, q11 = q11, cf = cf)
}

#' Quantile difference-in-differences (Athey and Imbens 2006, Theorem 3.2):
#' QDID(tau) = F^{-1}_{Y11}(tau) - F^{-1}_{Y10}(tau) - [F^{-1}_{Y01}(tau) - F^{-1}_{Y00}(tau)].
#' Arguments as in cic_qtt (cells Y_st, group first, period second).
qdid <- function(y00, y01, y10, y11, taus) {
  data.frame(tau = taus, qdid = emp_quantile(y11, taus) - emp_quantile(y10, taus) -
               (emp_quantile(y01, taus) - emp_quantile(y00, taus)))
}

## internal: validate a two-group, two-period data set (as chapter 20's .fz_prepare)
.cic_prepare <- function(df, y, d, group, time) {
  df <- as.data.frame(df)
  vars <- c(y, d, group, time)
  if (!is.character(vars) || length(vars) != 4) stop("y, d, group, time must be single column names (strings)")
  miss <- setdiff(vars, names(df))
  if (length(miss)) stop("column(s) not found: ", paste(miss, collapse = ", "))
  Y <- as.numeric(df[[y]]); D <- as.numeric(df[[d]]); G <- as.numeric(df[[group]]); Tt <- as.numeric(df[[time]])
  if (anyNA(Y) || anyNA(D) || anyNA(G) || anyNA(Tt)) stop("missing values in y, d, group or time")
  if (!all(D %in% c(0, 1)) || !all(G %in% c(0, 1)) || !all(Tt %in% c(0, 1))) stop("d, group and time must take values in {0, 1}")
  for (g in 0:1) for (t in 0:1) if (!any(G == g & Tt == t)) stop(sprintf("empty cell: group %d, time %d", g, t))
  list(Y = Y, D = D, G = G, T = Tt)
}

## internal: fuzzy CIC on prepared data x (list Y, D, G, T) at one level tau.
##   F_CIC,d(y) = [P(D11 = d) F_{Y_d11}(y) - P(D10 = d) F_{Q_d(Y_d10)}(y)] / [P(D11 = d) - P(D10 = d)]
## is the period-1 cdf of Y(d) for the switchers (dCDH 2018, Theorem 2.3: Assumptions 1, 2, 6, 7, 8); the
## switchers' LQTE is F^{-1}_{CIC,1}(tau) - F^{-1}_{CIC,0}(tau). num = E(Y11) - E(Q_{D10}(Y10)) is the
## Wald-CIC numerator (chapter 20's wald_cic: num / first_stage).
.cic_fuzzy_core <- function(x, tau) {
  in11 <- x$G == 1 & x$T == 1; in10 <- x$G == 1 & x$T == 0
  fs <- mean(x$D[in11]) - mean(x$D[in10])
  if (abs(fs) < 1e-8) stop("cic_lqte: treated-group first stage E(D11) - E(D10) is (near) zero")
  mapped <- list(numeric(0), numeric(0))                    # Q_d(Y) for the units of cell (1, 0) with D = d
  for (dd in 0:1) {
    k10 <- in10 & x$D == dd
    if (!any(k10)) next
    c00 <- x$D == dd & x$G == 0 & x$T == 0; c01 <- x$D == dd & x$G == 0 & x$T == 1
    if (!any(c00) || !any(c01)) stop(sprintf("cic_lqte: Q_%d undefined: no control-group unit with D = %d in both periods", dd, dd))
    mapped[[dd + 1]] <- qq_map(x$Y[k10], x$Y[c00], x$Y[c01])
  }
  num <- mean(x$Y[in11]) - (sum(mapped[[1]]) + sum(mapped[[2]])) / sum(in10)
  qcic <- function(dd) {
    a <- x$Y[in11 & x$D == dd]; b <- mapped[[dd + 1]]
    w11 <- mean(x$D[in11] == dd); w10 <- mean(x$D[in10] == dd)
    supp <- sort(unique(c(a, b)))
    Fa <- if (length(a)) findInterval(supp, sort(a)) / length(a) else 0
    Fb <- if (length(b)) findInterval(supp, sort(b)) / length(b) else 0
    gen_inverse(supp, (w11 * Fa - w10 * Fb) / (w11 - w10), tau)
  }
  list(est = qcic(1) - qcic(0), num = num, first_stage = fs)
}

#' Switchers' LQTE in a fuzzy two-group, two-period design (de Chaisemartin and D'Haultfoeuille 2018,
#' Theorem 2.3): LQTE_sw(tau) = F^{-1}_{CIC,1}(tau) - F^{-1}_{CIC,0}(tau), with the maps Q_d estimated in the
#' control group among units with D = d. df: data frame; y, d, group, time: column names of Y_i, D_i,
#' S_i, T_i (values 0 and 1, as in chapter 20's wald_cic); tau: one level. Inference: bootstrap that
#' resamples within each (group, period) cell (repeated cross sections; B draws, set.seed before the call),
#' interval est -/+ z se; B = 0 skips it (ci = NULL). Also returns the Wald-CIC numerator num =
#' E(Y11) - E(Q_{D10}(Y10)), first_stage = E(D11) - E(D10) and wald = num / first_stage (the switchers'
#' LATE, chapter 20's W_CIC), with its bootstrap se_wald and ci_wald from the same draws.
cic_lqte <- function(df, y, d, group, time, tau, B = 499, level = 0.95) {
  if (length(tau) != 1L) stop("cic_lqte: one tau at a time")
  x <- .cic_prepare(df, y, d, group, time)
  core <- .cic_fuzzy_core(x, tau)
  wald <- core$num / core$first_stage
  ci <- NULL; se <- NA_real_; boot <- NULL; ci_wald <- NULL; se_wald <- NA_real_
  if (B > 0) {
    cells <- split(seq_along(x$Y), paste(x$G, x$T))
    boot <- vapply(seq_len(B), function(b) {
      idx <- unlist(lapply(cells, function(ix) ix[sample.int(length(ix), replace = TRUE)]), use.names = FALSE)
      xb <- list(Y = x$Y[idx], D = x$D[idx], G = x$G[idx], T = x$T[idx])
      tryCatch({ r <- .cic_fuzzy_core(xb, tau); c(r$est, r$num / r$first_stage) }, error = function(e) c(NA_real_, NA_real_))
    }, numeric(2))
    if (anyNA(boot[1, ])) warning("cic_lqte: ", sum(is.na(boot[1, ])), " bootstrap draws failed and were dropped")
    z <- qnorm(1 - (1 - level) / 2)
    se <- sd(boot[1, ], na.rm = TRUE); se_wald <- sd(boot[2, ], na.rm = TRUE)
    ci <- c(core$est - z * se, core$est + z * se); ci_wald <- c(wald - z * se_wald, wald + z * se_wald)
  }
  list(est = core$est, ci = ci, se = se, num = core$num, first_stage = core$first_stage,
       wald = wald, se_wald = se_wald, ci_wald = ci_wald, tau = tau, B = B, boot = boot)
}

## =================================================================
## Results file shared by 03 to 06
## =================================================================

#' res_extra.rds = list(Q2, Q3, Q4, Q5, cross): each estimation script adds its block (01_dgp.R deletes the
#' file, so a full run never mixes old and new blocks).
save_extra <- function(name, value, file = "res_extra.rds") {
  x <- if (file.exists(file)) readRDS(file) else list()
  x[[name]] <- value
  saveRDS(x, file)
  invisible(file)
}

## =================================================================
## Marigold data generators (01_dgp.R calls each after its own set.seed; tests/test_dgp.R
## regenerates the saved data from the stored seeds and parameters)
## =================================================================

#' Q1 (AutoPrice experiment). P: n, n_treat, m, b, s, v = b^2 + s^2, k, logc. Monthly revenue in USD thousands:
#' log Y(0) = m + b X + s e, Y(1) = exp(logc) Y(0)^k (rank invariance: Y(1) an increasing function of Y(0)),
#' X ~ N(0, 1) the standardised log of last year's revenue, e ~ N(0, 1); D: exactly n_treat treated.
gen_Q1 <- function(P) {
  n <- P$n
  x <- rnorm(n); e <- rnorm(n)
  D <- sample(rep(0:1, c(n - P$n_treat, P$n_treat)))
  ly0 <- P$m + P$b * x + P$s * e
  y0 <- exp(ly0); y1 <- exp(P$logc + P$k * ly0)
  data.frame(merchant = seq_len(n), x = x, D = D, y = ifelse(D == 1, y1, y0), y0_true = y0, y1_true = y1)
}

#' Q1r: the same merchants, X, D and Y(0) as Q1; Y(1) has the same marginal distribution as in Q1 but is
#' re-paired with Y(0) through a Gaussian copula with correlation rho_r: with z0 = (log Y(0) - m) / sqrt(v)
#' and xi ~ N(0, 1), log Y(1) = logc + k (m + sqrt(v) (rho_r z0 + sqrt(1 - rho_r^2) xi)).
gen_Q1r <- function(P, Q1) {
  xi <- rnorm(nrow(Q1))
  z0 <- (log(Q1$y0_true) - P$m) / sqrt(P$v)
  y1 <- exp(P$logc + P$k * (P$m + sqrt(P$v) * (P$rho_r * z0 + sqrt(1 - P$rho_r^2) * xi)))
  data.frame(merchant = Q1$merchant, x = Q1$x, D = Q1$D, y = ifelse(Q1$D == 1, y1, Q1$y0_true),
             y0_true = Q1$y0_true, y1_true = y1)
}

#' Q2 (Fulfilled by Marigold). P: n, n_fbm, mu, g, bx, s0, s1. log monthly revenue (USD):
#' log Y = mu + g W + bx X + s_W e with W the FBM indicator (exactly n_fbm merchants, about the population
#' share p, placed at random), X ~ N(0, 1) the standardised log catalogue size (independent of W), e ~ N(0, 1).
gen_Q2 <- function(P) {
  n <- P$n
  w <- sample(rep(0:1, c(n - P$n_fbm, P$n_fbm))); x <- rnorm(n); e <- rnorm(n)
  ly <- P$mu + P$g * w + P$bx * x + ifelse(w == 1, P$s1, P$s0) * e
  data.frame(merchant = seq_len(n), fbm = w, x = x, logy = ly, y = exp(ly))
}

#' Q2 truth helpers: the cdf and density of log Y at q when the FBM share is p (X integrated out):
#' log Y | W = w ~ N(mu + g w, bx^2 + s_w^2).
mix_cdf_Q2 <- function(q, P, p = P$p) {
  p * pnorm(q, P$mu + P$g, sqrt(P$bx^2 + P$s1^2)) + (1 - p) * pnorm(q, P$mu, sqrt(P$bx^2 + P$s0^2))
}
mix_pdf_Q2 <- function(q, P, p = P$p) {
  p * dnorm(q, P$mu + P$g, sqrt(P$bx^2 + P$s1^2)) + (1 - p) * dnorm(q, P$mu, sqrt(P$bx^2 + P$s0^2))
}

#' Q3 (AdBoost adoption with a randomised ad-credit offer). P: n, n_z1, mu0, dmu, s0, s1, bx, rho, pi0, pi1.
#' Z: exactly n_z1 merchants offered the credit. (Z_U, Z_V) standard bivariate normal with correlation rho;
#' D(z) = 1{Phi(Z_V) > 1 - pi_z}; log Y(d) = mu0 + d dmu + s_d Z_U + bx X, X ~ N(0, 1) independent.
#' Always-takers Phi(Z_V) > 1 - pi0, compliers 1 - pi1 < Phi(Z_V) <= 1 - pi0, never-takers the rest.
gen_Q3 <- function(P) {
  n <- P$n
  Z <- sample(rep(0:1, c(n - P$n_z1, P$n_z1)))
  zv <- rnorm(n); eta <- rnorm(n); x <- rnorm(n)
  zu <- P$rho * zv + sqrt(1 - P$rho^2) * eta
  D0 <- as.integer(pnorm(zv) > 1 - P$pi0); D1 <- as.integer(pnorm(zv) > 1 - P$pi1)
  D <- ifelse(Z == 1, D1, D0)
  ly0 <- P$mu0 + P$s0 * zu + P$bx * x
  ly1 <- P$mu0 + P$dmu + P$s1 * zu + P$bx * x
  ly <- ifelse(D == 1, ly1, ly0)
  data.frame(merchant = seq_len(n), x = x, Z = Z, D = D, logy = ly, y = exp(ly), D0_true = D0, D1_true = D1,
             type_true = ifelse(D0 == 1, "always", ifelse(D1 == 1, "complier", "never")),
             logy0_true = ly0, logy1_true = ly1, zu_true = zu)
}

#' Q4 (merchant panel). P: N, mu, sa, cx, beta, sig, gam; Tn: periods (\mathcal{T}). log monthly revenue:
#' Y_it = mu + alpha_i + beta X_it + sig (1 + gam X_it) e_it, alpha_i ~ N(0, sa^2), X_it = plogis(cx alpha_i /
#' sa + nu_it) in (0, 1) the share of the catalogue on promotion, nu_it, e_it ~ N(0, 1). Hence
#' Quant_tau(Y_it | X_it, alpha_i) = (mu + alpha_i + sig z_tau) + (beta + sig gam z_tau) X_it: alpha_i is a
#' pure location shift (Canay 2011) and beta(tau) = beta + sig gam z_tau.
gen_Q4 <- function(P, Tn) {
  N <- P$N
  id <- rep(seq_len(N), each = Tn)
  a <- rnorm(N, 0, P$sa)[id]
  x <- plogis(P$cx * a / P$sa + rnorm(N * Tn))
  e <- rnorm(N * Tn)
  data.frame(id = id, t = rep(seq_len(Tn), N), x = x, y = P$mu + a + P$beta * x + P$sig * (1 + P$gam * x) * e,
             alpha_true = a)
}

#' Q5 structural function h(z, t) = A exp(g t + s (1 + kap t) z) (monthly revenue, USD thousands): increasing
#' in the rank z; the period-1 map y -> h(h^{-1}(y, 0), 1) is a power function, a time effect that is
#' nonlinear in the outcome (growth larger for larger merchants).
h_Q5 <- function(z, t, P) P$A * exp(P$g * t + P$s * (1 + P$kap * t) * z)

#' Q5 (sharp CIC: One-Click Checkout in the treated region). P: nc (merchants per cell), A, g, s, kap, mu1,
#' sd1, delta. Four repeated cross sections (S_i, T_i) in {0, 1}^2; rank z ~ N(0, 1) in the control region,
#' N(mu1, sd1^2) in the treated region, the same in both periods; Y(0) = h(z, t), Y(1) = (1 + delta) Y(0);
#' D = S T.
gen_Q5 <- function(P) {
  nc <- P$nc
  S <- rep(c(0, 0, 1, 1), each = nc); Tt <- rep(c(0, 1, 0, 1), each = nc)
  zt <- rnorm(4 * nc)
  z <- ifelse(S == 1, P$mu1 + P$sd1 * zt, zt)
  y0 <- h_Q5(z, Tt, P); y1 <- (1 + P$delta) * y0
  D <- S * Tt
  data.frame(merchant = seq_len(4 * nc), S = S, T = Tt, D = D, y = ifelse(D == 1, y1, y0),
             y0_true = y0, y1_true = y1, z_true = z)
}

#' Q5f (fuzzy CIC: voluntary adoption of Instant Payouts). P: as Q5 plus rho, p_at1, p_sw1, p_at0, tau.
#' Adoption type from v = Phi(Z_V): treated region always-takers v > 1 - p_at1, switchers
#' 1 - p_at1 - p_sw1 < v <= 1 - p_at1 (adopt in period 1 only), never-takers below; control region
#' always-takers v > 1 - p_at0 (a stable adoption rate), never-takers below. Rank z = mu_S + sd_S (rho Z_V +
#' sqrt(1 - rho^2) eta) (mu_0 = 0, sd_0 = 1; mu_1, sd_1 as in Q5): adopters have higher ranks.
#' Y(0) = h(z, t), Y(1) = (1 + delta) h(z, t); D = 1 for always-takers and for switchers in period 1.
gen_Q5f <- function(P) {
  nc <- P$nc
  S <- rep(c(0, 0, 1, 1), each = nc); Tt <- rep(c(0, 1, 0, 1), each = nc)
  zv <- rnorm(4 * nc); eta <- rnorm(4 * nc)
  zt <- P$rho * zv + sqrt(1 - P$rho^2) * eta
  z <- ifelse(S == 1, P$mu1 + P$sd1 * zt, zt)
  v <- pnorm(zv)
  type <- ifelse(S == 1,
                 ifelse(v > 1 - P$p_at1, "always", ifelse(v > 1 - P$p_at1 - P$p_sw1, "switcher", "never")),
                 ifelse(v > 1 - P$p_at0, "always", "never"))
  D <- as.integer(type == "always" | (type == "switcher" & Tt == 1))
  data.frame(merchant = seq_len(4 * nc), S = S, T = Tt, D = D, y = h_Q5(z, Tt, P) * (1 + P$delta)^D,
             type_true = type, z_true = z)
}
