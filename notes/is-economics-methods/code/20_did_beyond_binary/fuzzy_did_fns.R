# ------------------------------------------------------------------
# fuzzy_did_fns.R -- hand-coded fuzzy difference-in-differences estimators
#
#   wald_did(df, y, d, group, time) -> list(estimate, did_y, did_d)
#   wald_tc (df, y, d, group, time) -> list(estimate, delta, first_stage, did_d)
#   wald_cic(df, y, d, group, time) -> list(estimate, first_stage, did_d)
#
# Source: de Chaisemartin, C. and D'Haultfoeuille, X. (2018). Fuzzy
# differences-in-differences. Review of Economic Studies 85(2), 999-1028,
# doi:10.1093/restud/rdx049. Checked (Task 10, Step 0, 2026-10-20) against
# the preprint arXiv:1510.01757v4 (16 Feb 2017). The paper does not number
# these displays; references are to its sections/theorems. Assumption and
# theorem numbers below are those of the published Review of Economic Studies
# version (the preprint numbers them differently), as in the book chapter.
#
# Notation (paper, Section 2.1): G is the group dummy, T the period dummy;
# R_gt ~ R | G = g, T = t and R_dgt ~ R | D = d, G = g, T = t, so the FIRST
# index is the group and the SECOND the period. Here `group` plays G (the
# book's S_i, notation_v3.md 2.2) and `time` plays T.
#
#   Wald-DID (Section 2.2):  W_DID = DID_Y / DID_D,
#       DID_R = E(R11) - E(R10) - (E(R01) - E(R00)).
#   Wald-TC  (Section 2.3):  W_TC = [E(Y11) - E(Y10 + delta_{D10})] / [E(D11) - E(D10)],
#       delta_d = E(Y_d01) - E(Y_d00)   (control-group change among units with D = d).
#   Wald-CIC (Section 2.4):  W_CIC = [E(Y11) - E(Q_{D10}(Y10))] / [E(D11) - E(D10)],
#       Q_d(y) = F^{-1}_{Y_d01}(F_{Y_d00}(y)),  F^{-1}(q) = inf{y : F(y) >= q}.
#   Sample versions as in Section 4 (cell means; empirical cdf and quantile).
#
# Step 0 result: the brief's interface formulas coincide with the paper's
# definitions; no change to the estimators or to the toy values (15.0, 13.8).
# What the brief's interface leaves out is the identification condition:
#   Theorem 1: W_DID = Delta needs Assumptions 1, 3, 4, 5 (common trends; LATE of
#     units treated at T = 0 stable over time) and, unless the control-group
#     treatment rate is stable, Assumption 6 (Delta = Delta').
#   Theorem 2: W_TC = Delta needs Assumptions 1, 3, 4' (common trends within
#     subgroups sharing the period-0 treatment) and 2 (stable treatment rate in
#     the control group, 0 < E(D01) = E(D00) < 1).
#   Theorem 3: W_CIC = Delta needs Assumptions 1 to 3, 7, 8.
#   Section 3.1, Theorem 4: without Assumption 2 the paper only bounds Delta
#     (TC bounds, CIC bounds; the TC bounds are infinite when the outcome is
#     unbounded); W_TC and W_CIC are then not point-identifying in general.
# Tessera (01_dgp.R) follows the paper: block A has a stable control-group
# treatment rate, so W_TC = W_CIC = Delta, while W_DID is biased because
# Assumption 5 fails; variant A' lets the control rate rise to show the
# Theorem 4 case, where W_TC and W_CIC miss Delta.
#
# Edge cases (brief): `group`, `time`, `d` must take values in {0, 1}; when the
# first-stage DiD |DID_D| < 1e-8 every function returns estimate = NA_real_
# with warning "first-stage DiD is (near) zero". Wald-TC and Wald-CIC also
# return NA with a warning when their own denominator E(D11) - E(D10) is
# (near) zero or when delta_d / Q_d is undefined because a needed control cell
# (D = d, G = 0, T = 0 or 1) is empty.
#
# Cross-check (brief Step 5, end of tests/test_fuzzy_fns.R): Rfuzzydid 1.3.0
# (third-party port of Stata fuzzydid) gives the same Wald-DID and Wald-TC on
# tessera_A and tessera_A_rising (|difference| < 1e-6). Its Wald-CIC uses a
# floating-point type-1 quantile for F^{-1}(k/n) and so keeps the off-by-one
# fixed in .fz_qq below; it matches this file's Wald-CIC computed with that
# convention, and differs from the exact generalized inverse when a control
# cell has the same size in both periods (6.4e-4 on tessera_A). These
# functions are authoritative.
# ------------------------------------------------------------------

## internal: validate inputs and return plain vectors
.fz_prepare <- function(df, y, d, group, time) {
  df <- as.data.frame(df)
  vars <- c(y, d, group, time)
  if (!is.character(vars) || length(vars) != 4)
    stop("y, d, group, time must be single column names (strings)")
  miss <- setdiff(vars, names(df))
  if (length(miss)) stop("column(s) not found: ", paste(miss, collapse = ", "))
  Y <- as.numeric(df[[y]]); D <- as.numeric(df[[d]])
  G <- as.numeric(df[[group]]); Tt <- as.numeric(df[[time]])
  if (anyNA(Y) || anyNA(D) || anyNA(G) || anyNA(Tt)) stop("missing values in y, d, group or time")
  for (nm in c("d", "group", "time")) {
    v <- switch(nm, d = D, group = G, time = Tt)
    if (!all(v %in% c(0, 1))) stop("`", nm, "` must take values in {0, 1}")
  }
  for (g in 0:1) for (t in 0:1)
    if (!any(G == g & Tt == t)) stop(sprintf("empty cell: group %d, time %d", g, t))
  list(Y = Y, D = D, G = G, T = Tt)
}

## internal: cell mean E(R_gt)
.fz_cell <- function(R, x, g, t) mean(R[x$G == g & x$T == t])

## internal: first-stage DiD of the treatment
.fz_did_d <- function(x) {
  (.fz_cell(x$D, x, 1, 1) - .fz_cell(x$D, x, 1, 0)) -
    (.fz_cell(x$D, x, 0, 1) - .fz_cell(x$D, x, 0, 0))
}

## internal: shared NA-with-warning exits
.fz_na_first_stage <- function(extra = list()) {
  warning("first-stage DiD is (near) zero: the two groups' treatment rates change by the same amount",
          call. = FALSE)
  c(list(estimate = NA_real_), extra)
}

wald_did <- function(df, y, d, group, time) {
  x <- .fz_prepare(df, y, d, group, time)
  did <- function(R) (.fz_cell(R, x, 1, 1) - .fz_cell(R, x, 1, 0)) -
                     (.fz_cell(R, x, 0, 1) - .fz_cell(R, x, 0, 0))
  did_y <- did(x$Y); did_d <- did(x$D)
  if (abs(did_d) < 1e-8) return(.fz_na_first_stage(list(did_y = did_y, did_d = did_d)))
  list(estimate = did_y / did_d, did_y = did_y, did_d = did_d)
}

## internal: checks shared by Wald-TC and Wald-CIC; returns the treatment
## values present in cell (G = 1, T = 0) or a list with estimate = NA
.fz_tc_cic_checks <- function(x) {
  did_d <- .fz_did_d(x)
  if (abs(did_d) < 1e-8) return(.fz_na_first_stage(list(did_d = did_d)))
  fs <- .fz_cell(x$D, x, 1, 1) - .fz_cell(x$D, x, 1, 0)
  if (abs(fs) < 1e-8) {
    warning("treatment-group first stage E(D11) - E(D10) is (near) zero", call. = FALSE)
    return(list(estimate = NA_real_, first_stage = fs, did_d = did_d))
  }
  d_needed <- sort(unique(x$D[x$G == 1 & x$T == 0]))
  for (dd in d_needed) for (t in 0:1)
    if (!any(x$D == dd & x$G == 0 & x$T == t)) {
      warning(sprintf("delta_%d / Q_%d undefined: no control-group unit with D = %d at time %d",
                      dd, dd, dd, t), call. = FALSE)
      return(list(estimate = NA_real_, first_stage = fs, did_d = did_d))
    }
  list(ok = TRUE, d_needed = d_needed, first_stage = fs, did_d = did_d)
}

wald_tc <- function(df, y, d, group, time) {
  x <- .fz_prepare(df, y, d, group, time)
  ck <- .fz_tc_cic_checks(x)
  if (is.null(ck$ok)) return(ck)
  # delta_d = E(Y_d01) - E(Y_d00), for d in {0, 1} when defined
  delta <- c("0" = NA_real_, "1" = NA_real_)
  for (dd in 0:1) {
    s0 <- x$D == dd & x$G == 0 & x$T == 0; s1 <- x$D == dd & x$G == 0 & x$T == 1
    if (any(s0) && any(s1)) delta[as.character(dd)] <- mean(x$Y[s1]) - mean(x$Y[s0])
  }
  in10 <- x$G == 1 & x$T == 0
  counterfactual <- mean(x$Y[in10] + delta[as.character(x$D[in10])])   # E(Y10 + delta_{D10})
  num <- .fz_cell(x$Y, x, 1, 1) - counterfactual
  list(estimate = num / ck$first_stage, delta = delta, first_stage = ck$first_stage, did_d = ck$did_d)
}

## internal: quantile-quantile transform Q_d(y) = F^{-1}_{Y_d01}(F_{Y_d00}(y)),
## F the right-continuous empirical cdf, F^{-1}(q) = inf{y : F(y) >= q}.
## Exact integer arithmetic (fix round 1): with k0 = #{y_d00 <= y} = n0 F(y),
## the generalized inverse at q = k0 / n0 is the j-th order statistic of
## y_d01 with j = ceiling(n1 k0 / n0), the smallest j such that j n0 >= n1 k0.
## (The first version evaluated stats::quantile(type = 1) at probs = k0/n0;
## n1 * (k0/n0) can land just above an integer in floating point and returned
## the (j+1)-th order statistic, e.g. Q(7) = 108 instead of 107 for
## y_d00 = 1:25, y_d01 = 101:125.) Ties and unequal cell sizes are handled by
## the counts; y below min(y_d00) gives k0 = 0 (q = 0) and maps to min(y_d01).
## Products are formed in double precision, exact below 2^53, so large cells
## cannot overflow R's 32-bit integers; one correction step removes any
## rounding in the division.
.fz_qq <- function(yeval, y_d00, y_d01) {
  s00 <- sort(y_d00); s01 <- sort(y_d01)
  n0 <- length(s00); n1 <- length(s01)
  k0 <- findInterval(yeval, s00)                  # #{y_d00 <= y}, ties included
  num <- as.double(n1) * k0
  j <- ceiling(num / n0)
  j <- j + (j * n0 < num) - ((j - 1) * n0 >= num)
  s01[pmax(j, 1)]
}

wald_cic <- function(df, y, d, group, time) {
  x <- .fz_prepare(df, y, d, group, time)
  ck <- .fz_tc_cic_checks(x)
  if (is.null(ck$ok)) return(ck)
  in10 <- x$G == 1 & x$T == 0
  qy <- numeric(sum(in10)); y10 <- x$Y[in10]; d10 <- x$D[in10]
  for (dd in ck$d_needed) {
    k <- d10 == dd
    qy[k] <- .fz_qq(y10[k], x$Y[x$D == dd & x$G == 0 & x$T == 0], x$Y[x$D == dd & x$G == 0 & x$T == 1])
  }
  num <- .fz_cell(x$Y, x, 1, 1) - mean(qy)                                # E(Y11) - E(Q_{D10}(Y10))
  list(estimate = num / ck$first_stage, first_stage = ck$first_stage, did_d = ck$did_d)
}
