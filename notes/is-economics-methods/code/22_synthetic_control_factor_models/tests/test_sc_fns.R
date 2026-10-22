# tests/test_sc_fns.R
# Unit tests for the hand-written synthetic-control functions in sc_fns.R.
# Run from the code directory AFTER 01_dgp.R (the KKT checks read aster_S1.rds and
# aster_S2.rds):  Rscript tests/test_sc_fns.R   (a few seconds)
#
# Part 1 is the Task 9 brief's test, verbatim except that its Chinese comments are
# translated to English (Ruling S7; assertions unchanged). Part 2 adds checks on the
# parts of the interface the brief's toy problems do not reach (names, the weighting
# matrix V, the scaling of the penalty, the moving-block permutation grid, the
# placebo MSPE-ratio permutation and its guards, the convex-hull LP), and the KKT
# conditions of the QP solver on the real S1 and S2 pre-periods (24 x 20).
#
# Hand calculations:
#   sc_weights toy: X1 = (1, 2) = w1 (0, 2) + w2 (2, 2) + w3 (5, 5) with w on the
#     simplex. Row 1: 2 w2 + 5 w3 = 1; row 2: 2 (w1 + w2) + 5 w3 = 2, and
#     w1 + w2 = 1 - w3, so 3 w3 = 0: w = (0.5, 0.5, 0) is the unique exact fit.
#   penalized_sc, one predictor, donors at 0 and 3, X1 = 1, V = 1:
#     (1 - 3 w)^2 + pen [(1 - w) * 1 + w * 4]; the derivative in w (weight on the
#     second donor) is -6 (1 - 3 w) + 3 pen = 0, so w = 1/3 - pen / 6;
#     pen = 1 gives (5/6, 1/6), pen = 0 gives (2/3, 1/3) (exact fit).
#   V matters: donors (0, 0) and (2, 10), X1 = (1, 0). V = diag(1, 0) gives
#     w2 = 0.5; V = I gives (1 - 2 w2)^2 + (10 w2)^2, minimised at w2 = 4 / 208.
source("sc_fns.R")

## ---- Part 1: the brief ------------------------------------------------------------
X0 <- cbind(c(0, 2), c(2, 2), c(5, 5)); X1 <- c(1, 2)                   # unique solution: 0.5 d1 + 0.5 d2
stopifnot(all(abs(sc_weights(X1, X0, diag(2)) - c(0.5, 0.5, 0)) < 1e-6))
wp <- penalized_sc(c(4.9, 5.1), X0, pen = 1e6)$w                        # a huge penalty reduces to the nearest neighbour
stopifnot(which.max(wp) == 3, abs(max(wp) - 1) < 1e-6)
set.seed(3); ps <- replicate(500, conformal_p(rnorm(30), rep(0, 30), T0 = 24, scheme = "iid"))
stopifnot(abs(mean(ps <= 0.10) - 0.10) < 0.04)                          # p-values about uniform under the null

## ---- Part 2: sc_weights -------------------------------------------------------------
# names carried over from colnames(X0); simplex; exact fit when X1 is in the hull
# (8 donors, 6 predictors: the weights need not be unique, so test the fit, not w)
set.seed(11)
X0r <- matrix(rnorm(6 * 8), 6, 8, dimnames = list(NULL, paste0("d", 1:8)))
X1r <- drop(X0r %*% c(0.5, 0.3, 0.2, 0, 0, 0, 0, 0))
wr <- sc_weights(X1r, X0r, diag(6))
stopifnot(identical(names(wr), colnames(X0r)), all(wr >= 0), abs(sum(wr) - 1) < 1e-10,
          max(abs(X1r - X0r %*% wr)) < 1e-6)
# the weighting matrix V (the chapter's Omega) changes the answer
Xv0 <- cbind(c(0, 0), c(2, 10)); Xv1 <- c(1, 0)
stopifnot(abs(sc_weights(Xv1, Xv0, diag(c(1, 0)))[2] - 0.5) < 1e-6,
          abs(sc_weights(Xv1, Xv0, diag(2))[2] - 4 / 208) < 1e-6)
# V may be given as the vector of its diagonal
stopifnot(max(abs(sc_weights(X1, X0, c(1, 1)) - sc_weights(X1, X0, diag(2)))) < 1e-12)

## ---- Part 2: penalized_sc ----------------------------------------------------------
p1 <- penalized_sc(1, matrix(c(0, 3), 1), pen = 1)
p0 <- penalized_sc(1, matrix(c(0, 3), 1), pen = 0)
stopifnot(max(abs(p1$w - c(5/6, 1/6))) < 1e-6, max(abs(p0$w - c(2/3, 1/3))) < 1e-6,
          p1$n_nonzero == 2)
pt <- penalized_sc(X1, X0, pen = 0)                                      # pen = 0 is plain SC
stopifnot(max(abs(pt$w - c(0.5, 0.5, 0))) < 1e-6, pt$n_nonzero == 2)
stopifnot(all(wp >= 0), abs(sum(wp) - 1) < 1e-10)
# outside the hull: (3, 1.5) projects onto the edge d2-d3 at (2.25, 2.25); along that
# edge (w3 = t) the objective is (1 - 3t)^2 + (0.5 + 3t)^2 + pen (1.25 + 15 t), so
# t = (3 - 15 pen) / 36: two donors for pen < 0.2, the nearest donor d2 alone after
nn <- sapply(c(0, 0.01, 0.1, 1, 10), function(p) penalized_sc(c(3, 1.5), X0, pen = p)$n_nonzero)
stopifnot(identical(nn, c(2L, 2L, 2L, 1L, 1L)),
          max(abs(penalized_sc(c(3, 1.5), X0, pen = 0.1)$w - c(0, 23/24, 1/24))) < 1e-6,
          max(abs(penalized_sc(c(3, 1.5), X0, pen = 1)$w - c(0, 1, 0))) < 1e-9)

## ---- Part 2: conformal_p ------------------------------------------------------------
set.seed(5)
u <- rnorm(30)
pm <- conformal_p(u, rep(0, 30), T0 = 24, scheme = "moving_block")
stopifnot(pm > 0, pm <= 1, abs(pm * 30 - round(pm * 30)) < 1e-9)          # 30 cyclic shifts
ybig <- c(rnorm(24, 0, 0.1), rnorm(6, 5, 0.1))                          # a large post-period effect
stopifnot(abs(conformal_p(ybig, rep(0, 30), T0 = 24, scheme = "moving_block") - 1/30) < 1e-12,
          conformal_p(ybig, rep(0, 30), T0 = 24, scheme = "iid") <= 0.01)
# the p-value depends on y1 - y0hat only
stopifnot(abs(conformal_p(u + 7, rep(7, 30), T0 = 24, scheme = "moving_block") - pm) < 1e-12)
# moving-block p-values are exactly uniform on {1/30, ..., 1} under iid noise
set.seed(6); pmb <- replicate(600, conformal_p(rnorm(30), rep(0, 30), T0 = 24))
stopifnot(abs(mean(pmb <= 0.10) - 0.10) < 0.04)

## ---- Part 2: mspe_ratio_perm --------------------------------------------------------
# 6 units x 20 months, 15 pre-periods; unit 3 is treated (+3 in every post month);
# fit_fun = the plain donor mean, so every number can be recomputed by hand
set.seed(8)
Ym <- matrix(rnorm(6 * 20, 0, 0.1), 6, 20, dimnames = list(paste0("m", 1:6), NULL))
Ym[3, 16:20] <- Ym[3, 16:20] + 3
mean_fit <- function(y1, Y0, T0) rowMeans(Y0)
mr <- mspe_ratio_perm(Ym, treated = 3, T0 = 15, fit_fun = mean_fit, placebo_donors = "exclude_treated")
g3 <- Ym[3, ] - colMeans(Ym[-3, ])
stopifnot(abs(mr$ratio - mean(g3[16:20]^2) / mean(g3[1:15]^2)) < 1e-12,
          mr$rank == 1, abs(mr$p - 1/6) < 1e-12,
          identical(names(mr$all_ratios), rownames(Ym)))
# placebo for unit 1 under "exclude_treated": its donors leave out the treated unit 3
g1 <- Ym[1, ] - colMeans(Ym[c(2, 4, 5, 6), ])
stopifnot(abs(mr$all_ratios[["m1"]] - mean(g1[16:20]^2) / mean(g1[1:15]^2)) < 1e-12)
# the default is ADH (2010), who shift the treated unit into every placebo's donor
# pool (also tidysynth 0.2.1): unit 1's donors are all the other units, unit 3 included
mr2 <- mspe_ratio_perm(Ym, treated = 3, T0 = 15, fit_fun = mean_fit)
g1b <- Ym[1, ] - colMeans(Ym[-1, ])
stopifnot(abs(mr2$all_ratios[["m1"]] - mean(g1b[16:20]^2) / mean(g1b[1:15]^2)) < 1e-12,
          abs(mr2$ratio - mr$ratio) < 1e-12,                                   # the treated unit's own fit is the same
          identical(mr2, mspe_ratio_perm(Ym, treated = 3, T0 = 15, fit_fun = mean_fit, placebo_donors = "all_others")))
# sc_fit: the SC counterfactual on pre-period outcomes, usable as fit_fun
ms <- mspe_ratio_perm(Ym, treated = 3, T0 = 15, fit_fun = sc_fit)
stopifnot(ms$rank == 1, abs(ms$p - 1/6) < 1e-12)
# guards: T0 must be an integer; an exact pre-period fit (zero pre-MSPE) is refused
err <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
stopifnot(grepl("T0 must be an integer", err(mspe_ratio_perm(Ym, 3, 15.5, mean_fit))),
          grepl("T0 must be an integer", err(mspe_ratio_perm(Ym, 3, 20, mean_fit))))
Yz <- Ym; Yz[1, ] <- colMeans(Ym[-1, ])                                      # unit 1 = mean of the others
stopifnot(grepl("pre-period MSPE of unit 1 \\(m1\\) is zero", err(mspe_ratio_perm(Yz, 3, 15, mean_fit))))
# NA-safe guards (Task 10): NA T0 or treated and a fit returning NA stop with a message,
# not with R's "missing value where TRUE/FALSE needed"
na_fit <- function(y1, Y0, T0) rep(NA_real_, length(y1))
stopifnot(grepl("T0 must be an integer", err(mspe_ratio_perm(Ym, 3, NA, mean_fit))),
          grepl("treated must name one row", err(mspe_ratio_perm(Ym, NA, 15, mean_fit))),
          grepl("fit_fun returned missing", err(mspe_ratio_perm(Ym, 3, 15, na_fit))),
          grepl("T0 must be an integer", err(conformal_p(1:10, rep(0, 10), NA))))

## ---- Part 2: KKT conditions of the QP solver on the S1 and S2 pre-periods -----------------
# F(w) = ||X1 - X0 w||^2 + pen sum_j w_j ||X1 - X0[, j]||^2 on the simplex, X0 = 24 pre-months
# x 20 donors. At the optimum the gradient is equal (= mu) on the support and >= mu off it.
# Gradients are scaled by mean(diag(2 X0'X0)), the solver's own scaling. Weights must not
# depend on the order of the donors.
kkt <- function(X1, X0, pen) {
  w <- penalized_sc(X1, X0, pen = pen)$w
  d2 <- colSums((X1 - X0)^2)
  g <- drop(2 * crossprod(X0, X0 %*% w - X1) + pen * d2) / mean(diag(2 * crossprod(X0)))
  S <- w > 0
  mu <- mean(g[S])
  c(spread = diff(range(g[S])), min_mult = if (all(S)) 0 else min(g[!S] - mu),
    simplex = abs(sum(w) - 1) + max(0, -min(w)))
}
if (!all(file.exists(c("aster_S1.rds", "aster_S2.rds"))))
  stop("test_sc_fns.R: aster_S1.rds and aster_S2.rds are missing; run 01_dgp.R first (Rscript 01_dgp.R)", call. = FALSE)
set.seed(12)
for (blk in c("aster_S1.rds", "aster_S2.rds")) {
  Yb <- readRDS(blk)$Y
  X0b <- t(Yb[-1, 1:24]); X1b <- Yb[1, 1:24]
  for (pen in c(0, 0.01, 0.1, 1)) {
    k <- kkt(X1b, X0b, pen)
    stopifnot(k[["spread"]] < 1e-9, k[["min_mult"]] > -1e-9, k[["simplex"]] < 1e-12)
    w0 <- penalized_sc(X1b, X0b, pen = pen)$w
    for (r in 1:3) {
      o <- sample(ncol(X0b))
      wo <- penalized_sc(X1b, X0b[, o], pen = pen)$w
      stopifnot(max(abs(wo[order(o)] - w0)) < 1e-10)
    }
  }
}

## ---- Part 2: convex hull by LP feasibility -------------------------------------------
stopifnot(isTRUE(in_hull_lp(c(1, 2), X0)), isFALSE(in_hull_lp(c(6, 7), X0)),
          isTRUE(in_hull_lp(c(0, 2), X0)),                                  # a vertex is in the hull
          isTRUE(in_hull_lp(X1r, X0r)), isFALSE(in_hull_lp(X1r + c(50, rep(0, 5)), X0r)))

cat("test_sc_fns.R: all checks passed\n")
