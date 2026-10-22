# tests/test_rd_fns.R
# Unit tests for the hand-coded RD helpers in rd_fns.R (ll_jump, rk_ratio,
# lr_diffmeans, pooled_cutoff_weights) and the block generators gen_*().
# Run from the code directory:  Rscript tests/test_rd_fns.R   (a few seconds)
#
# Part 1 is the Task 6 brief's test, verbatim (toy data, exact values).
# Part 2 adds hand-checkable cases and non-vacuous guards (ruling R2).
source("rd_fns.R"); source("tests/helpers.R")      # helpers.R: is_block()

## ---- Part 1: the brief ------------------------------------------------------------
x <- seq(-1, 1, by = 0.01); y <- 1 + 0.5 * x + 2 * (x >= 0)
stopifnot(abs(ll_jump(y, x, c = 0, h = 0.5)$est - 2) < 1e-10)
b <- 0.1 * x + 0.2 * pmax(x, 0)
stopifnot(abs(rk_ratio(y = 3 * b, b = b, x = x, c = 0, h = 0.5)$est - 3) < 1e-10)
set.seed(1); z <- runif(400, -1, 1); yz <- rnorm(400) + 1 * (z >= 0)
lr <- lr_diffmeans(yz, z, c = 0, w0 = 0.5, reps = 999, seed = 2); stopifnot(lr$p < 0.01, abs(lr$est - 1) < 0.3)
w <- pooled_cutoff_weights(rep(c(75, 80, 85), each = 100), runif(300, 60, 100), h = 5); stopifnot(abs(sum(w) - 1) < 1e-12, all(w >= 0))

## ---- Part 2: ll_jump ----------------------------------------------------------------
j <- ll_jump(y, x, c = 0, h = 0.5)
stopifnot(identical(names(j), c("est", "n_left", "n_right")))
# window counts follow rdrobust's N_h: |x - c| <= h on each side (x = c is on the right;
# seq(-1, 1, by = 0.01) hits -0.5, 0 and 0.5 exactly)
stopifnot(j$n_left == sum(x < 0 & abs(x) <= 0.5), j$n_right == sum(x >= 0 & abs(x) <= 0.5),
          j$n_left == 50, j$n_right == 51)
# deriv = 1 returns the slope change: zero for the brief's y, 0.2 for b, under every kernel
for (k in c("triangular", "epanechnikov", "uniform")) {
  stopifnot(abs(ll_jump(y, x, c = 0, h = 0.5, kernel = k)$est - 2) < 1e-10,
            abs(ll_jump(y, x, c = 0, h = 0.5, kernel = k, deriv = 1)$est) < 1e-10,
            abs(ll_jump(b, x, c = 0, h = 0.5, kernel = k, deriv = 1)$est - 0.2) < 1e-10)
}
# a cutoff away from zero and a different intercept / slope on each side
x2 <- seq(60, 100, by = 0.25); y2 <- ifelse(x2 >= 80, 7 + 0.3 * (x2 - 80), 6.5 - 0.1 * (x2 - 80))
stopifnot(abs(ll_jump(y2, x2, c = 80, h = 6)$est - 0.5) < 1e-10,
          abs(ll_jump(y2, x2, c = 80, h = 6, deriv = 1)$est - 0.4) < 1e-10)
# the weighted fit uses only the window: points outside |x - c| <= h do not move it
y3 <- y2; y3[abs(x2 - 80) > 6] <- 1e6
stopifnot(abs(ll_jump(y3, x2, c = 80, h = 6)$est - 0.5) < 1e-8)
# triangular weights by hand on noisy data: matches lm() with weights = 1 - |u|
set.seed(3); xn <- runif(500, -2, 2); yn <- 1 + xn + 0.4 * (xn >= 0) + rnorm(500, 0, 0.2)
wl <- pmax(1 - abs(xn) / 1.2, 0)
fitR <- lm(yn ~ xn, weights = wl, subset = xn >= 0 & wl > 0)
fitL <- lm(yn ~ xn, weights = wl, subset = xn < 0 & wl > 0)
stopifnot(abs(ll_jump(yn, xn, c = 0, h = 1.2)$est - (coef(fitR)[[1]] - coef(fitL)[[1]])) < 1e-10,
          abs(ll_jump(yn, xn, c = 0, h = 1.2, deriv = 1)$est - (coef(fitR)[[2]] - coef(fitL)[[2]])) < 1e-10)
# bad input is refused
stopifnot(inherits(try(ll_jump(y, x, c = 0, h = -1), silent = TRUE), "try-error"),
          inherits(try(ll_jump(y, x, c = 0, h = 0.5, kernel = "gaussian"), silent = TRUE), "try-error"),
          inherits(try(ll_jump(y, x, c = 0, h = 0.5, deriv = 2), silent = TRUE), "try-error"),
          inherits(try(ll_jump(y, x, c = 0, h = 0.005), silent = TRUE), "try-error"))   # one point per side

## ---- Part 2: rk_ratio -------------------------------------------------------------------
rk <- rk_ratio(y = 3 * b, b = b, x = x, c = 0, h = 0.5)
stopifnot(identical(names(rk), c("est", "dslope_y", "dslope_b")),
          abs(rk$dslope_b - 0.2) < 1e-10, abs(rk$dslope_y - 0.6) < 1e-10)
# a level jump in y does not enter the kink ratio
stopifnot(abs(rk_ratio(y = 3 * b + 5 * (x >= 0), b = b, x = x, c = 0, h = 0.5)$est - 3) < 1e-10)
# no kink in the policy: refused rather than dividing by zero
stopifnot(inherits(try(rk_ratio(y = x, b = 0.1 * x, x = x, c = 0, h = 0.5), silent = TRUE), "try-error"))
# a bad policy vector is blamed on b, not on y
err_b <- try(rk_ratio(y = x, b = c(NA, 0.1 * x[-1]), x = x, c = 0, h = 0.5), silent = TRUE)
stopifnot(inherits(err_b, "try-error"), grepl("b and x", conditionMessage(attr(err_b, "condition")), fixed = TRUE))

## ---- Part 2: lr_diffmeans ------------------------------------------------------------------
stopifnot(identical(names(lr), c("est", "p")), lr$p >= 1 / (999 + 1), lr$p <= 1)   # valid p is never 0
inw <- abs(z) <= 0.5
stopifnot(abs(lr$est - (mean(yz[inw & z >= 0]) - mean(yz[inw & z < 0]))) < 1e-12)
# a strong effect: no permutation reaches the observed statistic, so p is exactly the
# smallest valid Monte Carlo value 1 / (reps + 1), not 0
lr_big <- lr_diffmeans(yz + 3 * (z >= 0), z, c = 0, w0 = 0.5, reps = 999, seed = 2)
stopifnot(lr_big$p >= 1 / (999 + 1), abs(lr_big$p - 1 / 1000) < 1e-15)
# the p formula by hand: (1 + #{|T*| >= |T_obs|}) / (1 + reps), same permutations
# (fixed margins, sample() after set.seed(seed)) on null data
set.seed(4); y0 <- rnorm(400)
dw <- as.integer(z[inw] >= 0); yw <- y0[inw]
dm <- function(d) mean(yw[d == 1]) - mean(yw[d == 0])
set.seed(5); perm <- vapply(1:999, function(r) dm(sample(dw)), 0)
p_hand <- (1 + sum(abs(perm) >= abs(dm(dw)))) / (1 + 999)
stopifnot(abs(lr_diffmeans(y0, z, c = 0, w0 = 0.5, reps = 999, seed = 5)$p - p_hand) < 1e-12)
# same seed, same p; the caller's random stream is left where it was
set.seed(9); before <- runif(1); set.seed(9)
lr2 <- lr_diffmeans(yz, z, c = 0, w0 = 0.5, reps = 999, seed = 2)
stopifnot(identical(lr2, lr), identical(runif(1), before))
# a null effect gives a large p
stopifnot(p_hand > 0.05)
# reps must be a positive whole number
stopifnot(inherits(try(lr_diffmeans(yz, z, c = 0, w0 = 0.5, reps = 99.5, seed = 2), silent = TRUE), "try-error"),
          inherits(try(lr_diffmeans(yz, z, c = 0, w0 = 0.5, reps = 0, seed = 2), silent = TRUE), "try-error"))

## ---- Part 2: pooled_cutoff_weights ------------------------------------------------------
# hand case: within h = 2 of their own cutoff are 2 units at 70 and 1 unit at 80 (the
# unit at 70 with x = 73 is outside), so the weights are 2/3 and 1/3; a cutoff with
# no unit in its window gets weight 0
cu <- c(70, 70, 70, 80, 80, 90); xu <- c(69, 71.5, 73, 81, 85, 95)
wu <- pooled_cutoff_weights(cu, xu, h = 2)
stopifnot(identical(names(wu), c("70", "80", "90")), all(abs(wu - c(2, 1, 0) / 3) < 1e-12))
# with many units uniform on [60, 100] and equal group sizes the weights are near 1/3
stopifnot(identical(names(w), c("75", "80", "85")), all(abs(w - 1 / 3) < 0.1))
# no unit within h of its cutoff: refused
stopifnot(inherits(try(pooled_cutoff_weights(c(70, 80), c(50, 50), h = 2), silent = TRUE), "try-error"))

## ---- Part 2: block generators (ruling R4) -------------------------------------------------
# every generator returns the ruling-R2 layout (equal-length column vectors plus meta),
# is reproducible from its seed, and changes with the seed
for (g in c("gen_R1c", "gen_R1h", "gen_R2", "gen_R3", "gen_R4", "gen_R4g")) {
  f <- get(g)
  a1 <- f(n = 2000, seed = 11); a2 <- f(n = 2000, seed = 11); a3 <- f(n = 2000, seed = 12)
  stopifnot(is_block(a1), length(a1$y) > 0, identical(a1, a2), !identical(a1$y, a3$y))
}
# gen_R2 with seed = NULL continues the caller's stream (Task 7 Monte Carlo)
set.seed(21); r_a <- gen_R2(n = 1500, seed = NULL); r_b <- gen_R2(n = 1500, seed = NULL)
stopifnot(!identical(r_a$y, r_b$y))
set.seed(21); stopifnot(identical(gen_R2(n = 1500, seed = NULL), r_a))

cat("test_rd_fns: OK\n")
