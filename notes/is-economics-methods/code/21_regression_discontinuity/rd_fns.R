# ------------------------------------------------------------------
# rd_fns.R -- hand-coded RD helpers and the Aster block generators for chapter 21
#
# Helpers (Task 6 brief interfaces):
#   ll_jump(y, x, c, h, kernel, deriv)        local-linear jump (deriv = 0) or slope change
#                                             (deriv = 1) at c: weighted least squares on each
#                                             side, kernel weights K((x - c) / h)
#   rk_ratio(y, b, x, c, h)                   regression kink ratio: slope change of the outcome
#                                             over slope change of the policy b(x)
#   lr_diffmeans(y, x, c, w0, reps, seed)     local randomization: difference in means inside
#                                             [c - w0, c + w0] and its Fisher permutation p-value
#   pooled_cutoff_weights(cutoff_of_unit, x, h)
#                                             weights the pooled (normalized-score) multi-cutoff
#                                             RD puts on each cutoff's effect: sample analog of
#                                             Cattaneo, Keele, Titiunik and Vazquez-Bare (2016,
#                                             Journal of Politics) Lemma 1
# Generators (ruling R4): gen_R1c, gen_R1h, gen_R2, gen_R3, gen_R4, gen_R4g, each
#   gen_<block>(n, seed, ...) -> list(<equal-length column vectors>, meta = list(...))
#   (ruling R2 layout). 01_dgp.R calls them with the seeds stored in each block's meta;
#   Task 7's Monte Carlo redraws R2 with gen_R2(n, seed = NULL) inside one set.seed().
#   Columns ending in _true are simulation truth, not observable data.
#
# Base R only. Code and comments in English; money in USD.
# ------------------------------------------------------------------

.rd_stop <- function(...) stop(..., call. = FALSE)

# `yname` is how the first vector is called in the error text (rk_ratio checks the policy b)
.rd_check_xy <- function(y, x, what, yname = "y") {
  if (!is.numeric(y) || !is.numeric(x) || length(y) != length(x) || length(x) == 0L)
    .rd_stop(what, ": ", yname, " and x must be numeric vectors of the same positive length")
  if (any(!is.finite(y)) || any(!is.finite(x)))
    .rd_stop(what, ": ", yname, " and x must be finite (no NA, NaN, Inf)")
}

.rd_check_scalar <- function(v, name, what, positive = FALSE) {
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || (positive && v <= 0))
    .rd_stop(what, ": ", name, " must be a finite numeric scalar", if (positive) " > 0" else "")
}

.rd_kernel <- function(u, kernel) {
  switch(kernel,
         triangular   = pmax(1 - abs(u), 0),
         epanechnikov = 0.75 * pmax(1 - u^2, 0),
         uniform      = 0.5 * (abs(u) <= 1),
         .rd_stop("ll_jump: kernel must be one of triangular, epanechnikov, uniform"))
}

#' Local-linear RD at cutoff c with bandwidth h.
#' On each side (left: x < c; right: x >= c, so a unit at the cutoff is treated) fit
#' y = a + b (x - c) by weighted least squares with weights K((x - c) / h); est is
#' a_right - a_left (deriv = 0, the jump in the conditional mean) or b_right - b_left
#' (deriv = 1, the slope change used by regression kink designs). n_left and n_right
#' count the units with |x - c| <= h on each side (rdrobust's N_h convention; with the
#' triangular kernel the units exactly at distance h get weight zero).
ll_jump <- function(y, x, c, h, kernel = "triangular", deriv = 0) {
  .rd_check_xy(y, x, "ll_jump")
  .rd_check_scalar(c, "c", "ll_jump")
  .rd_check_scalar(h, "h", "ll_jump", positive = TRUE)
  if (!is.character(kernel) || length(kernel) != 1L) .rd_stop("ll_jump: kernel must be a single string")
  if (!(identical(deriv, 0) || identical(deriv, 1) || identical(deriv, 0L) || identical(deriv, 1L)))
    .rd_stop("ll_jump: deriv must be 0 (jump) or 1 (slope change)")
  w <- .rd_kernel((x - c) / h, kernel)
  in_h <- abs(x - c) <= h
  right <- x >= c
  fit_side <- function(side, label) {
    keep <- side & w > 0
    if (length(unique(x[keep])) < 2L)
      .rd_stop("ll_jump: fewer than two distinct x values with positive weight on the ", label, " side")
    fit <- stats::lm.wfit(cbind(1, x[keep] - c), y[keep], w[keep])
    as.numeric(fit$coefficients)
  }
  bR <- fit_side(right, "right")
  bL <- fit_side(!right, "left")
  list(est = bR[deriv + 1] - bL[deriv + 1], n_left = sum(!right & in_h), n_right = sum(right & in_h))
}

#' Regression kink ratio at c (Card, Lee, Pei and Weber 2015): the change in the slope of
#' E[y | x] at c divided by the change in the slope of the policy b(x) at c, both from
#' ll_jump(deriv = 1) with the triangular kernel and bandwidth h. When b is the known
#' deterministic schedule (sharp kink) dslope_b is its exact slope change up to rounding.
rk_ratio <- function(y, b, x, c, h) {
  .rd_check_xy(b, x, "rk_ratio", yname = "b")
  dy <- ll_jump(y, x, c, h, deriv = 1)$est
  db <- ll_jump(b, x, c, h, deriv = 1)$est
  if (abs(db) < 1e-12) .rd_stop("rk_ratio: the policy b has no slope change at c (dslope_b is zero)")
  list(est = dy / db, dslope_y = dy, dslope_b = db)
}

#' Local randomization inference (Cattaneo, Frandsen and Titiunik 2015) in the window
#' [c - w0, c + w0] (closed, as in rdlocrand::rdrandinf with wl = c - w0, wr = c + w0):
#' est is the difference in means (x >= c minus x < c); p is the two-sided Fisher
#' randomization p-value of the sharp null of no effect for any unit, computed from
#' `reps` draws of complete randomization with the observed numbers of treated and
#' control units (fixed margins). The observed assignment counts as one draw of the
#' reference distribution, so the Monte Carlo p-value is valid in finite samples and
#' never 0: p = (1 + #{r : |T*_r| >= |T_obs|}) / (1 + reps), at least 1 / (1 + reps).
#' (rdlocrand::rdrandinf reports mean(|T*| >= |T_obs|), which can be 0.)
#' The permutations use set.seed(seed); the caller's random stream is restored on exit.
lr_diffmeans <- function(y, x, c, w0, reps = 2000, seed) {
  .rd_check_xy(y, x, "lr_diffmeans")
  .rd_check_scalar(c, "c", "lr_diffmeans")
  .rd_check_scalar(w0, "w0", "lr_diffmeans", positive = TRUE)
  .rd_check_scalar(reps, "reps", "lr_diffmeans", positive = TRUE)
  if (reps != round(reps)) .rd_stop("lr_diffmeans: reps must be a whole number")
  if (missing(seed)) .rd_stop("lr_diffmeans: seed is required")
  .rd_check_scalar(seed, "seed", "lr_diffmeans")
  in_w <- x >= c - w0 & x <= c + w0
  yw <- y[in_w]
  dw <- as.integer(x[in_w] >= c)
  n1 <- sum(dw); n0 <- sum(1L - dw)
  if (n1 < 1L || n0 < 1L) .rd_stop("lr_diffmeans: the window needs units on both sides of c")
  tot <- sum(yw)
  diffmeans <- function(d) { s1 <- sum(yw[d == 1L]); s1 / n1 - (tot - s1) / n0 }
  est <- diffmeans(dw)

  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = genv)
          else if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv),
          add = TRUE)
  set.seed(seed)
  perm <- vapply(seq_len(reps), function(r) diffmeans(sample(dw)), numeric(1))
  list(est = est, p = (1 + sum(abs(perm) >= abs(est))) / (1 + reps))
}

#' Weights of the pooled multi-cutoff RD estimand on the cutoff-specific effects.
#' Cattaneo, Keele, Titiunik and Vazquez-Bare (2016, Lemma 1): the RD on the normalized
#' score X - C pools to sum_c omega(c) tau(c) with omega(c) proportional to
#' f_{X|C}(c | c) P(C = c). Sample analog at one common bandwidth h: the share of all
#' units within h of their own cutoff that face cutoff c,
#'   omega_hat(c) = #{i : C_i = c, |X_i - c| <= h} / #{i : |X_i - C_i| <= h}.
#' Returns a numeric vector named by cutoff (sorted ascending, names as.character), >= 0,
#' summing to one; a cutoff with no unit in its window gets weight 0.
pooled_cutoff_weights <- function(cutoff_of_unit, x, h) {
  .rd_check_xy(cutoff_of_unit, x, "pooled_cutoff_weights")
  .rd_check_scalar(h, "h", "pooled_cutoff_weights", positive = TRUE)
  in_h <- abs(x - cutoff_of_unit) <= h
  if (!any(in_h)) .rd_stop("pooled_cutoff_weights: no unit lies within h of its own cutoff")
  cuts <- sort(unique(cutoff_of_unit))
  cnt <- vapply(cuts, function(cc) sum(in_h & cutoff_of_unit == cc), numeric(1))
  stats::setNames(cnt / sum(cnt), as.character(cuts))
}

## =================================================================================
## Block generators (ruling R4). Each gen_<block>(n, seed) draws one sample of the
## block's world; seed = NULL draws from the caller's current random stream (used by
## Task 7's Monte Carlo), any number calls set.seed(seed) first. The outcome y is the
## log of monthly GMV in USD throughout (exp(5.6), the R1 level at the cutoff, is about
## USD 270 a month). Block-level design notes are in 01_dgp.R; the numbers are here.
## =================================================================================

.rd_seed <- function(seed) {
  if (!is.null(seed)) {
    .rd_check_scalar(seed, "seed", "gen_*")
    set.seed(seed)
  }
}

.rd_check_n <- function(n, what) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 100 || n != round(n))
    .rd_stop(what, ": n must be a whole number >= 100")
}

# Y(0) mean of the inherited R1 world (01_dgp_rdd.R, unchanged): smooth, steep, wiggly
.r1_mu0 <- function(score) {
  s <- (score - 80) / 10
  5.60 + 0.95 * s + 0.20 * s^2 + 0.50 * sin(4 * s)
}

#' R1c: the R1 world (same score distribution, same Y(0) mean, same 0.20 badge jump at 80)
#' plus two predetermined covariates measured before the badge launched:
#'   gmv_prev  log monthly GMV (USD) in the year before launch = mu0(x) - 0.10 + a + e_prev
#'   tenure    years on Aster, rising mildly with the score (better sellers stay)
#' a ~ N(0, 0.20^2) is a persistent seller demand component that also drives today's GMV,
#' so gmv_prev predicts y; tenure adds repeat buyers (0.04 log points a year). The outcome
#' deviation u = a + 0.04 (tenure - E[tenure | x]) + eps keeps E[y | x] = mu0(x) + 0.20 d,
#' exactly the R1 conditional mean, and sd(u) is about 0.31 (R1: 0.30).
gen_R1c <- function(n = 12000, seed = NULL) {
  .rd_check_n(n, "gen_R1c")
  .rd_seed(seed)
  x <- pmin(pmax(rnorm(n, mean = 74, sd = 11), 45), 100)   # quality score, as in R1
  d <- as.integer(x >= 80)
  a <- rnorm(n, 0, 0.20)                                    # persistent demand component
  g <- exp(0.02 * (x - 74))
  tenure <- round(g * rgamma(n, shape = 2.5, scale = 1.2), 2)   # E[tenure | x] = 3 g
  gmv_prev <- .r1_mu0(x) - 0.10 + a + rnorm(n, 0, 0.12)
  y <- .r1_mu0(x) + 0.20 * d + a + 0.04 * (tenure - 3 * g) + rnorm(n, 0, 0.22)
  list(x = x, d = d, y = y, gmv_prev = gmv_prev, tenure = tenure,
       meta = list(block = "R1c", N = n, seed = seed, cutoff = 80, jump = 0.20,
                   sd_a = 0.20, sd_prev = 0.12, sd_eps = 0.22, tenure_coef = 0.04))
}

#' R1h: the R1 world a year later. Aster now reports the quality score with one decimal
#' and has published the formula; a seller whose true score is a little below 80 can make
#' a one-off push in the scoring week (answer every message within the hour, ship same
#' day) that lifts the reported score to exactly 80.0 and stops there. The push does not
#' change the seller's underlying quality, so y follows the true score:
#'   y = mu0(x_true) + 0.20 d + eps,  d = 1{x >= 80},  x = 80 for heapers, x_true otherwise.
#' The probability of pushing falls linearly from 0.5 at 79.9 to about zero at 76.0
#' (gap = 80 - x_true in (0, 4]); about 3.5% of sellers heap at 80.0.
gen_R1h <- function(n = 12000, seed = NULL) {
  .rd_check_n(n, "gen_R1h")
  .rd_seed(seed)
  x_true <- round(pmin(pmax(rnorm(n, mean = 74, sd = 11), 45), 100), 1)
  gap <- 80 - x_true
  p_push <- ifelse(gap > 0 & gap <= 4, 0.5 * (1 - (gap - 0.1) / 4), 0)
  heaped <- as.integer(runif(n) < p_push)
  x <- ifelse(heaped == 1L, 80, x_true)
  d <- as.integer(x >= 80)
  y <- .r1_mu0(x_true) + 0.20 * d + rnorm(n, 0, 0.30)
  list(x = x, d = d, y = y, x_true = x_true, heaped_true = heaped,
       meta = list(block = "R1h", N = n, seed = seed, cutoff = 80, jump = 0.20,
                   push_max = 0.5, push_width = 4, share_heaped = mean(heaped)))
}

# R2 rating model: latent seller quality q ~ N(0, 1); each review's stars are an ordered
# probit in 0.9 q + 0.3 + e, e ~ N(0, 0.8^2), with thresholds -2.3, -1.7, -1.0, 0.0
# (1|2, 2|3, 3|4, 4|5 stars). The running variable is the average of the 10 most recent
# reviews, on the grid 1.0, 1.1, ..., 5.0.
.r2_par <- list(K = 10L, th = c(-2.3, -1.7, -1.0, 0.0), load = 0.9, shift = 0.3, sd_e = 0.8,
                xmin = 3.0, cutoff = 4.25, base = 5.3, beta_q = 0.04, tau0 = 0.15,
                page_gain = 0.25, page_mid = 0.25, page_scale = 0.06, sd_eps = 0.30)

# the R2 filter effect: entering "4.5 stars & up" (tau0) plus the first-page gain, which
# rises with the underlying average because the filtered list is sorted by it
.r2_tau <- function(x, P = .r2_par) {
  page <- stats::plogis((x - P$cutoff - P$page_mid) / P$page_scale) - stats::plogis(-P$page_mid / P$page_scale)
  ifelse(x >= P$cutoff, P$tau0 + P$page_gain * page, 0)
}

# exact E[q | average = a] and P(average = a) under the rating model (grid quadrature)
.r2_cef_q <- function(P = .r2_par) {
  qg <- seq(-6, 6, by = 0.005)
  cdf <- vapply(P$th, function(t) stats::pnorm((t - P$load * qg - P$shift) / P$sd_e), numeric(length(qg)))
  pr <- cbind(cdf[, 1], cdf[, 2] - cdf[, 1], cdf[, 3] - cdf[, 2], cdf[, 4] - cdf[, 3], 1 - cdf[, 4])
  K <- P$K
  S <- matrix(0, length(qg), 4 * K + 1); S[, 1] <- 1          # column j: sum of (stars - 1) = j - 1
  for (j in seq_len(K)) {
    Sn <- matrix(0, length(qg), 4 * K + 1)
    for (k in 1:5) Sn[, k:(4 * K + 1)] <- Sn[, k:(4 * K + 1)] + S[, 1:(4 * K + 2 - k)] * pr[, k]
    S <- Sn
  }
  w <- stats::dnorm(qg)
  mass <- colSums(S * w)
  data.frame(x = round(1 + (0:(4 * K)) / K, 10), mass = mass / sum(mass), Eq = colSums(S * w * qg) / mass)
}

#' R2: discrete running variable from rating averages (Anderson and Magruder 2012 style).
#' Aster shows the average of a seller's 10 most recent reviews rounded to the nearest half
#' star, and its search has a "4.5 stars & up" filter: a seller is in it once the underlying
#' average reaches 4.25 (d = 1{x >= 4.25}; x = 4.2 shows 4.0 stars, x = 4.3 shows 4.5).
#' Sellers below 3.0 are suspended, so the sample has the 21 support points 3.0, ..., 5.0.
#'   y = 5.3 + 0.04 q + tau(x) + eps,  eps ~ N(0, 0.30^2)
#'   tau(x) = 1{x >= 4.25} [0.15 + 0.25 (L((x - 4.50) / 0.06) - L(-0.25 / 0.06))],  L = plogis
#' 0.15 is the filter effect at the threshold (the RD jump, truth21$R2$jump). The filtered
#' list is sorted by underlying average, so from about 4.4 up sellers reach the first page,
#' a further gain of up to 0.25 that levels off once on it. Quality q moves GMV only a
#' little (0.04 per sd) beyond the displayed stars. E[y(0) | x] = 5.3 + 0.04 E[q | x] is
#' smooth and nearly flat between 4.2 and 4.3; the curvature sits just right of the cutoff,
#' where local-linear extrapolation from 4.3, 4.4, 4.5 undershoots the jump.
#' tenure (years on Aster, 3 + 1.2 q + noise, floored at 0.2) is predetermined: it rises
#' with the underlying average through q, so covariate balance holds only in a narrow
#' window around 4.25 (what rdlocrand::rdwinselect looks for).
#' meta$M_bound = 0.04 * max|second difference of E[q | x]| / 0.1^2 + max|tau''| (about
#' 7.79; the first term, 1.11, comes from the top of the scale): a bound on |f''| of the
#' conditional mean on either side, the conservative M for honest (Armstrong-Kolesar) CIs.
#' meta$cef holds E[y(0) | x] (m0), tau(x) and the population mass at each support point.
gen_R2 <- function(n = 8000, seed = NULL) {
  .rd_check_n(n, "gen_R2")
  .rd_seed(seed)
  P <- .r2_par
  x <- numeric(0); q <- numeric(0)
  while (length(x) < n) {                       # draw sellers until n pass the 3.0 floor
    m <- n
    qq <- rnorm(m)
    stars <- 0
    for (j in seq_len(P$K)) {
      l <- P$load * qq + P$shift + rnorm(m, 0, P$sd_e)
      stars <- stars + 1L + (l > P$th[1]) + (l > P$th[2]) + (l > P$th[3]) + (l > P$th[4])
    }
    avg <- round(stars / P$K, 10)
    keep <- avg >= P$xmin
    x <- c(x, avg[keep]); q <- c(q, qq[keep])
  }
  x <- x[1:n]; q <- q[1:n]
  d <- as.integer(x >= P$cutoff)
  y <- P$base + P$beta_q * q + .r2_tau(x, P) + rnorm(n, 0, P$sd_eps)
  tenure <- round(pmax(0.2, 3 + 1.2 * q + rnorm(n, 0, 1.5)), 1)   # predetermined: years on Aster
  cef <- .r2_cef_q(P)
  cef <- cef[cef$x >= P$xmin, ]
  step <- 1 / P$K
  d2_q <- max(abs(diff(cef$Eq, differences = 2))) / step^2
  tau2 <- P$page_gain / (6 * sqrt(3)) / P$page_scale^2           # max |plogis''| = 1 / (6 sqrt 3)
  cef$m0 <- P$base + P$beta_q * cef$Eq                           # E[y(0) | x]
  cef$tau <- .r2_tau(cef$x, P)
  cef$mass <- cef$mass / sum(cef$mass)
  list(x = x, d = d, y = y, tenure = tenure, q_true = q,
       meta = list(block = "R2", N = n, seed = seed, cutoff = P$cutoff, jump = P$tau0,
                   support = cef$x, cef = cef[, c("x", "mass", "m0", "tau")],
                   M_bound = P$beta_q * d2_q + tau2, par = P))
}

#' R3: regression kink in Aster's commission schedule. The commission rate a seller pays
#' this quarter is set by last quarter's average monthly GMV x (USD, predetermined): 12%
#' up to USD 1,500, then a volume discount of 0.6 percentage points per USD 100 above it,
#' down to a 6% floor from USD 2,500 (a second kink, 1,000 USD away):
#'   b(x) = 12 - 0.006 * min(max(x - 1500, 0), 1000)    (percent; continuous, kinked at 1,500)
#' Outcome y = log monthly GMV this quarter = 0.85 log(x) + 1.1 + kink * (b(x) - 12) + eps,
#' eps ~ N(0, 0.15^2), kink = -0.05: one percentage point of commission lowers GMV by 5%
#' (sellers pass part of the fee into prices). E[y | x] is continuous at 1,500 (no level
#' jump) and its slope rises by kink * (-0.006) = 0.0003 per USD there. x is lognormal
#' around USD 1,400 and restricted to sellers with USD 500 to 3,000 a month.
#' meta stores every design constant the draw uses (the number registry reads them there).
gen_R3 <- function(n = 60000, seed = NULL) {
  .rd_check_n(n, "gen_R3")
  .rd_seed(seed)
  k <- 1500; s_r <- 0.006; r0 <- 12; w_disc <- 1000; kink <- -0.05
  x_median <- 1400; sdlog <- 0.40; x_lo <- 500; x_hi <- 3000      # running variable
  elast <- 0.85; const <- 1.1; sd_eps <- 0.15                     # outcome equation
  b_fun <- eval(bquote(function(x) .(r0) - .(s_r) * pmin(pmax(x - .(k), 0), .(w_disc))), baseenv())
  x <- numeric(0)
  while (length(x) < n) {
    z <- round(exp(rnorm(n, log(x_median), sdlog)), 2)
    x <- c(x, z[z >= x_lo & z <= x_hi])
  }
  x <- x[1:n]
  b <- b_fun(x)
  y <- elast * log(x) + const + kink * (b - r0) + rnorm(n, 0, sd_eps)
  list(x = x, b = b, y = y,
       meta = list(block = "R3", N = n, seed = seed, cutoff = k, b_fun = b_fun,
                   slope_left = 0, slope_right = -s_r, kink_effect = kink, floor_from = k + w_disc,
                   level_jump = 0, rate_base = r0, rate_floor = r0 - s_r * w_disc,
                   x_median = x_median, sdlog = sdlog, x_lo = x_lo, x_hi = x_hi,
                   elasticity_x = elast, intercept = const, sd_eps = sd_eps))
}

# R4 categories: cutoff, share of sellers, score distribution (normal, clamped to [45, 100]
# with pmin / pmax, so the few draws outside pile up at 45 and 100 and the density at the
# interior cutoffs is the plain normal density), category level of Y(0) in log points
.r4_cats <- data.frame(category = c("vintage", "home", "electronics"), cutoff = c(70, 80, 90),
                       share = c(8, 5, 2) / 15, mu = c(70, 77, 85), sd = c(9, 10, 9),
                       alpha = c(-0.25, 0, 0.30), stringsAsFactors = FALSE)
.r4_m0 <- function(x) 5.2 + 0.06 * (x - 75) + 0.0008 * (x - 75)^2      # common Y(0) shape
.r4_tau <- function(x) 0.40 - 0.015 * (x - 70)                         # badge effect at score x

#' R4: category-specific badge cutoffs (multi-cutoff RD; Cattaneo, Keele, Titiunik and
#' Vazquez-Bare 2016, 2021). Aster awards the badge at 70 in vintage (slow-shipping norms
#' pull scores down), 80 in home & living and 90 in electronics (buyers expect fast,
#' defect-free shipping). For a seller with score x in category g,
#'   y = m0(x) + alpha_g + d * tau(x) + eps,  d = 1{x >= cutoff_g},  eps ~ N(0, 0.30^2)
#'   m0(x) = 5.2 + 0.06 (x - 75) + 0.0008 (x - 75)^2,  tau(x) = 0.40 - 0.015 (x - 70)
#' The badge is worth less at a higher score (buyers already see strong signals), so the
#' cutoff effects are 0.40, 0.25 and 0.10. Y(0) differs across categories only by the
#' level alpha_g (-0.25, 0, +0.30): the "parallel" Y(0) curves of CKTV (2021) hold on
#' [70, 80], so the vintage effect at the interior score 75 (truth 0.325) is identified
#' from the home sellers below 80. Vintage is the biggest category and has most mass at
#' its cutoff, so pooling over cutoffs weights 0.40 most (Lemma 1 weights 0.571, 0.307,
#' 0.122; pooled target 0.317 against a simple mean of 0.25).
gen_R4 <- function(n = 24000, seed = NULL) {
  .rd_check_n(n, "gen_R4")
  .rd_seed(seed)
  G <- .r4_cats
  ng <- round(n * G$share); ng[length(ng)] <- n - sum(ng[-length(ng)])
  parts <- lapply(seq_len(nrow(G)), function(k) {
    x <- pmin(pmax(rnorm(ng[k], G$mu[k], G$sd[k]), 45), 100)
    d <- as.integer(x >= G$cutoff[k])
    y <- .r4_m0(x) + G$alpha[k] + d * .r4_tau(x) + rnorm(ng[k], 0, 0.30)
    list(x = x, cutoff = rep(G$cutoff[k], ng[k]), category = rep(G$category[k], ng[k]), d = d, y = y)
  })
  col <- function(nm) unlist(lapply(parts, `[[`, nm), use.names = FALSE)
  # population Lemma 1 weights: f_{X|C}(c | c) P(C = c); the clamping at 45 and 100 leaves
  # the density at the interior cutoffs equal to dnorm (no truncation rescaling)
  f_cc <- stats::dnorm(G$cutoff, G$mu, G$sd)
  omega <- stats::setNames(f_cc * G$share / sum(f_cc * G$share), as.character(G$cutoff))
  list(x = col("x"), cutoff = col("cutoff"), category = col("category"), d = col("d"), y = col("y"),
       meta = list(block = "R4", N = n, seed = seed, cutoffs = G$cutoff, categories = G,
                   tau_c = stats::setNames(.r4_tau(G$cutoff), as.character(G$cutoff)),
                   x_extrap = 75, extrap_cutoff = 70, ref_cutoff = 80, tau_extrap = .r4_tau(75),
                   omega = omega))
}

#' R4g: geographic RD at the edge of Aster's same-day delivery zone. Couriers can run a
#' same-day loop for sellers whose stock sits within 25 km of the metro fulfillment hub
#' (at the origin); those sellers get the same-day badge in search (d = 1). Coordinates
#' x1, x2 are km east and north of the hub. 45% of sellers are in the urban core (distance
#' half-normal with scale 14 km), 55% in a suburban ring (distance uniform on 12 to 40 km),
#' directions uniform. With r = distance to the hub,
#'   y = 5.4 + 0.35 exp(-r / 20) + 0.004 x2 + d * tau(x1, x2) + eps,  eps ~ N(0, 0.30^2)
#'   tau(x1, x2) = 0.20 + 0.10 * x1 / 25
#' (local demand falls with distance from the center; buyers on the east side value
#' same-day delivery more, so the effect rises with x1; seller directions are uniform, so
#' seller density does not differ by side). The six boundary points at 0, 60, ..., 300
#' degrees have effects 0.30, 0.25, 0.15, 0.10, 0.15, 0.25 (truth21$R4g$tau_b).
#' meta stores every design constant the draw uses (the number registry reads them there).
gen_R4g <- function(n = 40000, seed = NULL) {
  .rd_check_n(n, "gen_R4g")
  .rd_seed(seed)
  R0 <- 25
  share_urban <- 0.45; urban_scale <- 14; suburb_lo <- 12; suburb_hi <- 40        # seller locations
  base <- 5.4; demand_amp <- 0.35; demand_scale <- 20; north <- 0.004; sd_eps <- 0.30  # Y(0)
  tau_center <- 0.20; tau_east <- 0.10                                            # effect
  nu <- round(share_urban * n)
  r <- c(abs(rnorm(nu, 0, urban_scale)), runif(n - nu, suburb_lo, suburb_hi))
  th <- runif(n, 0, 2 * pi)
  x1 <- r * cos(th); x2 <- r * sin(th)
  rr <- sqrt(x1^2 + x2^2)
  d <- as.integer(rr <= R0)
  tau <- tau_center + tau_east * x1 / R0
  y <- base + demand_amp * exp(-rr / demand_scale) + north * x2 + d * tau + rnorm(n, 0, sd_eps)
  ang <- seq(0, 300, by = 60)
  bp <- data.frame(angle_deg = ang, b1 = R0 * cos(ang * pi / 180), b2 = R0 * sin(ang * pi / 180))
  list(x1 = x1, x2 = x2, d = d, y = y,
       meta = list(block = "R4g", N = n, seed = seed, hub = c(0, 0), radius = R0, boundary = bp,
                   tau_b = tau_center + tau_east * bp$b1 / R0,
                   share_urban = share_urban, urban_scale = urban_scale, suburb_lo = suburb_lo,
                   suburb_hi = suburb_hi, base = base, demand_amp = demand_amp, demand_scale = demand_scale,
                   north = north, sd_eps = sd_eps, tau_center = tau_center, tau_east = tau_east))
}
