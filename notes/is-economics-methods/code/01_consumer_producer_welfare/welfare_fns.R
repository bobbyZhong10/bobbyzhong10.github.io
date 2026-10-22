# ------------------------------------------------------------------
# welfare_fns.R -- functions for chapter 1 (consumer and producer theory, welfare, uncertainty).
#
# Sign convention (notation_v3.md 5.2): CV, EV and dCS are LOSSES, positive for a price rise
# (the opposite of MWG's sign). For a price change p0 -> p1 of good 1 at income w:
#   CV  = e(p1, u0) - w        (compensation at the new prices that restores u0)
#   EV  = w - e(p0, u1)        (income taken away at the old prices that gives u1)
#   dCS = integral from p0 to p1 of the Marshallian demand x1(p, w) dp
# Good 2 is the numeraire (p2 = 1) unless a function says otherwise.
#
# Chapter card interface:
#   cv_ev_exact(e_fun, p0, p1, w)               -> c(cv = , ev = ) from an expenditure function e_fun(p1, u)
#   vartia_cv(x1_fun, p0, p1, w, steps = 2000L) -> CV from Marshallian demand x1_fun(p1, m): fixed-step RK4
#                                                  for dm/dp1 = x1(p1, m), m(p0) = w (Vartia 1983)
#   willig_cv(dcs, eta, w)                      -> dCS + eta dCS^2 / (2 w) (Willig 1976, first order; see below)
#   garp_check(P, X)                            -> list(pass, n_viol, R): direct revealed preference plus
#                                                  Warshall transitive closure
#   ccei(P, X, tol = 1e-8)                      -> Afriat efficiency index by bisection (Afriat 1972; Varian 1990)
#   cara_reject_threshold(loss, gain, prob = 0.5) -> the CARA coefficient at which the bet "lose `loss` with
#                                                  probability prob, win `gain` otherwise" is just rejected
# Helpers used by 01 to 05:
#   dcs_numeric(x1_fun, p0, p1, w, steps)       area under x1_fun(p, w) from p0 to p1 (composite Simpson)
#   sg_demand, sg_indirect, sg_expenditure      Stone-Geary (linear expenditure system) demand for good 1,
#                                               indirect utility (log form) and expenditure function
#   sg_welfare(p0, p1, w, beta, g1, g2)         closed-form CV, EV, dCS and the income elasticity of good 1 at p0
#   ces_share1(p1, p2, alpha, rho)              optimal budget share of good 1 under CES utility (block W2)
#   fit_demand_forms(d)                         the analyst's three demand forms for block W1 (OLS)
#   gen_W1, gen_W2, gen_W2_bronars              the Plover data generators (parameters and seeds in 01_dgp.R)
#
# The Willig correction (Ruling P13; derivation recorded here and in the Task 01.1 report). Let
# L(p) = e(p, u0) - w be the compensating loss along the path of p1. By Shephard's lemma and the duality
# identity h1(p, u0) = x1(p, e(p, u0)), dL/dp = x1(p, w + L): the same ODE vartia_cv integrates. Expanding
# the right side to first order in L, x1(p, w + L) ~ x1(p, w) (1 + eta L / w) with eta the income elasticity
# of good 1 (held fixed along the path), and replacing L inside the correction term by its zeroth-order value
# CS(p) = integral from p0 to p of x1(s, w) ds:
#   CV = L(p1) ~ dCS + (eta / w) integral CS(p) x1(p, w) dp = dCS + eta dCS^2 / (2 w).
# So, under the loss-positive convention, the sign is + and the coefficient 1/2, as the card states; the same
# argument from u1 gives EV ~ dCS - eta dCS^2 / (2 w). Checked numerically in 02_welfare.R against the exact
# CV of block W0 (eta = 1) and of every W1 subscriber (tests/test_conclusions.R).
# Packages: base R (stats) only.
# ------------------------------------------------------------------

## =================================================================
## Card functions
## =================================================================

# Utility level u at which e_fun(p, u) = w (e_fun increasing in u). Two passes: an expanding bracket
# finds the root, a second uniroot on a narrow bracket refines it to rounding level.
.solve_u <- function(e_fun, p, w) {
  f <- function(u) e_fun(p, u) - w
  r <- stats::uniroot(f, c(-1, 1), extendInt = "upX", tol = 1e-10)$root
  h <- 1e-6 * max(1, abs(r))
  lo <- r - h; hi <- r + h
  while (f(lo) > 0) lo <- lo - 2 * h
  while (f(hi) < 0) hi <- hi + 2 * h
  stats::uniroot(f, c(lo, hi), tol = 1e-15 * max(1, abs(r)), maxiter = 500)$root
}

#' Exact CV and EV (loss positive) of a change in p1 from p0 to p1 at income w, from an expenditure
#' function e_fun(p1, u) (good 2's price fixed inside e_fun). u0 and u1 solve e_fun(p0, u0) = w and
#' e_fun(p1, u1) = w; then CV = e_fun(p1, u0) - w and EV = w - e_fun(p0, u1). Scalar w.
cv_ev_exact <- function(e_fun, p0, p1, w) {
  u0 <- .solve_u(e_fun, p0, w)
  u1 <- .solve_u(e_fun, p1, w)
  c(cv = e_fun(p1, u0) - w, ev = w - e_fun(p0, u1))
}

#' Vartia (1983): the exact CV from Marshallian demand alone. The compensated income m(p) = e(p, u0) solves
#' dm/dp1 = x1_fun(p1, m) with m(p0) = w (Shephard's lemma plus duality); CV = m(p1) - w. Classical fourth-order
#' Runge-Kutta with `steps` equal steps; p1 < p0 is allowed (a price fall gives a negative CV, and the path
#' from p1 back to p0 gives -EV). Vectorised over w: x1_fun must accept a vector m.
vartia_cv <- function(x1_fun, p0, p1, w, steps = 2000L) {
  steps <- as.integer(steps)
  if (steps < 1L) stop("vartia_cv: steps must be a positive integer")
  h <- (p1 - p0) / steps
  m <- w
  for (k in seq_len(steps)) {
    p <- p0 + (k - 1L) * h
    k1 <- x1_fun(p, m)
    k2 <- x1_fun(p + h / 2, m + h / 2 * k1)
    k3 <- x1_fun(p + h / 2, m + h / 2 * k2)
    k4 <- x1_fun(p + h, m + h * k3)
    m <- m + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  m - w
}

#' Willig (1976), first order, loss-positive convention: CV ~ dCS + eta dCS^2 / (2 w), with eta the income
#' elasticity of good 1 at the initial point (derivation in the header). willig_cv(dcs, -eta, w) is the
#' matching first-order EV. Vectorised.
willig_cv <- function(dcs, eta, w) dcs + eta * dcs^2 / (2 * w)

# GARP at Afriat efficiency level e (e = 1: GARP itself). Rows of P and X are observations t (prices p^t,
# chosen bundle x^t), columns goods. E[t, s] = p^t . x^s. Direct relation x^t R0 x^s: e p^t.x^t >= p^t.x^s;
# strict x^t P0 x^s: e p^t.x^t > p^t.x^s. R is the transitive closure of R0 (Warshall). A violation is a
# pair with x^t R x^s and x^s P0 x^t.
.garp_e <- function(P, X, e = 1) {
  P <- as.matrix(P); X <- as.matrix(X)
  if (!all(dim(P) == dim(X))) stop("garp_check: P and X must have the same dimensions")
  E <- P %*% t(X)
  own <- diag(E)
  R <- e * own >= E
  Pd <- e * own > E
  for (k in seq_len(nrow(R))) R <- R | outer(R[, k], R[k, ], "&")
  dimnames(R) <- NULL
  n_viol <- sum(R & t(Pd))
  list(pass = n_viol == 0L, n_viol = n_viol, R = R)
}

#' GARP test (Varian 1982): pass is TRUE when no x^t R x^s has x^s strictly directly revealed preferred to x^t;
#' n_viol counts the violating ordered pairs (t, s); R is the revealed-preference relation (transitive closure,
#' logical T x T).
garp_check <- function(P, X) .garp_e(P, X, 1)

#' Critical cost efficiency index (Afriat 1972; Varian 1990): the supremum of the efficiency levels e in [0, 1]
#' at which the data satisfy GARP when every budget is shrunk to e p^t.x^t. Exactly 1 when GARP holds; otherwise
#' bisection until the bracket is narrower than tol, returning its lower end (a level at which GARP holds).
ccei <- function(P, X, tol = 1e-8) {
  if (.garp_e(P, X, 1)$pass) return(1)
  lo <- 0; hi <- 1
  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (.garp_e(P, X, mid)$pass) lo <- mid else hi <- mid
  }
  lo
}

#' The CARA coefficient A > 0 at which a CARA expected-utility maximiser u(x) = -exp(-A x) is indifferent to
#' the bet "lose `loss` with probability prob, win `gain` with probability 1 - prob": the positive root of
#' prob exp(A loss) + (1 - prob) exp(-A gain) = 1. Every A above it rejects the bet, every A below accepts it
#' (wealth plays no role under CARA). Returns 0 when the bet's expected value is not positive (every risk
#' averter rejects it). The log of the left side is convex in A, zero at A = 0 and decreasing there, so the
#' positive root is unique; it is found on the log scale (log-sum-exp) by uniroot.
cara_reject_threshold <- function(loss, gain, prob = 0.5) {
  if (!(loss > 0 && gain > 0 && prob > 0 && prob < 1)) stop("cara_reject_threshold: need loss, gain > 0 and 0 < prob < 1")
  if ((1 - prob) * gain <= prob * loss) return(0)
  g <- function(A) {
    a <- c(log(prob) + A * loss, log(1 - prob) - A * gain)
    mx <- max(a)
    mx + log(sum(exp(a - mx)))
  }
  hi <- 1 / max(loss, gain)
  while (g(hi) <= 0) hi <- 2 * hi
  lo <- hi / 2
  while (g(lo) >= 0) lo <- lo / 2
  stats::uniroot(g, c(lo, hi), tol = 1e-15 * hi, maxiter = 1000)$root
}

## =================================================================
## Helpers
## =================================================================

#' Area under x1_fun(p, w) from p0 to p1 (composite Simpson, `steps` intervals; vectorised over w).
#' For p1 > p0 this is dCS of a price rise, loss positive.
dcs_numeric <- function(x1_fun, p0, p1, w, steps = 2000L) {
  steps <- 2L * ((as.integer(steps) + 1L) %/% 2L)
  h <- (p1 - p0) / steps
  acc <- x1_fun(p0, w) + x1_fun(p1, w)
  for (k in seq_len(steps - 1L)) acc <- acc + (if (k %% 2L == 1L) 4 else 2) * x1_fun(p0 + k * h, w)
  acc * h / 3
}

#' Stone-Geary preferences u = beta log(x1 - g1) + (1 - beta) log(x2 - g2) (g1, g2: subsistence quantities,
#' beta: marginal budget share of good 1). Demand x1 = g1 + beta (w - p1 g1 - p2 g2) / p1. Vectorised.
sg_demand <- function(p1, w, beta, g1, g2, p2 = 1) g1 + beta * (w - p1 * g1 - p2 * g2) / p1

#' Stone-Geary indirect utility (in the log form above) at prices (p1, p2) and income w.
sg_indirect <- function(p1, w, beta, g1, g2, p2 = 1) {
  log(w - p1 * g1 - p2 * g2) - beta * log(p1) - (1 - beta) * log(p2) + beta * log(beta) + (1 - beta) * log(1 - beta)
}

#' Stone-Geary expenditure function e(p, u) = p1 g1 + p2 g2 + exp(u) p1^beta p2^(1 - beta) / (beta^beta (1 - beta)^(1 - beta)).
sg_expenditure <- function(p1, u, beta, g1, g2, p2 = 1) {
  p1 * g1 + p2 * g2 + exp(u) * p1^beta * p2^(1 - beta) / (beta^beta * (1 - beta)^(1 - beta))
}

#' Closed-form welfare of a change in p1 from p0 to p1 for Stone-Geary consumers (vectorised). With y0 and y1
#' the supernumerary incomes w - p g1 - p2 g2 at p0 and p1 and r = p1 / p0:
#'   CV = (p1 - p0) g1 + y0 (r^beta - 1),  EV = y0 - y1 r^(-beta),
#'   dCS = (p1 - p0) g1 (1 - beta) + beta (w - p2 g2) log r,
#'   eta = income elasticity of good 1 at (p0, w) = (beta / p0) w / x1(p0, w).
#' Cobb-Douglas is g1 = g2 = 0.
sg_welfare <- function(p0, p1, w, beta, g1, g2, p2 = 1) {
  y0 <- w - p0 * g1 - p2 * g2
  y1 <- w - p1 * g1 - p2 * g2
  r <- p1 / p0
  list(cv = (p1 - p0) * g1 + y0 * (r^beta - 1),
       ev = y0 - y1 * r^(-beta),
       dcs = (p1 - p0) * g1 * (1 - beta) + beta * (w - p2 * g2) * log(r),
       eta = (beta / p0) * w / sg_demand(p0, w, beta, g1, g2, p2))
}

#' CES utility (alpha x1^rho + (1 - alpha) x2^rho)^(1 / rho), rho < 1: the optimal budget share of good 1,
#' p1 x1 / m = p1 / (p1 + p2 k) with k = x2 / x1 = ((p1 / p2) (1 - alpha) / alpha)^(1 / (1 - rho)). Vectorised.
ces_share1 <- function(p1, p2, alpha, rho) {
  k <- ((p1 / p2) * (1 - alpha) / alpha)^(1 / (1 - rho))
  p1 / (p1 + p2 * k)
}

#' Block W1: the analyst's three demand forms for monthly premium usage x1, fitted by OLS on the experiment
#' data d (columns p1, w, x1). Each form lets the income effect vary with the price, because the CV needs
#' x1(p1, m) along the whole price path:
#'   linear     x1 = a + b p1 + c w + d p1 w                  (truncated at zero above its choke price)
#'   loglinear  log x1 = a + b log p1 + c log w               (constant elasticities; level by Duan's smearing
#'                                                             factor mean(exp(residual)))
#'   income     x1 = t0 + t1 / p1 + t2 w / p1                 (the reciprocal-price form with an income term;
#'                                                             a linear expenditure system with marginal budget
#'                                                             share t2, so a Gorman form)
#' Returns list(linear, loglinear, income); each has coef, r2 (of its own dependent variable), x1_fun(p, m)
#' and eta_fun(p, m) (the income elasticity of x1_fun).
fit_demand_forms <- function(d) {
  r2 <- function(f) {                                            # R^2 of the form's own dependent variable
    y <- stats::model.response(stats::model.frame(f))
    1 - sum(stats::resid(f)^2) / sum((y - mean(y))^2)
  }
  f_lin <- stats::lm(x1 ~ p1 * w, data = d)
  cl <- unname(stats::coef(f_lin))
  lin_raw <- function(p, m) cl[1] + cl[2] * p + cl[3] * m + cl[4] * p * m
  linear <- list(coef = setNames(cl, c("a", "b_p1", "c_w", "d_p1w")), r2 = r2(f_lin),
                 x1_fun = function(p, m) pmax(0, lin_raw(p, m)),
                 eta_fun = function(p, m) (cl[3] + cl[4] * p) * m / lin_raw(p, m))
  f_log <- stats::lm(log(x1) ~ log(p1) + log(w), data = d)
  cg <- unname(stats::coef(f_log))
  smear <- mean(exp(stats::resid(f_log)))
  loglinear <- list(coef = setNames(cg, c("a", "b_logp1", "c_logw")), smear = smear, r2 = r2(f_log),
                    x1_fun = function(p, m) smear * exp(cg[1]) * p^cg[2] * m^cg[3],
                    eta_fun = function(p, m) rep(cg[3], max(length(p), length(m))))
  f_inc <- stats::lm(x1 ~ I(1 / p1) + I(w / p1), data = d)
  ci <- unname(stats::coef(f_inc))
  inc_raw <- function(p, m) ci[1] + (ci[2] + ci[3] * m) / p
  income <- list(coef = setNames(ci, c("t0", "t1_inv_p1", "t2_w_over_p1")), r2 = r2(f_inc),
                 x1_fun = function(p, m) pmax(0, inc_raw(p, m)),
                 eta_fun = function(p, m) (ci[3] * m / p) / inc_raw(p, m))
  list(linear = linear, loglinear = loglinear, income = income)
}

## =================================================================
## Data generators (parameters and seeds in 01_dgp.R)
## =================================================================

#' Block W1: Plover's premium price experiment. N subscribers with lognormal monthly budget w (median
#' w_median USD, sdlog w_sdlog), Stone-Geary preferences with marginal budget share beta_i ~ Beta(beta_phi
#' beta_mean, beta_phi (1 - beta_mean)) (independent of w) and subsistence quantities g1_i ~ U(g1_range),
#' g2_i ~ U(g2_range) (independent of w and beta). The premium unit price p1 is assigned by complete
#' randomisation, N / 5 subscribers at each price of p_grid; observed usage is the subscriber's demand.
#' common_beta = TRUE (variant W1g) sets every beta_i to beta_mean after drawing it, so that with the same seed
#' W1g has the same subscribers, incomes, subsistence quantities and prices as W1: Gorman holds in W1g.
#' Columns ending in _true are simulation truth, not observable data.
gen_W1 <- function(P, common_beta = FALSE) {
  N <- P$N
  w <- exp(log(P$w_median) + P$w_sdlog * stats::rnorm(N))
  beta <- stats::rbeta(N, P$beta_phi * P$beta_mean, P$beta_phi * (1 - P$beta_mean))
  if (common_beta) beta <- rep(P$beta_mean, N)
  g1 <- stats::runif(N, P$g1_range[1], P$g1_range[2])
  g2 <- stats::runif(N, P$g2_range[1], P$g2_range[2])
  p1 <- sample(rep(P$p_grid, each = N / length(P$p_grid)))
  data.frame(subscriber = seq_len(N), p1 = p1, w = w, x1 = sg_demand(p1, w, beta, g1, g2, P$p2),
             beta_true = beta, g1_true = g1, g2_true = g2)
}

#' Block W2: Plover's credit-bundle panel. n_users users each face n_budgets random budgets: a monthly credit
#' of m USD and per-hour prices p1 (premium) and p2 (standard) drawn independently from U(price_range).
#' Users have CES utility with alpha ~ U(alpha_range) and rho ~ U(rho_range). A share share_rational of the
#' users (assigned at random) choose the CES optimum exactly; the others stay on the budget line with
#' multiplicative noise on the odds of the budget share: logit(s1) = logit(s1*) + noise_sd eps, eps ~ N(0, 1).
#' Long format, one row per user and budget.
gen_W2 <- function(P) {
  nu <- P$n_users; Tn <- P$n_budgets
  n_rat <- round(P$share_rational * nu)
  type <- sample(rep(c("rational", "noisy"), c(n_rat, nu - n_rat)))
  alpha <- stats::runif(nu, P$alpha_range[1], P$alpha_range[2])
  rho <- stats::runif(nu, P$rho_range[1], P$rho_range[2])
  p1 <- stats::runif(nu * Tn, P$price_range[1], P$price_range[2])
  p2 <- stats::runif(nu * Tn, P$price_range[1], P$price_range[2])
  eps <- stats::rnorm(nu * Tn)
  user <- rep(seq_len(nu), each = Tn)
  s_opt <- ces_share1(p1, p2, alpha[user], rho[user])
  s <- ifelse(type[user] == "noisy", stats::plogis(stats::qlogis(s_opt) + P$noise_sd * eps), s_opt)
  m <- rep(P$m, nu * Tn)
  data.frame(user = user, t = rep(seq_len(Tn), times = nu), type_true = type[user],
             alpha_true = alpha[user], rho_true = rho[user], p1 = p1, p2 = p2, m = m,
             x1 = s * m / p1, x2 = (1 - s) * m / p2,
             x1_opt_true = s_opt * m / p1, x2_opt_true = (1 - s_opt) * m / p2, stringsAsFactors = FALSE)
}

#' Block W2 power benchmark (Bronars 1987): n_per_user uniform random choosers on each user's own budgets.
#' Each chooser spends a U(0, 1) share of the budget on good 1 in every period (uniform on the budget line).
gen_W2_bronars <- function(W2, n_per_user) {
  rows <- unlist(lapply(split(seq_len(nrow(W2)), W2$user), function(r) rep(r, times = n_per_user)), use.names = FALSE)
  Tn <- max(W2$t)
  d <- W2[rows, c("user", "t", "p1", "p2", "m")]
  d$chooser <- rep(seq_len(length(rows) / Tn), each = Tn)
  s <- stats::runif(nrow(d))
  d$x1 <- s * d$m / d$p1
  d$x2 <- (1 - s) * d$m / d$p2
  rownames(d) <- NULL
  d[c("chooser", "user", "t", "p1", "p2", "m", "x1", "x2")]
}
