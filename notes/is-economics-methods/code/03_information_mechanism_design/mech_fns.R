# ------------------------------------------------------------------
# mech_fns.R -- functions for chapter 3 (information economics, contracts, mechanism design, auctions,
# matching). Sourced by every script and test of this directory.
#
# Models, in the text's notation (notation_v3.md 5.4):
#   Auctions (assumptions I.1, I.2). n bidders with values v_i iid with distribution F and density f > 0 on
#     [lower, upper]. Virtual value Psi(v) = v - (1 - F(v)) / f(v); optimal reserve r* solves Psi(r*) = 0.
#     First-price symmetric equilibrium with reserve r (r = 0: no reserve):
#       b(v) = v - int_r^v F(z)^(n-1) dz / F(v)^(n-1) = E[ max(v_{-i,(1)}, r) | v_{-i,(1)} < v ],  v >= r.
#     GPV in quantile form: F^{-1}(tau) = F_b^{-1}(tau) + tau / (n - 1) * d/dtau F_b^{-1}(tau), with the
#     empirical quantiles of bids and a central finite difference of half-width h (no kernel; chapter 14).
#   Screening (Mussa and Rosen 1978). Types theta with density f on [theta_lo, theta_hi]; utility
#     theta q - t; cost cost2 q^2 / 2. Optimal q(theta) = max(0, Psi_bar(theta)) / cost2, where Psi_bar is
#     Psi ironed: the isotonic regression of Psi weighted by the type mass (the slope of the convex hull of
#     int Psi dF in quantile space), computed by pool adjacent violators. Rent U(theta) = int q (envelope,
#     U(theta_lo) = 0), payment t = theta q - U. Two types: share mu of high types (MWG 14.C style).
#   LEN contract (Holmstrom and Milgrom 1987). Output x = e + eps, eps ~ N(0, sigma2); agent CARA with
#     coefficient A (len_bstar's r), cost c(e) = c2 e^2 / 2; linear pay w = w0 + beta x; beta* = 1 / (1 + A sigma2 c2).
#     With signals z ~ N(mu e, Sigma) the best linear pay loads on the sufficient statistic with variance
#     sigma_s^2 = 1 / (mu' Sigma^{-1} mu): weights b = beta* sigma_s^2 Sigma^{-1} mu, beta* = 1 / (1 + A sigma_s^2 c2).
#   Matching (MT.). One-to-one, strict preferences given as lists (most preferred first; agents not listed
#     are unacceptable). deferred_acceptance returns, for each proposer, the receiver it is matched to (NA
#     if unmatched): the proposer-optimal stable matching (Gale and Shapley 1962).
#
# Card functions: fpa_bid, gpv_quantile, optimal_menu, len_bstar, deferred_acceptance, is_stable.
# Money is in US dollars (USD per click in the auction blocks).
# ------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

## =================================================================
## Distributions
## =================================================================

#' Value distribution from a parameter list: family "uniform" (lower, upper) or "truncnorm" (mean, sd of the
#' untruncated normal, lower, upper). Returns list(params, lower, upper, cdf, pdf, quantile).
value_dist <- function(par) {
  fam <- par$family
  if (identical(fam, "uniform")) {
    lo <- par$lower; hi <- par$upper
    cdf <- function(v) pmin(1, pmax(0, (v - lo) / (hi - lo)))
    pdf <- function(v) ifelse(v < lo | v > hi, 0, 1 / (hi - lo))
    quantile <- function(tau) lo + tau * (hi - lo)
  } else if (identical(fam, "truncnorm")) {
    m <- par$mean; s <- par$sd; lo <- par$lower; hi <- par$upper
    pa <- stats::pnorm((lo - m) / s); Z <- stats::pnorm((hi - m) / s) - pa
    cdf <- function(v) pmin(1, pmax(0, (stats::pnorm((v - m) / s) - pa) / Z))
    pdf <- function(v) ifelse(v < lo | v > hi, 0, stats::dnorm((v - m) / s) / (s * Z))
    quantile <- function(tau) pmin(hi, pmax(lo, m + s * stats::qnorm(pa + tau * Z)))
  } else stop("value_dist: unknown family ", fam)
  list(params = par, lower = lo, upper = hi, cdf = cdf, pdf = pdf, quantile = quantile)
}

#' Type density (not normalised) on a grid theta, for block V2: family "truncnorm" (mean, sd) or
#' "truncnorm_mix" (weights, means, sds), both truncated to [lower, upper].
type_density <- function(par, theta) {
  inside <- theta >= par$lower & theta <= par$upper
  d <- if (identical(par$family, "truncnorm")) {
    stats::dnorm(theta, par$mean, par$sd)
  } else if (identical(par$family, "truncnorm_mix")) {
    Reduce(`+`, Map(function(w, m, s) w * stats::dnorm(theta, m, s), par$weights, par$means, par$sds))
  } else stop("type_density: unknown family ", par$family)
  ifelse(inside, d, 0)
}

#' Virtual value Psi(v) = v - (1 - F(v)) / f(v).
virtual_value <- function(v, dist) v - (1 - dist$cdf(v)) / dist$pdf(v)

#' Optimal reserve r*: the root of Psi on the support (the lower end when Psi >= 0 there).
optimal_reserve <- function(dist) {
  lo <- dist$lower; hi <- dist$upper
  if (virtual_value(lo, dist) >= 0) return(lo)
  stats::uniroot(function(v) virtual_value(v, dist), c(lo, hi), tol = 1e-14)$root
}

## =================================================================
## Auctions
## =================================================================

#' First-price symmetric equilibrium bid b(v) = v - int_r^v F^(n-1) / F(v)^(n-1) (card function).
#' v: values; n >= 2 bidders; Fcdf: the value CDF; reserve: r (0 = no reserve). NA for v < r (no bid);
#' b(r) = r; where F(v) = 0 the bid is v.
fpa_bid <- function(v, n, Fcdf, reserve = 0) {
  stopifnot(n >= 2)
  vapply(v, function(x) {
    if (is.na(x) || x < reserve) return(NA_real_)
    if (x == reserve) return(reserve)
    Fx <- Fcdf(x)^(n - 1)
    if (Fx <= 0) return(x)
    I <- stats::integrate(function(z) Fcdf(z)^(n - 1), lower = reserve, upper = x,
                          rel.tol = 1e-11, abs.tol = 0, subdivisions = 1000L)$value
    x - I / Fx
  }, numeric(1))
}

#' A fast version of fpa_bid for the Monte Carlo: the exact bid on a grid of n_grid points over
#' [max(lower, reserve), upper] (the integral accumulated piece by piece), interpolated by a monotone
#' cubic spline. Returns a function of v (NA below the reserve).
fpa_bid_fun <- function(n, dist, reserve = 0, n_grid = 2001L) {
  lo <- max(dist$lower, reserve); hi <- dist$upper
  g <- seq(lo, hi, length.out = n_grid)
  Fn <- function(z) dist$cdf(z)^(n - 1)
  pieces <- vapply(seq_len(n_grid - 1L), function(k)
    stats::integrate(Fn, g[k], g[k + 1L], rel.tol = 1e-12, abs.tol = 0)$value, numeric(1))
  I <- c(0, cumsum(pieces))                          # F = 0 below the support: nothing to add there
  Fg <- Fn(g)
  b <- g - I / ifelse(Fg > 0, Fg, 1)                 # where F(v) = 0 the bid is v (then I = 0)
  sf <- stats::splinefun(g, b, method = "monoH.FC")
  function(v) {
    out <- rep(NA_real_, length(v))
    ok <- !is.na(v) & v >= reserve
    out[ok] <- sf(pmin(pmax(v[ok], lo), hi))
    out
  }
}

#' Expected revenue per auction with n bidders: E[v_(2)] without a reserve (reserve <= lower), and
#' n int_r^upper Psi(v) F(v)^(n-1) f(v) dv with a reserve r (Myerson's revenue formula).
expected_revenue <- function(n, dist, reserve = 0) {
  F <- dist$cdf; f <- dist$pdf
  if (reserve <= dist$lower)
    stats::integrate(function(v) v * n * (n - 1) * F(v)^(n - 2) * (1 - F(v)) * f(v),
                     dist$lower, dist$upper, rel.tol = 1e-12)$value
  else
    stats::integrate(function(v) virtual_value(v, dist) * n * F(v)^(n - 1) * f(v),
                     reserve, dist$upper, rel.tol = 1e-12)$value
}

#' Revenue of each auction (rows of the value matrix V, one column per bidder) under four formats: second
#' price (truthful bids), first price (bid_fp, the equilibrium bid function), and both with reserve r
#' (bid_fp_res). Unsold: revenue 0. Returns a data.table with v1, v2 and the four revenues.
auction_outcomes <- function(V, reserve, bid_fp, bid_fp_res) {
  V <- as.matrix(V)
  cols <- lapply(seq_len(ncol(V)), function(j) V[, j])
  v1 <- do.call(pmax, cols)
  v2 <- do.call(pmax, lapply(cols, function(x) ifelse(x < v1, x, -Inf)))
  sold <- v1 >= reserve
  fp_res <- rep(0, length(v1)); fp_res[sold] <- bid_fp_res(v1[sold])
  data.table(v1 = v1, v2 = v2, rev_sp = v2, rev_fp = bid_fp(v1),
             rev_sp_res = ifelse(sold, pmax(v2, reserve), 0), rev_fp_res = fp_res)
}

#' GPV inversion in quantile form (card function): v_hat(tau) = b(tau) + tau b'(tau) / (n - 1), where b(tau)
#' is the empirical tau-quantile of the bids (type 7) and b'(tau) a finite difference over
#' [tau - h, tau + h] cut to [0, 1]: central when h <= tau <= 1 - h, off-centre (shorter on the boundary
#' side) within h of 0 or 1. Returns data.table(tau, b, v_hat).
gpv_quantile <- function(bids, n, taus, h = 0.02) {
  q <- function(p) stats::quantile(bids, p, type = 7, names = FALSE)
  lo <- pmax(taus - h, 0); hi <- pmin(taus + h, 1)
  b <- q(taus)
  slope <- (q(hi) - q(lo)) / (hi - lo)
  data.table(tau = taus, b = b, v_hat = b + taus * slope / (n - 1))
}

## =================================================================
## Screening
## =================================================================

#' Weighted isotonic (non-decreasing) regression by pool adjacent violators. Returns the fitted vector with
#' attribute "block": the id of the pooled block of each element.
pava <- function(y, w) {
  n <- length(y); val <- numeric(n); wt <- numeric(n); len <- integer(n); k <- 0L
  for (i in seq_len(n)) {
    k <- k + 1L; val[k] <- y[i]; wt[k] <- w[i]; len[k] <- 1L
    while (k > 1L && val[k - 1L] > val[k]) {
      nw <- wt[k - 1L] + wt[k]
      val[k - 1L] <- (wt[k - 1L] * val[k - 1L] + wt[k] * val[k]) / nw
      wt[k - 1L] <- nw; len[k - 1L] <- len[k - 1L] + len[k]; k <- k - 1L
    }
  }
  out <- rep(val[seq_len(k)], len[seq_len(k)])
  attr(out, "block") <- rep(seq_len(k), len[seq_len(k)])
  out
}

#' Optimal nonlinear menu for types on the grid theta with density values dens (card function). F is the
#' cumulative trapezoid of dens, normalised to F(theta_hi) = 1; w are the trapezoid masses of the grid
#' points (they sum to 1). iron = TRUE irons Psi (pava); iron = FALSE uses Psi itself (the menu is then not
#' incentive compatible where Psi falls). Returns data.table(theta, f, F, w, psi, psi_bar, q, q_eff, U, t,
#' block, bunched).
optimal_menu <- function(theta, dens, cost2, iron = TRUE) {
  K <- length(theta); d <- diff(theta)
  stopifnot(K >= 3L, all(d > 0), all(dens > 0))
  cum <- c(0, cumsum((dens[-1] + dens[-K]) / 2 * d)); tot <- cum[K]
  Fth <- cum / tot; f <- dens / tot
  w <- f * c(d[1], d[-1] + d[-(K - 1L)], d[K - 1L]) / 2
  psi <- theta - (1 - Fth) / f
  if (iron) {
    pb <- pava(psi, w); block <- attr(pb, "block"); psi_bar <- as.numeric(pb)
  } else {
    psi_bar <- psi; block <- seq_len(K)
  }
  q <- pmax(psi_bar, 0) / cost2
  U <- c(0, cumsum((q[-1] + q[-K]) / 2 * d))          # envelope: U(theta) = int q, U(theta_lo) = 0
  size <- tabulate(block)[block]
  data.table(theta = theta, f = f, F = Fth, w = w, psi = psi, psi_bar = psi_bar, q = q, q_eff = theta / cost2,
             U = U, t = theta * q - U, block = block, bunched = size > 1L)
}

#' Largest gain any type on the menu's grid gets by choosing another type's item (0 when the menu is
#' incentive compatible on the grid).
menu_ic_gain <- function(menu) {
  pay <- outer(menu$theta, menu$q) - matrix(menu$t, nrow(menu), nrow(menu), byrow = TRUE)
  max(apply(pay, 1, max) - (menu$theta * menu$q - menu$t))
}

#' Summary of a menu: expected profit, surplus, rent, the efficient surplus, exclusion, and the widest
#' bunch with positive quantity (an ironed block of more than one grid point).
menu_summary <- function(menu, cost2) {
  w <- menu$w
  surplus <- sum(w * (menu$theta * menu$q - cost2 * menu$q^2 / 2))
  rent <- sum(w * menu$U)
  profit <- sum(w * (menu$t - cost2 * menu$q^2 / 2))
  surplus_eff <- sum(w * menu$theta^2 / (2 * cost2))
  served <- menu$q > 0
  bl <- unique(menu$block[menu$bunched & served])
  if (length(bl)) {
    wd <- vapply(bl, function(b) diff(range(menu$theta[menu$block == b])), numeric(1))
    b0 <- bl[which.max(wd)]; ib <- menu$block == b0
    bunch <- list(n_bunches = length(bl), bunch_lo = min(menu$theta[ib]), bunch_hi = max(menu$theta[ib]),
                  bunch_width = max(wd), bunch_q = menu$q[ib][1], bunch_mass = sum(w[ib]),
                  bunch_psi_bar = menu$psi_bar[ib][1])
  } else {
    bunch <- list(n_bunches = 0L, bunch_lo = NA_real_, bunch_hi = NA_real_, bunch_width = 0, bunch_q = 0,
                  bunch_mass = 0, bunch_psi_bar = NA_real_)
  }
  c(list(profit = profit, surplus = surplus, rent = rent, surplus_efficient = surplus_eff,
         dwl = surplus_eff - surplus, theta_exclusion = if (any(served)) min(menu$theta[served]) else NA_real_,
         share_served = sum(w[served]), q_top = tail(menu$q, 1), q_eff_top = tail(menu$q_eff, 1),
         ic_gain = menu_ic_gain(menu)), bunch)
}

#' Two-type screening: utility theta q - t, cost cost2 q^2 / 2, share mu of high types. IR of the low type
#' and IC of the high type bind; the high type gets the efficient q, the low type's q is pushed down
#' (to 0 when mu / (1 - mu) (theta_H - theta_L) >= theta_L). Also the first best (types observed).
screening_two_type <- function(theta_L, theta_H, mu, cost2) {
  q_H <- theta_H / cost2
  q_L <- max(0, theta_L - mu / (1 - mu) * (theta_H - theta_L)) / cost2
  t_L <- theta_L * q_L; rent_H <- (theta_H - theta_L) * q_L; t_H <- theta_H * q_H - rent_H
  q_L_fb <- theta_L / cost2
  prof <- function(qL, tL, qH, tH) (1 - mu) * (tL - cost2 * qL^2 / 2) + mu * (tH - cost2 * qH^2 / 2)
  list(q_L = q_L, q_H = q_H, t_L = t_L, t_H = t_H, rent_H = rent_H,
       profit_L = t_L - cost2 * q_L^2 / 2, profit_H = t_H - cost2 * q_H^2 / 2, profit = prof(q_L, t_L, q_H, t_H),
       q_L_fb = q_L_fb, q_H_fb = q_H, t_L_fb = theta_L * q_L_fb, t_H_fb = theta_H * q_H,
       profit_fb = prof(q_L_fb, theta_L * q_L_fb, q_H, theta_H * q_H))
}

## =================================================================
## Moral hazard: the LEN contract
## =================================================================

#' Optimal slope of the linear contract (card function): beta* = 1 / (1 + r sigma2 c2), r = CARA coefficient
#' (A in the text), sigma2 = output noise variance, c2 = c''.
len_bstar <- function(r, sigma2, c2) 1 / (1 + r * sigma2 * c2)

#' Optimal linear contract on several signals z ~ N(mu e, Sigma) (the first is output x = e + eps). Returns
#' the variance of the sufficient statistic sigma_s2, the incentive power beta (= b' mu), the weights b, the
#' effort beta / c2, the risk premium A b' Sigma b / 2, the effort cost, the principal's expected profit
#' e - c(e) - risk premium (the agent's certainty equivalent held at 0 by the fixed salary w0), the sd of pay.
len_contract <- function(A, c2, Sigma, mu) {
  Sigma <- as.matrix(Sigma); mu <- as.numeric(mu)
  Si_mu <- solve(Sigma, mu)
  sigma_s2 <- 1 / sum(mu * Si_mu)
  beta <- len_bstar(A, sigma_s2, c2)
  b <- beta * sigma_s2 * Si_mu
  effort <- beta / c2
  rp <- A * sum(b * (Sigma %*% b)) / 2
  cost <- c2 * effort^2 / 2
  list(sigma_s2 = sigma_s2, beta = beta, b = b, effort = effort, risk_premium = rp, effort_cost = cost,
       profit = effort - cost - rp, pay_sd = sqrt(sum(b * (Sigma %*% b))), w0 = cost + rp - beta * effort)
}

## =================================================================
## Matching
## =================================================================

# rank[i, j] = position of j on i's list (Inf: j is unacceptable to i)
.rank_matrix <- function(prefs, n_other) {
  R <- matrix(Inf, length(prefs), n_other)
  for (i in seq_along(prefs)) { l <- prefs[[i]]; if (length(l)) R[i, l] <- seq_along(l) }
  R
}

#' Turn a matching from one side's view into the other side's: match[i] = partner of i (NA if single).
invert_match <- function(match, n_other) {
  out <- rep(NA_integer_, n_other)
  ok <- !is.na(match)
  out[match[ok]] <- which(ok)
  out
}

#' Deferred acceptance with pref_prop proposing (card function). Lists are vectors of partner indices, most
#' preferred first; unlisted partners are unacceptable. Proposers propose one at a time (McVitie-Wilson
#' order; the outcome does not depend on the order). Returns an integer vector: the receiver of each
#' proposer, NA if unmatched.
deferred_acceptance <- function(pref_prop, pref_recv) {
  np <- length(pref_prop); nr <- length(pref_recv)
  pref_prop <- lapply(pref_prop, as.integer)
  rank_r <- .rank_matrix(pref_recv, np)
  nxt <- rep(1L, np); held <- rep(NA_integer_, nr)
  for (p0 in seq_len(np)) {
    cur <- p0
    while (!is.na(cur)) {
      l <- pref_prop[[cur]]
      if (nxt[cur] > length(l)) break                       # list exhausted: cur stays single
      r <- l[nxt[cur]]; nxt[cur] <- nxt[cur] + 1L
      rk <- rank_r[r, cur]
      if (is.infinite(rk)) next                             # r finds cur unacceptable
      h <- held[r]
      if (is.na(h)) {
        held[r] <- cur; cur <- NA_integer_
      } else if (rk < rank_r[r, h]) {
        held[r] <- cur; cur <- h                            # r trades up; h is free again
      }
    }
  }
  invert_match(held, np)
}

#' Is the matching (match[i] = receiver of proposer i, NA if single) stable for the lists? Checks that it is
#' one to one, individually rational (matched partners find each other acceptable) and has no blocking pair
#' (a proposer and a receiver who are acceptable to each other and each prefer the other to its partner,
#' being single counting as worse than any acceptable partner) (card function).
is_stable <- function(match, pref_prop, pref_recv) {
  np <- length(pref_prop); nr <- length(pref_recv)
  match <- as.integer(match)
  if (length(match) != np) return(FALSE)
  mm <- match[!is.na(match)]
  if (any(mm < 1L | mm > nr) || anyDuplicated(mm)) return(FALSE)
  rp <- .rank_matrix(pref_prop, nr); rr <- .rank_matrix(pref_recv, np)
  partner <- invert_match(match, nr)
  cur_p <- rep(Inf, np); i <- which(!is.na(match)); cur_p[i] <- rp[cbind(i, match[i])]
  cur_r <- rep(Inf, nr); j <- which(!is.na(partner)); cur_r[j] <- rr[cbind(j, partner[j])]
  if (any(is.infinite(cur_p[i])) || any(is.infinite(cur_r[j]))) return(FALSE)
  block <- (rp < cur_p) & (t(rr) < matrix(cur_r, np, nr, byrow = TRUE))
  !any(block)
}

#' Position of each agent's partner on its own list (Inf if single or unacceptable).
partner_rank <- function(match, prefs) {
  vapply(seq_along(prefs), function(i) {
    if (is.na(match[i])) return(Inf)
    k <- match(match[i], prefs[[i]])
    if (is.na(k)) Inf else as.numeric(k)
  }, numeric(1))
}

#' Preference lists from a utility matrix (row i: i's utility from each partner), most preferred first.
prefs_from_utility <- function(U) lapply(seq_len(nrow(U)), function(i) order(U[i, ], decreasing = TRUE))

#' Misreports: keep the first k entries of a list; swap the entries at two positions.
truncate_pref <- function(l, k) l[seq_len(k)]
swap_pref <- function(l, pos) { l[pos] <- l[rev(pos)]; l }

## =================================================================
## Block V0 helpers (hand calculations)
## =================================================================

#' VCG in the position auction: values v (sorted, decreasing) per click, slot click rates ctr (decreasing).
#' The bidder in slot k pays sum_{j >= k} (ctr_j - ctr_{j+1}) v_{j+1}; bidders without a slot pay 0.
vcg_position <- function(v, ctr) {
  K <- length(ctr); nv <- length(v); ctr0 <- c(ctr, 0); v0 <- c(v, 0)   # v0: value 0 beyond the last bidder
  pay <- numeric(nv); per_click <- numeric(nv)
  for (k in seq_len(min(K, nv))) {
    j <- k:K
    pay[k] <- sum((ctr0[j] - ctr0[j + 1L]) * v0[pmin(j + 1L, nv + 1L)])
    per_click[k] <- pay[k] / ctr[k]
  }
  list(pay = pay, per_click = per_click, revenue = sum(pay))
}

#' GSP revenue under truthful bids: the bidder in slot k pays the next value per click.
gsp_revenue <- function(v, ctr) {
  K <- min(length(ctr), length(v) - 1L)
  sum(v[seq_len(K) + 1L] * ctr[seq_len(K)])
}

#' Bayesian persuasion with two states: prior P(good), the receiver acts iff the posterior is at least thresh.
#' persuasion_x: probability a bad state still sends the "act" signal; persuasion_value: the sender's
#' optimal probability of action, min(1, prior / thresh) (the concave closure at the prior).
persuasion_x <- function(prior, thresh) prior * (1 - thresh) / (thresh * (1 - prior))
persuasion_value <- function(prior, thresh) if (prior >= thresh) 1 else prior + (1 - prior) * persuasion_x(prior, thresh)

#' Cheap talk (Crawford and Sobel 1982), uniform state and quadratic loss, sender bias b: the largest number of
#' intervals N with N (N - 1) < 1 / (2 b); and the partition a_i = i / N + 2 b i (i - N), i = 0, ..., N.
cheap_talk_N <- function(b) vapply(b, function(bb) {
  N <- 1L
  while ((N + 1) * N < 1 / (2 * bb)) N <- N + 1L
  N
}, integer(1))
cheap_talk_partition <- function(b, N) { i <- 0:N; i / N + 2 * b * i * (i - N) }

## =================================================================
## Generators of the simulated blocks (01_dgp.R; tests/test_dgp.R re-runs them from the stored seeds)
## =================================================================

#' Block V1: n_auc_per_n auctions for each n in n_set; values iid from F; exact first-price equilibrium bids
#' without a reserve and with the optimal reserve (NA: no bid below r*). Call after set.seed().
gen_V1 <- function(P1) {
  D <- value_dist(P1$F); rs <- optimal_reserve(D)
  out <- vector("list", length(P1$n_set))
  for (k in seq_along(P1$n_set)) {
    n <- P1$n_set[k]; N <- P1$n_auc_per_n
    v <- D$quantile(stats::runif(N * n))
    out[[k]] <- data.table(n = as.integer(n), bidder = rep(seq_len(n), N), v_true = v,
                           b_fp = fpa_bid(v, n, D$cdf), b_fp_res = fpa_bid(v, n, D$cdf, reserve = rs))
  }
  dt <- rbindlist(out)
  dt[, auction := cumsum(bidder == 1L)]
  dt[, .(auction, n, bidder, v_true, b_fp, b_fp_res)]
}

#' Block V4: workers (proposers) and firms (receivers). Worker i's utility from firm j is a_j + e_ij, firm
#' j's from worker i is s_i + e'_ji: a common component (sd_common) and an idiosyncratic one (sd_idio).
#' Call after set.seed().
gen_V4 <- function(P4) {
  nw <- P4$n_workers; nf <- P4$n_firms
  a <- stats::rnorm(nf, 0, P4$sd_common)
  s <- stats::rnorm(nw, 0, P4$sd_common)
  Uw <- matrix(a, nw, nf, byrow = TRUE) + matrix(stats::rnorm(nw * nf, 0, P4$sd_idio), nw, nf)
  Uf <- matrix(s, nf, nw, byrow = TRUE) + matrix(stats::rnorm(nf * nw, 0, P4$sd_idio), nf, nw)
  list(prefs = list(prop = prefs_from_utility(Uw), recv = prefs_from_utility(Uf)),
       common = list(firm = a, worker = s))
}
