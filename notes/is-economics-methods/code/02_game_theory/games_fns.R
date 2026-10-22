# ------------------------------------------------------------------
# games_fns.R -- functions for chapter 2 (game theory).
#
# Conventions. A 2x2 game is a pair of payoff matrices indexed [Osprey's action, Boreal's action]:
# A holds Osprey's payoffs, B holds Boreal's (so a symmetric game has B = t(A)). Action 1 is Osprey's first
# row and Boreal's first column (High, S1, Yield ...). Mixed strategies are p = Pr(Osprey plays action 1)
# and q = Pr(Boreal plays action 1).
#
# Chapter card interface:
#   nash_2x2(A, B)                      -> list(pure = list(c(i, j), ...), mixed = c(p, q) or numeric(0))
#   logit_qre(A, B, lambda, start)      -> c(p = , q = ): a logit quantal response equilibrium (nleqslv);
#                                          lambda is the precision chi of the text
#   qre_logodds_ols(df)                 -> list(coef, ci, ...): payoffs from group log-odds by OLS; with
#                                          several replications (a rep column), coef is the mean estimate and
#                                          ci the 2.5% to 97.5% range of the estimates (Ruling P23: no
#                                          standard-error formulas in Part I)
#   gp_simulate(par, T, K, seed)        -> data.table(market, t, war, price, trigger_hit): Green-Porter panel
#   osd_check(par)                      -> logical: no profitable one-shot deviation from the trigger strategy
#   rubinstein_share(delta1, delta2)    -> the first proposer's share in Rubinstein's alternating offers
#   regret_matching(A, B, T, seed)      -> 2x2 matrix: empirical joint distribution of play (Hart and
#                                          Mas-Colell 2000), with the action histories as attributes
#   ce_violation(mu, A, B)              -> the largest expected gain from disobeying a recommendation of mu
#   gg_threshold(sigma, ...)            -> threshold signal of the global game's unique equilibrium
# Helpers used by 01 to 06:
#   qre_payoff_diff, qre_all, qre_bifurcation            logit QRE: payoff differences, all QRE, the
#                                                          precision at which the principal branch splits
#   g1_game, g1_equilibria, gen_G1, g1_group_freq        block G1 (QRE markets)
#   gp_values, gp_equilibrium, gp_optimal_trigger        block G2 (Green and Porter 1984)
#   regret_matching_reps, fictitious_play               block G3 (learning)
#   gg_roots, gg_min_slope, gg_sigma_max, gen_G4         block G4 (global game)
#   minmax_2x2, folk_sets                               the minmax payoff of a 2x2 game; the feasible payoff
#                                                       set and its individually rational part (folk theorem)
#
# Models, in the notation the text uses (notation_v3.md 5.3):
#   Logit QRE: p = Lambda(chi [u1(1, q) - u1(2, q)]), q = Lambda(chi [u2(1, p) - u2(2, p)]), Lambda the
#     logistic cdf; so log(p / (1 - p)) is chi times Osprey's expected payoff difference (McKelvey and
#     Palfrey 1995).
#   Green-Porter: Cournot duopoly with inverse demand p_t = a - Q_t + sigma z_t, z_t ~ N(0, 1) unseen by the
#     firms, unit cost `cost`, discount factor delta. Firms produce q_c each while the price stays at or above
#     the trigger price; a price below the trigger starts a war: n_pun periods of Cournot quantities, then
#     collusion again. beta = Pr(war starts | collusion) = Phi((trigger - (a - 2 q_c)) / sigma). Values:
#       Vc = pi_c + delta [(1 - beta) Vc + beta Vp],   Vp = pi_N (1 - delta^n_pun) / (1 - delta) + delta^n_pun Vc.
#     The one-shot deviation principle (finite punishment, discounting: continuity at infinity holds) makes
#     the strategy an equilibrium iff no single-period deviation pays in either phase. In a war period the
#     continuation does not depend on today's quantity, so Cournot is optimal; in collusion a firm choosing
#     q gets pi(q, q_c) + delta [(1 - beta(q)) Vc + beta(q) Vp] with beta(q) = Phi((trigger - (a - q_c - q)) / sigma).
#     With sigma = 0 the price is seen exactly and the trigger is read as "price below the trigger" (strict).
#   Global game: adopting pays theta + l - 1 (l the share adopting), not adopting pays 0. Prior (public
#     signal) theta ~ N(prior_mean, prior_sd^2), private signals x_i = theta + sigma e_i. In a threshold
#     equilibrium x*, l(theta) = Phi((theta - x*) / sigma), and the marginal player is indifferent:
#       m(x*) + Phi((m(x*) - x*) / sqrt(sigma^2 + v)) - 1 = 0,
#     m(x) = (alpha y + beta x) / (alpha + beta) the posterior mean, v = 1 / (alpha + beta), alpha = 1/prior_sd^2,
#     beta = 1/sigma^2. With a flat prior (no public signal), m(x) = x and E[l | x*] = 1/2: x* = 1/2.
# Packages: stats, nleqslv (logit_qre), data.table (gp_simulate).
# ------------------------------------------------------------------

## =================================================================
## 2x2 games: Nash equilibria, correlated equilibria, minmax
## =================================================================

.check_2x2 <- function(...) {
  for (M in list(...)) if (!is.matrix(M) || !identical(dim(M), c(2L, 2L)) || !is.numeric(M) || anyNA(M))
    stop("expected a numeric 2x2 matrix", call. = FALSE)
  invisible(TRUE)
}

#' Pure-strategy Nash equilibria (cells where neither player gains by switching; weak inequalities) and
#' the interior mixed equilibrium from the two indifference conditions: p makes Boreal indifferent,
#' q makes Osprey indifferent. `pure` lists c(i, j) in row-major order; `mixed` is numeric(0) when the
#' game has no interior mixed equilibrium (a dominant action, or a degenerate game). In a degenerate game
#' (a player indifferent between its actions against some or all opponent strategies) the equilibrium set
#' can contain a continuum of mixed equilibria; nash_2x2 does not report them (no game of this chapter is
#' degenerate).
nash_2x2 <- function(A, B) {
  .check_2x2(A, B)
  pure <- list()
  for (i in 1:2) for (j in 1:2)
    if (A[i, j] >= A[3L - i, j] && B[i, j] >= B[i, 3L - j]) pure[[length(pure) + 1L]] <- c(i, j)
  dA <- A[1, 1] - A[1, 2] - A[2, 1] + A[2, 2]
  dB <- B[1, 1] - B[1, 2] - B[2, 1] + B[2, 2]
  mixed <- numeric(0)
  if (dA != 0 && dB != 0) {
    p <- (B[2, 2] - B[2, 1]) / dB          # Boreal: p B[1,1] + (1-p) B[2,1] = p B[1,2] + (1-p) B[2,2]
    q <- (A[2, 2] - A[1, 2]) / dA          # Osprey: q A[1,1] + (1-q) A[1,2] = q A[2,1] + (1-q) A[2,2]
    if (p > 0 && p < 1 && q > 0 && q < 1) mixed <- c(p, q)
  }
  list(pure = pure, mixed = mixed)
}

#' Largest expected gain, over both players and every recommendation, from disobeying the joint
#' distribution mu (mu[i, j] = Pr(Osprey told i, Boreal told j)); 0 when mu is a correlated equilibrium.
#' For player 1 told i: sum_j mu[i, j] (A[3-i, j] - A[i, j]) (Aumann 1974; the unconditional form of the
#' constraint, as in Hart and Mas-Colell 2000).
ce_violation <- function(mu, A, B) {
  .check_2x2(mu, A, B)
  v <- c(sum(mu[1, ] * (A[2, ] - A[1, ])), sum(mu[2, ] * (A[1, ] - A[2, ])),
         sum(mu[, 1] * (B[, 2] - B[, 1])), sum(mu[, 2] * (B[, 1] - B[, 2])))
  max(0, v)
}

#' Minmax payoff of the row player of payoff matrix M (rows: own actions, columns: the opponent's):
#' min over the opponent's mixed strategy q of max over own pure actions. Use t(B) for Boreal.
minmax_2x2 <- function(M) {
  .check_2x2(M)
  val <- function(q) max(M %*% c(q, 1 - q))
  cand <- c(0, 1)
  den <- (M[1, 1] - M[1, 2]) - (M[2, 1] - M[2, 2])
  if (den != 0) {
    qx <- (M[2, 2] - M[1, 2]) / den        # where the two rows' expected payoffs cross
    if (qx > 0 && qx < 1) cand <- c(cand, qx)
  }
  min(vapply(cand, val, numeric(1)))
}

#' Folk-theorem sets of a 2x2 game. hull: the feasible payoff set, the convex hull of the four payoff pairs
#' (u1, u2) = (A[i, j], B[i, j]), vertices counter-clockwise; ir: its individually rational part
#' {u in hull : u1 >= minmax1, u2 >= minmax2}, the hull clipped by the two half-planes (Sutherland-Hodgman),
#' counter-clockwise from the vertex nearest the minmax point; minmax = c(minmax_2x2(A), minmax_2x2(t(B)));
#' area_hull, area_ir by the shoelace formula.
folk_sets <- function(A, B) {
  .check_2x2(A, B)
  pts <- cbind(u1 = as.vector(A), u2 = as.vector(B))
  pts <- pts[!duplicated(pts), , drop = FALSE]
  hull <- pts[rev(grDevices::chull(pts)), , drop = FALSE]                    # chull is clockwise
  mm <- c(minmax_2x2(A), minmax_2x2(t(B)))
  cut <- function(p, q, k, m) p + (m - p[k]) / (q[k] - p[k]) * (q - p)
  clip <- function(P, k, m) {                                                # keep the part with P[, k] >= m
    out <- matrix(numeric(0), 0, 2, dimnames = list(NULL, colnames(P))); n <- nrow(P)
    for (i in seq_len(n)) {
      cur <- P[i, ]; prv <- P[if (i == 1L) n else i - 1L, ]
      if (cur[k] >= m) {
        if (prv[k] < m) out <- rbind(out, cut(prv, cur, k, m))
        out <- rbind(out, cur)
      } else if (prv[k] >= m) out <- rbind(out, cut(prv, cur, k, m))
    }
    out
  }
  ir <- clip(clip(hull, 1L, mm[1]), 2L, mm[2])
  ir <- ir[!duplicated(round(ir, 12)), , drop = FALSE]
  start <- which.min((ir[, 1] - mm[1])^2 + (ir[, 2] - mm[2])^2)
  ir <- ir[c(start:nrow(ir), seq_len(start - 1L)), , drop = FALSE]
  rownames(hull) <- rownames(ir) <- NULL
  area <- function(P) 0.5 * abs(sum(P[, 1] * c(P[-1, 2], P[1, 2]) - c(P[-1, 1], P[1, 1]) * P[, 2]))
  list(hull = hull, ir = ir, minmax = mm, area_hull = area(hull), area_ir = area(ir))
}

## =================================================================
## Logit quantal response equilibrium
## =================================================================

.qre_d1 <- function(A, q) (A[1, 1] - A[2, 1]) * q + (A[1, 2] - A[2, 2]) * (1 - q)   # u1(1, q) - u1(2, q)
.qre_d2 <- function(B, p) (B[1, 1] - B[1, 2]) * p + (B[2, 1] - B[2, 2]) * (1 - p)   # u2(1, p) - u2(2, p)

#' Expected payoff differences c(d1 = u1(1, q) - u1(2, q), d2 = u2(1, p) - u2(2, p)) at (p, q).
qre_payoff_diff <- function(A, B, p, q) c(d1 = .qre_d1(A, q), d2 = .qre_d2(B, p))

#' A logit QRE: solves log(p/(1-p)) = lambda d1(q), log(q/(1-q)) = lambda d2(p) in log-odds with nleqslv,
#' from `start` (a pair of probabilities). With several QRE, the one found depends on `start`; qre_all
#' finds them all.
logit_qre <- function(A, B, lambda, start = c(.5, .5)) {
  .check_2x2(A, B)
  if (!is.numeric(lambda) || length(lambda) != 1L || !(lambda >= 0)) stop("logit_qre: lambda must be >= 0")
  f <- function(l) {
    p <- stats::plogis(l[1]); q <- stats::plogis(l[2])
    c(l[1] - lambda * .qre_d1(A, q), l[2] - lambda * .qre_d2(B, p))
  }
  st <- stats::qlogis(pmin(pmax(as.numeric(start), 1e-12), 1 - 1e-12))
  s <- nleqslv::nleqslv(st, f, control = list(xtol = 1e-15, ftol = 1e-13, maxit = 1000))
  if (max(abs(s$fvec)) > 1e-9) stop("logit_qre: no convergence (", s$message, ")")
  c(p = stats::plogis(s$x[1]), q = stats::plogis(s$x[2]))
}

#' Every logit QRE of a 2x2 game: the roots of g(l) = l - lambda d1(Lambda(lambda d2(Lambda(l)))) in Osprey's
#' log-odds l, found by a sign scan over [-(lambda M + 1), lambda M + 1] (M = max |d1|, so no root lies
#' outside) and uniroot. Returns a matrix with columns p, q, one row per QRE, sorted by p.
qre_all <- function(A, B, lambda, n_grid = 8001L) {
  .check_2x2(A, B)
  g <- function(l) l - lambda * .qre_d1(A, stats::plogis(lambda * .qre_d2(B, stats::plogis(l))))
  lim <- lambda * max(abs(c(A[1, 1] - A[2, 1], A[1, 2] - A[2, 2]))) + 1
  ls <- seq(-lim, lim, length.out = n_grid)
  gv <- g(ls)
  roots <- ls[gv == 0]
  idx <- which(gv[-n_grid] * gv[-1] < 0)
  roots <- sort(c(roots, vapply(idx, function(i)
    stats::uniroot(g, ls[c(i, i + 1)], tol = 1e-14 * max(1, abs(ls[i])), maxiter = 500)$root, numeric(1))))
  p <- stats::plogis(roots)
  q <- stats::plogis(lambda * .qre_d2(B, p))
  cbind(p = p, q = q)
}

#' The precision at which the principal QRE branch (the one that starts at (1/2, 1/2) as lambda -> 0)
#' stops being the only QRE: the Jacobian of the log-odds system is singular there,
#' 1 - lambda^2 b1 b2 p (1 - p) q (1 - q) = 0, with b1 = d d1/dq and b2 = d d2/dp. Searched in (lower, upper),
#' where the principal branch is the unique QRE below the root and the middle QRE above it.
qre_bifurcation <- function(A, B, lower, upper) {
  b1 <- (A[1, 1] - A[2, 1]) - (A[1, 2] - A[2, 2]); b2 <- (B[1, 1] - B[1, 2]) - (B[2, 1] - B[2, 2])
  principal <- function(lam) { r <- qre_all(A, B, lam); r[ceiling(nrow(r) / 2), ] }
  h <- function(lam) { r <- principal(lam); 1 - lam^2 * b1 * b2 * r[["p"]] * (1 - r[["p"]]) * r[["q"]] * (1 - r[["q"]]) }
  lam <- stats::uniroot(h, c(lower, upper), tol = 1e-12)$root
  c(lambda = lam, principal(lam))
}

#' Grouped-data inversion of the logit QRE by OLS (block G1). `df` has one row per group: p and q, the
#' groups' observed frequencies of action 1 for Osprey and Boreal, x the payoff shifter that adds to
#' Osprey's action-1 payoffs (USD, known), chi the known precision and, optionally, rep, a replication
#' index. Under the QRE,
#'   log(p / (1 - p)) / chi - x = u11 q - u22 (1 - q),
#' u11 = A[1,1] - A[2,1] (Osprey's gain from S1 over S2 when Boreal is on S1) and u22 = A[2,2] - A[1,2]
#' (the gain from S2 when Boreal is on S2). OLS without intercept on (q, -(1 - q)), one fit per replication
#' (closed-form 2x2 normal equations; at least 3 groups each). est: the per-replication estimates; coef:
#' their mean (the estimate itself without a rep column). With several replications, ci and ci_u22 are the
#' central `level` ranges of the estimates across replications (quantiles (1 - level)/2 and (1 + level)/2,
#' type 7) and sd their standard deviations: a statement about the estimator over simulated samples, not a
#' standard-error formula from one sample (Ruling P23). With one replication, ci, ci_u22 and sd are NULL.
qre_logodds_ols <- function(df, level = 0.95) {
  need <- c("p", "q", "x", "chi")
  if (!all(need %in% names(df))) stop("qre_logodds_ols: df needs columns p, q, x, chi")
  if (any(df$p <= 0 | df$p >= 1)) stop("qre_logodds_ols: a choice frequency of 0 or 1 has no log-odds")
  rep <- if (is.null(df$rep)) rep(1L, nrow(df)) else df$rep
  y <- stats::qlogis(df$p) / df$chi - df$x
  a <- df$q; b <- -(1 - df$q)
  S <- rowsum(cbind(aa = a * a, ab = a * b, bb = b * b, ay = a * y, by = b * y, n = 1), rep, reorder = TRUE)
  if (any(S[, "n"] < 3)) stop("qre_logodds_ols: needs at least 3 groups per replication")
  det <- S[, "aa"] * S[, "bb"] - S[, "ab"]^2
  est <- cbind(u11 = (S[, "bb"] * S[, "ay"] - S[, "ab"] * S[, "by"]) / det,
               u22 = (S[, "aa"] * S[, "by"] - S[, "ab"] * S[, "ay"]) / det)
  out <- list(coef = colMeans(est), est = est, n_rep = nrow(est), n_groups = unname(S[, "n"]))
  if (nrow(est) > 1L) {
    pr <- c((1 - level) / 2, (1 + level) / 2)
    rng <- function(v) stats::setNames(stats::quantile(v, pr, type = 7, names = FALSE), c("lo", "hi"))
    out$ci <- rng(est[, "u11"]); out$ci_u22 <- rng(est[, "u22"]); out$sd <- apply(est, 2, stats::sd)
  }
  out
}

## =================================================================
## Block G1: QRE markets
## =================================================================

#' Market m's game: the base matrices plus the shifter x on Osprey's S1 row and Boreal's S1 column.
g1_game <- function(P1, x) list(A = P1$A0 + rbind(c(x, x), 0), B = P1$B0 + cbind(c(x, x), 0))

#' Every QRE at every shifter value of both blocks: data.frame(block, x, role, p, q, n_qre). role is
#' "unique" (one QRE), or "S2" (lowest p), "middle", "S1" (highest p) when there are three.
g1_equilibria <- function(P1) {
  one <- function(block, x) {
    g <- g1_game(P1, x); r <- qre_all(g$A, g$B, P1$lambda[[block]])
    role <- if (nrow(r) == 1L) "unique" else if (nrow(r) == 3L) c("S2", "middle", "S1") else paste0("root", seq_len(nrow(r)))
    data.frame(block = block, x = x, role = role, p = r[, "p"], q = r[, "q"], n_qre = nrow(r), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, c(lapply(P1$x_unique, function(x) one("unique", x)), lapply(P1$x_high, function(x) one("high", x))))
  rownames(out) <- NULL
  out
}

#' Draw block G1 (call set.seed first). Unique block: n_markets_unique markets spread equally over
#' x_unique, each playing the unique QRE of its game at precision lambda["unique"]. High block:
#' n_markets_high markets spread equally over x_high, precision lambda["high"] (three QRE at every x);
#' history: S1 came first in the market with probability plogis(kappa x), and the market plays the QRE
#' that favours its first standard (S1: highest p; S2: lowest p). Each market makes n_dec integration
#' decisions; k_osprey and k_boreal count the decisions on S1. Draw order: unique-block Osprey counts,
#' Boreal counts, high-block histories, Osprey counts, Boreal counts. Columns ending in _true are truth.
gen_G1 <- function(P1, eq) {
  pick <- function(block, x, role) {
    i <- vapply(seq_along(x), function(k) which(eq$block == block & abs(eq$x - x[k]) < 1e-12 & eq$role == role[k]), integer(1))
    eq[i, c("p", "q")]
  }
  xu <- rep(P1$x_unique, each = P1$n_markets_unique / length(P1$x_unique))
  eu <- pick("unique", xu, rep("unique", length(xu)))
  ku_o <- stats::rbinom(length(xu), P1$n_dec, eu$p)
  ku_b <- stats::rbinom(length(xu), P1$n_dec, eu$q)
  xh <- rep(P1$x_high, each = P1$n_markets_high / length(P1$x_high))
  hist <- ifelse(stats::runif(length(xh)) < stats::plogis(P1$kappa * xh), "S1", "S2")
  eh <- pick("high", xh, hist)
  kh_o <- stats::rbinom(length(xh), P1$n_dec, eh$p)
  kh_b <- stats::rbinom(length(xh), P1$n_dec, eh$q)
  data.frame(block = rep(c("unique", "high"), c(length(xu), length(xh))),
             market = seq_len(length(xu) + length(xh)), x = c(xu, xh),
             history = c(rep(NA_character_, length(xu)), hist), n_dec = as.integer(P1$n_dec),
             k_osprey = c(ku_o, kh_o), k_boreal = c(ku_b, kh_b),
             p_true = c(eu$p, eh$p), q_true = c(eu$q, eh$q), stringsAsFactors = FALSE)
}

#' Group frequencies for the inversion: one row per value of the `by` columns (x, or x and history),
#' with n (decisions), p and q (frequencies of S1), sorted by the grouping columns.
g1_group_freq <- function(d, by) {
  key <- interaction(d[by], drop = TRUE, lex.order = TRUE)
  n <- rowsum(d$n_dec, key); ko <- rowsum(d$k_osprey, key); kb <- rowsum(d$k_boreal, key)
  first <- !duplicated(key)
  out <- d[first, by, drop = FALSE][order(key[first]), , drop = FALSE]
  out$n_markets <- as.vector(table(key))
  out$n <- as.vector(n); out$p <- as.vector(ko) / out$n; out$q <- as.vector(kb) / out$n
  rownames(out) <- NULL
  out
}

## =================================================================
## Block G2: Green and Porter (1984), imperfect public monitoring
## =================================================================
## par: list(a, cost, delta, n_pun, sigma, trigger, q_c). Prices in USD, quantities in thousand units,
## profits in thousand USD per period.

.gp_profit <- function(q, q_other, par) (par$a - q_other - q - par$cost) * q
.gp_cournot <- function(par) { q <- (par$a - par$cost) / 3; list(q = q, profit = q^2) }
.gp_detect <- function(price_mean, par) {                 # Pr(observed price < trigger | expected price)
  if (par$sigma == 0) as.numeric(price_mean < par$trigger) else stats::pnorm((par$trigger - price_mean) / par$sigma)
}
.gp_annuity <- function(par) {
  d <- par$delta; N <- par$n_pun
  dN <- d^N
  list(AN = (1 - dN) / (1 - d), dN = dN)
}

#' Values of the trigger strategy: beta (Pr(war starts | collusion)), Vc (value in collusion), Vp (value at
#' the start of a war), pi_c, and the Cournot quantity and profit q_N, pi_N.
gp_values <- function(par) {
  cn <- .gp_cournot(par); an <- .gp_annuity(par); d <- par$delta
  beta <- .gp_detect(par$a - 2 * par$q_c, par)
  pi_c <- .gp_profit(par$q_c, par$q_c, par)
  Vc <- (pi_c + d * beta * cn$profit * an$AN) / (1 - d * (1 - beta) - d * beta * an$dN)
  list(beta = beta, Vc = Vc, Vp = cn$profit * an$AN + an$dN * Vc, pi_c = pi_c, q_N = cn$q, pi_N = cn$profit)
}

#' One-shot deviation check of the trigger strategy in par. Collusive phase: the best single-period
#' quantity q (all later periods on the strategy) must not beat Vc; searched on a grid of 4001 quantities
#' in [0, a - cost] refined by optimize (sigma > 0), or exactly on the two pieces "price stays at or above
#' the trigger" and "price falls below it" (sigma = 0). War phase: the continuation does not depend on
#' today's quantity, so the gain is max_q pi(q, q_N) - pi_N (optimize). TRUE when both gains are at most
#' tol * Vc; the gains and the best collusive-phase deviation are attributes.
osd_check <- function(par, tol = 1e-8) {
  v <- gp_values(par)
  qmax <- par$a - par$cost
  dev <- function(q) {
    b <- .gp_detect(par$a - par$q_c - q, par)
    .gp_profit(q, par$q_c, par) + par$delta * ((1 - b) * v$Vc + b * v$Vp)
  }
  br <- (par$a - par$cost - par$q_c) / 2                    # static best response to q_c
  if (par$sigma == 0) {
    qs <- max(0, par$a - par$q_c - par$trigger)              # largest quantity that keeps the price at the trigger
    q1 <- min(max(br, 0), qs)                                # best undetected quantity
    q2 <- max(br, qs)                                        # best detected quantity (a supremum if br <= qs)
    vals <- c(.gp_profit(q1, par$q_c, par) + par$delta * v$Vc, .gp_profit(q2, par$q_c, par) + par$delta * v$Vp)
    q_dev <- c(q1, q2)[which.max(vals)]; best <- max(vals)
  } else {
    qs <- sort(unique(c(seq(0, qmax, length.out = 4001), par$q_c, br)))
    w <- dev(qs); i <- which.max(w)
    o <- stats::optimize(dev, qs[c(max(1L, i - 1L), min(length(qs), i + 1L))], maximum = TRUE, tol = 1e-12)
    if (o$objective >= w[i]) { best <- o$objective; q_dev <- o$maximum } else { best <- w[i]; q_dev <- qs[i] }
  }
  gain_c <- best - v$Vc
  op <- stats::optimize(function(q) .gp_profit(q, v$q_N, par), c(0, qmax), maximum = TRUE, tol = 1e-12)
  gain_p <- max(0, op$objective - v$pi_N)
  ok <- gain_c <= tol * max(1, abs(v$Vc)) && gain_p <= tol * max(1, abs(v$Vc))
  structure(ok, gain_collusive = gain_c, gain_punish = gain_p, q_dev = q_dev)
}

#' The most collusive symmetric equilibrium quantity for the trigger price in par (sigma > 0): roots of the
#' first-order condition  a - cost - 3 q = delta (Vc(q) - Vp(q)) phi(z) / sigma,  z = (trigger - (a - 2q)) / sigma
#' (marginal profit = marginal expected punishment), found by a sign scan on [q_M / 2, q_N) and uniroot;
#' each root is kept only if osd_check passes (a global best response), and the one with the largest Vc is
#' returned in par$q_c. Cournot (q_N, no collusion) if none passes. sigma = 0: q_c = (a - trigger) / 2, the
#' quantity that puts the collusive price exactly at the trigger (the caller checks it with osd_check).
gp_equilibrium <- function(par) {
  qN <- (par$a - par$cost) / 3; qM <- (par$a - par$cost) / 4
  if (par$sigma == 0) return(modifyList(par, list(q_c = (par$a - par$trigger) / 2)))
  foc <- function(qc) {
    pp <- par; pp$q_c <- qc; v <- gp_values(pp)
    z <- (par$trigger - (par$a - 2 * qc)) / par$sigma
    (par$a - par$cost - 3 * qc) - par$delta * (v$Vc - v$Vp) * stats::dnorm(z) / par$sigma
  }
  qs <- seq(qM / 2, qN - 1e-6, length.out = 1201)
  fv <- vapply(qs, foc, numeric(1))
  idx <- which(fv[-length(fv)] * fv[-1] < 0)
  roots <- vapply(idx, function(i) stats::uniroot(foc, qs[c(i, i + 1)], tol = 1e-14)$root, numeric(1))
  roots <- roots[roots < qN - 1e-5]
  best <- NULL; best_v <- -Inf
  for (r in roots) {
    pp <- modifyList(par, list(q_c = r))
    if (isTRUE(as.vector(osd_check(pp)))) { vc <- gp_values(pp)$Vc; if (vc > best_v) { best <- r; best_v <- vc } }
  }
  modifyList(par, list(q_c = if (is.null(best)) qN else best))
}

#' The trigger price that maximises the collusive value Vc, with q_c the most collusive equilibrium quantity
#' for each trigger (gp_equilibrium): a grid of 81 triggers on [a - 2 q_N - 3 sigma, a - 2 q_M + sigma] and
#' optimize around the best. sigma = 0: the smallest q_c >= q_M whose trigger at the collusive price
#' a - 2 q_c deters the only tempting deviation, the static best response (a higher quantity, always
#' detected): slack(q_c) = Vc - (pi(br(q_c), q_c) + delta Vp) >= 0; the monopoly split when patience allows
#' it, otherwise the root of the slack found by a scan of [q_M, q_N) and uniroot. When no trigger supports
#' any collusion the result is Cournot forever: q_c = q_N and trigger = -Inf (no war is ever started).
gp_optimal_trigger <- function(par) {
  qN <- (par$a - par$cost) / 3; qM <- (par$a - par$cost) / 4
  no_collusion <- modifyList(par, list(q_c = qN, trigger = -Inf))
  if (par$sigma == 0) {
    at_q <- function(qc) modifyList(par, list(q_c = qc, trigger = par$a - 2 * qc))
    slack <- function(qc) {
      v <- gp_values(at_q(qc)); br <- (par$a - par$cost - qc) / 2
      v$Vc - (.gp_profit(br, qc, par) + par$delta * v$Vp)
    }
    if (slack(qM) >= 0) return(at_q(qM))
    qs <- seq(qM, qN, length.out = 2001)[-2001]
    sv <- vapply(qs, slack, numeric(1)); i <- which(sv > 0)[1]
    if (is.na(i)) return(no_collusion)
    return(at_q(stats::uniroot(slack, qs[c(i - 1L, i)], tol = 1e-14)$root))
  }
  vc_of <- function(trig) gp_values(gp_equilibrium(modifyList(par, list(trigger = trig))))$Vc
  trs <- seq(par$a - 2 * qN - 3 * par$sigma, par$a - 2 * qM + par$sigma, length.out = 81)
  vs <- vapply(trs, vc_of, numeric(1)); i <- which.max(vs)
  if (max(vs) <= gp_values(no_collusion)$Vc * (1 + 1e-12)) return(no_collusion)
  if (i == 1L || i == length(trs)) stop("gp_optimal_trigger: the best trigger is at the edge of the search range")
  o <- stats::optimize(vc_of, trs[c(i - 1L, i + 1L)], maximum = TRUE, tol = 1e-10)
  gp_equilibrium(modifyList(par, list(trigger = o$maximum)))
}

#' Simulate K markets for T periods under the trigger strategy in par (the equilibrium path: no firm ever
#' deviates). Draw order after set.seed(seed): the K x T standard normal demand shocks z (market by row),
#' then the starting state of each market from the stationary distribution of the phase chain (war with
#' probability beta n_pun / (1 + beta n_pun), periods left uniform on 1..n_pun), so the same seed gives the
#' same shocks at every sigma. war: the period is a war (Cournot) period; trigger_hit: a collusive period
#' whose price fell below the trigger (a war starts next period). Rows ordered by market, then t.
gp_simulate <- function(par, T, K, seed) {
  set.seed(seed)
  z <- matrix(stats::rnorm(K * T), K, T)
  v <- gp_values(par); N <- as.integer(par$n_pun)
  p_war <- v$beta * N / (1 + v$beta * N)
  u0 <- stats::runif(K); left0 <- sample.int(N, K, replace = TRUE)
  left <- ifelse(u0 < p_war, left0, 0L)
  war <- hit <- matrix(FALSE, K, T); price <- matrix(NA_real_, K, T)
  p_coll <- par$a - 2 * par$q_c; p_cour <- par$a - 2 * v$q_N
  for (t in seq_len(T)) {
    w <- left > 0L
    pr <- ifelse(w, p_cour, p_coll) + par$sigma * z[, t]
    h <- !w & pr < par$trigger
    war[, t] <- w; price[, t] <- pr; hit[, t] <- h
    left <- ifelse(w, left - 1L, ifelse(h, N, 0L))
  }
  data.table::data.table(market = rep(seq_len(K), each = T), t = rep(seq_len(T), times = K),
                         war = as.vector(t(war)), price = as.vector(t(price)), trigger_hit = as.vector(t(hit)))
}

## =================================================================
## Bargaining
## =================================================================

#' Rubinstein (1982): the unique SPE share of the first proposer when the players discount by delta1
#' (proposer) and delta2 (responder): (1 - delta2) / (1 - delta1 delta2); 1 / (1 + delta) when equal.
#' Vectorised.
rubinstein_share <- function(delta1, delta2) (1 - delta2) / (1 - delta1 * delta2)

## =================================================================
## Block G3: learning
## =================================================================

# Regret matching for R independent runs at once (vectors over runs). Player i, having played j last
# period, switches to the other action with probability max(D_i(j), 0) / (t mu_c), where D_i(j) is the sum,
# over past periods in which i played j, of the payoff gain had i played the other action instead
# (Hart and Mas-Colell 2000, the procedure of their Section 2; mu_c > 2 M (m - 1) with M = max |u|, m = 2).
# Draw order: initial actions of player 1 (R draws), of player 2 (R draws), then per period 2R uniforms.
# At the periods in `checkpoints` the ce_violation of each run's empirical distribution so far is recorded;
# with keep = TRUE the action histories are returned as T x R matrices a1, a2.
.rm_core <- function(A, B, T, R, keep = FALSE, checkpoints = integer(0)) {
  mu_c <- 2 * max(abs(c(A, B))) + 1
  s1 <- sample.int(2L, R, replace = TRUE); s2 <- sample.int(2L, R, replace = TRUE)
  D1 <- D2 <- matrix(0, R, 2)
  cnt <- matrix(0L, R, 4)
  ri <- seq_len(R)
  if (keep) { a1 <- a2 <- matrix(0L, T, R) }
  path <- matrix(NA_real_, R, length(checkpoints)); k <- 1L
  for (t in seq_len(T)) {
    if (keep) { a1[t, ] <- s1; a2[t, ] <- s2 }
    cell <- cbind(ri, (s1 - 1L) * 2L + s2)
    cnt[cell] <- cnt[cell] + 1L
    if (k <= length(checkpoints) && t == checkpoints[k]) { path[, k] <- .ce_violation_rows(cnt / t, A, B); k <- k + 1L }
    D1[cbind(ri, s1)] <- D1[cbind(ri, s1)] + A[cbind(3L - s1, s2)] - A[cbind(s1, s2)]
    D2[cbind(ri, s2)] <- D2[cbind(ri, s2)] + B[cbind(s1, 3L - s2)] - B[cbind(s1, s2)]
    sw1 <- pmax(D1[cbind(ri, s1)], 0) / (t * mu_c); sw2 <- pmax(D2[cbind(ri, s2)], 0) / (t * mu_c)
    u <- stats::runif(2L * R)
    s1 <- ifelse(u[ri] < sw1, 3L - s1, s1); s2 <- ifelse(u[R + ri] < sw2, 3L - s2, s2)
  }
  out <- list(counts = cnt, mu_c = mu_c, path = path)
  if (keep) { out$a1 <- a1; out$a2 <- a2 }
  out
}

# ce_violation for many joint distributions at once: M has one row per distribution, columns mu11, mu12,
# mu21, mu22.
.ce_violation_rows <- function(M, A, B) {
  pmax(0, M[, 1] * (A[2, 1] - A[1, 1]) + M[, 2] * (A[2, 2] - A[1, 2]),
          M[, 3] * (A[1, 1] - A[2, 1]) + M[, 4] * (A[1, 2] - A[2, 2]),
          M[, 1] * (B[1, 2] - B[1, 1]) + M[, 3] * (B[2, 2] - B[2, 1]),
          M[, 2] * (B[1, 1] - B[1, 2]) + M[, 4] * (B[2, 1] - B[2, 2]))
}

#' One regret-matching run of T periods: the empirical joint distribution of play (mu[i, j] = share of
#' periods with Osprey on i and Boreal on j), with attributes a1, a2 (the action histories) and mu_c.
regret_matching <- function(A, B, T, seed) {
  .check_2x2(A, B)
  set.seed(seed)
  r <- .rm_core(A, B, T, R = 1L, keep = TRUE)
  mu <- matrix(r$counts[1, ] / T, 2, 2, byrow = TRUE)
  attr(mu, "a1") <- r$a1[, 1]; attr(mu, "a2") <- r$a2[, 1]; attr(mu, "mu_c") <- r$mu_c
  mu
}

#' R independent regret-matching runs of T periods (call set.seed first): data.frame(run, mu11, mu12, mu21,
#' mu22, viol) with each run's empirical joint distribution at T and its ce_violation. With checkpoints, the
#' attribute "path" is the R x length(checkpoints) matrix of each run's ce_violation after t periods.
regret_matching_reps <- function(A, B, T, R, checkpoints = integer(0)) {
  .check_2x2(A, B)
  r <- .rm_core(A, B, T, R, checkpoints = checkpoints)
  mu <- r$counts / T
  viol <- vapply(seq_len(R), function(i) ce_violation(matrix(mu[i, ], 2, 2, byrow = TRUE), A, B), numeric(1))
  out <- data.frame(run = seq_len(R), mu11 = mu[, 1], mu12 = mu[, 2], mu21 = mu[, 3], mu22 = mu[, 4], viol = viol)
  if (length(checkpoints)) { attr(out, "path") <- r$path; attr(out, "checkpoints") <- checkpoints }
  out
}

#' Fictitious play for T periods: each player best-responds to the empirical frequency of the opponent's
#' past actions, seeded with one fictitious observation of a random action; exact ties broken at random.
#' Draw order after set.seed(seed): the two fictitious observations, then two uniforms per period.
#' Returns list(freq = c(p, q) (shares of periods on action 1), a1, a2).
fictitious_play <- function(A, B, T, seed) {
  .check_2x2(A, B)
  set.seed(seed)
  c1 <- c2 <- c(0, 0)                                       # c1: Osprey's counts of Boreal's actions; c2: Boreal's of Osprey's
  c1[sample.int(2L, 1L)] <- 1; c2[sample.int(2L, 1L)] <- 1
  a1 <- a2 <- integer(T)
  for (t in seq_len(T)) {
    u <- stats::runif(2L)
    e1 <- drop(A %*% c1); e2 <- drop(crossprod(B, c2))       # expected payoffs (unnormalised beliefs)
    a1[t] <- if (e1[1] > e1[2]) 1L else if (e1[1] < e1[2]) 2L else if (u[1] < 0.5) 1L else 2L
    a2[t] <- if (e2[1] > e2[2]) 1L else if (e2[1] < e2[2]) 2L else if (u[2] < 0.5) 1L else 2L
    c1[a2[t]] <- c1[a2[t]] + 1; c2[a1[t]] <- c2[a1[t]] + 1
  }
  list(freq = c(p = mean(a1 == 1L), q = mean(a2 == 1L)), a1 = a1, a2 = a2)
}

## =================================================================
## Block G4: global game
## =================================================================

.gg_F <- function(x, sigma, y, tau) {
  a <- 1 / tau^2; b <- 1 / sigma^2
  m <- (a * y + b * x) / (a + b)
  m + stats::pnorm((m - x) / sqrt(sigma^2 + 1 / (a + b))) - 1
}

#' Every threshold equilibrium: the roots of the marginal player's indifference condition in x (sign scan on
#' [prior_mean - 20, prior_mean + 20] and uniroot).
gg_roots <- function(sigma, prior_mean, prior_sd, n_grid = 40001L) {
  xs <- seq(prior_mean - 20, prior_mean + 20, length.out = n_grid)
  f <- function(x) .gg_F(x, sigma, prior_mean, prior_sd)
  fv <- f(xs); idx <- which(fv[-n_grid] * fv[-1] < 0)
  sort(c(xs[fv == 0], vapply(idx, function(i) stats::uniroot(f, xs[c(i, i + 1)], tol = 1e-15)$root, numeric(1))))
}

#' Threshold signal x* of the global game's equilibrium (the state at which half the players adopt), for
#' each sigma. Without prior_mean and prior_sd the prior is flat: the marginal player's posterior over the
#' share adopting is uniform (Laplacian belief), E[l | x*] = Phi(0) = 1/2, so x* = 1 - Phi(0) = 1/2 for every
#' sigma. With a normal public prior the indifference condition is solved numerically; an error if it has
#' more than one root (multiple equilibria).
gg_threshold <- function(sigma, prior_mean = NULL, prior_sd = NULL) {
  if (is.null(prior_mean) || is.null(prior_sd)) return(rep(1 - stats::pnorm(0), length(sigma)))
  vapply(sigma, function(s) {
    r <- gg_roots(s, prior_mean, prior_sd)
    if (length(r) != 1L) stop("gg_threshold: ", length(r), " threshold equilibria at sigma = ", s)
    r
  }, numeric(1))
}

#' Minimum over x of the slope of the indifference condition, (beta - alpha / (s sqrt(2 pi))) / (alpha + beta),
#' s = sqrt(sigma^2 + 1 / (alpha + beta)); positive means at most one threshold equilibrium for every prior mean.
gg_min_slope <- function(sigma, prior_sd) {
  a <- 1 / prior_sd^2; b <- 1 / sigma^2; s <- sqrt(sigma^2 + 1 / (a + b))
  (b - a / (s * sqrt(2 * pi))) / (a + b)
}

#' The largest private noise sigma for which the threshold equilibrium is unique whatever the prior mean
#' (gg_min_slope = 0), given the public prior's sd.
gg_sigma_max <- function(prior_sd) stats::uniroot(function(s) gg_min_slope(s, prior_sd), c(1e-4, 50), tol = 1e-15)$root

#' Draw block G4 (call set.seed first). pop: for each sigma (in P4$sigma order) and each theta in
#' P4$theta_grid, n_pop players draw signals theta + sigma e and adopt when the signal is at least x*;
#' n_adopt counts them. mc: for each sigma, n_mc draws of theta from the marginal player's posterior at x*
#' and of one other player's signal; El_sim is the share of draws in which that player adopts (E[l | x*]).
gen_G4 <- function(P4, x_star) {
  pop <- do.call(rbind, lapply(seq_along(P4$sigma), function(k) {
    s <- P4$sigma[k]
    n_adopt <- vapply(P4$theta_grid, function(th) sum(th + s * stats::rnorm(P4$n_pop) >= x_star[k]), integer(1))
    data.frame(sigma = s, theta = P4$theta_grid, n_pop = as.integer(P4$n_pop), n_adopt = n_adopt)
  }))
  a <- 1 / P4$prior_sd^2
  mc <- do.call(rbind, lapply(seq_along(P4$sigma), function(k) {
    s <- P4$sigma[k]; b <- 1 / s^2
    m <- (a * P4$prior_mean + b * x_star[k]) / (a + b); v <- 1 / (a + b)
    th <- m + sqrt(v) * stats::rnorm(P4$n_mc)
    xj <- th + s * stats::rnorm(P4$n_mc)
    data.frame(sigma = s, x_star = x_star[k], m_post = m, El_sim = mean(xj >= x_star[k]), n_mc = as.integer(P4$n_mc))
  }))
  list(pop = pop, mc = mc)
}
