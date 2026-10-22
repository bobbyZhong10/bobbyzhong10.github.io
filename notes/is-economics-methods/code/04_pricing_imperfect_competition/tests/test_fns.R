# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of pricing_fns.R against closed forms and brute force.
# The first block is the chapter card, verbatim (comment translated, Ruling P6).
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("pricing_fns.R")

## ---- chapter card (verbatim) -------------------------------------------------------------
stopifnot(abs(cournot_n(12, 1, 2, 2)$p - 16/3) < 1e-12, abs(passthrough_cournot(3) - 3/4) < 1e-12)
pr <- bertrand_logit(delta = c(1, 1), alpha = 1, c = c(1, 1)); s <- logit_share(pr, c(1, 1), 1)
stopifnot(max(abs(pr - 1 - 1 / (1 - s))) < 1e-10)                                   # single-product logit markup
bb <- bbm_segment(values = c(1, 2, 3), probs = rep(1/3, 3))
stopifnot(abs(bb$profit - 4/3) < 1e-12, abs(bb$profit - bb$uniform_profit) < 1e-12, abs(bb$total - 2) < 1e-12)

## ---- cournot_n and passthrough_cournot: closed forms, best responses, limits -------------------
for (n in c(1, 2, 3, 5, 10)) {
  r <- cournot_n(12, 1, 2, n)
  stopifnot(abs(r$q - 10 / (n + 1)) < 1e-12, abs(r$p - (12 - n * r$q)) < 1e-12,
            abs(r$profit - (r$p - 2) * r$q) < 1e-12,
            abs((r$p - 2) - 1 * r$q) < 1e-12)                       # each firm's FOC: p - c + q p'(Q) = 0
  h <- 1e-4                                                          # pass-through = dp/dc (exact for linear demand)
  stopifnot(abs((cournot_n(12, 1, 2 + h, n)$p - cournot_n(12, 1, 2 - h, n)$p) / (2 * h) - passthrough_cournot(n)) < 1e-9)
}
m1 <- cournot_n(12, 1, 2, 1)
stopifnot(m1$q == 5, m1$p == 7, m1$profit == 25)                    # n = 1 is the monopoly
stopifnot(abs(cournot_n(12, 1, 2, 1e6)$p - 2) < 1e-4)               # n to infinity: price to marginal cost
stopifnot(all(diff(sapply(1:20, function(n) cournot_n(12, 1, 2, n)$p)) < 0),
          all(diff(passthrough_cournot(1:20)) > 0), abs(passthrough_cournot(1) - 1/2) < 1e-15)
# a different slope b: quantities scale by 1/b, the price does not depend on b
stopifnot(abs(cournot_n(12, 2, 2, 3)$q - 10 / (2 * 4)) < 1e-12, abs(cournot_n(12, 2, 2, 3)$p - cournot_n(12, 1, 2, 3)$p) < 1e-12)

## ---- passthrough_monopoly: 1 / (2 + Q p''(Q) / p'(Q)) ------------------------------------------
# linear p = a - bQ: 1/2; constant elasticity p = (Q/A)^(-1/e): e / (e - 1); exponential p = lambda log(A/Q): 1
stopifnot(abs(passthrough_monopoly(Q = 5, dp = -1, d2p = 0) - 1/2) < 1e-15)
for (e in c(1.5, 2, 3)) {
  Q <- 4; A <- 10; p <- (Q / A)^(-1 / e)
  stopifnot(abs(passthrough_monopoly(Q, dp = -p / (e * Q), d2p = (1 / e) * (1 / e + 1) * p / Q^2) - e / (e - 1)) < 1e-12)
}
stopifnot(abs(passthrough_monopoly(Q = 3, dp = -0.5 / 3, d2p = 0.5 / 9) - 1) < 1e-12)
# the formula against a brute-force monopoly price under constant elasticity e = 2: p = c e / (e - 1)
mono_ce <- function(cc) optimize(function(p) (p - cc) * 10 * p^(-2), c(cc, 50 * cc), maximum = TRUE, tol = 1e-12)$maximum
stopifnot(abs((mono_ce(2.01) - mono_ce(1.99)) / 0.02 - 2) < 1e-4)

## ---- logit_share ---------------------------------------------------------------------------
s3 <- logit_share(c(1, 2, 3), c(0.5, 1, 2), 0.8)
e3 <- exp(c(0.5, 1, 2) - 0.8 * c(1, 2, 3))
stopifnot(max(abs(s3 - e3 / (1 + sum(e3)))) < 1e-15, all(s3 > 0), sum(s3) < 1)
stopifnot(abs(logit_share(0, 0, 1) - 1/2) < 1e-15)                  # one product with zero utility: half the market

## ---- bertrand_logit, joint_max_logit, logit_foc: first-order conditions and closed forms ----------------
cases <- list(list(delta = c(1, 1), alpha = 1, c = c(1, 1)),
              list(delta = c(4, 3), alpha = 1, c = c(2, 2.5)),
              list(delta = c(2, 2, 1.5), alpha = 0.7, c = c(1, 1.2, 0.8)))
for (k in cases) {
  pb <- bertrand_logit(k$delta, k$alpha, k$c); sb <- logit_share(pb, k$delta, k$alpha)
  pj <- joint_max_logit(k$delta, k$alpha, k$c); sj <- logit_share(pj, k$delta, k$alpha)
  stopifnot(max(abs(pb - k$c - 1 / (k$alpha * (1 - sb)))) < 1e-10,          # Bertrand: p - c = 1 / (alpha (1 - s))
            max(abs(pj - k$c - 1 / (k$alpha * (1 - sum(sj))))) < 1e-10,     # joint: equal markups 1 / (alpha s0)
            all(pj > pb),
            max(abs(logit_foc(pb, k$delta, k$alpha, k$c, kappa = 0))) < 1e-10,
            max(abs(logit_foc(pj, k$delta, k$alpha, k$c, kappa = 1))) < 1e-10)
  # each Bertrand price is a best response: no unilateral deviation on a fine grid around it pays
  for (j in seq_along(pb)) {
    prof_j <- function(x) { p <- pb; p[j] <- x; (x - k$c[j]) * logit_share(p, k$delta, k$alpha)[j] }
    grid <- pb[j] + seq(-1, 1, by = 0.01)
    stopifnot(prof_j(pb[j]) >= max(sapply(grid, prof_j)) - 1e-12)
  }
  # the joint price vector maximises total profit (gradient zero, and no better point near it)
  tot <- function(p) sum((p - k$c) * logit_share(p, k$delta, k$alpha))
  g <- numDeriv::grad(tot, pj)
  stopifnot(max(abs(g)) < 1e-7, tot(pj) > tot(pb))
  set.seed(1); for (r in 1:200) stopifnot(tot(pj) >= tot(pj + rnorm(length(pj), 0, 0.05)))
}
# one product: the Bertrand and joint prices are both the monopoly price
pm1 <- bertrand_logit(2, 1, 1); stopifnot(abs(pm1 - joint_max_logit(2, 1, 1)) < 1e-10)
# kappa between 0 and 1 gives prices between the Bertrand and the joint prices
pk <- solve_logit(c(4, 3), 1, c(2, 2.5), kappa = 0.5)
stopifnot(all(pk > bertrand_logit(c(4, 3), 1, c(2, 2.5))), all(pk < joint_max_logit(c(4, 3), 1, c(2, 2.5))))

## ---- passthrough_logit: implicit function theorem against re-solving; closed forms ----------------------
# three duopolies, two of them with alpha != 1 (so a misplaced alpha in the Jacobians cannot pass)
duos <- list(list(delta = c(4, 3), alpha = 1, c = c(2, 2.5)), list(delta = c(2.5, 3.2), alpha = 0.6, c = c(1.5, 2)),
             list(delta = c(6, 5), alpha = 1.8, c = c(1, 1.4)))
for (du in duos) for (kap in c(0, 1)) {
  d <- du$delta; a <- du$alpha; cc <- du$c
  PT <- passthrough_logit(d, a, cc, kappa = kap)
  h <- 1e-5; fd <- matrix(NA_real_, 2, 2)
  for (j in 1:2) {
    up <- cc; up[j] <- up[j] + h; dn <- cc; dn[j] <- dn[j] - h
    fd[, j] <- (solve_logit(d, a, up, kap) - solve_logit(d, a, dn, kap)) / (2 * h)
  }
  stopifnot(max(abs(PT - fd)) < 1e-7, all(diag(PT) > 0), all(diag(PT) < 1))
  if (kap == 0) stopifnot(all(PT[row(PT) != col(PT)] > 0))        # Bertrand: prices are strategic complements
  # joint maximisation: every markup equals m with alpha m = 1 + sum_j exp(delta_j - alpha c_j - alpha m), so
  # dm/dc_j = -s_j and PT[k, j] = 1{k = j} - s_j: the cross pass-through is negative
  if (kap == 1) {
    sj <- logit_share(joint_max_logit(d, a, cc), d, a)
    stopifnot(max(abs(PT - (diag(2) - matrix(sj, 2, 2, byrow = TRUE)))) < 1e-10, all(PT[row(PT) != col(PT)] < 0))
  }
}
# single-product logit monopoly: PT = 1 - s
for (d in c(-1, 0, 2, 4)) {
  p1 <- bertrand_logit(d, 1.3, 1); s1 <- logit_share(p1, d, 1.3)
  stopifnot(abs(passthrough_logit(d, 1.3, 1, kappa = 0) - (1 - s1)) < 1e-9)
}
# symmetric products under joint ownership: a common cost shock passes through at the outside share s0
dS <- c(2, 2, 2); pS <- joint_max_logit(dS, 1, c(1, 1, 1)); s0 <- 1 - sum(logit_share(pS, dS, 1))
stopifnot(max(abs(rowSums(passthrough_logit(dS, 1, c(1, 1, 1), kappa = 1)) - s0)) < 1e-9)

## ---- price_leader: Stackelberg in prices, product 1 leads ----------------------------------------------
for (k in c(cases[1:2], duos[2:3])) {                                # alpha = 1, 1, 0.6, 1.8
  pl <- price_leader(k$delta, k$alpha, k$c); pb <- bertrand_logit(k$delta, k$alpha, k$c)
  # the follower's price is a best response to the leader's (Bertrand FOC of product 2 at the leader's price)
  stopifnot(abs(logit_foc(pl, k$delta, k$alpha, k$c, kappa = 0)[2]) < 1e-10)
  # the leader's price maximises its profit along the follower's best response (brute force with optimize)
  br2 <- function(p1) uniroot(function(p2) logit_foc(c(p1, p2), k$delta, k$alpha, k$c, kappa = 0)[2],
                              c(k$c[2], k$c[2] + 50 / k$alpha), tol = 1e-13)$root
  lead_profit <- function(p1) (p1 - k$c[1]) * logit_share(c(p1, br2(p1)), k$delta, k$alpha)[1]
  opt <- optimize(lead_profit, c(k$c[1], k$c[1] + 20 / k$alpha), maximum = TRUE, tol = 1e-10)
  stopifnot(abs(opt$maximum - pl[1]) < 1e-5, lead_profit(pl[1]) >= opt$objective - 1e-12)
  # both prices exceed the simultaneous (Bertrand-Nash) prices; both firms earn more than at Bertrand
  prof <- function(p) (p - k$c) * logit_share(p, k$delta, k$alpha)
  stopifnot(all(pl > pb), all(prof(pl) > prof(pb)))
}

## ---- opt_price and bundle_opt: empirical optimal price, revenue, deadweight loss ----------------------
stopifnot(identical(opt_price(c(1, 2, 3)), 2), opt_price(c(1, 2, 3), cost = 1.5) == 3, opt_price(c(5, 5, 1)) == 5)
b0 <- bundle_opt(c(1, 2, 3))
stopifnot(b0$price == 2, abs(b0$rev - 4/3) < 1e-15, abs(b0$dwl - 1/3) < 1e-15)
# legacy 6.5: two merchant types valuing (A, B) at (10, 2) and (2, 10); the pure bundle sells at 12 to both
bl <- bundle_opt(matrix(c(10, 2, 2, 10), nrow = 2, byrow = TRUE))
stopifnot(bl$price == 12, bl$rev == 12, bl$dwl == 0)
# separate selling of good A to the same two types: price 10 to one of two (revenue 5 per consumer), 2 lost
sa <- bundle_opt(c(10, 2)); stopifnot(sa$price == 10, sa$rev == 5, sa$dwl == 1)
# U[0, 1] values on a fine grid: price 1/2, revenue 1/4, deadweight loss 1/8
ug <- bundle_opt((seq_len(1e5) - 0.5) / 1e5)
stopifnot(abs(ug$price - 0.5) < 1e-4, abs(ug$rev - 0.25) < 1e-5, abs(ug$dwl - 0.125) < 1e-4)
# ties: the price counts every consumer at or above it
bt <- bundle_opt(c(2, 2, 2, 1)); stopifnot(bt$price == 2, bt$rev == 1.5, bt$dwl == 0.25)

## ---- bbm_segment: the greedy construction of Bergemann, Brooks and Morris (2015) -------------------------
check_bbm <- function(values, probs) {
  b <- bbm_segment(values, probs); M <- b$segments$mass
  unif <- max(sapply(values, function(p) p * sum(probs[values >= p])))
  pstar <- values[which.max(sapply(values, function(p) p * sum(probs[values >= p])))]
  stopifnot(max(abs(colSums(M) - probs)) < 1e-12, all(M >= 0), abs(sum(b$segments$weight) - sum(probs)) < 1e-12,
            abs(b$uniform_profit - unif) < 1e-12, abs(b$profit - b$uniform_profit) < 1e-10,
            abs(b$total - sum(values * probs)) < 1e-12, abs(b$efficient_total - sum(values * probs)) < 1e-12)
  for (r in seq_len(nrow(M))) {
    S <- which(M[r, ] > 0); w <- sum(M[r, ])
    rev <- sapply(values[S], function(p) p * sum(M[r, values >= p]))
    stopifnot(max(abs(rev - rev[1])) < 1e-10 * max(1, rev[1]),                     # extremal: indifferent on its support
              all(sapply(values, function(p) p * sum(M[r, values >= p])) <= rev[1] + 1e-10 * max(1, rev[1])),
              b$segments$price[r] == values[S[1]],                                    # priced at its lowest value
              pstar %in% values[S])                                                   # the uniform price stays optimal
  }
  b
}
check_bbm(c(1, 2, 3), rep(1/3, 3))
check_bbm(c(1, 2, 3, 4, 5), c(0.1, 0.3, 0.25, 0.2, 0.15))
check_bbm(c(2, 3, 7, 10), c(0.4, 0.1, 0.3, 0.2))
b1 <- check_bbm(5, 1); stopifnot(nrow(b1$segments$mass) == 1, b1$profit == 5)        # one value: one segment
# the test example by hand: segments of weight 2/3, 1/6, 1/6 priced at 1, 2, 2
bh <- bbm_segment(c(1, 2, 3), rep(1/3, 3))
stopifnot(max(abs(bh$segments$weight - c(2/3, 1/6, 1/6))) < 1e-12, identical(bh$segments$price, c(1, 2, 2)))

## ---- linear_pd: third-degree discrimination with linear segment demands (output test) ------------------
lp <- linear_pd(a = c(12, 8), b = c(1, 1), cost = 2)                   # legacy 6.2 segments
stopifnot(abs(lp$p_uniform - 6) < 1e-12, max(abs(lp$p_pd - c(7, 5))) < 1e-12,
          abs(lp$Q_uniform - 8) < 1e-12, abs(lp$Q_pd - 8) < 1e-12,
          abs(lp$profit_uniform - 32) < 1e-12, abs(lp$profit_pd - 34) < 1e-12,
          abs(lp$CS_uniform - 20) < 1e-12, abs(lp$CS_pd - 17) < 1e-12,
          abs(lp$W_uniform - 52) < 1e-12, abs(lp$W_pd - 51) < 1e-12, lp$n_served_uniform == 2)
lo <- linear_pd(a = c(12, 4), b = c(1, 1), cost = 2)                   # the weak segment is not served at one price
stopifnot(abs(lo$p_uniform - 7) < 1e-12, lo$n_served_uniform == 1, abs(lo$Q_uniform - 5) < 1e-12,
          abs(lo$Q_pd - 6) < 1e-12, abs(lo$W_uniform - 37.5) < 1e-12, abs(lo$W_pd - 39) < 1e-12)
# brute force on a price grid for an uneven three-segment market
l3 <- linear_pd(a = c(20, 14, 9), b = c(2, 1, 0.5), cost = 3)
grid <- seq(3, 20, by = 1e-4)
prof <- sapply(grid, function(p) (p - 3) * sum(pmax(0, (c(20, 14, 9) - p) / c(2, 1, 0.5))))
stopifnot(abs(l3$p_uniform - grid[which.max(prof)]) < 2e-4, abs(l3$profit_uniform - max(prof)) < 1e-6)
