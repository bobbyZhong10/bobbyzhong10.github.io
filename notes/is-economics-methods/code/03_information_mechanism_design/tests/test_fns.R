# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of mech_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# The first block is the chapter card's assertions, verbatim (comments translated, Ruling P6);
# the rest checks each function against closed forms (uniform first-price bids with and without a
# reserve, the uniform virtual value and its optimal menu, the two-type screening solution, the LEN
# slope with one and two signals) and against small cases whose answer is known (stable matchings
# enumerated by brute force, the pool-adjacent-violators algorithm on hand examples), and the ironing
# against its defining conditions (weighted means on bunches, the convex minorant, the virtual surplus).
# ------------------------------------------------------------------
source("mech_fns.R")

## ---- chapter card (verbatim) -------------------------------------------------------
stopifnot(abs(fpa_bid(0.6, 2, punif) - 0.3) < 1e-8)                          # U[0,1], n = 2: b(v) = v/2
stopifnot(abs(fpa_bid(0.6, 3, punif) - 0.4) < 1e-8)                          # n = 3: b(v) = 2v/3
pp <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3)); pr <- list(c(2, 1, 3), c(1, 2, 3), c(3, 2, 1))
m <- deferred_acceptance(pp, pr); stopifnot(is_stable(m, pp, pr))
stopifnot(abs(len_bstar(r = 1, sigma2 = 0, c2 = 1) - 1) < 1e-12)              # no noise: first-best slope b = 1

## ---- fpa_bid: b(v) = v - int_r^v F^(n-1) / F(v)^(n-1), uniform closed forms -----------------
vs <- c(0.05, 0.3, 0.6, 0.95, 1)
for (n in 2:6) stopifnot(max(abs(fpa_bid(vs, n, punif) - (n - 1) / n * vs)) < 1e-9)
# with a reserve r on U[0,1]: b(v) = v - (v^n - r^n) / (n v^(n-1)) for v >= r; no bid below r
for (n in 2:4) {
  r <- 0.5; v_hi <- c(0.5, 0.7, 0.9)
  stopifnot(max(abs(fpa_bid(v_hi, n, punif, reserve = r) - (v_hi - (v_hi^n - r^n) / (n * v_hi^(n - 1))))) < 1e-9,
            all(is.na(fpa_bid(c(0.1, 0.49), n, punif, reserve = r))), fpa_bid(r, n, punif, reserve = r) == r)
}
stopifnot(fpa_bid(0, 2, punif) == 0, length(fpa_bid(numeric(0), 2, punif)) == 0L)
# b(v) = E[v_{-i,(1)} | v_{-i,(1)} < v] for the truncated normal of block V1 (direct integral of the rival maximum)
D1 <- value_dist(list(family = "truncnorm", mean = 1.3, sd = 0.6, lower = 0, upper = 3))
for (n in c(2, 4, 6)) for (v in c(0.4, 1.3, 2.6)) {
  num <- integrate(function(y) y * (n - 1) * D1$cdf(y)^(n - 2) * D1$pdf(y), 0, v, rel.tol = 1e-12)$value
  stopifnot(abs(fpa_bid(v, n, D1$cdf) - num / D1$cdf(v)^(n - 1)) < 1e-8)
}
b_grid <- fpa_bid(seq(0.05, 2.95, by = 0.05), 4, D1$cdf)
stopifnot(all(diff(b_grid) > 0), all(b_grid < seq(0.05, 2.95, by = 0.05)))   # increasing, shaded below value
# fpa_bid_fun (spline of the exact bid on a grid, used by the Monte Carlo) agrees with fpa_bid
for (n in c(2, 6)) for (r in c(0, 1)) {
  bf <- fpa_bid_fun(n, D1, reserve = r)
  vv <- c(r, r + 0.003, seq(r + 0.01, 2.99, length.out = 37), 3)
  stopifnot(max(abs(bf(vv) - fpa_bid(vv, n, D1$cdf, reserve = r))) < 1e-7)
  if (r > 0) stopifnot(all(is.na(bf(c(0.2, r - 1e-9)))))
}

## ---- value_dist: the truncated normal and its virtual value -------------------------------------
stopifnot(abs(D1$cdf(0)) < 1e-15, abs(D1$cdf(3) - 1) < 1e-15, D1$pdf(0) > 0, D1$pdf(3) > 0,       # I.1: interval support, positive density
          abs(integrate(D1$pdf, 0, 3)$value - 1) < 1e-9,
          max(abs(D1$cdf(D1$quantile(c(0.01, 0.3, 0.5, 0.99))) - c(0.01, 0.3, 0.5, 0.99))) < 1e-12)
ps <- virtual_value(seq(0.01, 3, by = 0.01), D1)
stopifnot(all(diff(ps) > 0))                                                                 # regular (log-concave density)
rs <- optimal_reserve(D1)
stopifnot(abs(virtual_value(rs, D1)) < 1e-9, abs(rs - (1 - D1$cdf(rs)) / D1$pdf(rs)) < 1e-9)
DU <- value_dist(list(family = "uniform", lower = 0, upper = 1))
stopifnot(abs(optimal_reserve(DU) - 0.5) < 1e-9, abs(virtual_value(0.8, DU) - 0.6) < 1e-12)
# expected revenue: second-highest value (no reserve) and the virtual-value integral (reserve r)
for (n in 2:4) stopifnot(abs(expected_revenue(n, DU) - (n - 1) / (n + 1)) < 1e-9)
stopifnot(abs(expected_revenue(2, DU, reserve = 0.5) - 5 / 12) < 1e-9)
# auction_outcomes: revenues per auction for the four formats on a value matrix
Vm <- rbind(c(0.9, 0.2, 0.6), c(0.3, 0.4, 0.1), c(0.45, 0.2, 0.3))
ao <- auction_outcomes(Vm, reserve = 0.5, bid_fp = function(v) 2 * v / 3,
                       bid_fp_res = function(v) ifelse(v < 0.5, NA, v - (v^3 - 0.125) / (3 * v^2)))
stopifnot(isTRUE(all.equal(ao$rev_sp, c(0.6, 0.3, 0.3))), isTRUE(all.equal(ao$rev_fp, 2 * c(0.9, 0.4, 0.45) / 3)),
          isTRUE(all.equal(ao$rev_sp_res, c(0.6, 0, 0))),
          isTRUE(all.equal(ao$rev_fp_res, c(0.9 - (0.9^3 - 0.125) / (3 * 0.9^2), 0, 0))))

## ---- gpv_quantile: v(tau) = b(tau) + tau b'(tau) / (n - 1) ------------------------------------
# a deterministic "sample": uniform values at the midpoints of N cells, equilibrium bids 2v/3 (n = 3)
N <- 20000; v_mid <- (seq_len(N) - 0.5) / N
g <- gpv_quantile(2 * v_mid / 3, 3, taus = c(0.05, 0.25, 0.5, 0.75, 0.95))
stopifnot(inherits(g, "data.table"), identical(names(g), c("tau", "b", "v_hat")),
          max(abs(g$b - 2 * g$tau / 3)) < 1e-4, max(abs(g$v_hat - g$tau)) < 1e-4)
g6 <- gpv_quantile(5 * v_mid / 6, 6, taus = seq(0.1, 0.9, by = 0.1), h = 0.05)
stopifnot(max(abs(g6$v_hat - g6$tau)) < 1e-4)
# the finite difference stays inside (0, 1) at the ends
ge <- gpv_quantile(2 * v_mid / 3, 3, taus = c(0.01, 0.99), h = 0.02)
stopifnot(all(is.finite(ge$v_hat)), max(abs(ge$v_hat - ge$tau)) < 1e-3)

## ---- pava and optimal_menu -------------------------------------------------------------------------
stopifnot(isTRUE(all.equal(pava(c(1, 3, 2, 4), c(1, 1, 1, 1)), c(1, 2.5, 2.5, 4), check.attributes = FALSE)),
          isTRUE(all.equal(pava(c(3, 1), c(1, 3)), c(1.5, 1.5), check.attributes = FALSE)),
          identical(attr(pava(c(3, 1, 5), c(1, 1, 1)), "block"), c(1L, 1L, 2L)),
          isTRUE(all.equal(pava(1:5, rep(1, 5)), 1:5, check.attributes = FALSE)))
# uniform types on [2, 4], cost q^2 / 2 * cost2: Psi = 2 theta - 4, q = max(0, 2 theta - 4) / cost2
th <- seq(2, 4, length.out = 401)
mu_ <- optimal_menu(th, rep(1, length(th)), cost2 = 2)
stopifnot(inherits(mu_, "data.table"),
          all(c("theta", "f", "F", "w", "psi", "psi_bar", "q", "q_eff", "U", "t", "bunched") %in% names(mu_)),
          max(abs(mu_$psi - (2 * th - 4))) < 1e-9, max(abs(mu_$q - pmax(0, 2 * th - 4) / 2)) < 1e-9,
          abs(sum(mu_$w) - 1) < 1e-12, !any(mu_$bunched),
          abs(tail(mu_$q, 1) - tail(mu_$q_eff, 1)) < 1e-12)                               # no distortion at the top
# a bimodal type density needs ironing: the ironed allocation is monotone, the raw one is not,
# and the ironed menu is incentive compatible on the grid while the raw one is not
th2 <- seq(2, 10, length.out = 1601)
d2 <- 0.7 * dnorm(th2, 4, 0.7) + 0.3 * dnorm(th2, 8, 0.7)
mi <- optimal_menu(th2, d2, cost2 = 1); mn <- optimal_menu(th2, d2, cost2 = 1, iron = FALSE)
stopifnot(any(diff(mi$psi) < 0), all(diff(mi$q) >= 0), any(diff(mn$q) < 0), any(mi$bunched),
          menu_ic_gain(mi) < 1e-9, menu_ic_gain(mn) > 1,
          abs(tail(mi$q, 1) - tail(mi$q_eff, 1)) < 1e-12)
sm <- menu_summary(mi, cost2 = 1)
stopifnot(sm$bunch_width > 0, sm$bunch_q > 0, sm$theta_exclusion > 2,
          abs(sm$profit - sum(mi$w * (mi$t - mi$q^2 / 2))) < 1e-12, sm$profit < sm$surplus_efficient,
          abs(sm$surplus - sm$profit - sm$rent) < 1e-12)
stopifnot(all(abs(diff(mi$q[mi$bunched & mi$q > 0])) < 1e-12))                              # flat on the bunch
# the ironing is the right one, not just some monotone fix of Psi (a monotone, IC menu with a bunch passes all
# the checks above). Bunches are read off Psi_bar itself (maximal runs of one value), not off the block column.
# On each bunch: Psi_bar is the density-weighted (w) mean of Psi; every prefix of the bunch has
# sum w (Psi - Psi_bar) >= 0 and the whole bunch 0 (int Psi dF lies above its convex minorant inside the bunch
# and touches it at both ends); so at both ends Psi = Psi_bar up to one grid step: Psi_bar lies between Psi at
# the last point outside and the first point inside the bunch.
K2 <- nrow(mi)
bu <- Filter(function(ii) length(ii) > 1L, split(seq_len(K2), cumsum(c(TRUE, diff(mi$psi_bar) != 0))))
stopifnot(length(bu) >= 1L, length(bu) == length(unique(mi$block[mi$bunched])))             # the block column marks the same bunches
for (ii in bu) {
  lo <- min(ii); hi <- max(ii); pb <- mi$psi_bar[lo]; S <- cumsum(mi$w[ii] * (mi$psi[ii] - pb))
  stopifnot(all(mi$block[ii] == mi$block[lo]),
            abs(pb - sum(mi$w[ii] * mi$psi[ii]) / sum(mi$w[ii])) < 1e-10,                 # Psi_bar = w-weighted mean of Psi
            all(S > -1e-12), abs(S[length(S)]) < 1e-12,                                    # convex-minorant condition
            lo == 1L || mi$psi[lo - 1L] <= pb, pb <= mi$psi[lo],                            # Psi ~ Psi_bar at the lower end
            mi$psi[hi] <= pb, hi == K2 || pb <= mi$psi[hi + 1L])                            # and at the upper end
}
# Psi_bar is the slope of the greatest convex minorant of the points (sum w, sum w Psi): a lower hull computed
# here without pava
X <- c(0, cumsum(mi$w)); Y <- c(0, cumsum(mi$w * mi$psi)); h <- 1L
for (i in 2:length(X)) {
  while (length(h) >= 2L) {
    a <- h[length(h) - 1L]; b <- h[length(h)]
    if ((Y[b] - Y[a]) * (X[i] - X[a]) >= (Y[i] - Y[a]) * (X[b] - X[a])) h <- h[-length(h)] else break
  }
  h <- c(h, i)
}
stopifnot(max(abs((diff(Y[h]) / diff(X[h]))[findInterval(seq_len(K2), h)] - mi$psi_bar)) < 1e-8)
# the expected virtual surplus with Psi_bar equals the seller's profit from the transfers (on a bunch q is flat
# and Psi_bar averages Psi with the weights w, so ironing moves no virtual surplus), and the ironed menu beats
# another monotone, IC menu: Psi flattened by its running maximum, q = max(0, cummax(Psi)) / cost2, U = int q
cm <- pmax(cummax(mi$psi), 0); Ucm <- c(0, cumsum((cm[-1] + cm[-K2]) / 2 * diff(mi$theta)))
stopifnot(abs(sum(mi$w * (mi$psi_bar * mi$q - mi$q^2 / 2)) - sm$profit) < 1e-5,
          sm$profit > sum(mi$w * (mi$theta * cm - Ucm - cm^2 / 2)) + 1e-8)

## ---- screening_two_type: utility theta q - t, cost cost2 q^2 / 2, share mu of high types ------------
s2t <- screening_two_type(theta_L = 4, theta_H = 7, mu = 0.4, cost2 = 1)
stopifnot(abs(s2t$q_H - 7) < 1e-12, abs(s2t$q_L - 2) < 1e-12, abs(s2t$t_L - 8) < 1e-12, abs(s2t$rent_H - 6) < 1e-12,
          abs(s2t$t_H - 43) < 1e-12, abs(s2t$profit - 11) < 1e-12, abs(s2t$profit_fb - 14.6) < 1e-12,
          abs(s2t$q_L_fb - 4) < 1e-12,
          abs((7 * s2t$q_H - s2t$t_H) - (7 * s2t$q_L - s2t$t_L)) < 1e-12,                     # IC of the high type binds
          abs(4 * s2t$q_L - s2t$t_L) < 1e-12, 4 * s2t$q_H - s2t$t_H < 0)                        # IR of the low type binds; low type does not envy
stopifnot(screening_two_type(4, 7, mu = 0.6, cost2 = 1)$q_L == 0)                                # mu / (1 - mu) (7 - 4) > 4: the low type is excluded

## ---- LEN: len_bstar and len_contract with several signals ---------------------------------------------
stopifnot(abs(len_bstar(0.5, 4, 1) - 1 / 3) < 1e-15, abs(len_bstar(2, 1, 0.5) - 0.5) < 1e-15,
          isTRUE(all.equal(len_bstar(1, c(0, 1, 3), 1), c(1, 1 / 2, 1 / 4))))
k1 <- len_contract(A = 0.5, c2 = 1, Sigma = matrix(4), mu = 1)
stopifnot(abs(k1$beta - 1 / 3) < 1e-15, abs(k1$effort - 1 / 3) < 1e-15, abs(k1$profit - 1 / 6) < 1e-15,
          abs(k1$risk_premium - 0.5 * (1 / 9) * 4 / 2) < 1e-15, abs(k1$b - 1 / 3) < 1e-15)
# an independent signal y = -0.5 e + eta, Var(eta) = 1: sigma_s^2 = 1 / (1/4 + 0.25) = 2, beta = 1/2
k2 <- len_contract(A = 0.5, c2 = 1, Sigma = diag(c(4, 1)), mu = c(1, -0.5))
stopifnot(abs(k2$sigma_s2 - 2) < 1e-12, abs(k2$beta - 0.5) < 1e-12, max(abs(k2$b - c(0.25, -0.5))) < 1e-12,
          abs(sum(k2$b * c(1, -0.5)) - k2$beta) < 1e-12, abs(k2$profit - 0.25) < 1e-12)
# pure noise (no effort sensitivity, independent): weight 0, nothing changes
k0 <- len_contract(A = 0.5, c2 = 1, Sigma = diag(c(4, 1)), mu = c(1, 0))
stopifnot(abs(k0$beta - 1 / 3) < 1e-12, abs(k0$b[2]) < 1e-12, abs(k0$profit - 1 / 6) < 1e-12)
# a peer benchmark: no effort sensitivity but correlated noise (rho = 0.6): relative performance evaluation
kp <- len_contract(A = 0.5, c2 = 1, Sigma = matrix(c(4, 2.4, 2.4, 4), 2), mu = c(1, 0))
stopifnot(abs(kp$sigma_s2 - 4 * (1 - 0.36)) < 1e-12, abs(kp$b[2] / kp$b[1] + 0.6) < 1e-12, kp$beta > 1 / 3)

## ---- deferred_acceptance and is_stable --------------------------------------------------------------
# the card example: proposers 1 and 2 both rank receiver 1 first; receiver 1 keeps proposer 2, proposer 1
# moves on to receiver 2, which rejects proposer 3, who ends with receiver 3
stopifnot(identical(m, c(2L, 1L, 3L)), is.integer(m))
# legacy 3 x 3 (section 6.5): workers propose -> (1, 2, 3); companies propose -> worker view (2, 1, 3)
wp <- list(c(1, 2, 3), c(2, 1, 3), c(1, 2, 3)); fp <- list(c(2, 1, 3), c(1, 2, 3), c(1, 2, 3))
stopifnot(identical(deferred_acceptance(wp, fp), c(1L, 2L, 3L)),
          identical(invert_match(deferred_acceptance(fp, wp), 3L), c(2L, 1L, 3L)))
# an unstable matching is detected: in (3, 2, 1), worker 1 and company 1 prefer each other to their partners
stopifnot(!is_stable(c(3L, 2L, 1L), wp, fp))
# brute force on random 4 x 4 markets: DA is stable and proposer optimal, receiver-proposing DA is receiver optimal
perm <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perm(v[-i]), function(p) c(v[i], p))))
all_m <- perm(1:4)
set.seed(99)
for (rep in 1:30) {
  P <- replicate(4, sample(4), simplify = FALSE); R <- replicate(4, sample(4), simplify = FALSE)
  st <- Filter(function(x) is_stable(x, P, R), all_m)
  mP <- deferred_acceptance(P, R); mR <- invert_match(deferred_acceptance(R, P), 4L)
  stopifnot(length(st) >= 1L, is_stable(mP, P, R), is_stable(mR, P, R))
  rkP <- partner_rank(mP, P); rkR_recv <- partner_rank(invert_match(mR, 4L), R)
  for (x in st) stopifnot(all(rkP <= partner_rank(x, P)), all(rkR_recv <= partner_rank(invert_match(x, 4L), R)))
}
# incomplete lists: unacceptable partners are never matched, unmatched agents get NA
P3 <- list(c(1L), c(1L, 2L), integer(0)); R3 <- list(c(2L, 1L), c(1L))
m3 <- deferred_acceptance(P3, R3)
stopifnot(identical(m3, c(NA_integer_, 1L, NA_integer_)), is_stable(m3, P3, R3),
          !is_stable(c(2L, 1L, NA), P3, R3),                    # proposer 1 does not find receiver 2 acceptable (IR fails)
          !is_stable(c(NA, NA, NA), P3, R3))                    # proposer 2 and receiver 1 would block
stopifnot(identical(partner_rank(m3, P3), c(Inf, 1, Inf)))
# preference helpers
U <- matrix(c(3, 1, 2, 0, 5, 4), nrow = 2, byrow = TRUE)
stopifnot(identical(prefs_from_utility(U), list(c(1L, 3L, 2L), c(2L, 3L, 1L))),
          identical(truncate_pref(c(4L, 2L, 3L, 1L), 2L), c(4L, 2L)),
          identical(swap_pref(c(4L, 2L, 3L, 1L), c(1L, 3L)), c(3L, 2L, 4L, 1L)),
          identical(invert_match(c(2L, NA, 1L), 3L), c(3L, 1L, NA)))

## ---- V0 helpers: VCG and GSP in the position auction, persuasion, cheap talk, bilateral trade --------
pa <- vcg_position(v = c(10, 6, 4), ctr = c(100, 40))
stopifnot(isTRUE(all.equal(pa$pay, c(520, 160, 0))), isTRUE(all.equal(pa$per_click, c(5.2, 4, 0))),
          pa$revenue == 680, gsp_revenue(c(10, 6, 4), c(100, 40)) == 760)
stopifnot(abs(persuasion_x(0.3, 0.5) - 3 / 7) < 1e-15, abs(persuasion_value(0.3, 0.5) - 0.6) < 1e-15,
          persuasion_value(0.6, 0.5) == 1)
stopifnot(identical(cheap_talk_N(c(0.01, 0.05, 0.1, 0.2, 0.25, 0.4)), c(7L, 3L, 2L, 2L, 1L, 1L)),
          max(abs(cheap_talk_partition(0.05, 3) - c(0, 2 / 15, 7 / 15, 1))) < 1e-15,
          max(abs(diff(diff(cheap_talk_partition(0.05, 3))) - 4 * 0.05)) < 1e-15)
