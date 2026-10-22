# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of games_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# The first block is the chapter card's assertions, verbatim (comments translated, Ruling P6);
# the rest checks each function against closed forms (2x2 indifference conditions, the logit
# fixed point, lm on each replication, Rubinstein's formula, the regret-matching switching rule step by
# step, the grim-trigger threshold, the sigma = 0 collusive boundary, the folk-theorem polygons, the
# global-game indifference condition integrated numerically) and against small cases whose answer is known.
# ------------------------------------------------------------------
source("games_fns.R")

## ---- chapter card (verbatim) -------------------------------------------------------
A <- matrix(c(10, 2, 14, 6), 2, byrow = TRUE)
stopifnot(identical(nash_2x2(A, t(A))$pure, list(c(2L, 2L))))                     # (Low, Low)
B1 <- matrix(c(3, 0, 0, 1), 2, byrow = TRUE); B2 <- matrix(c(1, 0, 0, 3), 2, byrow = TRUE)
stopifnot(max(abs(nash_2x2(B1, B2)$mixed - c(3/4, 1/4))) < 1e-12)
stopifnot(max(abs(logit_qre(B1, B2, lambda = 1e-8) - 0.5)) < 1e-6)                # precision -> 0: uniform
stopifnot(abs(rubinstein_share(0.9, 0.9) - 1 / 1.9) < 1e-12)
stopifnot(ce_violation(matrix(c(.5, 0, 0, .5), 2), B1, B2) <= 1e-12)              # a public coin is a correlated equilibrium

## ---- nash_2x2: pure equilibria by best-response check, mixed by indifference ------------
stopifnot(identical(nash_2x2(B1, B2)$pure, list(c(1L, 1L), c(2L, 2L))),          # compatibility game: (S1, S1), (S2, S2)
          length(nash_2x2(A, t(A))$mixed) == 0L)                                    # dominance: no interior mixed equilibrium
MP <- matrix(c(1, -1, -1, 1), 2, byrow = TRUE)                                      # matching pennies (Osprey matches)
ne_mp <- nash_2x2(MP, -MP)
stopifnot(length(ne_mp$pure) == 0L, max(abs(ne_mp$mixed - 0.5)) < 1e-15)
C1 <- matrix(c(5, 2, 6, 0), 2, byrow = TRUE)                                        # Chicken: rows Yield, Enter
ne_ch <- nash_2x2(C1, t(C1))
stopifnot(identical(ne_ch$pure, list(c(1L, 2L), c(2L, 1L))), max(abs(ne_ch$mixed - 2/3)) < 1e-15)

## ---- logit_qre and qre_all: the logit fixed point -------------------------------------------
# log(p / (1 - p)) = lambda (u1(S1, q) - u1(S2, q)); here u1(S1, q) - u1(S2, q) = 3q - (1 - q) and
# u2(S1, p) - u2(S2, p) = p - 3(1 - p)
for (lam in c(0.5, 2)) for (st in list(c(.5, .5), c(.95, .9), c(.05, .05))) {
  r <- logit_qre(B1, B2, lambda = lam, start = st)
  stopifnot(identical(names(r), c("p", "q")),
            abs(qlogis(r[["p"]]) - lam * (3 * r[["q"]] - (1 - r[["q"]]))) < 1e-9,
            abs(qlogis(r[["q"]]) - lam * (r[["p"]] - 3 * (1 - r[["p"]]))) < 1e-9)
}
stopifnot(nrow(qre_all(B1, B2, lambda = 0.8)) == 1L)                              # low precision: unique QRE
qa <- qre_all(B1, B2, lambda = 2)
stopifnot(nrow(qa) == 3L, identical(colnames(qa), c("p", "q")), all(diff(qa[, "p"]) > 0))
for (i in 1:3) stopifnot(max(abs(logit_qre(B1, B2, lambda = 2, start = qa[i, ]) - qa[i, ])) < 1e-8)
qa200 <- qre_all(B1, B2, lambda = 200)                                             # high precision: the Nash limits
stopifnot(nrow(qa200) == 3L, qa200[1, "p"] < 0.01, qa200[3, "p"] > 0.99,
          max(abs(qa200[2, ] - c(3/4, 1/4))) < 0.01)

## ---- qre_logodds_ols: noise-free frequencies give back the payoffs exactly -------------------
xs <- seq(-1, 1, by = 0.25)
eqs <- t(vapply(xs, function(x) logit_qre(B1 + rbind(c(x, x), 0), B2 + cbind(c(x, x), 0), lambda = 0.8), numeric(2)))
o <- qre_logodds_ols(data.frame(p = eqs[, 1], q = eqs[, 2], x = xs, chi = 0.8))
stopifnot(identical(names(o$coef), c("u11", "u22")), max(abs(o$coef - c(3, 1))) < 1e-8,
          is.null(o$ci), is.null(o$ci_u22), o$n_rep == 1L, o$n_groups == length(xs))     # one sample: no interval
stopifnot(inherits(try(qre_logodds_ols(data.frame(p = c(0, .5, .7), q = c(.2, .5, .6), x = 0, chi = 1)), silent = TRUE),
                   "try-error"))                                                   # a frequency of 0 or 1 has no log-odds
# several replications: one OLS fit each (checked against lm), coef = their mean, ci = the 2.5% to 97.5% range
set.seed(15)
rd <- do.call(rbind, lapply(1:40, function(r) {
  q <- runif(6, 0.1, 0.9); x <- runif(6, -1, 1); y <- 3 * q - (1 - q) + rnorm(6, sd = 0.05)
  data.frame(rep = r, p = plogis(0.8 * (y + x)), q = q, x = x, chi = 0.8)
}))
orr <- qre_logodds_ols(rd)
ref <- t(vapply(split(rd, rd$rep), function(g) {
  b <- coef(lm(I(qlogis(p) / chi - x) ~ 0 + q + I(-(1 - q)), data = g)); c(u11 = b[[1]], u22 = b[[2]])
}, numeric(2)))
stopifnot(orr$n_rep == 40L, all(orr$n_groups == 6), max(abs(orr$est - ref)) < 1e-10,
          max(abs(orr$coef - colMeans(ref))) < 1e-10,
          max(abs(orr$ci - quantile(ref[, "u11"], c(0.025, 0.975), names = FALSE))) < 1e-12,
          max(abs(orr$ci_u22 - quantile(ref[, "u22"], c(0.025, 0.975), names = FALSE))) < 1e-12,
          identical(names(orr$ci), c("lo", "hi")), max(abs(orr$sd - apply(ref, 2, sd))) < 1e-10)

## ---- rubinstein_share: proposer's share (1 - delta2) / (1 - delta1 delta2) -------------------
stopifnot(abs(rubinstein_share(0.8, 0.8) - 1 / 1.8) < 1e-14,
          abs(rubinstein_share(0.9, 0.5) - 0.5 / 0.55) < 1e-14,
          max(abs(rubinstein_share(c(0.8, 0.99), c(0.8, 0.99)) - 1 / (1 + c(0.8, 0.99)))) < 1e-14)
# Binmore, Rubinstein and Wolinsky: delta_i = exp(-r_i dt), dt -> 0 gives r2 / (r1 + r2)
stopifnot(abs(rubinstein_share(exp(-0.1 * 1e-6), exp(-0.05 * 1e-6)) - 1 / 3) < 1e-6)

## ---- ce_violation: the largest gain from disobeying a recommendation ----------------------------
stopifnot(ce_violation(outer(c(3/4, 1/4), c(1/4, 3/4)), B1, B2) < 1e-12,       # the mixed Nash equilibrium is a CE
          abs(ce_violation(matrix(c(0, 0, 1, 0), 2), B1, B2) - 1) < 1e-12)       # all mass on (S1, S2): gain 1
mu3 <- matrix(c(1, 1, 1, 0) / 3, 2, byrow = TRUE)                                    # Chicken: 1/3 on (Y,Y), (Y,E), (E,Y)
stopifnot(ce_violation(mu3, C1, t(C1)) < 1e-12)

## ---- regret_matching and fictitious_play ----------------------------------------------------------
mu <- regret_matching(B1, B2, T = 5000, seed = 11)
stopifnot(is.matrix(mu), identical(dim(mu), c(2L, 2L)), abs(sum(mu) - 1) < 1e-12, all(mu >= 0),
          identical(mu, regret_matching(B1, B2, T = 5000, seed = 11)))
a1 <- attr(mu, "a1"); a2 <- attr(mu, "a2")
stopifnot(length(a1) == 5000L, all(a1 %in% 1:2), all(a2 %in% 1:2),
          abs(mu[1, 2] - mean(a1 == 1L & a2 == 2L)) < 1e-15, abs(mu[2, 1] - mean(a1 == 2L & a2 == 1L)) < 1e-15)
set.seed(13); rr <- regret_matching_reps(C1, t(C1), T = 2000, R = 5, checkpoints = c(10L, 500L, 2000L))
stopifnot(nrow(rr) == 5L, identical(dim(attr(rr, "path")), c(5L, 3L)),
          max(abs(attr(rr, "path")[, 3] - rr$viol)) < 1e-12,                         # the last checkpoint is T
          all(abs(rowSums(rr[, c("mu11", "mu12", "mu21", "mu22")]) - 1) < 1e-12))
# The switching rule itself, one step at a time (compatibility game; mu_c = 2 max|u| + 1 = 7). After period t,
# a player that played j switches with probability max(D(j), 0) / (t mu_c), D(j) the summed gain, over the
# periods in which it played j, from having played the other action:
#   t = 1 from (S1, S2): Osprey's gain from S2 is A[2,2] - A[1,2] = 1, Boreal's from S1 is B[1,1] - B[1,2] = 1:
#     1/7 each; from (S2, S1): A[1,1] - A[2,1] = 3 and B[2,2] - B[2,1] = 3: 3/7 each; on (S1, S1) or (S2, S2): 0.
#   t = 2 after (S1, S2) twice: D = 2, so 2 / (2 x 7) = 1/7 (without the 1/t averaging it would be 2/7);
#     after (S2, S1) twice: 6 / 14 = 3/7 (6/7 without it);
#     after (S2, S1) then (S2, S2): Osprey's D(S2) = 3 + (A[1,2] - A[2,2]) = 2, so 2/14 = 1/7; Boreal's D(S2) =
#     B[2,1] - B[2,2] = -3, so 0.
set.seed(16); h <- .rm_core(B1, B2, T = 3L, R = 400000L, keep = TRUE)
stopifnot(h$mu_c == 7, identical(dim(h$a1), c(3L, 400000L)))
sw <- function(t, path1, path2) {                    # switch shares at t + 1 among runs whose first t profiles are given
  s <- rep(TRUE, ncol(h$a1))
  for (k in seq_len(t)) s <- s & h$a1[k, ] == path1[k] & h$a2[k, ] == path2[k]
  c(osprey = mean(h$a1[t + 1, s] != path1[t]), boreal = mean(h$a2[t + 1, s] != path2[t]), n = sum(s))
}
chk <- function(r, want) all(r[["n"]] > 10000, abs(r[c("osprey", "boreal")] - want) < 0.01)
stopifnot(chk(sw(1, 1, 2), c(1, 1) / 7), chk(sw(1, 2, 1), c(3, 3) / 7),
          all(sw(1, 1, 1)[1:2] == 0), all(sw(1, 2, 2)[1:2] == 0),
          chk(sw(2, c(1, 1), c(2, 2)), c(1, 1) / 7), chk(sw(2, c(2, 2), c(1, 1)), c(3, 3) / 7),
          chk(sw(2, c(2, 2), c(1, 2)), c(1, 0) / 7))
D <- matrix(c(2, 2, 0, 0), 2, byrow = TRUE)                                          # action 1 strictly dominant for both
muD <- regret_matching(D, t(D), T = 5000, seed = 12)
stopifnot(muD[1, 1] > 0.99, ce_violation(muD, D, t(D)) < 0.01)
fp <- fictitious_play(MP, -MP, T = 20000, seed = 3)
stopifnot(length(fp$a1) == 20000L, max(abs(fp$freq - 0.5)) < 0.02,
          identical(fp, fictitious_play(MP, -MP, T = 20000, seed = 3)))

## ---- gg_threshold: the threshold signal of the global game ---------------------------------------
stopifnot(all(gg_threshold(c(0.3, 0.1, 0.01)) == 0.5))                             # flat prior: Laplacian belief, 1/2
th <- gg_threshold(c(0.1, 0.05, 0.01), prior_mean = 0.25, prior_sd = 0.5)
# indifference of the marginal player, E[l | x*] by numerical integration over the posterior of theta
indiff <- function(x, s, y = 0.25, tau = 0.5) {
  a <- 1 / tau^2; b <- 1 / s^2; m <- (a * y + b * x) / (a + b); v <- 1 / (a + b)
  El <- integrate(function(t) pnorm((t - x) / s) * dnorm(t, m, sqrt(v)), m - 12 * sqrt(v), m + 12 * sqrt(v),
                  rel.tol = 1e-12)$value
  m + El - 1
}
stopifnot(max(abs(mapply(indiff, th, c(0.1, 0.05, 0.01)))) < 1e-9, all(diff(abs(th - 0.5)) < 0), all(th > 0.5))
stopifnot(identical(gg_roots(0.1, prior_mean = 0.5, prior_sd = 0.5), 0.5))           # symmetric prior: exactly 1/2 (a grid point)
r3 <- gg_roots(0.1, prior_mean = 0.5, prior_sd = 0.1)                                # public as precise as private: three
stopifnot(length(r3) == 3L, abs(r3[2] - 0.5) < 1e-12,
          inherits(try(gg_threshold(0.1, prior_mean = 0.5, prior_sd = 0.1), silent = TRUE), "try-error"))
sm <- gg_sigma_max(prior_sd = 0.5)                                                   # uniqueness bound on the private noise
stopifnot(sm > 0.1, abs(gg_min_slope(sm, prior_sd = 0.5)) < 1e-9,
          gg_min_slope(0.9 * sm, prior_sd = 0.5) > 0, gg_min_slope(1.1 * sm, prior_sd = 0.5) < 0)

## ---- Green-Porter: values, equilibrium quantity, one-shot deviation check, simulation -------------
# sigma = 0 and an (almost) endless punishment is grim trigger with Cournot reversion: pi_m = 8, pi_d = 9,
# pi_N = 64/9, so the check switches at delta* = (9 - 8) / (9 - 64/9) = 9/17
gpar <- list(a = 10, cost = 2, n_pun = 1e6, sigma = 0, q_c = 2, trigger = 6)
dstar <- (9 - 8) / (9 - 64 / 9)
stopifnot(abs(dstar - 9 / 17) < 1e-12, osd_check(modifyList(gpar, list(delta = dstar + 1e-4))),
          !osd_check(modifyList(gpar, list(delta = dstar - 1e-4))))
# the Bellman equations behind gp_values
p2 <- list(a = 10, cost = 2, delta = 0.95, n_pun = 6, sigma = 0.2, trigger = 5.4, q_c = 2.1)
v <- gp_values(p2)
stopifnot(abs(v$beta - pnorm((5.4 - (10 - 4.2)) / 0.2)) < 1e-15,
          abs(v$Vc - ((10 - 4.2 - 2) * 2.1 + 0.95 * ((1 - v$beta) * v$Vc + v$beta * v$Vp))) < 1e-9,
          abs(v$Vp - (sum(0.95^(0:5)) * 64 / 9 + 0.95^6 * v$Vc)) < 1e-9)
# the collusive quantity for a given trigger: the first-order condition holds and no one-shot deviation pays
e <- gp_equilibrium(modifyList(p2, list(q_c = NULL)))
ve <- gp_values(e); z <- (e$trigger - (10 - 2 * e$q_c)) / 0.2
stopifnot(e$q_c > 2, e$q_c < 8 / 3, osd_check(e),
          abs((8 - 3 * e$q_c) - 0.95 * (ve$Vc - ve$Vp) * dnorm(z) / 0.2) < 1e-8,
          attr(osd_check(e), "gain_punish") <= 1e-9)
stopifnot(!osd_check(modifyList(e, list(delta = 0.5))))                             # impatient firms deviate
# sigma = 0: the most collusive quantity is the monopoly split, the trigger is the collusive price itself
e0 <- gp_optimal_trigger(list(a = 10, cost = 2, delta = 0.95, n_pun = 6, sigma = 0))
stopifnot(abs(e0$q_c - 2) < 1e-12, abs(e0$trigger - 6) < 1e-12, osd_check(e0))
# sigma = 0 when the monopoly split fails (delta = 0.5): the binding deviation is the static best response
# (8 - q)/2, detected for sure. With pi(q) = (8 - 2q) q, sum(0.5^(0:5)) = 63/32 and 0.5^6 = 1/64, the
# constraint 2 pi(q) = (8 - q)^2 / 4 + 0.5 (64/9) (63/32) + (1/64) pi(q) is 135 q^2 - 636 q + 736 = 0, whose
# smaller root is 92/45 (the larger is Cournot's 8/3)
e05 <- gp_optimal_trigger(list(a = 10, cost = 2, delta = 0.5, n_pun = 6, sigma = 0))
stopifnot(abs(e05$q_c - 92 / 45) < 1e-10, abs(e05$trigger - (10 - 2 * 92 / 45)) < 1e-10, osd_check(e05),
          !osd_check(modifyList(e05, list(q_c = 92 / 45 - 1e-4, trigger = 10 - 2 * (92 / 45 - 1e-4)))),
          abs((1 - 0.5) * gp_values(e05)$Vc - (8 - 2 * 92 / 45) * 92 / 45) < 1e-10)
# no trigger supports collusion (sigma = 0.4, delta = 0.5): Cournot forever, no trigger, no wars
e_none <- gp_optimal_trigger(list(a = 10, cost = 2, delta = 0.5, n_pun = 6, sigma = 0.4))
stopifnot(abs(e_none$q_c - 8 / 3) < 1e-12, e_none$trigger == -Inf, gp_values(e_none)$beta == 0,
          abs((1 - 0.5) * gp_values(e_none)$Vc - 64 / 9) < 1e-12, osd_check(e_none),
          !any(gp_simulate(e_none, T = 20, K = 3, seed = 2)$war))

## ---- folk_sets: the pricing game's feasible set and its individually rational part -------------
fs <- folk_sets(A, t(A))
same_set <- function(P, Q) nrow(P) == nrow(Q) && all(apply(Q, 1, function(v) any(abs(P[, 1] - v[1]) + abs(P[, 2] - v[2]) < 1e-12)))
stopifnot(identical(fs$minmax, c(6, 6)), same_set(fs$hull, rbind(c(10, 10), c(2, 14), c(14, 2), c(6, 6))),
          same_set(fs$ir, rbind(c(6, 6), c(12, 6), c(10, 10), c(6, 12))), all(fs$ir[1, ] == c(6, 6)),
          abs(fs$area_hull - 48) < 1e-12, abs(fs$area_ir - 24) < 1e-12)
s0 <- gp_simulate(e0, T = 50, K = 3, seed = 1)
stopifnot(data.table::is.data.table(s0), nrow(s0) == 150L, identical(names(s0), c("market", "t", "war", "price", "trigger_hit")),
          !any(s0$war), all(s0$price == 6))
s1 <- gp_simulate(e, T = 400, K = 30, seed = 5)
stopifnot(identical(s1, gp_simulate(e, T = 400, K = 30, seed = 5)), any(s1$war), any(!s1$war))
# every war that starts inside the window follows a collusive period whose price fell below the trigger,
# and lasts n_pun periods unless the window ends first
for (k in unique(s1$market)) {
  w <- s1$war[s1$market == k]; hit <- s1$trigger_hit[s1$market == k]
  starts <- which(w & !c(TRUE, head(w, -1)))
  for (st in starts[starts > 1]) {
    stopifnot(hit[st - 1], !w[st - 1])
    run <- rle(w[st:length(w)])$lengths[1]
    stopifnot(run == e$n_pun || st + run - 1 == length(w))
  }
  stopifnot(!any(hit & w))                                                          # the trigger is only read in collusion
}
