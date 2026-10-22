# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Osprey-Boreal data and truth02.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("games_fns.R")
tr <- readRDS("truth02.rds"); P <- tr$params
G1 <- readRDS("osprey_G1.rds"); G2 <- readRDS("osprey_G2.rds"); G3 <- readRDS("osprey_G3.rds"); G4 <- readRDS("osprey_G4.rds")
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))

## ---- truth02 fields named on the chapter card --------------------------------------------
stopifnot(identical(tr$G1$payoff, 3), identical(names(tr$G1$lambda), c("unique", "high")),
          identical(unname(tr$G1$lambda), c(0.8, 2)), identical(tr$G2$sigma, c(0, 0.1, 0.2, 0.4)),
          identical(tr$G4$theta_star_limit, 0.5))
# seeds: integers, all different
sd <- unlist(P$seed)
stopifnot(length(sd) == 8L, !anyDuplicated(sd), all(sd == round(sd)))

## ---- the generators reproduce the saved data from the stored seeds ------------------------
set.seed(P$seed$G1); stopifnot(identical(gen_G1(P$G1, tr$G1$eq), G1))
for (l in seq_along(tr$G2$sigma)) {
  sim <- gp_simulate(tr$G2$eq[[l]], T = P$G2$T, K = P$G2$K, seed = P$seed$G2)
  sub <- G2[G2$level == l, ]
  for (v in names(sim)) stopifnot(identical(sim[[v]], sub[[v]]))
}
stopifnot(identical(regret_matching(P$G3$compat_A, P$G3$compat_B, T = P$G3$T, seed = P$seed$G3_rm_compat), G3$rm_compat),
          identical(regret_matching(P$G3$chicken_A, P$G3$chicken_B, T = P$G3$T, seed = P$seed$G3_rm_chicken), G3$rm_chicken),
          identical(fictitious_play(P$G3$mp_A, P$G3$mp_B, T = P$G3$T, seed = P$seed$G3_fp), G3$fp_mp))
set.seed(P$seed$G3_rm_reps)
stopifnot(identical(regret_matching_reps(P$G3$compat_A, P$G3$compat_B, T = P$G3$T, R = P$G3$R, checkpoints = P$G3$checkpoints),
                    G3$reps_compat),
          identical(regret_matching_reps(P$G3$chicken_A, P$G3$chicken_B, T = P$G3$T, R = P$G3$R, checkpoints = P$G3$checkpoints),
                    G3$reps_chicken))
set.seed(P$seed$G4); stopifnot(identical(gen_G4(P$G4, tr$G4$x_star), G4))

## ---- G0: hand block, no random numbers -----------------------------------------------------
stopifnot(is.list(tr$G0), fin(tr$G0$compat_p), fin(tr$G0$grim_delta_star), fin(tr$G0$brw_limit))
# folk-theorem sets of the pricing game (fig_02_folk.svg): feasible area 48, individually rational area 24
stopifnot(identical(tr$G0$folk_minmax, c(6, 6)), nrow(tr$G0$folk_hull) == 4L, nrow(tr$G0$folk_ir) == 4L,
          abs(tr$G0$folk_area_hull - 48) < 1e-12, abs(tr$G0$folk_area_ir - 24) < 1e-12)

## ---- G1: QRE markets -------------------------------------------------------------------------------
P1 <- P$G1
stopifnot(nrow(G1) == P1$n_markets_unique + P1$n_markets_high, P1$n_markets_unique >= 2000,
          identical(names(G1), c("block", "market", "x", "history", "n_dec", "k_osprey", "k_boreal", "p_true", "q_true")),
          all(table(G1$x[G1$block == "unique"]) == P1$n_markets_unique / length(P1$x_unique)),
          all(table(G1$x[G1$block == "high"]) == P1$n_markets_high / length(P1$x_high)),
          all(G1$n_dec == P1$n_dec), all(G1$k_osprey >= 0 & G1$k_osprey <= G1$n_dec),
          all(G1$k_boreal >= 0 & G1$k_boreal <= G1$n_dec),
          all(is.na(G1$history[G1$block == "unique"])), all(G1$history[G1$block == "high"] %in% c("S1", "S2")))
eq <- tr$G1$eq
stopifnot(all(eq$n_qre[eq$block == "unique"] == 1L), all(eq$n_qre[eq$block == "high"] == 3L))
# every market plays a logit QRE of its own game; in the high block, S1-first markets play the QRE with the
# highest p and S2-first markets the one with the lowest p
for (i in seq_len(nrow(eq))) {
  r <- eq[i, ]; lam <- tr$G1$lambda[[r$block]]
  g <- g1_game(P1, r$x)
  stopifnot(abs(qlogis(r$p) - lam * qre_payoff_diff(g$A, g$B, r$p, r$q)[1]) < 1e-9,
            abs(qlogis(r$q) - lam * qre_payoff_diff(g$A, g$B, r$p, r$q)[2]) < 1e-9)
}
hi <- G1[G1$block == "high", ]
for (x in P1$x_high) {
  e <- eq[eq$block == "high" & abs(eq$x - x) < 1e-12, ]
  stopifnot(all(hi$p_true[abs(hi$x - x) < 1e-12 & hi$history == "S1"] == max(e$p)),
            all(hi$p_true[abs(hi$x - x) < 1e-12 & hi$history == "S2"] == min(e$p)))
}

## ---- G2: Green-Porter equilibria and panels ----------------------------------------------------------
P2 <- P$G2
stopifnot(nrow(G2) == length(tr$G2$sigma) * P2$K * P2$T, length(tr$G2$eq) == length(tr$G2$sigma))
e0 <- tr$G2$eq[[1]]
stopifnot(e0$sigma == 0, abs(e0$q_c - (P2$a - P2$cost) / 4) < 1e-12, abs(e0$trigger - (P2$a - 2 * e0$q_c)) < 1e-12,
          !any(G2$war[G2$level == 1]), all(G2$price[G2$level == 1] == e0$trigger))
eqv <- lapply(tr$G2$eq, gp_values)
stopifnot(all(vapply(tr$G2$eq, function(e) e$delta == P2$delta && e$n_pun == P2$n_pun, logical(1))),
          eqv[[1]]$beta == 0, all(diff(vapply(eqv, `[[`, numeric(1), "beta")) > 0),
          all(diff(vapply(tr$G2$eq, `[[`, numeric(1), "q_c")) > 0))
# the price during a war is the Cournot price plus the shock; in collusion the collusive price plus the shock
qN <- (P2$a - P2$cost) / 3
for (l in 2:length(tr$G2$sigma)) {
  e <- tr$G2$eq[[l]]; s <- G2[G2$level == l, ]
  set.seed(P$seed$G2); z <- matrix(rnorm(P2$K * P2$T), P2$K, P2$T)
  zz <- z[cbind(s$market, s$t)]
  stopifnot(max(abs(s$price - ifelse(s$war, P2$a - 2 * qN, P2$a - 2 * e$q_c) - e$sigma * zz)) < 1e-12,
            identical(s$trigger_hit, !s$war & s$price < e$trigger))
}

## ---- G3: learning histories ------------------------------------------------------------------------------
stopifnot(P$G3$T >= 50000, length(attr(G3$rm_compat, "a1")) == P$G3$T, length(attr(G3$rm_chicken, "a2")) == P$G3$T,
          length(G3$fp_mp$a1) == P$G3$T, nrow(G3$reps_compat) == P$G3$R, nrow(G3$reps_chicken) == P$G3$R,
          all(abs(rowSums(G3$reps_compat[, c("mu11", "mu12", "mu21", "mu22")]) - 1) < 1e-12))

## ---- G4: global-game populations -------------------------------------------------------------------------
P4 <- P$G4
stopifnot(identical(tr$G4$sigma, c(0.1, 0.05, 0.01)),
          max(abs(tr$G4$x_star - gg_threshold(tr$G4$sigma, prior_mean = P4$prior_mean, prior_sd = P4$prior_sd))) < 1e-12,
          nrow(G4$pop) == length(tr$G4$sigma) * length(P4$theta_grid),
          all(G4$pop$n_adopt >= 0 & G4$pop$n_adopt <= P4$n_pop))
