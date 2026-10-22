# ------------------------------------------------------------------
# 01_dgp.R -- the Osprey-Boreal world for chapter 2 (game theory)
#
# Osprey and Boreal are two (fictional) digital platforms (Osprey replaces the legacy Aster, Ruling P12).
# Money is in US dollars. Five blocks:
#   G0  the legacy hand calculations (no random numbers): pricing prisoner's dilemma (10, 14, 6, 2), its
#       minmax and its feasible and individually rational payoff sets (new), the compatibility game (pure and mixed equilibria, 3/4, 1/4, 3/8), the sequential version,
#       the BNE with beta = 2/5 (7.2 against 6.6), grim trigger (delta >= 1/2), Nash bargaining (50/50, 60/40,
#       68), Rubinstein 1/(1 + delta), Binmore-Rubinstein-Wolinsky shares as the time between offers shrinks
#       (new), the public-coin correlated equilibrium (2 against 3/4), the global-game threshold 1/2.
#   G1  QRE markets (new): 4,000 local markets in which the two platforms play the compatibility game with a
#       market-level payoff shifter x_m (thousand USD a month) on S1. Unique block: 2,000 markets at
#       precision chi = 0.8, one QRE at every x_m. High-precision block: 2,000 markets at chi = 2, three QRE at
#       every x_m; a market plays the QRE that favours the standard that came first in it (history), and S1
#       came first with probability plogis(8 x_m). Each market makes 50 integration decisions per platform.
#   G2  imperfect public monitoring (new; Green and Porter 1984): Cournot duopoly in K = 200 local markets
#       for T = 300 months, demand shocks the firms never see, trigger price and a 6-month punishment; for
#       each shock sd sigma in {0, 0.1, 0.2, 0.4} USD the trigger price that maximises the collusive value and
#       the collusive quantity are solved numerically from the incentive constraint.
#   G3  learning (new): regret matching (Hart and Mas-Colell 2000) for 50,000 periods on the compatibility
#       game and on a Chicken game, plus 100 independent runs of each; fictitious play for 50,000 periods
#       on matching pennies.
#   G4  global game (new): adoption game with payoff theta + l - 1, a public prior theta ~ N(0.25, 0.5^2)
#       and private noise sigma in {0.1, 0.05, 0.01}; thresholds solved numerically; simulated populations
#       of 2,000 players on a grid of theta, and a Monte Carlo check of the marginal player's indifference.
# Every number in chapter 2 comes from these data, truth02.rds and the res_*.rds files.
#
# Outputs (in this directory):
#   osprey_G1.rds   data.frame: block market x history n_dec k_osprey k_boreal p_true q_true
#   osprey_G2.rds   data.table: level sigma market t war price trigger_hit
#   osprey_G3.rds   list: rm_compat, rm_chicken (regret_matching output), fp_mp (fictitious_play output),
#                   reps_compat, reps_chicken (regret_matching_reps output)
#   osprey_G4.rds   list: pop (sigma theta n_pop n_adopt), mc (sigma x_star m_post El_sim n_mc)
#   truth02.rds     list(G0, G1, G2, G3, G4, params)
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 5.3): QRE precision chi (code: lambda), Green-Porter trigger price underline-p
# (code: trigger), punishment length N_pun (code: n_pun), number of markets N_mkt (code: K), periods
# script-T (code: T), private noise sigma of the global game, threshold theta*.
#
# Seeds (every simulated block its own seed, set at its start): G1 2101, G1 Monte Carlo 2102 (used in
# 02_qre.R), G2 2103 (the same seed at every sigma: common random numbers, so the four panels share their
# demand shocks), regret matching 2104 (compatibility game) and 2105 (Chicken), fictitious play 2106, the
# 100 regret-matching runs 2107, G4 2108. The block parameters were calibrated on exploratory seeds outside
# this directory before these seeds were fixed; the seeds were not tuned afterwards.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("games_fns.R")

seed <- list(G1 = 2101, G1_mc = 2102, G2 = 2103, G3_rm_compat = 2104, G3_rm_chicken = 2105, G3_fp = 2106,
             G3_rm_reps = 2107, G4 = 2108)

## =================================================================
## G0: the legacy hand calculations (no random numbers)
## =================================================================
## Payoffs in million USD a month (legacy section 1). Matrices are [Osprey's action, Boreal's action].
lab <- function(r, c) list(Osprey = r, Boreal = c)
A_pd <- matrix(c(10, 2, 14, 6), 2, byrow = TRUE, dimnames = lab(c("High", "Low"), c("High", "Low")))
B_pd <- matrix(t(A_pd), 2, dimnames = dimnames(A_pd))                       # symmetric game: B = t(A)
A_cp <- matrix(c(3, 0, 0, 1), 2, byrow = TRUE, dimnames = lab(c("S1", "S2"), c("S1", "S2")))
B_cp <- matrix(c(1, 0, 0, 3), 2, byrow = TRUE, dimnames = lab(c("S1", "S2"), c("S1", "S2")))
G0 <- local({
  out <- list()
  ## 6.1 pricing game: Low strictly dominates High; (Low, Low) is the unique equilibrium
  ne <- nash_2x2(A_pd, B_pd)
  out$pd_A <- A_pd; out$pd_B <- B_pd
  out$pd_high_dominated <- all(A_pd["Low", ] > A_pd["High", ])
  out$pd_ne <- ne$pure
  out$pd_ne_payoff <- A_pd["Low", "Low"]; out$pd_coop_payoff <- A_pd["High", "High"]
  out$pd_temptation <- A_pd["Low", "High"]; out$pd_sucker <- A_pd["High", "Low"]
  out$pd_industry_coop <- A_pd["High", "High"] + B_pd["High", "High"]
  out$pd_industry_ne <- A_pd["Low", "Low"] + B_pd["Low", "Low"]
  out$pd_forecast_excess <- out$pd_industry_coop / out$pd_industry_ne - 1      # the analyst's 20 is 2/3 above 12
  out$pd_minmax_osprey <- minmax_2x2(A_pd); out$pd_minmax_boreal <- minmax_2x2(t(B_pd))   # folk theorem
  fs <- folk_sets(unname(A_pd), unname(B_pd))                                  # fig_02_folk.svg
  out$folk_hull <- fs$hull; out$folk_ir <- fs$ir; out$folk_minmax <- fs$minmax
  out$folk_area_hull <- fs$area_hull; out$folk_area_ir <- fs$area_ir
  ## 6.2 compatibility game: two pure equilibria and the mixed one
  nc <- nash_2x2(A_cp, B_cp)
  p <- nc$mixed[1]; q <- nc$mixed[2]
  joint <- outer(c(p, 1 - p), c(q, 1 - q))
  out$compat_A <- A_cp; out$compat_B <- B_cp; out$compat_pure <- nc$pure
  out$compat_p <- p; out$compat_q <- q
  out$compat_coord_prob <- joint[1, 1] + joint[2, 2]; out$compat_miscoord_prob <- 1 - out$compat_coord_prob
  out$compat_eu_osprey <- sum(joint * A_cp); out$compat_eu_boreal <- sum(joint * B_cp)
  ## 6.3 Osprey moves first: Boreal matches, Osprey picks S1 (backward induction)
  follow <- c(unname(which.max(B_cp[1, ])), unname(which.max(B_cp[2, ])))   # Boreal's reply to S1, to S2
  pay_if <- c(A_cp[1, follow[1]], A_cp[2, follow[2]])
  out$seq_boreal_follow <- follow; out$seq_osprey_choice <- unname(which.max(pay_if))
  out$seq_payoff_osprey <- A_cp[out$seq_osprey_choice, follow[out$seq_osprey_choice]]
  out$seq_payoff_boreal <- B_cp[out$seq_osprey_choice, follow[out$seq_osprey_choice]]
  ## 6.4 BNE: Boreal is Normal (plays Low) with probability beta, Committed (plays High) otherwise
  A_bne <- matrix(c(10, 3, 7, 6), 2, byrow = TRUE, dimnames = lab(c("High", "Low"), c("High", "Low")))
  beta <- 2 / 5
  eu <- drop(A_bne %*% c(1 - beta, beta))
  out$bne_A <- A_bne; out$bne_beta <- beta
  out$bne_eu_high <- eu[["High"]]; out$bne_eu_low <- eu[["Low"]]
  out$bne_threshold <- (A_bne[1, 1] - A_bne[2, 1]) / ((A_bne[1, 1] - A_bne[2, 1]) + (A_bne[2, 2] - A_bne[1, 2]))
  out$bne_osprey_high <- eu[["High"]] >= eu[["Low"]]
  ## 6.5 and 4.1 grim trigger: cooperate on High while the rival does, Low forever after a cut
  R <- A_pd["High", "High"]; Tm <- A_pd["Low", "High"]; P <- A_pd["Low", "Low"]
  coop <- function(d) R / (1 - d); cheat <- function(d) Tm + d * P / (1 - d)
  out$R <- R; out$T <- Tm; out$P <- P
  out$grim_delta_star <- (Tm - R) / (Tm - P)
  out$grim_coop_at_star <- coop(out$grim_delta_star); out$grim_cheat_at_star <- cheat(out$grim_delta_star)
  out$grim_deltas <- c(0.4, 0.5, 0.7, 0.9)
  out$grim_slack <- coop(out$grim_deltas) - cheat(out$grim_deltas)
  ## 6.6 bargaining over V = 100: Nash solution d1 + alpha (V - d1 - d2); Rubinstein 1 / (1 + delta)
  V <- 100
  nb <- function(d1, d2, alpha = 0.5) d1 + alpha * (V - d1 - d2)
  out$V <- V; out$nb_d1 <- 20; out$nb_alpha <- 0.6
  out$nb_sym <- nb(0, 0); out$nb_asym <- nb(20, 0); out$nb_asym_boreal <- V - nb(20, 0); out$nb_weighted <- nb(20, 0, 0.6)
  out$rub_delta <- c(0.8, 0.99)
  out$rub_mover <- rubinstein_share(out$rub_delta, out$rub_delta); out$rub_responder <- 1 - out$rub_mover
  ## (new) Binmore, Rubinstein and Wolinsky (1986): offers every dt years, delta_i = exp(-r_i dt). Osprey's
  ## yearly discount rate 0.10, Boreal's 0.05 (Boreal is more patient). As dt -> 0 Osprey's share tends to
  ## r_B / (r_O + r_B) = 1/3 whoever proposes first: generalized Nash bargaining with weight 1/3.
  out$brw_r <- c(osprey = 0.10, boreal = 0.05)
  out$brw_dt <- c(1, 1 / 12, 1 / 52, 1 / 365)
  dO <- exp(-out$brw_r[["osprey"]] * out$brw_dt); dB <- exp(-out$brw_r[["boreal"]] * out$brw_dt)
  out$brw_share_osprey_first <- rubinstein_share(dO, dB)                   # Osprey proposes first
  out$brw_share_osprey_second <- 1 - rubinstein_share(dB, dO)              # Boreal proposes first
  out$brw_limit <- out$brw_r[["boreal"]] / sum(out$brw_r)
  out$brw_nash_share <- nb(0, 0, out$brw_limit)                            # generalized Nash, weight 1/3, of V = 100
  ## 6.7 a public fair coin over (S1, S1) and (S2, S2)
  coin <- matrix(c(0.5, 0, 0, 0.5), 2)
  out$ce_osprey <- sum(coin * A_cp); out$ce_boreal <- sum(coin * B_cp)
  out$ce_mixed_ne <- out$compat_eu_osprey
  out$ce_violation_coin <- ce_violation(coin, A_cp, B_cp)
  ## 6.8 global game with a flat prior: E[l | x*] = 1/2, threshold 1/2 (for any private noise)
  out$gg_E_l <- stats::pnorm(0); out$gg_theta_star <- gg_threshold(0.1)
  out
})
stopifnot(G0$compat_coord_prob == 3 / 8, G0$grim_delta_star == 0.5, G0$nb_asym == 60, G0$ce_osprey == 2,
          abs(G0$brw_share_osprey_first[4] - G0$brw_limit) < 1e-3)

## =================================================================
## G1: QRE markets
## =================================================================
## Payoffs in thousand USD a month in each local market: the compatibility game of G0 (Osprey 3 on (S1, S1),
## 1 on (S2, S2); Boreal 1 and 3; 0 on a mismatch) plus x_m on Osprey's S1 row and on Boreal's S1 column
## (the local S1 installed base, worth the same to both). The analyst knows chi and x_m and inverts Osprey's
## log-odds for u11 = 3 (truth G1$payoff) and u22 = 1.
P1 <- list(A0 = unname(A_cp), B0 = unname(B_cp), lambda = c(unique = 0.8, high = 2),
           x_unique = round(seq(-1.2, 1.2, by = 0.1), 1), x_high = round(seq(-0.2, 0.2, by = 0.1), 1),
           n_markets_unique = 2000L, n_markets_high = 2000L, n_dec = 50L, kappa = 8)
eq1 <- g1_equilibria(P1)
stopifnot(all(eq1$n_qre[eq1$block == "unique"] == 1L), all(eq1$n_qre[eq1$block == "high"] == 3L))
set.seed(seed$G1)
G1 <- gen_G1(P1, eq1)
# the QRE of the G0 game (x = 0) at the two precisions, and where the principal branch splits
qre_x0 <- list(unique = qre_all(P1$A0, P1$B0, P1$lambda[["unique"]]), high = qre_all(P1$A0, P1$B0, P1$lambda[["high"]]))
bif <- qre_bifurcation(P1$A0, P1$B0, 1, 1.2)
truth_G1 <- list(payoff = P1$A0[1, 1] - P1$A0[2, 1], payoff_u22 = P1$A0[2, 2] - P1$A0[1, 2], lambda = P1$lambda,
                 eq = eq1, qre_x0 = qre_x0, bifurcation = bif,
                 share_S1_first = mean(G1$history[G1$block == "high"] == "S1"))

## =================================================================
## G2: imperfect public monitoring (Green and Porter 1984)
## =================================================================
## Each local market: inverse demand p = 10 - Q + sigma z (USD per unit, Q in thousand units), unit cost
## 2 USD, monthly discount factor 0.95. Monopoly split 2 + 2 (price 6, profit 8 thousand USD a month each),
## Cournot 8/3 each (price 14/3, profit 64/9). Wars last 6 months. For each sigma, gp_optimal_trigger
## picks the trigger price that maximises the collusive value given that the collusive quantity must satisfy
## the incentive constraint (the first-order condition plus a global one-shot deviation check).
## delta_low = 0.5 is the lower patience at which 03_green_porter.R (a) re-checks the delta = 0.95 strategies
## held fixed and (b) re-optimises the trigger to find the best collusion impatient firms can sustain.
P2 <- list(a = 10, cost = 2, delta = 0.95, n_pun = 6L, sigma = c(0, 0.1, 0.2, 0.4), K = 200L, T = 300L,
           delta_low = 0.5)
eq2 <- lapply(P2$sigma, function(s) gp_optimal_trigger(list(a = P2$a, cost = P2$cost, delta = P2$delta,
                                                             n_pun = P2$n_pun, sigma = s)))
G2 <- data.table::rbindlist(lapply(seq_along(eq2), function(l) {
  s <- gp_simulate(eq2[[l]], T = P2$T, K = P2$K, seed = seed$G2)
  cbind(data.table::data.table(level = l, sigma = P2$sigma[l]), s)
}))
truth_G2 <- list(sigma = P2$sigma, eq = eq2, values = lapply(eq2, gp_values),
                 monopoly = list(q = (P2$a - P2$cost) / 4, price = (P2$a + P2$cost) / 2, profit = (P2$a - P2$cost)^2 / 8),
                 cournot = list(q = (P2$a - P2$cost) / 3, price = P2$a - 2 * (P2$a - P2$cost) / 3, profit = ((P2$a - P2$cost) / 3)^2),
                 deviation = list(q = 3 * (P2$a - P2$cost) / 8, profit = (3 * (P2$a - P2$cost) / 8)^2))

## =================================================================
## G3: learning
## =================================================================
## Payoffs in million USD a month, as G0. Chicken: each platform Yields or Enters a new city; both Yield
## (5, 5), one Enters alone (6 to the entrant, 2 to the other), both Enter (0, 0). Matching pennies: Osprey
## gains 1 when it matches Boreal's feature choice, Boreal gains 1 when it does not.
## The 100 independent runs record their violation at about 200 log-spaced checkpoints (10 to 50,000
## periods, plus 100, 1,000 and 10,000).
P3 <- list(T = 50000L, R = 100L, compat_A = P1$A0, compat_B = P1$B0,
           chicken_A = matrix(c(5, 2, 6, 0), 2, byrow = TRUE), chicken_B = t(matrix(c(5, 2, 6, 0), 2, byrow = TRUE)),
           mp_A = matrix(c(1, -1, -1, 1), 2, byrow = TRUE), mp_B = -matrix(c(1, -1, -1, 1), 2, byrow = TRUE))
P3$checkpoints <- as.integer(sort(unique(c(round(10^seq(1, log10(P3$T), length.out = 200)), 100, 1000, 10000))))
rm_compat <- regret_matching(P3$compat_A, P3$compat_B, T = P3$T, seed = seed$G3_rm_compat)
rm_chicken <- regret_matching(P3$chicken_A, P3$chicken_B, T = P3$T, seed = seed$G3_rm_chicken)
fp_mp <- fictitious_play(P3$mp_A, P3$mp_B, T = P3$T, seed = seed$G3_fp)
set.seed(seed$G3_rm_reps)
reps_compat <- regret_matching_reps(P3$compat_A, P3$compat_B, T = P3$T, R = P3$R, checkpoints = P3$checkpoints)
reps_chicken <- regret_matching_reps(P3$chicken_A, P3$chicken_B, T = P3$T, R = P3$R, checkpoints = P3$checkpoints)
G3 <- list(rm_compat = rm_compat, rm_chicken = rm_chicken, fp_mp = fp_mp, reps_compat = reps_compat, reps_chicken = reps_chicken)
# Chicken's equilibria: the two pure ones, the mixed one, and the correlated equilibrium with 1/3 on each of
# (Y, Y), (Y, E), (E, Y), whose payoff lies outside the convex hull of the Nash payoffs
nch <- nash_2x2(P3$chicken_A, P3$chicken_B)
jm <- outer(c(nch$mixed[1], 1 - nch$mixed[1]), c(nch$mixed[2], 1 - nch$mixed[2]))
mu3 <- matrix(c(1, 1, 1, 0) / 3, 2, byrow = TRUE)
truth_G3 <- list(chicken_pure = nch$pure, chicken_mixed = nch$mixed,
                 chicken_mixed_payoff = sum(jm * P3$chicken_A), chicken_ce3_payoff = sum(mu3 * P3$chicken_A),
                 chicken_ce3_violation = ce_violation(mu3, P3$chicken_A, P3$chicken_B),
                 chicken_pure_payoffs = rbind(c(P3$chicken_A[1, 2], P3$chicken_B[1, 2]), c(P3$chicken_A[2, 1], P3$chicken_B[2, 1])),
                 mp_mixed = nash_2x2(P3$mp_A, P3$mp_B)$mixed)

## =================================================================
## G4: global game
## =================================================================
## Adoption of a new Osprey-Boreal interconnection standard by merchants: adopting pays theta + l - 1, not
## adopting 0. Public prior (an analysts' consensus) theta ~ N(0.25, 0.5^2), pessimistic relative to the
## complete-information midpoint 1/2. Private noise sigma = 0.1, 0.05, 0.01. Populations of 2,000 merchants
## at each theta from 0.25 to 0.85 in steps of 0.01; 200,000 Monte Carlo draws for the indifference check.
P4 <- list(sigma = c(0.1, 0.05, 0.01), prior_mean = 0.25, prior_sd = 0.5,
           theta_grid = round(seq(0.25, 0.85, by = 0.01), 2), n_pop = 2000L, n_mc = 200000L)
x_star <- gg_threshold(P4$sigma, prior_mean = P4$prior_mean, prior_sd = P4$prior_sd)
set.seed(seed$G4)
G4 <- gen_G4(P4, x_star)
truth_G4 <- list(theta_star_limit = 0.5, sigma = P4$sigma, x_star = x_star,
                 sigma_max = gg_sigma_max(P4$prior_sd), min_slope = gg_min_slope(P4$sigma, P4$prior_sd))

## =================================================================
## Save
## =================================================================
saveRDS(G1, "osprey_G1.rds")
saveRDS(G2, "osprey_G2.rds")
saveRDS(G3, "osprey_G3.rds")
saveRDS(G4, "osprey_G4.rds")
saveRDS(list(G0 = G0, G1 = truth_G1, G2 = truth_G2, G3 = truth_G3, G4 = truth_G4,
             params = list(seed = seed, G1 = P1, G2 = P2, G3 = P3, G4 = P4)), "truth02.rds")

cat(sprintf("G0: mixed (p, q) = (%.2f, %.2f), coordination %.3f; BNE %.1f vs %.1f; grim delta* %.2f; Rubinstein %.4f\n",
            G0$compat_p, G0$compat_q, G0$compat_coord_prob, G0$bne_eu_high, G0$bne_eu_low, G0$grim_delta_star, G0$rub_mover[1]))
cat(sprintf("G1: %d markets; principal QRE branch splits at chi = %.4f; S1 first in %.3f of the high block\n",
            nrow(G1), bif[["lambda"]], truth_G1$share_S1_first))
for (l in seq_along(eq2)) {
  v <- gp_values(eq2[[l]])
  cat(sprintf("G2: sigma %.2f  trigger %.4f  q_c %.4f  beta %.4f  value per month %.4f\n",
              eq2[[l]]$sigma, eq2[[l]]$trigger, eq2[[l]]$q_c, v$beta, v$Vc * (1 - P2$delta)))
}
cat(sprintf("G3: %d periods; regret matching joint (compat) %s; FP frequencies %.4f %.4f\n", P3$T,
            paste(round(as.vector(t(rm_compat)), 4), collapse = " "), fp_mp$freq[1], fp_mp$freq[2]))
cat(sprintf("G4: thresholds %s; uniqueness for sigma < %.4f\n", paste(sprintf("%.5f", x_star), collapse = " "), truth_G4$sigma_max))
cat("Saved osprey_G1.rds, osprey_G2.rds, osprey_G3.rds, osprey_G4.rds, truth02.rds\n")
