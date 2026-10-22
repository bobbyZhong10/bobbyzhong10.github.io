# ------------------------------------------------------------------
# tests/test_inherit.R -- block G0 reproduces every number of the legacy chapter-2 script
# (snapshot: dev markdown/upgrade_2026-10-19/sdd_partI/legacy/code/02_game_theory/01_games.R, sections
# 6.1 to 6.8 and 4.1) and of the legacy text (introduction, 3.2, 3.3, 3.7, 3.8, 3.9, 4, 6), to the precision
# printed there; the card's G0 list: pricing prisoner's dilemma (10, 14, 6, 2); compatibility game 3/4, 1/4,
# coordination probability 3/8; BNE beta = 2/5, 7.2 against 6.6; grim trigger delta >= 1/2; Rubinstein
# 1/(1 + delta); Nash solution 50/50 and 60/40; correlated equilibrium 2 against 3/4; global-game threshold 1/2.
# (The card's last G0 item, the cheap-talk interval count, has no legacy number in chapter 2: cheap talk is
# chapter 3's section 3.8. See the Task 02.1 report.)
# Aster is renamed Osprey (Ruling P12).
# Run from the chapter code directory after 01_dgp.R (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("games_fns.R")
G0 <- readRDS("truth02.rds")$G0
at <- function(x, v, d) length(x) == length(v) && all(abs(round(x, d) - v) < 1e-9)   # equal at d printed decimals

## ---- legacy 6.1 and introduction: pricing game, dominance, Nash equilibrium -----------------
A <- matrix(c(10, 2, 14, 6), nrow = 2, byrow = TRUE)
stopifnot(identical(unname(G0$pd_A), A), identical(unname(G0$pd_B), t(A)),
          isTRUE(G0$pd_high_dominated), identical(G0$pd_ne, list(c(2L, 2L))),
          at(G0$pd_ne_payoff, 6, 9), at(G0$pd_coop_payoff, 10, 9), at(G0$pd_temptation, 14, 9), at(G0$pd_sucker, 2, 9),
          at(G0$pd_industry_coop, 20, 9), at(G0$pd_industry_ne, 12, 9), at(G0$pd_forecast_excess, 2 / 3, 12))

## ---- legacy 6.2 and 3.2: compatibility game -----------------------------------------------------
stopifnot(identical(unname(G0$compat_A), matrix(c(3, 0, 0, 1), 2, byrow = TRUE)),
          identical(unname(G0$compat_B), matrix(c(1, 0, 0, 3), 2, byrow = TRUE)),
          identical(G0$compat_pure, list(c(1L, 1L), c(2L, 2L))),
          at(G0$compat_p, 0.75, 3), at(G0$compat_q, 0.25, 3), at(G0$compat_coord_prob, 0.375, 3),
          at(G0$compat_eu_osprey, 0.75, 3), at(G0$compat_eu_boreal, 0.75, 3), at(G0$compat_miscoord_prob, 0.625, 3),
          abs(G0$compat_coord_prob - 3 / 8) < 1e-12)

## ---- legacy 6.3: sequential compatibility game, backward induction ----------------------------
stopifnot(identical(G0$seq_boreal_follow, c(1L, 2L)), identical(G0$seq_osprey_choice, 1L),
          at(G0$seq_payoff_osprey, 3, 9), at(G0$seq_payoff_boreal, 1, 9))

## ---- legacy 6.4: Bayesian Nash equilibrium ----------------------------------------------------------
stopifnot(at(G0$bne_beta, 0.4, 9), at(G0$bne_eu_high, 7.2, 1), at(G0$bne_eu_low, 6.6, 1),
          at(G0$bne_threshold, 0.5, 9), isTRUE(G0$bne_osprey_high),
          abs(G0$bne_eu_high - (10 - 7 * 2 / 5)) < 1e-12, abs(G0$bne_eu_low - (7 - 2 / 5)) < 1e-12)

## ---- legacy 6.5 and 4.1: grim trigger -----------------------------------------------------------------
stopifnot(at(c(G0$R, G0$T, G0$P), c(10, 14, 6), 9), at(G0$grim_delta_star, 0.5, 9),
          at(G0$grim_coop_at_star, 20, 1), at(G0$grim_cheat_at_star, 20, 1),
          abs(G0$grim_coop_at_star - G0$grim_cheat_at_star) < 1e-9)
stopifnot(at(G0$grim_deltas, c(0.4, 0.5, 0.7, 0.9), 9), at(G0$grim_slack, c(-1.33, 0, 5.33, 32), 2),
          all(diff(G0$grim_slack) > 0), abs(G0$grim_slack[G0$grim_deltas == G0$grim_delta_star]) < 1e-9)

## ---- legacy 6.6 and 3.8: Nash bargaining and Rubinstein ----------------------------------------------
stopifnot(at(G0$V, 100, 9), at(G0$nb_sym, 50, 9), at(G0$nb_asym, 60, 9), at(G0$nb_asym_boreal, 40, 9),
          at(G0$nb_weighted, 68, 9), at(G0$nb_alpha, 0.6, 9), at(G0$nb_d1, 20, 9))
stopifnot(at(G0$rub_delta, c(0.8, 0.99), 9),
          at(100 * G0$rub_mover, c(55.5556, 50.2513), 4), at(100 * G0$rub_responder, c(44.4444, 49.7487), 4),
          at(100 * G0$rub_mover, c(55.56, 50.25), 2), at(100 * G0$rub_responder, c(44.44, 49.75), 2),
          max(abs(G0$rub_mover - 1 / (1 + G0$rub_delta))) < 1e-14,
          max(abs(G0$rub_mover + G0$rub_responder - 1)) < 1e-14)

## ---- legacy 6.7 and 3.7: correlated equilibrium with a public coin ------------------------------------
stopifnot(at(G0$ce_osprey, 2, 2), at(G0$ce_boreal, 2, 2), at(G0$ce_mixed_ne, 0.75, 2),
          G0$ce_violation_coin <= 1e-12, G0$ce_osprey > G0$ce_mixed_ne)

## ---- legacy 6.8 and 3.9: global-game threshold ------------------------------------------------------------
stopifnot(at(G0$gg_theta_star, 0.5, 9), at(G0$gg_E_l, 0.5, 9))
