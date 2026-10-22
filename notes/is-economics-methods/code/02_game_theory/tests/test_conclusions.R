# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 2.
# The first block is the chapter card, verbatim (comments translated, Ruling P6). Under Ruling P23 the
# intervals g1$<block>$ci are simulation ranges (the 2.5% to 97.5% range of the estimates over 1,000 Monte
# Carlo replications), not standard-error intervals. The blocks after the card lock the conclusions the
# text will draw from the same results: the by-market inversion in the high-precision block (Task 02.1
# brief: g1$by_market must cover the truth), the stationary price-war share, what patience does to the
# cartel (fix round 1), the learning runs and the global-game threshold.
# Run from the chapter code directory after 01 to 05 (tests/run_tests.R does this).
# ------------------------------------------------------------------

## ---- chapter card (verbatim) --------------------------------------------------------------
source("../../tools/R/chapter_kit.R"); tr <- readRDS("truth02.rds")
g1 <- readRDS("res_qre.rds"); g2 <- readRDS("res_gp.rds"); g3 <- readRDS("res_learn.rds"); g4 <- readRDS("res_gg.rds")
stopifnot(covers(g1$unique$ci, tr$G1$payoff), !covers(g1$pooled$ci, tr$G1$payoff))  # a unique equilibrium inverts; mixing in multiple equilibria does not
stopifnot(g2$war_freq[1] == 0, all(diff(g2$war_freq) > 0), g2$osd_ok, !g2$osd_low_delta)
stopifnot(g3$ce_violation_T < 0.01)
stopifnot(all(diff(abs(g4$theta_star - 0.5)) < 0), abs(tail(g4$theta_star, 1) - 0.5) < 0.02)

## ---- G1: simulation statements (Ruling P23). g1$<block>$ci is the 2.5% to 97.5% range of the estimates over the
## Monte Carlo replications and coef their mean. Multiplicity alone does not break the inversion (grouping by
## equilibrium covers the truth); pooling across equilibria does -----------------------------------------------
for (k in c("unique", "pooled", "by_market")) stopifnot(g1[[k]]$n_rep >= 1000, length(g1[[k]]$ci) == 2L)
stopifnot(covers(g1$by_market$ci, tr$G1$payoff),                                   # selection rule known: consistent
          g1$pooled$ci[["hi"]] < tr$G1$payoff,                                       # pooled: below the truth (this design)
          covers(g1$unique$ci_u22, tr$G1$payoff_u22), covers(g1$by_market$ci_u22, tr$G1$payoff_u22),
          !covers(g1$pooled$ci_u22, tr$G1$payoff_u22),
          abs(g1$unique$coef[["u11"]] - tr$G1$payoff) < 0.01, abs(g1$by_market$coef[["u11"]] - tr$G1$payoff) < 0.01,
          g1$pooled$coef[["u11"]] < 0.5 * tr$G1$payoff)
# the displayed draw's plain OLS estimates lie inside their block's simulation range
for (k in c("unique", "pooled", "by_market")) stopifnot(covers(g1[[k]]$ci, g1[[k]]$draw[["u11"]]))

## ---- G2: price wars are on the equilibrium path, no firm ever deviates ---------------------------------
# simulated share of market-months at war matches the stationary share beta N / (1 + beta N) of the
# equilibrium's Markov chain; collusion still beats Cournot at every noise level
stopifnot(max(abs(g2$war_freq - g2$war_share_theory)) < 0.02,
          all(abs(g2$war_freq[-1] / g2$war_share_theory[-1] - 1) < 0.1),
          all(g2$value_per_period > g2$cournot_profit), all(diff(g2$value_per_period) < 0),
          all(diff(g2$trigger) < 0), all(g2$osd_each))

## ---- G2: the comparative statics of section 4.2 (the "这说明" paragraph of Theorem 18). With more noise the loss from a price
## war Vc - Vp falls, the marginal effect of one's own quantity on the war probability (dzeta, defined for sigma > 0) rises, and the
## collusive quantity rises: deterrence weakens, so colluders produce more ---------------------------------------------
stopifnot(all(diff(g2$gap_Vc_Vw) < 0), all(diff(g2$dzeta[-1]) > 0), all(diff(g2$q_c) > 0))

## ---- G2: patience. The card's g2$osd_low_delta concerns the delta = 0.95 strategies held fixed: at delta_low = 0.5
## they fail the one-shot deviation check at every noise level. That does not mean impatient firms cannot collude:
## with the trigger re-optimised at delta_low, collusion still pays at sigma = 0, 0.1, 0.2 (less than at
## delta = 0.95, more than Cournot) and collapses to Cournot at sigma = 0.4 ------------------------------------------
stopifnot(!any(g2$fixed_strategy_pass_low), identical(g2$osd_low_delta, any(g2$fixed_strategy_pass_low)),
          all(g2$low_osd_each), all(g2$value_low_per_period < g2$value_per_period),
          identical(g2$collusive_low, c(TRUE, TRUE, TRUE, FALSE)),
          all(g2$value_low_per_period[g2$collusive_low] > g2$cournot_profit),
          all(abs(g2$value_low_per_period[!g2$collusive_low] - g2$cournot_profit) < 1e-9))

## ---- G3: learning ----------------------------------------------------------------------------------------
stopifnot(length(g3$ce_violation_T) == 2L, max(g3$reps_max_violation) < 0.01,
          max(abs(g3$fp_freq - 0.5)) < 0.01)

## ---- G4: the threshold falls to 1/2 and the simulated populations tip there -------------------------------
stopifnot(identical(g4$sigma, c(0.1, 0.05, 0.01)), all(g4$theta_star > 0.5),           # pessimistic public signal: above 1/2
          max(abs(g4$theta_star_sim - g4$theta_star)) < 0.01, max(abs(g4$mc_indiff)) < 0.005,
          all(g4$min_slope > 0), all(g4$sigma < g4$sigma_max))
