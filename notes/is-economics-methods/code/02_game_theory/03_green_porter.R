# ------------------------------------------------------------------
# 03_green_porter.R -- block G2: price wars on the equilibrium path under imperfect public monitoring.
#
# For each demand-noise level sigma (0, 0.1, 0.2, 0.4 USD) the equilibrium of 01_dgp.R (trigger price and
# collusive quantity that maximise the collusive value subject to the incentive constraint) is checked
# with the one-shot deviation principle at delta = 0.95. The simulated panels (200 markets x 300 months, no
# firm ever deviates) give the share of market-months at war, compared with the stationary share
# beta N_pun / (1 + beta N_pun) of the phase chain.
#
# Patience, two different questions at delta_low = 0.5:
#   fixed_strategy_pass_low  do the delta = 0.95 strategies, held fixed, still pass the one-shot deviation
#                            check? (No at every sigma: their collusive quantity solves the delta = 0.95
#                            first-order condition.) The card's field osd_low_delta = any() of these, so the
#                            card's !osd_low_delta reads "none of the delta = 0.95 strategies survives".
#   value_low_per_period     the best collusion impatient firms can sustain: the trigger re-optimised at
#                            delta_low (gp_optimal_trigger), its value as a per-month annuity and whether it
#                            colludes at all (collusive_low); Cournot forever when no trigger works.
# Output: res_gp.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("games_fns.R")
tr <- readRDS("truth02.rds"); P2 <- tr$params$G2
G2 <- readRDS("osprey_G2.rds")
eqs <- tr$G2$eq; L <- length(eqs)
vals <- lapply(eqs, gp_values)
num <- function(f) vapply(seq_len(L), f, numeric(1))

osd_main <- lapply(eqs, osd_check)
osd_fixed_low <- lapply(eqs, function(e) osd_check(modifyList(e, list(delta = P2$delta_low))))
eq_low <- lapply(P2$sigma, function(s) gp_optimal_trigger(list(a = P2$a, cost = P2$cost, delta = P2$delta_low,
                                                                n_pun = P2$n_pun, sigma = s)))
vals_low <- lapply(eq_low, gp_values)
lvl <- function(l) G2[G2$level == l, ]

res <- list(
  sigma = P2$sigma,
  trigger = num(function(l) eqs[[l]]$trigger),
  q_c = num(function(l) eqs[[l]]$q_c),
  price_collusive = num(function(l) P2$a - 2 * eqs[[l]]$q_c),
  beta = num(function(l) vals[[l]]$beta),
  # the L345 comparative statics (section 4.2): the loss from a price war Vc - Vp, and the marginal effect of
  # one's own quantity on the war probability dzeta = phi(z) / sigma (undefined at sigma = 0, where the war probability is a step)
  gap_Vc_Vw = num(function(l) vals[[l]]$Vc - vals[[l]]$Vp),
  dzeta = num(function(l) if (P2$sigma[l] > 0) stats::dnorm((eqs[[l]]$trigger - (eqs[[l]]$a - 2 * eqs[[l]]$q_c)) / P2$sigma[l]) / P2$sigma[l] else NA_real_),
  value_per_period = num(function(l) (1 - P2$delta) * vals[[l]]$Vc),
  value_war_start_per_period = num(function(l) (1 - P2$delta) * vals[[l]]$Vp),
  cournot_profit = vals[[1]]$pi_N, cournot_price = P2$a - 2 * vals[[1]]$q_N, cournot_q = vals[[1]]$q_N,
  monopoly_profit = tr$G2$monopoly$profit, monopoly_price = tr$G2$monopoly$price,
  war_share_theory = num(function(l) vals[[l]]$beta * P2$n_pun / (1 + vals[[l]]$beta * P2$n_pun)),
  war_freq = num(function(l) mean(lvl(l)$war)),
  wars_per_market_year = num(function(l) 12 * sum(lvl(l)$trigger_hit) / (P2$K * P2$T)),
  mean_price_collusion = num(function(l) mean(lvl(l)$price[!lvl(l)$war])),
  mean_price_war = num(function(l) if (any(lvl(l)$war)) mean(lvl(l)$price[lvl(l)$war]) else NA_real_),
  share_markets_with_war = num(function(l) mean(tapply(lvl(l)$war, lvl(l)$market, any))),
  osd_each = vapply(osd_main, function(o) isTRUE(as.vector(o)), logical(1)),
  # (a) the delta = 0.95 strategies held fixed at delta_low
  fixed_strategy_pass_low = vapply(osd_fixed_low, function(o) isTRUE(as.vector(o)), logical(1)),
  fixed_strategy_gain_low = num(function(l) attr(osd_fixed_low[[l]], "gain_collusive")),
  fixed_strategy_q_dev_low = num(function(l) attr(osd_fixed_low[[l]], "q_dev")),
  # (b) the best collusion at delta_low, trigger re-optimised
  collusive_low = vapply(eq_low, function(e) e$q_c < vals[[1]]$q_N - 1e-9, logical(1)),
  value_low_per_period = num(function(l) (1 - P2$delta_low) * vals_low[[l]]$Vc),
  q_c_low = num(function(l) eq_low[[l]]$q_c),
  trigger_low = num(function(l) eq_low[[l]]$trigger),
  beta_low = num(function(l) vals_low[[l]]$beta),
  war_share_theory_low = num(function(l) vals_low[[l]]$beta * P2$n_pun / (1 + vals_low[[l]]$beta * P2$n_pun)),
  low_osd_each = vapply(eq_low, function(e) isTRUE(as.vector(osd_check(e))), logical(1)),
  eq_low = eq_low,
  delta = P2$delta, delta_low = P2$delta_low, n_pun = P2$n_pun, K = P2$K, T = P2$T)
res$osd_ok <- all(res$osd_each)
res$osd_low_delta <- any(res$fixed_strategy_pass_low)     # card field: see the header, (a)
# standard error of the simulated war share, from the market-level shares (markets are independent)
res$war_freq_se <- num(function(l) { s <- tapply(lvl(l)$war, lvl(l)$market, mean); stats::sd(s) / sqrt(length(s)) })
# one market's path at sigma = 0.2 for the text (first 120 months of market 1)
res$path_example <- as.data.frame(lvl(3)[lvl(3)$market == 1L & lvl(3)$t <= 120L, c("t", "war", "price", "trigger_hit")])
saveRDS(res, "res_gp.rds")

for (l in seq_len(L))
  cat(sprintf("sigma %.2f  trigger %.4f  price %.4f  beta %.4f  war share: theory %.4f simulated %.4f (se %.4f)  value %.4f  OSD %s\n",
              res$sigma[l], res$trigger[l], res$price_collusive[l], res$beta[l], res$war_share_theory[l], res$war_freq[l],
              res$war_freq_se[l], res$value_per_period[l], res$osd_each[l]))
for (l in seq_len(L))
  cat(sprintf("delta %.2f, sigma %.2f: fixed delta-%.2f strategies pass %s (gain %.4f); re-optimised: collusive %s, q_c %.6f, trigger %.4f, beta %.4f, value %.6f (Cournot %.6f)\n",
              P2$delta_low, res$sigma[l], P2$delta, res$fixed_strategy_pass_low[l], res$fixed_strategy_gain_low[l],
              res$collusive_low[l], res$q_c_low[l], res$trigger_low[l], res$beta_low[l], res$value_low_per_period[l],
              res$cournot_profit))
cat("Saved res_gp.rds\n")
