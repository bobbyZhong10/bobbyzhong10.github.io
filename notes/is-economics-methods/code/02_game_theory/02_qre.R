# ------------------------------------------------------------------
# 02_qre.R -- block G1: inverting the logit QRE by OLS on group log-odds.
#
# The analyst sees, in each local market, the shifter x_m, how many of the 50 integration decisions each
# platform made on S1 and (in the high-precision block) which standard came first. She knows the precision
# chi of each block. Under the QRE, Osprey's log-odds of S1 in a group of markets that play the same QRE is
# chi times its expected payoff difference at Boreal's frequency q:
#     log(p / (1 - p)) / chi - x = u11 q - u22 (1 - q)
# and qre_logodds_ols fits u11 (truth 3) and u22 (truth 1) by OLS. Three groupings:
#   unique      unique-QRE block, groups by x (25 groups of 80 markets)
#   pooled      high-precision block, groups by x only (5 groups of 400 markets): each group mixes markets
#               playing the S1-favouring and the S2-favouring QRE
#   by_market   high-precision block, groups by x and history (10 groups): markets in a group play the
#               same QRE, i.e. the selection rule is known; the market-by-market inversion with each market's
#               own observed opponent frequency, pooled only over markets that play the same equilibrium
# Ruling P23: Part I uses no standard-error formulas. The displayed draw (osprey_G1.rds) gives plain OLS point
# estimates (`draw`); the statements about the estimator come from simulation: 1,000 fresh draws of block G1
# (seed G1_mc = 2102: same markets and equilibria, new histories and decisions), the three inversions on each,
# and for each grouping the mean estimate (`coef`) and the 2.5% to 97.5% range of the estimates (`ci`,
# `ci_u22`), with their standard deviation (`sd`).
# Output: res_qre.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("games_fns.R")
tr <- readRDS("truth02.rds"); P1 <- tr$params$G1
G1 <- readRDS("osprey_G1.rds")

designs <- list(unique = list(block = "unique", by = "x"), pooled = list(block = "high", by = "x"),
                by_market = list(block = "high", by = c("x", "history")))
groups_of <- function(d, k) {
  g <- g1_group_freq(d[d$block == designs[[k]]$block, ], designs[[k]]$by)
  g$chi <- P1$lambda[[designs[[k]]$block]]
  g
}

## ---- the displayed draw: plain OLS point estimates --------------------------------------------------------
shown <- lapply(setNames(names(designs), names(designs)), function(k) groups_of(G1, k))

## ---- Monte Carlo: 1,000 fresh draws, the same three inversions on each ----------------------------------------
R <- 1000L
sims <- setNames(vector("list", length(designs)), names(designs))
for (k in names(designs)) sims[[k]] <- vector("list", R)
set.seed(tr$params$seed$G1_mc)
for (r in seq_len(R)) {
  d <- gen_G1(P1, tr$G1$eq)
  for (k in names(designs)) { g <- groups_of(d, k); g$rep <- r; sims[[k]][[r]] <- g }
}
res <- lapply(setNames(names(designs), names(designs)), function(k) {
  mc <- qre_logodds_ols(do.call(rbind, sims[[k]]))
  list(draw = qre_logodds_ols(shown[[k]])$coef, coef = mc$coef, ci = mc$ci, ci_u22 = mc$ci_u22, sd = mc$sd,
       est = mc$est, n_rep = mc$n_rep, n_groups = nrow(shown[[k]]), groups = shown[[k]])
})
res$R <- R

## ---- what the pooled groups look like: the mixture of the two QRE in each x-group ------------------------
eqh <- tr$G1$eq[tr$G1$eq$block == "high", ]
mix <- res$pooled$groups
mix$share_S1 <- vapply(mix$x, function(x) mean(G1$history[G1$block == "high" & abs(G1$x - x) < 1e-12] == "S1"), numeric(1))
mix$p_S1 <- vapply(mix$x, function(x) eqh$p[abs(eqh$x - x) < 1e-12 & eqh$role == "S1"], numeric(1))
mix$p_S2 <- vapply(mix$x, function(x) eqh$p[abs(eqh$x - x) < 1e-12 & eqh$role == "S2"], numeric(1))
mix$q_S1 <- vapply(mix$x, function(x) eqh$q[abs(eqh$x - x) < 1e-12 & eqh$role == "S1"], numeric(1))
mix$q_S2 <- vapply(mix$x, function(x) eqh$q[abs(eqh$x - x) < 1e-12 & eqh$role == "S2"], numeric(1))
res$pooled_mixture <- mix
res$payoff <- tr$G1$payoff; res$payoff_u22 <- tr$G1$payoff_u22
saveRDS(res, "res_qre.rds")

for (k in names(designs)) {
  z <- res[[k]]
  cat(sprintf("%-9s draw: u11 %.4f u22 %.4f (%d groups) | %d replications: u11 mean %.4f, 2.5%%-97.5%% range [%.4f, %.4f], sd %.4f; u22 mean %.4f range [%.4f, %.4f]\n",
              k, z$draw[["u11"]], z$draw[["u22"]], z$n_groups, z$n_rep, z$coef[["u11"]], z$ci[["lo"]], z$ci[["hi"]],
              z$sd[["u11"]], z$coef[["u22"]], z$ci_u22[["lo"]], z$ci_u22[["hi"]]))
}
cat("Saved res_qre.rds\n")
