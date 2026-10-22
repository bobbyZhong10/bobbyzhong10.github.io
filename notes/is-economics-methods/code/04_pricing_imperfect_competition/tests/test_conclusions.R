# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 4.
# The first block is the chapter card, verbatim (comment translated, Ruling P6). The blocks after it lock the
# further conclusions the text states from C1 to C4.
# Run from the chapter code directory after 01 to 05 (tests/run_tests.R does this).
# ------------------------------------------------------------------

## ---- chapter card (verbatim) ---------------------------------------------------------------------------
source("../../tools/R/chapter_kit.R")
c1 <- readRDS("res_logit.rds"); c2 <- readRDS("res_pd.rds"); c3 <- readRDS("res_bundle.rds"); c4 <- readRDS("res_leader.rds")
stopifnot(c1$pt_own_bertrand > 0, c1$pt_own_bertrand < 1, c1$pt_cross_bertrand > 0, all(c1$p_joint > c1$p_bertrand))
stopifnot(c2$profit[["fine"]] >= c2$profit[["coarse"]] * 0.995, c2$profit[["coarse"]] >= c2$profit[["uniform"]] * 0.995)
stopifnot(abs(c2$linear$Q_pd - c2$linear$Q_uniform) < 1e-8, c2$linear$W_pd < c2$linear$W_uniform)  # Schmalensee
stopifnot(tail(c3$rev_per_good, 1) > 0.40, tail(c3$dwl_share, 1) < 0.05, c3$rev_per_good[1] < 0.26)
stopifnot(all(c4$p_seq > c4$p_sim))

source("pricing_fns.R")
tr <- readRDS("truth04.rds")

## ---- C1: pass-through by OLS under two conducts -------------------------------------------------------------
# OLS recovers the average pass-through: every coefficient within 4 conventional standard errors of the truth
o <- c1$ols
stopifnot(nrow(o) == 8L, all(o$se > 0),
          all(mapply(function(e, s, v) covers(c(e - 4 * s, e + 4 * s), v), o$estimate, o$se, o$truth_mean)))
# the conduct shows in the cross pass-through, not in the own one: joint maximisation makes the cross pass-through
# negative (PT_kj = 1{k = j} - s_j), Bertrand makes it positive; the own pass-through hardly differs
stopifnot(c1$pt_cross_joint < 0, c1$pt_cross_bertrand - c1$pt_cross_joint > 4 * (c1$se_cross_bertrand + c1$se_cross_joint),
          abs(tr$C1$pt_bertrand_mean[1, 1] - tr$C1$pt_joint_mean[1, 1]) < 0.05,
          abs(c1$pt_own_bertrand - c1$pt_own_joint) < 2 * (c1$se_own_bertrand + c1$se_own_joint))
# homogeneous products: Cournot passes a firm's cost through at 1/3 (n / (n + 1) = 2/3 for a common shock);
# Bertrand prices at the higher cost, so the two cost coefficients sum to about 1
h <- c1$hom
stopifnot(all(mapply(function(e, s) covers(c(e - 4 * s, e + 4 * s), 1/3), h$cournot_coef, h$cournot_se)),
          covers(c(0.9, 1.1), sum(h$bertrand_coef)),
          sum(h$bertrand_coef) > sum(h$cournot_coef), h$markup_bertrand_mean < h$markup_cournot_mean)
# the single-product logit monopoly passes through 1 - s
stopifnot(abs(c1$mono_example$pt - (1 - c1$mono_example$s)) < 1e-9)
# pass-through depends on curvature and on conduct (figure): Cournot with linear demand below 1 and rising in n,
# constant elasticity above 1 and falling, exponential exactly 1; logit with a common cost shock: Bertrand rises
# with n, joint ownership equals the outside share s0 and falls with n
pc <- c1$pt_curves
g <- function(m) pc$pt[pc$model == m][order(pc$n[pc$model == m])]
stopifnot(all(g("cournot_linear") < 1), all(diff(g("cournot_linear")) > 0), all(g("cournot_ce2") > 1),
          all(diff(g("cournot_ce2")) < 0), all(abs(g("cournot_exp") - 1) < 1e-12),
          all(diff(g("logit_bertrand")) > 0), all(g("logit_bertrand") < 1), all(diff(g("logit_joint")) < 0),
          max(abs(g("logit_joint") - c1$pt_curves_s0)) < 1e-9, abs(g("logit_bertrand")[1] - g("logit_joint")[1]) < 1e-9)

## ---- C2: personalised pricing ---------------------------------------------------------------------------------
# finer segments: more profit, less consumer surplus, more merchants served
stopifnot(c2$profit[["coarse"]] > c2$profit[["uniform"]], c2$profit[["fine"]] > c2$profit[["coarse"]],
          c2$cs[["fine"]] < c2$cs[["uniform"]], c2$served[["fine"]] > c2$served[["uniform"]])
# BBM triangle of the evaluation half: floor = the evaluation half's best single-price profit, CS >= 0, total <= efficient.
# The consumer-optimal segmentation reaches the floor with efficient total surplus. Coarse and fine segmentation lie
# inside the triangle: their profit is above the floor (the substantive claim). The one-price regime does NOT: its price
# (24 USD) was set on the training half, while the evaluation half's best single price is 23 USD, so its profit is
# 0.04 USD below the floor (pd_bbm_uniform_gap); the gap is small, under 1% of the floor.
b <- c2$bbm
stopifnot(abs(b$profit - b$uniform_profit) < 1e-9, abs(b$total - b$efficient) < 1e-9, b$cs_max > max(c2$cs),
          c2$profit[["coarse"]] > b$uniform_profit, c2$profit[["fine"]] > b$uniform_profit,
          c2$cs[["coarse"]] > 0, c2$cs[["fine"]] > 0,
          c2$profit[["coarse"]] + c2$cs[["coarse"]] < b$efficient, c2$profit[["fine"]] + c2$cs[["fine"]] < b$efficient,
          is.numeric(b$uniform_gap), length(b$uniform_gap) == 1L,
          abs(b$uniform_gap - (b$uniform_profit - c2$profit[["uniform"]])) < 1e-12,
          b$uniform_gap > 0, b$uniform_gap < 0.01 * b$uniform_profit, c2$prices$uniform != b$uniform_price)
# the output test: when discrimination opens a market the single price leaves unserved, output and welfare rise
lo <- c2$linear_open
stopifnot(lo$n_served_uniform < lo$n_served_pd, lo$Q_pd > lo$Q_uniform, lo$W_pd > lo$W_uniform,
          c2$linear$profit_pd > c2$linear$profit_uniform, c2$linear$CS_pd < c2$linear$CS_uniform)
# Monte Carlo over fresh populations and splits: the profit ranking holds in (almost) every replication, and the
# consumer-surplus loss and the gain in merchants served are not artefacts of one draw
mc <- c2$mc$summary
q <- function(stat, col) mc[[col]][mc$stat == stat]
stopifnot(c2$mc$R == tr$params$C2$R_mc, tr$params$C2$R_mc >= 200,
          q("profit_gain_coarse", "lo") > 0, q("profit_fine_over_coarse", "lo") > -0.005,
          q("cs_change_fine", "hi") < 0, q("served_change_fine", "lo") > 0)

## ---- C3: bundling ----------------------------------------------------------------------------------------------
stopifnot(identical(c3$J, tr$C3$J), all(diff(c3$rev_per_good) > 0), all(diff(c3$dwl_share) < 0),
          abs(c3$rev_per_good[1] - 1/4) < 0.005, abs(c3$dwl_share[1] - 1/4) < 0.01, all(c3$rev_per_good < 1/2),
          all(diff(1/2 - c3$rev_per_good) < 0), all(diff(c3$share_buy) > 0))

## ---- C4: price leadership ----------------------------------------------------------------------------------------
# whoever leads, both prices exceed the simultaneous prices and both firms earn more (the leader can always
# repeat its Bertrand price; prices are strategic complements, so it raises it and the follower follows)
stopifnot(all(c4$p_seq_garnet > c4$p_sim), all(c4$profit_seq > c4$profit_sim), all(c4$profit_seq_garnet > c4$profit_sim))
# calibration-specific, not a general result: at this baseline both sequential price vectors stay below the joint prices
stopifnot(all(c4$p_seq < c4$p_joint), all(c4$p_seq_garnet < c4$p_joint))
