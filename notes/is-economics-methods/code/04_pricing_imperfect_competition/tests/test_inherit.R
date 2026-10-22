# ------------------------------------------------------------------
# tests/test_inherit.R -- block C0 reproduces every number of the legacy chapter-4 script
# (snapshot: dev markdown/upgrade_2026-10-19/sdd_partI/legacy/code/04_pricing_imperfect_competition/01_pricing.R)
# and of the legacy text's running case (section 6 and the figures of 3.1 and 4.1), to the precision they print,
# and the card's C0 list: monopoly 7 and 25, CS and DWL 12.5, tiered menu 26, perfect discrimination and
# two-part tariff 50, Cournot 5.33 and 22.2 (200/9), Bertrand 2, Hotelling price 3 and profit 0.5, bundling
# 20 against 24, double marginalisation 9.5 and 18.75, the Hotelling's Rule path.
# Run from the chapter code directory after 01_dgp.R (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("pricing_fns.R")
C0 <- readRDS("truth04.rds")$C0
at <- function(x, v, d) length(x) == length(v) && all(abs(round(x, d) - v) < 1e-9)   # equal at d printed decimals
ex <- function(x, v) length(x) == length(v) && all(abs(x - v) < 1e-12)                # exact up to rounding error

## ---- the Cirrus market: p = 12 - q, c = 2 ---------------------------------------------------
stopifnot(C0$a == 12, C0$b == 1, C0$cost == 2)

## ---- legacy 6.1 and 3.1: single-price monopoly, CS, DWL, Lerner ----------------------------------
stopifnot(ex(C0$qm, 5), ex(C0$pm, 7), ex(C0$profit_m, 25), ex(C0$qc, 10), ex(C0$cs_m, 12.5), ex(C0$dwl_m, 12.5),
          at(C0$lerner, 0.7142857, 7), at(C0$inv_elas, 0.7142857, 7), at(C0$lerner, 0.714, 3), ex(C0$lerner, 5/7),
          ex(C0$elas_abs, 1.4), abs(C0$lerner - C0$inv_elas) < 1e-12, ex(C0$mr_at_qm, 2),
          ex(C0$ts_competitive, 50), ex(C0$ts_m, 37.5), ex(C0$ts_m + C0$dwl_m, C0$ts_competitive))
# the same numbers from the chapter's functions: n = 1 Cournot is the monopoly
stopifnot(ex(cournot_n(12, 1, 2, 1)$p, C0$pm), ex(cournot_n(12, 1, 2, 1)$profit, C0$profit_m))

## ---- legacy 6.2: third-degree discrimination, the inverse elasticity rule ---------------------------
stopifnot(ex(C0$seg_a, c(12, 8)), ex(C0$seg_q, c(5, 3)), ex(C0$seg_p, c(7, 5)),
          at(C0$seg_elas, c(1.4, 1.6666667), 7), at(C0$seg_elas, c(1.4, 1.67), 2),
          at(C0$seg_markup, c(0.7142857, 0.6), 7), ex(C0$seg_markup, c(5/7, 3/5)), ex(C0$seg_elas, c(7/5, 5/3)),
          C0$seg_markup[1] > C0$seg_markup[2], C0$seg_elas[1] < C0$seg_elas[2])

## ---- legacy 6.3: the two-type tiered commission menu ------------------------------------------------
stopifnot(C0$thH == 10, C0$thL == 6, C0$mu == 0.5)
stopifnot(ex(C0$single_both, 18), ex(C0$single_high, 25), ex(C0$single_bench, 25),
          ex(C0$unit_p, 4), ex(C0$unit_profit, 16), C0$unit_profit < C0$single_bench)
stopifnot(ex(C0$qL, 2), ex(C0$qH, 10), ex(C0$qL_eff, 6), ex(C0$rentH, 8), ex(C0$TL, 10), ex(C0$TH, 42),
          ex(C0$profit_menu, 26), ex(C0$profit_fb, 34), ex(C0$info_cost, 8), ex(C0$gross_H, 50), ex(C0$fb_L, 18),
          ex(C0$uH_own, 8), ex(C0$uH_mimic, 8), ex(C0$uL_own, 0), C0$profit_menu > C0$single_bench)
# the legacy script's optimize() benchmarks agree with the closed forms
unit_rev <- function(p) p * (0.5 * pmax(10 - p, 0) + 0.5 * pmax(6 - p, 0))
uo <- optimize(unit_rev, c(0, 10), maximum = TRUE)
stopifnot(abs(uo$maximum - C0$unit_p) < 1e-3, abs(uo$objective - C0$unit_profit) < 1e-6,
          abs(optimize(function(q) 6 * q - q^2 / 2, c(0, 12), maximum = TRUE)$objective - C0$single_both) < 1e-6)

## ---- legacy 6.4 and 4.1: Cournot, Bertrand, Hotelling ------------------------------------------------
stopifnot(ex(C0$cournot_q, 10/3), ex(C0$cournot_Q, 20/3), ex(C0$cournot_p, 16/3), ex(C0$cournot_profit_each, 100/9),
          ex(C0$cournot_profit_total, 200/9), at(C0$cournot_p, 5.333333, 6), at(C0$cournot_p, 5.33, 2),
          at(C0$cournot_profit_each, 11.111111, 6), at(C0$cournot_profit_each, 11.11, 2),
          at(C0$cournot_profit_total, 22.2, 1), C0$pm > C0$cournot_p, C0$cournot_p > C0$bertrand_p)
stopifnot(ex(C0$br_intercept, 5), ex(C0$br_slope, -1/2), ex(C0$collusion_q_each, 5/2), ex(C0$competitive_q_each, 5),
          ex(C0$br_intercept + C0$br_slope * C0$cournot_q, C0$cournot_q))            # the equilibrium is on both best responses
stopifnot(ex(C0$bertrand_p, 2), ex(C0$bertrand_profit, 0), ex(C0$hotelling_t, 1), ex(C0$hotelling_p, 3),
          ex(C0$hotelling_profit, 0.5))
stopifnot(ex(cournot_n(12, 1, 2, 2)$p, C0$cournot_p), ex(cournot_n(12, 1, 2, 2)$profit, C0$cournot_profit_each))

## ---- legacy 6.5 and 3.4: two-part tariff and bundling -----------------------------------------------
stopifnot(ex(C0$tpt_q, 10), ex(C0$tpt_F, 50), ex(C0$tpt_profit, 50), ex(C0$tpt_profit_single, 25),
          ex(C0$tpt_profit / C0$tpt_profit_single, 2))
stopifnot(ex(C0$vA, c(10, 2)), ex(C0$vB, c(2, 10)), ex(C0$rev_A, 10), ex(C0$rev_B, 10), ex(C0$rev_separate, 20),
          ex(C0$bundle_price, 12), ex(C0$rev_bundle, 24), C0$rev_bundle > C0$rev_separate)
# the same numbers from bundle_opt (two consumers: revenue per consumer times 2)
V <- cbind(C0$vA, C0$vB)
stopifnot(ex(2 * bundle_opt(V)$rev, C0$rev_bundle), ex(2 * (bundle_opt(V[, 1])$rev + bundle_opt(V[, 2])$rev), C0$rev_separate))

## ---- legacy 6.6: the rent path of Hotelling's Rule ------------------------------------------------------
stopifnot(ex(C0$rate, 0.05), ex(C0$rent0, 4), ex(C0$cext, 2), ex(C0$hot_t, 0:2),
          at(C0$hot_rent, c(4, 4.2, 4.41), 4), at(C0$hot_price, c(6, 6.2, 6.41), 4), at(C0$hot_pv, c(4, 4, 4), 4),
          ex(C0$hot_pv, c(4, 4, 4)), ex(C0$hot_rent[-1] / C0$hot_rent[-3], c(1.05, 1.05)))

## ---- legacy 6.7: double marginalisation -------------------------------------------------------------------
stopifnot(ex(C0$vert_pm, 7), ex(C0$vert_profit_int, 25), ex(C0$vert_w, 7), ex(C0$vert_qd, 2.5), ex(C0$vert_pd, 9.5),
          ex(C0$vert_profit_up, 12.5), ex(C0$vert_profit_down, 6.25), ex(C0$vert_profit_joint, 18.75),
          ex(C0$vert_cs, 3.125), ex(C0$vert_cs_int, 12.5), C0$vert_pd > C0$vert_pm, C0$vert_profit_joint < C0$vert_profit_int)
stopifnot(abs(optimize(function(w) (w - 2) * (12 - w) / 2, c(2, 12), maximum = TRUE)$maximum - C0$vert_w) < 1e-3)

## ---- legacy 4.8: linear-demand pass-through, monopoly to Cournot n = 8 -------------------------------------
stopifnot(ex(C0$pt_n, 1:8), ex(C0$pt_cournot, (1:8) / (2:9)), at(C0$pt_cournot[c(1, 8)], c(0.5, 0.889), 3),
          all(diff(C0$pt_cournot) > 0), ex(C0$pt_cournot, passthrough_cournot(1:8)))

## ---- section 1 summary: the profit range of the Cirrus market -------------------------------------------------
# single price 25; menu 26; full information 34 (two types); two-part tariff 50; Cournot 22.2; Bertrand 0
stopifnot(C0$profit_m == 25, C0$profit_menu == 26, C0$profit_fb == 34, C0$tpt_profit == 50,
          C0$bertrand_profit == 0, C0$tpt_profit == C0$cs_m + C0$profit_m + C0$dwl_m)
