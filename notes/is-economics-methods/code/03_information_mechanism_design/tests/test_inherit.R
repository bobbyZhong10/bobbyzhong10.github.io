# ------------------------------------------------------------------
# tests/test_inherit.R -- block V0 reproduces every number of the legacy chapter-3 script
# (snapshot: dev markdown/upgrade_2026-10-19/sdd_partI/legacy/code/03_information_mechanism_design/01_mechanisms.R)
# and of the legacy chapter-3 text, to the precision they print, and the card's list for V0: Vireo's
# 1240 and 680; the GSP locally envy-free equilibrium; 1/3 and 1/2 for n = 2, 3; 5/12 under the reserve 1/2;
# persuasion x = 3/7 and 0.6; Myerson-Satterthwaite's 9/32; moral hazard's 25/3 and 5/3; the 3 x 3 DA.
# Additions (Ruling P21 moves the cheap-talk interval count here): the cheap-talk N(b) and one partition;
# the bilateral-trade gains behind 9/32; the moral-hazard principal's payoffs.
# Run from the chapter code directory after 01_dgp.R (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("mech_fns.R")
V0 <- readRDS("truth03.rds")$V0
at <- function(x, v, d) length(x) == length(v) && all(abs(round(x, d) - v) < 1e-9)   # equal at d printed decimals

## ---- legacy section 1: the product manager's pay-your-bid estimate -----------------------
stopifnot(identical(V0$v, c(10, 6, 4)), identical(V0$ctr, c(100, 40, 0)),
          at(V0$pm_top, 1000, 9), at(V0$pm_second, 240, 9), at(V0$pm_naive, 1240, 9))

## ---- legacy 6.1: efficient allocation and VCG payments -------------------------------------
stopifnot(at(V0$vcg_others_without_1, 760, 9), at(V0$vcg_others_with_1, 240, 9), at(V0$vcg_pay1, 520, 9),
          at(V0$vcg_others_without_2, 1160, 9), at(V0$vcg_others_with_2, 1000, 9), at(V0$vcg_pay2, 160, 9),
          at(V0$vcg_pay3, 0, 9), at(V0$vcg_per_click_1, 5.2, 9), at(V0$vcg_per_click_2, 4, 9),
          at(V0$vcg_rev, 680, 9), at(V0$vcg_util1, 480, 9), at(V0$vcg_util2, 80, 9), V0$vcg_util1 >= 0, V0$vcg_util2 >= 0)
pa <- vcg_position(V0$v, V0$ctr[1:2])                                    # the same numbers from the chapter's function
stopifnot(isTRUE(all.equal(pa$pay, c(V0$vcg_pay1, V0$vcg_pay2, V0$vcg_pay3))), pa$revenue == V0$vcg_rev)

## ---- legacy 6.2: GSP, truthful bids are a locally envy-free equilibrium; the lowest one is VCG ----
stopifnot(at(V0$gsp_truth, 760, 9), gsp_revenue(V0$v, V0$ctr[1:2]) == V0$gsp_truth,
          at(V0$gsp_dev_adv2_up, 0, 9), at(V0$gsp_now_adv2, 80, 9), at(V0$gsp_dev_adv3_up, 0, 9), at(V0$gsp_now_adv3, 0, 9),
          at(V0$gsp_dev_adv1_down, 240, 9), at(V0$gsp_now_adv1, 400, 9), at(V0$gsp_dev_adv2_exit, 0, 9),
          V0$gsp_dev_adv2_up <= V0$gsp_now_adv2, V0$gsp_dev_adv3_up <= V0$gsp_now_adv3,
          V0$gsp_dev_adv1_down <= V0$gsp_now_adv1, V0$gsp_dev_adv2_exit <= V0$gsp_now_adv2, isTRUE(V0$gsp_truth_lef),
          at(V0$gsp_p1_low, 5.2, 9), at(V0$gsp_low, 680, 9), abs(V0$gsp_low - V0$vcg_rev) < 1e-9)

## ---- legacy 4.2 and 6.3: revenue equivalence, uniform values, n = 2 and 3 ---------------------
stopifnot(identical(V0$re_n, c(2, 3)), at(V0$re_spa, c(1 / 3, 1 / 2), 12), at(V0$re_fpa, c(1 / 3, 1 / 2), 12),
          at(V0$re_spa, c(0.3333, 0.5), 4), at(V0$re_spa, c(0.3333333, 0.5), 7), max(abs(V0$re_spa - V0$re_fpa)) < 1e-12)
for (n in 2:3) stopifnot(abs(fpa_bid(0.7, n, punif) - (n - 1) / n * 0.7) < 1e-9,
                         abs(expected_revenue(n, value_dist(list(family = "uniform", lower = 0, upper = 1))) - (n - 1) / (n + 1)) < 1e-9)

## ---- legacy 4.3 and 6.3: the optimal reserve 1/2 and 5/12 --------------------------------------
stopifnot(at(V0$res_r, 0.5, 12), at(V0$res_rev, 5 / 12, 12), at(V0$res_rev, 0.4166667, 7), at(V0$res_rev, 0.4167, 4),
          at(V0$res_rev_no_reserve, 1 / 3, 12), at(V0$res_unsold, 1 / 4, 12),
          abs(optimal_reserve(value_dist(list(family = "uniform", lower = 0, upper = 1))) - V0$res_r) < 1e-9)

## ---- legacy 4.6 (script): first-price shading 1/n at v = 1 falls with n ------------------------
stopifnot(identical(V0$shading_n, 2:8), at(V0$shading, c(0.5, 0.333, 0.25, 0.2, 0.167, 0.143, 0.125), 3),
          all(diff(V0$shading) < 0))

## ---- legacy 3.1: Akerlof, the market unravels ------------------------------------------------------
stopifnot(at(V0$akerlof_buyer_mult, 1.5, 9), at(V0$akerlof_wtp_slope, 0.75, 9), isTRUE(V0$akerlof_unravels),
          at(V0$akerlof_p, seq(0.1, 1, by = 0.1), 9), at(V0$akerlof_wtp, 0.75 * seq(0.1, 1, by = 0.1), 9))

## ---- legacy 3.4: moral hazard, two effort levels, limited liability ---------------------------------
stopifnot(at(V0$mh_cost, 5, 9), at(V0$mh_p_high, 0.8, 9), at(V0$mh_p_low, 0.2, 9), at(V0$mh_value, 100, 9),
          at(V0$mh_ic_coef, 0.6, 9), at(V0$mh_gap, 25 / 3, 12), at(V0$mh_wH, 8.33, 2), at(V0$mh_wH, 8.333, 3), at(V0$mh_wL, 0, 9),
          at(V0$mh_rent, 5 / 3, 12), at(V0$mh_rent, 1.67, 2), at(V0$mh_rent, 1.667, 3),
          abs(V0$mh_ic_coef * (V0$mh_wH - V0$mh_wL) - V0$mh_cost) < 1e-12)
# additions: the principal's payoffs (high effort induced 220/3, low effort 20, effort observable 75); the agency
# cost equals the limited-liability rent
stopifnot(at(V0$mh_principal_high, 220 / 3, 12), at(V0$mh_principal_low, 20, 12), at(V0$mh_principal_fb, 75, 12),
          at(V0$mh_agency_cost, 5 / 3, 12))

## ---- legacy 3.7 and 6.4: Bayesian persuasion ------------------------------------------------------------
stopifnot(at(V0$pers_prior, 0.3, 9), at(V0$pers_thresh, 0.5, 9), at(V0$pers_x, 3 / 7, 12), at(V0$pers_x, 0.4285714, 7),
          at(V0$pers_x, 0.429, 3), at(V0$pers_p_buy, 0.6, 12), at(V0$pers_cav, 0.6, 12), at(V0$pers_no_info, 0, 12),
          at(V0$pers_posterior_rec, 0.5, 12), at(V0$pers_full_disclosure, 0.3, 12),
          abs(persuasion_x(V0$pers_prior, V0$pers_thresh) - V0$pers_x) < 1e-15)

## ---- legacy 4.5: Myerson-Satterthwaite, uniform values and costs ----------------------------------------
stopifnot(at(V0$ms_p_trade_fb, 1 / 2, 12), at(V0$ms_p_trade_sb, 9 / 32, 12), at(V0$ms_p_trade_sb, 0.281, 3),
          at(V0$ms_wedge, 1 / 4, 12))
# additions: the linear double-auction equilibrium (seller asks 2c/3 + 1/4, buyer bids 2v/3 + 1/12) trades iff
# v >= c + 1/4; gains 9/64 against 1/6, a share of 27/32; the efficient VCG mechanism runs an expected deficit of 1/6
stopifnot(at(V0$ms_gains_fb, 1 / 6, 12), at(V0$ms_gains_sb, 9 / 64, 12), at(V0$ms_gains_share, 27 / 32, 12),
          at(V0$ms_vcg_deficit, 1 / 6, 12))
pt <- integrate(function(c) sapply(c, function(cc) max(0, 1 - (cc + 1 / 4))), 0, 1)$value
gs <- integrate(function(c) sapply(c, function(cc) if (cc + 1 / 4 >= 1) 0 else ((1 - cc)^2 - (1 / 4)^2) / 2), 0, 1)$value
stopifnot(abs(pt - V0$ms_p_trade_sb) < 1e-8, abs(gs - V0$ms_gains_sb) < 1e-8)
ask <- function(c) 2 * c / 3 + 1 / 4; bid <- function(v) 2 * v / 3 + 1 / 12
stopifnot(abs(bid(0.75) - ask(0.5)) < 1e-12, bid(0.76) > ask(0.5), bid(0.74) < ask(0.5))     # trade iff v - c >= 1/4

## ---- legacy 3.8 (Ruling P21): cheap talk, the largest N with N (N - 1) < 1 / (2 b) -----------------------------
stopifnot(at(V0$ct_b, c(0.01, 0.05, 0.1, 0.2, 0.25), 9), identical(V0$ct_N, c(7L, 3L, 2L, 2L, 1L)),
          identical(V0$ct_N, cheap_talk_N(V0$ct_b)),
          all(V0$ct_N * (V0$ct_N - 1) < 1 / (2 * V0$ct_b)), all((V0$ct_N + 1) * V0$ct_N >= 1 / (2 * V0$ct_b)),
          at(V0$ct_partition_b0p05, c(0, 2 / 15, 7 / 15, 1), 12))

## ---- legacy 6.5: deferred acceptance on the 3 x 3 talent market ---------------------------------------------
stopifnot(identical(V0$da_worker_pref, list(c(1, 2, 3), c(2, 1, 3), c(1, 2, 3))),
          identical(V0$da_firm_pref, list(c(2, 1, 3), c(1, 2, 3), c(1, 2, 3))),
          identical(V0$da_worker_opt, c(1L, 2L, 3L)), identical(V0$da_firm_opt, c(2L, 1L, 3L)),
          is_stable(V0$da_worker_opt, V0$da_worker_pref, V0$da_firm_pref),
          is_stable(V0$da_firm_opt, V0$da_worker_pref, V0$da_firm_pref),
          identical(partner_rank(V0$da_worker_opt, V0$da_worker_pref), c(1, 1, 3)),
          identical(partner_rank(V0$da_firm_opt, V0$da_worker_pref), c(2, 2, 3)))
