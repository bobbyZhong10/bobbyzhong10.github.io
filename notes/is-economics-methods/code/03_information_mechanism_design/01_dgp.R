# ------------------------------------------------------------------
# 01_dgp.R -- the Vireo world for chapter 3 (information economics, contracts, mechanism design)
#
# Vireo is a (fictional) search and ad-exchange platform. Money is in US dollars; auction values and bids
# are USD per click. Five blocks:
#   V0   by hand (inherited from the legacy chapter; no random numbers): the two-slot position auction
#        (values 10, 6, 4 USD per click; 100 and 40 clicks), VCG 680 against the pay-your-bid 1240, GSP
#        and its lowest locally envy-free equilibrium; uniform revenue equivalence (1/3, 1/2) and the
#        reserve 1/2 (5/12); first-price shading; Akerlof; two-action moral hazard (25/3, 5/3); two-state
#        persuasion (3/7, 0.6); Myerson-Satterthwaite (9/32); the 3 x 3 deferred acceptance. Additions:
#        the cheap-talk interval count (Ruling P21), the bilateral-trade gains and the VCG deficit behind
#        9/32, the moral-hazard principal's payoffs, Theorem 2's signaling example (1, sqrt 2, 1/2).
#   V1   ad-exchange auctions (new): 2,000 auctions for each n in {2, ..., 6} (10,000 auctions); values iid
#        from a normal truncated to [0, 3] USD per click (mean parameter 1.3, sd parameter 0.6: smooth,
#        log-concave, so regular, with positive density on the whole support, assumption I.1). The same
#        value draws are run under four formats: first price (exact equilibrium bids by numerical
#        integration), second price (truthful), and both with the optimal reserve r* (Psi(r*) = 0).
#        02_auction.R compares revenues and inverts the first-price bids by quantile GPV.
#   V2   SaaS menus (new, deterministic): a usage tier q (thousand API calls a month) for customers whose
#        type theta (USD per thousand calls) lies in [2, 10]; serving cost q^2 / 2 USD. A regular type
#        density (normal(6, 2) truncated) and a bimodal one that needs ironing (0.7 normal(4, 0.7) +
#        0.3 normal(8, 0.7), truncated); grid of 2,001 types. Two-type version for the IR/IC picture:
#        theta_L = 4, theta_H = 7, share of high types mu = 0.4.
#   V3   moderator contract (new, deterministic plus a small simulation in 04_contract.R): Vireo's
#        ad-review moderators; output x = e + eps (thousand USD a month, sd 2), CARA A = 0.5, c(e) = e^2 / 2.
#        Extra signals: the user-report rate on approved ads (falls with effort: y = -0.5 e + eta, sd 1),
#        a peer team's output (no effort content, noise correlated 0.6 with eps), a pure-noise metric.
#   V4   talent matching market (new): 200 workers (proposers) and 200 firms (receivers); utilities are a
#        common component plus an idiosyncratic one with equal variances; complete strict lists.
# Every number in chapter 3 comes from these data, truth03.rds and the res_*.rds files.
#
# Outputs (in this directory):
#   vireo_V1.rds   data.table: auction n bidder v_true b_fp b_fp_res (one row per bidder; b_fp_res is NA below r*)
#   truth03.rds    list(V0, V1, V2, V3, V4, params); V1$F (value distribution parameters), V1$v (values),
#                  V2$theta_dist, V4$prefs as the card asks
# Columns ending in _true are simulation truth, not observable data.
#
# Seeds (every simulated block its own set.seed at its start): V1 3101, V1 Monte Carlo 3102 (02_auction.R),
# V3 simulation 3103 (04_contract.R), V4 3104, V4 manipulation draws 3105 (05_matching.R). The block
# parameters were calibrated in a scratch session outside this directory (GPV accuracy against the sample
# size and h on seed 5; the size of the set of stable matchings against the common share on seed 11); a
# scratch check of the manipulation logic ran on 3104 and 3105. No seed was tuned.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("mech_fns.R")

seed <- list(V1 = 3101, V1_mc = 3102, V3_sim = 3103, V4 = 3104, V4_manip = 3105)

## =================================================================
## V0: by hand (legacy sections 1, 3.1, 3.4, 3.7, 3.8, 4.2 to 4.6, 6.1 to 6.5; no random numbers)
## =================================================================
V0 <- local({
  ## Vireo's two slots: 100 and 40 clicks a month; values 10, 6, 4 USD per click (sorted)
  v <- c(10, 6, 4); ctr <- c(100, 40, 0)                       # third entry: no slot
  o <- list(v = v, ctr = ctr)
  ## section 1: pay-your-bid with truthful reports (the product manager's estimate)
  o$pm_top <- v[1] * ctr[1]; o$pm_second <- v[2] * ctr[2]; o$pm_naive <- o$pm_top + o$pm_second
  ## 6.1: VCG payment = others' best total without i - others' total under the chosen allocation
  o$vcg_others_without_1 <- v[2] * ctr[1] + v[3] * ctr[2]; o$vcg_others_with_1 <- v[2] * ctr[2] + v[3] * ctr[3]
  o$vcg_pay1 <- o$vcg_others_without_1 - o$vcg_others_with_1
  o$vcg_others_without_2 <- v[1] * ctr[1] + v[3] * ctr[2]; o$vcg_others_with_2 <- v[1] * ctr[1] + v[3] * ctr[3]
  o$vcg_pay2 <- o$vcg_others_without_2 - o$vcg_others_with_2
  o$vcg_pay3 <- (v[1] * ctr[1] + v[2] * ctr[2]) - (v[1] * ctr[1] + v[2] * ctr[2])   # advertiser 3 changes nothing
  o$vcg_per_click_1 <- o$vcg_pay1 / ctr[1]; o$vcg_per_click_2 <- o$vcg_pay2 / ctr[2]
  o$vcg_rev <- o$vcg_pay1 + o$vcg_pay2 + o$vcg_pay3
  o$vcg_util1 <- v[1] * ctr[1] - o$vcg_pay1; o$vcg_util2 <- v[2] * ctr[2] - o$vcg_pay2
  ## 6.2: GSP with truthful bids: slot 1 pays 6 per click, slot 2 pays 4 per click
  o$gsp_truth <- v[2] * ctr[1] + v[3] * ctr[2]
  o$gsp_now_adv1 <- (v[1] - v[2]) * ctr[1]; o$gsp_now_adv2 <- (v[2] - v[3]) * ctr[2]; o$gsp_now_adv3 <- 0
  o$gsp_dev_adv2_up <- (v[2] - v[2]) * ctr[1]                  # advertiser 2 takes slot 1 at the price 6
  o$gsp_dev_adv3_up <- (v[3] - v[3]) * ctr[2]                  # advertiser 3 takes slot 2 at the price 4
  o$gsp_dev_adv1_down <- (v[1] - v[3]) * ctr[2]                # advertiser 1 drops to slot 2 at the price 4
  o$gsp_dev_adv2_exit <- 0                                     # advertiser 2 drops out
  o$gsp_truth_lef <- o$gsp_dev_adv2_up <= o$gsp_now_adv2 && o$gsp_dev_adv3_up <= o$gsp_now_adv3 &&
    o$gsp_dev_adv1_down <= o$gsp_now_adv1 && o$gsp_dev_adv2_exit <= o$gsp_now_adv2
  ## the lowest-revenue locally envy-free equilibrium: slot 2 at v3, slot 1 at the price that leaves advertiser 2 indifferent
  o$gsp_p1_low <- (v[2] * (ctr[1] - ctr[2]) + v[3] * ctr[2]) / ctr[1]
  o$gsp_low <- o$gsp_p1_low * ctr[1] + v[3] * ctr[2]
  ## 4.2 and 6.3: revenue equivalence with values uniform on [0, 1]
  n <- c(2, 3)
  o$re_n <- n; o$re_spa <- (n - 1) / (n + 1); o$re_fpa <- (n - 1) / n * n / (n + 1)
  ## 4.3 and 6.3: the optimal reserve r* = 1/2 solves r = (1 - F(r)) / f(r) = 1 - r; n = 2
  r <- 0.5
  o$res_r <- r
  o$res_rev <- 2 * r * (1 - r) * r + (1 - r)^2 * (r + (1 - r) / 3)    # one bidder above r pays r; both: E[min | both >= r]
  o$res_rev_no_reserve <- 1 / 3; o$res_unsold <- r^2
  ## 4.6 (legacy script): first-price shading v - b(v) = v / n at v = 1
  o$shading_n <- 2:8; o$shading <- 1 - (o$shading_n - 1) / o$shading_n
  ## 3.1: Akerlof; quality uniform on [0, 1], sellers value theta, buyers 1.5 theta
  o$akerlof_buyer_mult <- 1.5; o$akerlof_p <- seq(0.1, 1, by = 0.1)
  o$akerlof_wtp <- o$akerlof_buyer_mult * (o$akerlof_p / 2); o$akerlof_wtp_slope <- o$akerlof_buyer_mult / 2
  o$akerlof_unravels <- all(o$akerlof_wtp < o$akerlof_p)
  ## 3.2: signaling (Theorem 2's hand example); theta_L = 1, theta_H = 2, cost c(s, theta) = s^2 / theta; the ends
  ## s_low, s_bar of the separating interval solve theta_H - c(s, theta) = theta_L for theta_L and for theta_H
  o$sig_theta_L <- 1; o$sig_theta_H <- 2
  o$sig_s_low <- sqrt((o$sig_theta_H - o$sig_theta_L) * o$sig_theta_L)
  o$sig_s_bar <- sqrt((o$sig_theta_H - o$sig_theta_L) * o$sig_theta_H)
  o$sig_cost_H <- o$sig_s_low^2 / o$sig_theta_H                # the high type's cost at s_low
  ## 3.4: moral hazard; high effort costs 5, success 0.8 (0.2 under low effort), success worth 100; w >= 0
  o$mh_cost <- 5; o$mh_p_high <- 0.8; o$mh_p_low <- 0.2; o$mh_value <- 100
  o$mh_ic_coef <- o$mh_p_high - o$mh_p_low; o$mh_gap <- o$mh_cost / o$mh_ic_coef
  o$mh_wL <- 0; o$mh_wH <- o$mh_wL + o$mh_gap
  o$mh_rent <- o$mh_p_high * o$mh_wH + (1 - o$mh_p_high) * o$mh_wL - o$mh_cost
  o$mh_principal_high <- o$mh_p_high * (o$mh_value - o$mh_wH) - (1 - o$mh_p_high) * o$mh_wL
  o$mh_principal_low <- o$mh_p_low * o$mh_value                # low effort needs no incentive pay
  o$mh_principal_fb <- o$mh_p_high * o$mh_value - o$mh_cost     # effort observable: pay its cost
  o$mh_agency_cost <- o$mh_principal_fb - o$mh_principal_high
  ## 3.7 and 6.4: persuasion; prior P(good) = 0.3, the buyer buys iff the posterior is at least 0.5
  o$pers_prior <- 0.3; o$pers_thresh <- 0.5
  o$pers_x <- persuasion_x(o$pers_prior, o$pers_thresh)
  o$pers_p_buy <- o$pers_prior + (1 - o$pers_prior) * o$pers_x
  o$pers_posterior_rec <- o$pers_prior / (o$pers_prior + (1 - o$pers_prior) * o$pers_x)
  o$pers_cav <- o$pers_prior / o$pers_thresh                   # concave closure 2 mu on [0, 0.5], at the prior
  o$pers_no_info <- 0; o$pers_full_disclosure <- o$pers_prior
  ## 4.5: Myerson-Satterthwaite with v, c uniform on [0, 1]; the linear double-auction equilibrium trades iff v >= c + 1/4
  wd <- 1 / 4
  o$ms_p_trade_fb <- 1 / 2; o$ms_wedge <- wd; o$ms_p_trade_sb <- (1 - wd)^2 / 2
  o$ms_gains_fb <- 1 / 6                                       # E[(v - c)+]
  o$ms_gains_sb <- ((1 - wd^3) / 3 - wd^2 * (1 - wd)) / 2       # int_0^{3/4} ((1 - c)^2 - (1/4)^2) / 2 dc
  o$ms_gains_share <- o$ms_gains_sb / o$ms_gains_fb
  o$ms_vcg_deficit <- o$ms_gains_fb                            # efficient pivot mechanism: buyer pays c, seller gets v
  ## 3.8 (Ruling P21): cheap talk, uniform state, quadratic loss
  o$ct_b <- c(0.01, 0.05, 0.1, 0.2, 0.25); o$ct_N <- cheap_talk_N(o$ct_b)
  o$ct_partition_b0p05 <- cheap_talk_partition(0.05, cheap_talk_N(0.05))
  ## 6.5: deferred acceptance, 3 workers and 3 companies
  o$da_worker_pref <- list(c(1, 2, 3), c(2, 1, 3), c(1, 2, 3))
  o$da_firm_pref <- list(c(2, 1, 3), c(1, 2, 3), c(1, 2, 3))
  o$da_worker_opt <- deferred_acceptance(o$da_worker_pref, o$da_firm_pref)                   # worker -> company
  o$da_firm_opt <- invert_match(deferred_acceptance(o$da_firm_pref, o$da_worker_pref), 3L)    # worker view
  o
})
stopifnot(V0$vcg_rev == 680, V0$pm_naive == 1240, abs(V0$gsp_low - V0$vcg_rev) < 1e-9, abs(V0$res_rev - 5 / 12) < 1e-12,
          abs(V0$ms_p_trade_sb - 9 / 32) < 1e-12, abs(V0$pers_x - 3 / 7) < 1e-12)

## =================================================================
## V1: ad-exchange auctions (seed 3101)
## =================================================================
P1 <- list(F = list(family = "truncnorm", mean = 1.3, sd = 0.6, lower = 0, upper = 3),
           n_set = 2:6, n_auc_per_n = 2000L, taus = seq(0.05, 0.95, by = 0.01), h = 0.02, R_mc = 1000L)
set.seed(seed$V1)
V1 <- gen_V1(P1)
D1 <- value_dist(P1$F); rs <- optimal_reserve(D1)
rev_true <- data.frame(n = P1$n_set,
                       none = vapply(P1$n_set, expected_revenue, numeric(1), dist = D1),
                       reserve = vapply(P1$n_set, expected_revenue, numeric(1), dist = D1, reserve = rs))
T1 <- list(F = P1$F, v = V1$v_true, r_star = rs, rev_true = rev_true,
           rev_true_pooled = mean(rev_true$none), rev_reserve_true_pooled = mean(rev_true$reserve),
           p_unsold = D1$cdf(rs)^P1$n_set,
           mean_value = stats::integrate(function(v) v * D1$pdf(v), 0, 3, rel.tol = 1e-12)$value,
           median_value = D1$quantile(0.5), pdf_lower = D1$pdf(0), pdf_upper = D1$pdf(3),
           v_quantile = data.frame(tau = P1$taus, v = D1$quantile(P1$taus)))
cat(sprintf("V1: %d auctions, %d bids; r* = %.4f USD per click; expected revenue n = 2..6: %s (reserve: %s)\n",
            max(V1$auction), nrow(V1), rs, paste(sprintf("%.4f", rev_true$none), collapse = " "),
            paste(sprintf("%.4f", rev_true$reserve), collapse = " ")))

## =================================================================
## V2: SaaS menus (deterministic)
## =================================================================
P2 <- list(theta_lo = 2, theta_hi = 10, n_grid = 2001L, cost2 = 1,
           regular = list(family = "truncnorm", mean = 6, sd = 2, lower = 2, upper = 10),
           irregular = list(family = "truncnorm_mix", weights = c(0.7, 0.3), means = c(4, 8), sds = c(0.7, 0.7),
                            lower = 2, upper = 10),
           two_type = list(theta_L = 4, theta_H = 7, mu = 0.4, cost2 = 1))
T2 <- list(theta_dist = list(regular = P2$regular, irregular = P2$irregular),
           grid = seq(P2$theta_lo, P2$theta_hi, length.out = P2$n_grid), cost2 = P2$cost2, two_type = P2$two_type)

## =================================================================
## V3: the moderators' linear contract (parameters; 04_contract.R computes and simulates with seed 3103)
## =================================================================
P3 <- list(A = 0.5, c2 = 1, sigma2 = 4,
           signals = list(report = list(k = -0.5, var = 1, rho = 0),       # user-report rate: falls with effort
                          peer = list(k = 0, var = 4, rho = 0.6),          # a peer team's output: common shocks only
                          noise = list(k = 0, var = 1, rho = 0)),          # a metric unrelated to effort or to eps
           N_sim = 1000L, R_sim = 500L)
T3 <- P3

## =================================================================
## V4: the 200 x 200 talent market (seed 3104)
## =================================================================
P4 <- list(n_workers = 200L, n_firms = 200L, sd_common = 1, sd_idio = 1,
           n_trunc_prop = 1000L, n_swap_prop = 1000L, n_trunc_recv = 1000L)
set.seed(seed$V4)
g4 <- gen_V4(P4)
T4 <- list(prefs = g4$prefs, common = g4$common)
cat(sprintf("V4: %d workers x %d firms; common share of utility variance %.2f\n", P4$n_workers, P4$n_firms,
            P4$sd_common^2 / (P4$sd_common^2 + P4$sd_idio^2)))

## =================================================================
## save
## =================================================================
truth <- list(V0 = V0, V1 = T1, V2 = T2, V3 = T3, V4 = T4,
              params = list(seed = seed, V1 = P1, V2 = P2, V3 = P3, V4 = P4))
saveRDS(V1, "vireo_V1.rds")
saveRDS(truth, "truth03.rds")
cat("Saved vireo_V1.rds, truth03.rds\n")
