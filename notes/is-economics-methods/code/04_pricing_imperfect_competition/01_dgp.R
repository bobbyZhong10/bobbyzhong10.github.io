# ------------------------------------------------------------------
# 01_dgp.R -- the Cirrus world for chapter 4 (pricing and imperfect competition)
#
# Cirrus is the dominant e-commerce platform in a niche; Garnet is its rival (Ruling P12, names.md). Money is in
# US dollars. Five blocks:
#   C0   the Cirrus market by hand (inherited from the legacy chapter; no random numbers): demand p = 12 - q,
#        marginal cost c = 2; monopoly, third-degree discrimination, the two-type commission menu, Cournot,
#        Bertrand, Hotelling, two-part tariff, bundling, Hotelling's Rule, double marginalisation, the linear
#        pass-through n / (n + 1); plus (new) the Cournot path to competition, the Cournot Lerner check and the
#        monopoly pass-through under three demand curvatures.
#   C1   logit duopoly markets (new): 500 local markets in which Cirrus and Garnet sell differentiated operating
#        services; market-specific marginal costs (cost shifters) and demand shocks; prices under Bertrand-Nash
#        (kappa = 0) and under joint profit maximisation (kappa = 1), solved with nleqslv; a homogeneous-product
#        variant with the C0 demand and the same costs, under Cournot and Bertrand.
#   C2   personalised pricing of Cirrus Pro (new): 20,000 merchants with unit demand; willingness to pay depends
#        on observed features plus noise; a random half is the training half (WTP observed, as in an earlier
#        price-elicitation survey), the other half is priced.
#   C3   bundling (new): 50,000 consumers, independent U[0, 1] values for 100 information goods; bundles of
#        the first J goods, J in {1, 2, 5, 20, 100}.
#   C4   price leadership (new, deterministic): the C1 duopoly at its baseline, simultaneous against sequential
#        pricing; parameters only, 05_leader.R computes.
# Every number in chapter 4 comes from these data, truth04.rds and the res_*.rds files.
#
# Outputs (in this directory):
#   cirrus_C1.rds   data.frame: market c1 c2 xi1_true xi2_true a_hom_true p1_bertrand p2_bertrand s1_bertrand
#                   s2_bertrand p1_joint p2_joint s1_joint s2_joint p_cournot_hom q1_cournot_hom q2_cournot_hom
#                   p_bertrand_hom
#   cirrus_C2.rds   data.frame: customer sales_z tenure_z category wtp train index_true eps_true
#   cirrus_C3.rds   matrix N x 5 of bundle values (columns J1 J2 J5 J20 J100)
#   truth04.rds     list(C0, C1, C2, C3, params)
# Columns ending in _true are simulation truth, not observable data. In C2 the analyst observes wtp only in
# the training half.
#
# Notation (notation_v3.md 5.5): inverse demand p(q), marginal cost c, quantities q and Q, profit pi; logit
# u_ij = xi_j - alpha p_j + eps_ij (the code's delta is xi_j), shares s_j, s_0; conduct kappa; pass-through
# PT_jk = dp_k / dc_j (the code's matrices are indexed [k, j], prices by costs); values v_i, v_ij; bundle size J.
#
# Seeds (every simulated block its own set.seed at its start): C1 4101, C2 4102, C3 4103; the C2 Monte Carlo in
# 03_personalized.R uses 4104. The block parameters were calibrated on exploratory seeds (1 to 3 for C1, 1, 2 and
# 11 to 16 for C2, outside this directory) before these seeds were fixed; the seeds were not tuned afterwards.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("pricing_fns.R")

seed <- list(C1 = 4101, C2 = 4102, C3 = 4103, C2_mc = 4104)

## =================================================================
## C0: the Cirrus market by hand (legacy 6.1 to 6.7 and 4.8; no random numbers)
## =================================================================
## Demand for operating services p = 12 - q (price per unit of operating activity, USD), marginal cost 2 USD.
C0 <- local({
  a <- 12; b <- 1; cost <- 2
  o <- list(a = a, b = b, cost = cost)
  # 6.1 single-price monopoly: MR = 12 - 2q = 2
  qm <- (a - cost) / (2 * b); pm <- a - b * qm; qc <- (a - cost) / b
  o <- c(o, list(qm = qm, pm = pm, profit_m = (pm - cost) * qm, qc = qc, mr_at_qm = a - 2 * b * qm,
                 cs_m = 0.5 * (a - pm) * qm, dwl_m = 0.5 * (qc - qm) * (pm - cost),
                 lerner = (pm - cost) / pm, elas_abs = (1 / b) * pm / qm, inv_elas = qm / pm,
                 ts_competitive = 0.5 * (a - cost) * qc))
  o$ts_m <- o$cs_m + o$profit_m
  # 6.2 third-degree discrimination: segments p = 12 - q and p = 8 - q, MC 2 in both
  seg_a <- c(12, 8); seg_q <- (seg_a - cost) / 2; seg_p <- seg_a - seg_q
  o <- c(o, list(seg_a = seg_a, seg_q = seg_q, seg_p = seg_p, seg_elas = seg_p / seg_q, seg_markup = (seg_p - cost) / seg_p))
  # 6.3 two merchant types theta_H = 10, theta_L = 6, share of high types mu = 1/2, gross value theta q - q^2/2,
  # platform marginal cost 0. Single-package benchmarks, the uniform unit price, the screening menu.
  thH <- 10; thL <- 6; mu <- 0.5
  single_both <- thL * thL - thL^2 / 2                 # q = theta_L priced at the low type's gross value
  single_high <- mu * (thH * thH - thH^2 / 2)          # q = 10 at 50, bought by the high type only
  # unit price p <= theta_L: both buy theta - p, revenue p (mu (thH - p) + (1 - mu)(thL - p)), maximised at p = 4
  unit_p <- (mu * thH + (1 - mu) * thL) / 2
  qH <- thH; qL <- thL - (mu / (1 - mu)) * (thH - thL); rentH <- (thH - thL) * qL
  TL <- thL * qL - qL^2 / 2; TH <- (thH * qH - qH^2 / 2) - rentH
  profit_menu <- mu * TH + (1 - mu) * TL
  profit_fb <- mu * thH^2 / 2 + (1 - mu) * thL^2 / 2
  o <- c(o, list(thH = thH, thL = thL, mu = mu, single_both = single_both, single_high = single_high,
                 single_bench = max(single_both, single_high),
                 unit_p = unit_p, unit_profit = unit_p * (mu * (thH - unit_p) + (1 - mu) * (thL - unit_p)),
                 qL = qL, qH = qH, qL_eff = thL, rentH = rentH, TL = TL, TH = TH, profit_menu = profit_menu,
                 profit_fb = profit_fb, info_cost = profit_fb - profit_menu, gross_H = thH * qH - qH^2 / 2,
                 fb_L = thL^2 / 2, uH_own = thH * qH - qH^2 / 2 - TH, uH_mimic = thH * qL - qL^2 / 2 - TL,
                 uL_own = thL * qL - qL^2 / 2 - TL))
  # 6.4 Cournot duopoly, Bertrand, Hotelling (quadratic transport cost t = 1, firms at the ends)
  cn <- cournot_n(a, b, cost, 2)
  o <- c(o, list(cournot_q = cn$q, cournot_Q = 2 * cn$q, cournot_p = cn$p, cournot_profit_each = cn$profit,
                 cournot_profit_total = 2 * cn$profit,
                 br_intercept = (a - cost) / (2 * b), br_slope = -1 / 2,        # q1 = (10 - q2) / 2
                 collusion_q_each = qm / 2, competitive_q_each = qc / 2,
                 bertrand_p = cost, bertrand_profit = 0,
                 hotelling_t = 1, hotelling_p = cost + 1, hotelling_profit = 1 * 0.5))
  # Cournot and the Lerner index: (p - c) / p = s_i / |E| with s_i = 1/2, |E| = p / Q at the equilibrium
  o$cournot_lerner <- (cn$p - cost) / cn$p
  o$cournot_elas_abs <- (1 / b) * cn$p / (2 * cn$q)
  # 6.5 two-part tariff for homogeneous merchants with demand p = 12 - q; bundling with negatively correlated values
  o <- c(o, list(tpt_q = a - cost, tpt_F = 0.5 * (a - cost) * (a - cost), tpt_profit = 0.5 * (a - cost)^2,
                 tpt_profit_single = (pm - cost) * qm))
  vA <- c(10, 2); vB <- c(2, 10)
  revA <- max(2 * sum(vA >= 2), 10 * sum(vA >= 10)); revB <- max(2 * sum(vB >= 2), 10 * sum(vB >= 10))
  o <- c(o, list(vA = vA, vB = vB, rev_A = revA, rev_B = revB, rev_separate = revA + revB,
                 bundle_price = 12, rev_bundle = 12 * sum((vA + vB) >= 12)))
  # 6.6 Hotelling's Rule: interest rate 0.05, initial scarcity rent 4, extraction cost 2
  rate <- 0.05; rent0 <- 4; cext <- 2; tt <- 0:2
  rent <- rent0 * (1 + rate)^tt
  o <- c(o, list(rate = rate, rent0 = rent0, cext = cext, hot_t = tt, hot_rent = rent, hot_price = cext + rent,
                 hot_pv = rent / (1 + rate)^tt))
  # 6.7 double marginalisation: upstream cost 2, wholesale price w; downstream MR = w gives q = (12 - w) / 2;
  # upstream maximises (w - 2)(12 - w) / 2, so w = 7
  w <- (a + cost) / 2; qd <- (a - w) / 2; pd <- (a + w) / 2
  o <- c(o, list(vert_pm = pm, vert_profit_int = (pm - cost) * qm, vert_w = w, vert_qd = qd, vert_pd = pd,
                 vert_profit_up = (w - cost) * qd, vert_profit_down = (pd - w) * qd,
                 vert_profit_joint = (w - cost) * qd + (pd - w) * qd, vert_cs = 0.5 * (a - pd) * qd,
                 vert_cs_int = 0.5 * (a - pm) * qm))
  # 4.8 linear-demand pass-through, monopoly to Cournot n = 8
  o$pt_n <- 1:8
  o$pt_cournot <- passthrough_cournot(1:8)
  # new: the Cournot path to competition (theorem box), n = 1 to 8 and 100
  o$path_n <- c(1:8, 100)
  path <- lapply(o$path_n, function(n) cournot_n(a, b, cost, n))
  o$path_p <- vapply(path, `[[`, numeric(1), "p")
  o$path_Q <- o$path_n * vapply(path, `[[`, numeric(1), "q")
  o$path_profit_total <- o$path_n * vapply(path, `[[`, numeric(1), "profit")
  o$path_dwl <- 0.5 * (qc - o$path_Q) * (o$path_p - cost)
  # new: monopoly pass-through 1 / (2 + Q p''/p') under three curvatures; constant elasticity e = 2 and 3
  o$ptm_linear <- passthrough_monopoly(qm, -b, 0)
  o$ptm_ce_e <- c(2, 3)
  o$ptm_ce <- vapply(o$ptm_ce_e, function(e) passthrough_monopoly(1, -1 / e, (1 / e) * (1 / e + 1)), numeric(1))
  o$ptm_exp <- passthrough_monopoly(1, -1, 1)
  o
})

## =================================================================
## C1: logit duopoly markets (Cirrus and Garnet)
## =================================================================
## 500 local markets. Logit demand with quality xi_Cirrus = 4, xi_Garnet = 3 and price sensitivity alpha = 1 per
## USD; base marginal cost 2 USD (as C0) for both, plus a cost shifter U(-0.5, 0.5) per firm and market (the
## per-unit payment-processing fee the analyst observes), independent across firms; demand shocks xi_jm ~
## N(0, 0.25^2) added to the qualities, independent of the costs. At the baseline Cirrus holds about 43% of the
## market, Garnet 24%, the outside option 33%: a dominant platform and a smaller rival. Homogeneous-product
## variant: inverse demand p = a_m - Q with a_m = 12 + eta_m, eta_m ~ N(0, 1), the same costs.
P1 <- list(M = 500, delta = c(Cirrus = 4, Garnet = 3), alpha = 1, c_base = c(2, 2), cost_halfwidth = 0.5,
           xi_sd = 0.25, hom_a = 12, hom_b = 1, hom_eta_sd = 1)
set.seed(seed$C1)
C1 <- gen_C1(P1)
truth_C1 <- local({
  pt_m <- function(kappa, cols) {
    lst <- lapply(seq_len(P1$M), function(m)
      passthrough_logit(P1$delta + c(C1$xi1_true[m], C1$xi2_true[m]), P1$alpha, c(C1$c1[m], C1$c2[m]), kappa,
                        p = c(C1[[cols[1]]][m], C1[[cols[2]]][m])))
    list(mean = Reduce(`+`, lst) / P1$M, own1 = vapply(lst, function(x) x[1, 1], numeric(1)))
  }
  ptb <- pt_m(0, c("p1_bertrand", "p2_bertrand")); ptj <- pt_m(1, c("p1_joint", "p2_joint"))
  pb <- bertrand_logit(P1$delta, P1$alpha, P1$c_base); pj <- joint_max_logit(P1$delta, P1$alpha, P1$c_base)
  list(alpha = P1$alpha, delta = P1$delta, c_base = P1$c_base,
       p_bertrand_base = pb, s_bertrand_base = logit_share(pb, P1$delta, P1$alpha),
       p_joint_base = pj, s_joint_base = logit_share(pj, P1$delta, P1$alpha),
       pt_bertrand_base = unname(passthrough_logit(P1$delta, P1$alpha, P1$c_base, 0)),
       pt_joint_base = unname(passthrough_logit(P1$delta, P1$alpha, P1$c_base, 1)),
       pt_bertrand_mean = unname(ptb$mean), pt_joint_mean = unname(ptj$mean),
       pt_own1_bertrand_range = range(ptb$own1), pt_own1_joint_range = range(ptj$own1),
       pt_cournot_hom = 1 / 3)                                      # dp/dc_j in the homogeneous Cournot variant
})

## =================================================================
## C2: personalised pricing of Cirrus Pro
## =================================================================
## Cirrus Pro is a monthly analytics and promotion toolkit for merchants (software: marginal cost 0). 20,000
## merchants, each buys at most one month. Features: standardised log monthly sales and tenure (independent
## N(0, 1)), category apparel 50%, electronics 30%, home 20%. log WTP = 3.2 + 0.45 sales + 0.15 tenure + 0.2
## electronics - 0.15 home + eps, eps ~ N(0, 0.5^2): median WTP about exp(3.2) = 24.5 USD a month; the features
## explain about half of the variance of log WTP. WTP is recorded in whole USD (at least 1), so the value
## distribution is discrete, as BBM's construction needs. A random half (10,000) is the training half; prices
## are set on it and evaluated on the other half. Regimes: one price (K = 1), coarse segments (K = 4 quartiles of
## predicted WTP), fine segments (K = 50, "personalised"); the full ladder K in {1, 2, 3, 4, 5, 10, 20, 50, 100,
## 200} shows the gain and the overfitting of very small cells. Monte Carlo: 200 fresh populations and splits.
P2 <- list(N = 20000, beta = c(const = 3.2, sales = 0.45, tenure = 0.15, electronics = 0.2, home = -0.15), sigma = 0.5,
           cat_levels = c("apparel", "electronics", "home"), cat_probs = c(0.5, 0.3, 0.2), cost = 0,
           K = c(uniform = 1, coarse = 4, fine = 50), K_ladder = c(1, 2, 3, 4, 5, 10, 20, 50, 100, 200), R_mc = 200)
set.seed(seed$C2)
C2 <- gen_C2(P2)
truth_C2 <- list(wtp = C2$wtp, beta = P2$beta, sigma = P2$sigma,
                 efficient = mean(C2$wtp), efficient_eval = mean(C2$wtp[!C2$train]),
                 r2_latent = var(C2$index_true) / var(C2$index_true + C2$eps_true))

## =================================================================
## C3: bundling information goods
## =================================================================
## 50,000 consumers, values for 100 goods independent U[0, 1] (mean 1/2), zero marginal cost; bundles of the first J
## goods, J in {1, 2, 5, 20, 100}. Only the bundle values are saved (the 5 million values are regenerated from the
## seed by gen_C3).
P3 <- list(N = 50000, J = c(1, 2, 5, 20, 100))
set.seed(seed$C3)
C3 <- gen_C3(P3)
truth_C3 <- list(J = P3$J, mean_value = 1 / 2, rev_per_good_separate = 1 / 4, dwl_share_separate = 1 / 4)

## =================================================================
## C4: price leadership (deterministic): the C1 duopoly at its baseline
## =================================================================
P4 <- list(delta = P1$delta, alpha = P1$alpha, c = P1$c_base)

## =================================================================
## Save
## =================================================================
saveRDS(C1, "cirrus_C1.rds")
saveRDS(C2, "cirrus_C2.rds")
saveRDS(C3, "cirrus_C3.rds")
saveRDS(list(C0 = C0, C1 = truth_C1, C2 = truth_C2, C3 = truth_C3,
             params = list(seed = seed, C1 = P1, C2 = P2, C3 = P3, C4 = P4)), "truth04.rds")
cat(sprintf("C0: monopoly p = %g, profit %g; menu %g; Cournot p = %.4f\n", C0$pm, C0$profit_m, C0$profit_menu, C0$cournot_p))
cat(sprintf("C1: %d markets; baseline Bertrand prices %.3f, %.3f; joint %.3f, %.3f\n", P1$M,
            truth_C1$p_bertrand_base[1], truth_C1$p_bertrand_base[2], truth_C1$p_joint_base[1], truth_C1$p_joint_base[2]))
cat(sprintf("C2: %d merchants, median WTP %g USD, mean %.2f USD\n", P2$N, median(C2$wtp), mean(C2$wtp)))
cat(sprintf("C3: %d consumers, J = %s\n", P3$N, paste(P3$J, collapse = ", ")))
cat("Saved cirrus_C1.rds, cirrus_C2.rds, cirrus_C3.rds, truth04.rds\n")
