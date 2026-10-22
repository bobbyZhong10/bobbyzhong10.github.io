# ------------------------------------------------------------------
# 01_dgp.R -- the Tamarind simulation world for chapter 7 (GMM and minimum distance)
#
# Tamarind is a (fictional) US online tutoring marketplace: households buy lesson credits
# (one credit books one lesson; base prices USD 8.99 to 11.49 per credit). Money is in US dollars. One demand
# elasticity, theta = -1.6, is learned three ways (all blocks new in v3; nothing is inherited,
# the legacy Helix GMM and weak-identification numbers retire and are redone in G2 and G3):
#   G1   16 metro markets each run one price A/B test (base price against a price higher by
#        3% to 12%); each test gives a reduced-form price slope and its standard error; the
#        elasticity is the same in every market. Minimum distance combines the 16 slopes.
#   G1h  the same 16 tests when the elasticity differs across markets (less elastic where
#        household income is higher): the MD overidentification test (J) detects it.
#   G2   200 regions, one week: the price paid is partly the regional manager's discretionary
#        promotions (correlated with the demand shock), partly 10 concurrent randomised pricing
#        experiments, which are the instruments. Linear moments, r = 11, k = 2: two-step,
#        iterated and continuously updated GMM, the Windmeijer (2005) correction, identity weights.
#   G3   1,000 regions in which the only experiments are four penny-sized list-price
#        perturbations (0.8 to 2.0 cents per lesson credit): weak moments, Wald against the
#        Stock-Wright S confidence set.
# Every number in chapter 7 comes from these data and from truth07.rds.
#
# Outputs (in this directory):
#   tamarind_G1.rds, tamarind_G1h.rds   data.frame: market B log_price y (one row per household)
#   tamarind_G2.rds                     data.frame: y log_price z1 ... z10 v_true eps_true (one row per region)
#   tamarind_G3.rds                     data.frame: y log_price z1 ... z4 v_true eps_true
#   truth07.rds                         list(G1, G1h, G2, G3, params): truth and parameters
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 1.3, 1.6): market index m = 1, ..., r (r = 16 in G1, so the number of
# markets is the number of reduced-form slopes pi_m); theta_m the market elasticities of G1h; in G2
# and G3 the moments are g_i(theta) = z_i (y_i - beta0 - theta log p_i) with z_i = (1, z_i1, ..., z_iK)
# (r = K + 1); D_0, S_0 as in gmm_fns.R. Field names follow the chapter card (truth07: G1$theta,
# G2$theta, G3$theta).
#
# Seeds (every block its own set.seed, nothing inherited): G1 7101, G1h 7102, G2 7103, G3 7104; the
# Monte Carlo seeds used by 02 to 04 are stored here too (G1 7111, G2 7113 and 73100 + j for grid
# point j, G3 7114).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("gmm_fns.R")

seed <- list(G1 = 7101, G1h = 7102, G2 = 7103, G3 = 7104, G1_mc = 7111, G2_mc = 7113, G2_grid = 73100, G3_mc = 7114)
theta0 <- -1.6                     # the structural price elasticity of tutoring demand (every block)

## =================================================================
## G1: one price A/B test in each of 16 metro markets, homogeneous elasticity
## =================================================================
## Market m: n[m] households (2,000 to 8,000), half at random in arm A (base price p0[m] per
## lesson credit), half in arm B (p0[m] (1 + delta[m]), delta from 3% to 12%). Outcome y: log of the
## household's average weekly lessons booked over the 8-week test. Median weekly lessons at the base
## price q0[m] (3.1 to 4.1); sd of log weekly lessons sigma[m] (0.58 to 0.80) in arm A and 10% larger in
## arm B (het = 1.1). Elasticity theta = -1.6 in every market, so the reduced-form slope of market m,
## pi_m (the OLS slope of y on log price within the test), has pi_m = theta for every m:
## h(theta) = theta iota_r, r = 16, k = 1, 15 overidentifying restrictions.
## The variance of pi_hat_m is exactly (sigma_A^2 / n_A + sigma_B^2 / n_B) / log(1 + delta)^2 (normal
## errors, difference in arm means divided by the log price step): stored as avar.
## Monte Carlo (02_md.R): R_mc = 2000 replications of all 16 tests (seed G1_mc).
P1 <- list(markets = c("Atlanta", "Austin", "Boston", "Charlotte", "Chicago", "Columbus", "Dallas", "Denver", "Houston",
                       "Miami", "Minneapolis", "Nashville", "Philadelphia", "Phoenix", "Portland", "Seattle"),
           n = c(4000, 2400, 6000, 2000, 8000, 2000, 5200, 3600, 6400, 3000, 2800, 2200, 5600, 4400, 2600, 4800),
           delta = c(0.08, 0.05, 0.10, 0.12, 0.06, 0.03, 0.10, 0.08, 0.05, 0.12, 0.10, 0.06, 0.08, 0.04, 0.10, 0.12),
           p0 = c(9.99, 10.49, 11.49, 9.49, 10.99, 8.99, 9.99, 10.49, 9.49, 10.99, 10.49, 9.99, 11.49, 9.49, 10.99, 11.49),
           q0 = c(3.6, 3.8, 3.2, 3.9, 3.4, 4.1, 3.7, 3.5, 3.9, 3.3, 3.6, 3.8, 3.1, 4.0, 3.4, 3.2),
           sigma = c(0.62, 0.70, 0.58, 0.75, 0.60, 0.80, 0.66, 0.64, 0.72, 0.68, 0.62, 0.74, 0.60, 0.78, 0.64, 0.58),
           het = 1.1, theta = theta0, R_mc = 2000,
           # G1h: median household income (USD thousands) and Tamarind active households (thousands) per market
           income = c(72, 86, 94, 68, 78, 64, 70, 88, 62, 66, 84, 74, 76, 66, 82, 104),
           households = c(310, 180, 290, 150, 520, 120, 400, 210, 450, 260, 170, 140, 380, 230, 160, 270),
           kappa = 3)
r1 <- length(P1$markets)
set.seed(seed$G1)
G1 <- gen_G1(P1)
avar <- with(P1, sigma^2 * (1 / (n / 2) + het^2 / (n / 2)) / log(1 + delta)^2)
w_md <- (1 / avar) / sum(1 / avar)
truth_G1 <- list(theta = theta0, r = r1, avar = avar, se_market = sqrt(avar), weights = w_md,
                 se_md = sqrt(1 / sum(1 / avar)), log_step = log(1 + P1$delta))

## =================================================================
## G1h: the same 16 tests, elasticity heterogeneous across markets
## =================================================================
## theta_m = -1.6 + kappa (log income_m - mean_m log income_m), kappa = 3: from -2.22 (Houston, USD 62k)
## to -0.67 (Seattle, USD 104k); the simple mean of theta_m is exactly -1.6. Same design (n, delta,
## p0, q0, sigma) as G1, its own seed. Under heterogeneity the MD estimate with the optimal weight
## converges to theta_star = sum_m w_m theta_m with w_m proportional to 1 / avar_m: the weights are set
## by each test's sample size, price step and noise, not by the market's size (theta_hh weights by
## Tamarind households). The MD J statistic is noncentral chi^2_15 with noncentrality
## lambda = sum_m (theta_m - theta_star)^2 / avar_m.
theta_m <- P1$theta + P1$kappa * (log(P1$income) - mean(log(P1$income)))
set.seed(seed$G1h)
G1h <- gen_G1(P1, theta_m)
theta_star <- sum(w_md * theta_m)
lambda_J <- sum((theta_m - theta_star)^2 / avar)
truth_G1h <- list(theta_m = theta_m, theta_star = theta_star, theta_mean = mean(theta_m),
                  theta_hh = sum(P1$households * theta_m) / sum(P1$households), lambda_J = lambda_J,
                  power_J_01 = pchisq(qchisq(0.99, r1 - 1), r1 - 1, ncp = lambda_J, lower.tail = FALSE))

## =================================================================
## G2: regions, 10 randomised pricing experiments as instruments (linear moments)
## =================================================================
## 200 regions, one week. Experiment j puts n share[j] regions in its upper arm; the arms
## differ by step[j] in the experiment's own unit (USD, or percentage points for the two discount-depth
## tests); the instrument is the centred assignment z_j = step[j] (T_j - share[j]) (mean exactly 0).
## gamma[j] is the change in the region's average log price paid per unit (a pack is 8 lesson credits at
## USD 10.49, USD 83.92):
##    1 booking fee per pack        arms +/- USD 1.50   gamma = 1 / 83.92
##    2 deep-discount promotion     20% of regions get a discount 10 pp deeper; gamma = -0.006 (60% of purchases use it)
##    3 service fee per purchase    +/- USD 0.40        gamma = 1 / 83.92
##    4 small-pack surcharge        +/- USD 2.00        gamma = 0.4 / 83.92 (40% of purchases)
##    5 premium-tutor upcharge      +/- USD 0.75        gamma = 0.3 / 10.49 (per credit, 30% of lessons)
##    6 loyalty rebate per pack     +/- USD 3.00        gamma = -0.5 / 83.92 (half the households hold rebates)
##    7 fee-waiver threshold        +/- USD 10          gamma = 0.001
##    8 two-pack bundle discount    +/- USD 4.00        gamma = -0.004
##    9 weekend-lesson surcharge    +/- USD 1.00        gamma = 0.6 / 83.92 (60% of lessons on weekends)
##   10 household-size discount     +/- 2.5 pp          gamma = -0.003
## The two arms of an experiment differ in price by 1% to 6% (sd of its price effect 0.5% to 2.4%). The
## manager's discretionary promotions v have sd sigma_v = 0.03 and correlation rho = 0.8 with the demand
## shock (managers cut prices in weak weeks), so OLS is biased towards zero.
## Demand: y = log weekly lesson bookings per 1,000 households = log(60) - 1.6 (log price - log 10.49) + eps,
## eps = sigma_i e, e ~ N(0, 1): sigma_i = 0.10, four times larger (0.40) in the deep-discount regions (deal
## seekers book erratically; het = 3 on experiment 2). Because the noisy arm is the 20% arm, the
## discount moment is far noisier than E[z z'] suggests: S_0 is not proportional to E[z z'] and efficient
## GMM is asymptotically more precise than 2SLS (pop_moments, gmm_fns.R, gives D_0, S_0 and the
## asymptotic standard errors). With r = 11 moments and n = 200 the estimated weight matrix is noisy
## and the conventional two-step standard errors are too small (Windmeijer 2005).
## Monte Carlo (03_gmm_mc.R): R_mc = 2000 at n = 200 (seed G2_mc) and R_grid = 1000 at each other n on
## grid_n (seed G2_grid + j).
P2 <- list(n = 200, step = c(3, 10, 0.8, 4, 1.5, 6, 20, 8, 2, 5), share = c(0.5, 0.2, rep(0.5, 8)),
           gamma = c(1 / 83.92, -0.006, 1 / 83.92, 0.4 / 83.92, 0.3 / 10.49, -0.5 / 83.92, 0.001, -0.004, 0.6 / 83.92, -0.003),
           experiments = c("delivery_fee", "deep_coupon", "service_fee", "small_basket_fee", "premium_upcharge",  # first-draft ids, kept for a byte-identical truth07.rds (Ruling 32); labels above
                           "loyalty_credit", "free_delivery_threshold", "bundle_discount", "weekend_surcharge",
                           "household_discount"),
           units = c("USD", "pp", "USD", "USD", "USD", "USD", "USD", "USD", "USD", "pp"),
           sigma_v = 0.03, rho = 0.8, sigma_eps = 0.10, het = 3, het_on = 2, mix_p = 0, mix_scale = 1,
           p_base = 10.49, q_base = 60, theta = theta0, R_mc = 2000, grid_n = c(100, 200, 400, 800), R_grid = 1000)
set.seed(seed$G2)
G2 <- gen_G2(P2)
pp2 <- pop_moments(P2)
truth_G2 <- c(list(theta = theta0, par = c(beta0 = pp2$beta0, theta = theta0)), pp2)

## =================================================================
## G3: weak moments, penny-sized price perturbations
## =================================================================
## 1,000 regions; the only experiments are four list-price perturbations of +/- 0.8, 1.2, 1.6 and
## 2.0 cents per lesson credit (steps of 1.6 to 4.0 cents, each balanced: 500 regions up, 500 down), passing one
## for one into the price (gamma = 1 / 10.49 per USD: effects of +/- 0.08% to 0.19% of the price). Same
## demand, same discretionary promotions (sigma_v = 0.03, rho = 0.8) and a homoskedastic normal demand
## shock (sigma_eps = 0.10; there is no deep-discount promotion in these regions). The concentration parameter
## mu^2 = n sum_j Var(gamma_j z_j) / sigma_v^2 is about 8.7 (population first-stage F about
## 1 + mu^2 / 4 = 3.2): the moments barely identify theta. Moments: g_conc (gmm_fns.R), the K = 4 shock
## moments with the intercept partialled out, k = 1, so the S statistic has 4 degrees of freedom.
## Monte Carlo (04_weak.R): R_mc = 2000 (seed G3_mc).
P3 <- list(n = 1000, step = c(0.016, 0.024, 0.032, 0.040), share = rep(0.5, 4), gamma = rep(1 / 10.49, 4),
           experiments = c("menu_perturbation_0p8c", "menu_perturbation_1p2c", "menu_perturbation_1p6c", "menu_perturbation_2p0c"),  # first-draft ids, kept (Ruling 32)
           units = rep("USD", 4),
           sigma_v = 0.03, rho = 0.8, sigma_eps = 0.10, het = 0, het_on = 1, mix_p = 0, mix_scale = 1,
           p_base = 10.49, q_base = 60, theta = theta0, R_mc = 2000)
set.seed(seed$G3)
G3 <- gen_G3(P3)
pp3 <- pop_moments(P3)
truth_G3 <- list(theta = theta0, beta0 = pp3$beta0, eff_sd = pp3$eff_sd, mu2 = pp3$mu2, F_pop = 1 + pp3$mu2 / length(P3$step),
                 ols_plim = pp3$ols_plim, r2_first_stage = pp3$r2_first_stage, se_2sls = pp3$se_2sls)

## =================================================================
## Save
## =================================================================
saveRDS(G1, "tamarind_G1.rds"); saveRDS(G1h, "tamarind_G1h.rds"); saveRDS(G2, "tamarind_G2.rds"); saveRDS(G3, "tamarind_G3.rds")
tr <- list(G1 = truth_G1, G1h = truth_G1h, G2 = truth_G2, G3 = truth_G3,
           params = list(seed = seed, G1 = P1, G2 = P2, G3 = P3))
saveRDS(tr, "truth07.rds")

cat(sprintf("G1: %d markets, %d households; true MD standard error %.4f; largest weight %.3f (%s)\n",
            r1, nrow(G1), tr$G1$se_md, max(w_md), P1$markets[which.max(w_md)]))
cat(sprintf("G1h: theta_m from %.3f to %.3f; theta_star %.4f, mean %.4f, household-weighted %.4f; lambda %.1f, power %.4f\n",
            min(theta_m), max(theta_m), theta_star, mean(theta_m), tr$G1h$theta_hh, lambda_J, tr$G1h$power_J_01))
cat(sprintf("G2: n = %d, r = %d; mu^2 %.1f, first-stage R^2 %.3f, OLS plim %.4f; asymptotic SE 2SLS %.4f, efficient %.4f, identity %.4f\n",
            nrow(G2), length(P2$step) + 1, pp2$mu2, pp2$r2_first_stage, pp2$ols_plim, pp2$se_2sls, pp2$se_eff, pp2$se_identity))
cat(sprintf("G3: n = %d, r = %d; mu^2 %.2f, population F %.2f, OLS plim %.4f\n",
            nrow(G3), length(P3$step), tr$G3$mu2, tr$G3$F_pop, tr$G3$ols_plim))
cat("Saved tamarind_G1.rds, tamarind_G1h.rds, tamarind_G2.rds, tamarind_G3.rds, truth07.rds\n")
