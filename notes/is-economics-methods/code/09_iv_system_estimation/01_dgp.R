# ------------------------------------------------------------------
# 01_dgp.R -- the Bramble simulation world for chapter 9
#             (instrumental variables and system estimation)
#
# Bramble is a (fictional) two-sided marketplace for local food producers: households (buyers)
# order from home bakers, small farms and independent kitchens (sellers) in 1,500 local markets.
# Money is in US dollars. One world, six blocks (all new in v3; nothing is inherited):
#   S1   one month of the two-sided equilibrium: a buyer equation (active buyers respond to the
#        number of sellers) and a seller equation (active sellers respond to the number of buyers),
#        solved jointly. Buyer-side shifters that enter only the buyer equation: the ad budget shock
#        and the first-order coupon; seller-side shifter that enters only the seller equation: the
#        monthly onboarding fee. The buyer equation is just identified, the seller equation is
#        over-identified by one (ruling 14). The city plan pairs deeper coupons with fee waivers
#        (corr(coupon, fee) = -0.7); this link is what lets 3SLS improve on 2SLS for the cross-side
#        coefficient n_sellers (with independent shifters the gain is exactly zero, see pop_avar).
#   S1m  the same markets, but the true seller equation has a coupon x buyers term (coupon-driven
#        buyers order less, so each 1,000 extra buyers draw fewer sellers where coupons are deep);
#        the analyst's linear seller equation omits it.
#   S2   two outcome equations for the same kind of markets (average basket in USD, orders per buyer),
#        correlated errors: SUR with identical regressors (Theorem 7.6) and with different ones.
#   S3   the weak-instrument variant: a pilot that nudged the onboarding fee by about USD 0.35 in each
#        market, independently of the coupon: the fee barely moves seller counts (concentration
#        parameter mu2 = 1). Same structural equations and errors as S1.
#   S4   an over-identified buyer equation: three seller-side programs (fee, payout delay, referral
#        bonus) randomized in one quarter; in S4 the referral bonus also brings buyers directly
#        (it violates the exclusion restriction), in S4v it does not. Errors are heteroskedastic in
#        the referral bonus.
# Every number in chapter 9 comes from these data and from truth09.rds.
#
# Outputs (in this directory):
#   bramble_S1.rds, bramble_S1m.rds, bramble_S3.rds  data.frame: market buyers n_sellers ad_budget coupon fee log_pop
#                                                    eps_buyer_true eps_seller_true
#   bramble_S2.rds    data.frame: market free_ship_min premium_share loyalty basket_usd orders_pb eps_basket_true eps_orders_true
#   bramble_S4.rds, bramble_S4v.rds  data.frame: market buyers n_sellers ad_budget coupon fee payout_days referral_bonus
#                                    log_pop eps_buyer_true eps_seller_true h_true
#   truth09.rds       list(S1, S1m, S2, S3, S4, params): truth, population quantities and parameters
# Columns ending in _true are simulation truth, not observable data. Units: buyers in thousands of
# active buyers, n_sellers in hundreds of active sellers, ad_budget in USD thousands per month, coupon,
# fee, free_ship_min, basket_usd and referral_bonus in USD, payout_days in days, log_pop standardized.
#
# Notation (notation_v3.md 1.7): M = 2 equations; y_i = (buyers, n_sellers)'; z_i = (1, ad_budget,
# coupon, fee, log_pop)'; y_i' Gamma + z_i' Delta + eps_i' = 0, Pi = -Delta Gamma^-1, Pi_(m) the
# identification submatrix; Omega = Var(eps_i), Omega_v = (Gamma^-1)' Omega Gamma^-1 the reduced-form
# error variance; mu2 the concentration parameter of the weak fee pilot.
#
# Seeds (every block its own set.seed, nothing inherited): S1 9101, S1m 9102, S2 9103, S3 9104,
# S4 9105, S4v 9106; large-sample checks S1 9111, S1m 9112; the Monte Carlo seeds used by 02, 04 and
# 05 are stored here too (S1 9121, S1m 9122, S3 9124 and 94000 + j for grid point j, S4 9125 for the
# valid design and 9126 for the invalid one). The large-sample contagion step 02b_contagion_r0.R has its own
# seeds, set in that script (9113 for corr(coupon, fee) = 0, 9114 for -0.7).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("iv_fns.R")

# Seeds were fixed before any result was seen and were not tuned afterwards.
seed <- list(S1 = 9101, S1m = 9102, S2 = 9103, S3 = 9104, S4 = 9105, S4v = 9106,
             S1_large = 9111, S1m_large = 9112, S1_mc = 9121, S1m_mc = 9122, S3_mc = 9124, S3_grid = 94000,
             S4_mc = 9125, S4_power = 9126)
N_large <- 1e6

## =================================================================
## S1: the two-sided equilibrium
## =================================================================
## 1,500 markets, one month. Buyer-side shifters: ad_budget ~ N(40, 10^2) USD thousands (HQ
## reallocates the ad budget at random across markets), coupon = 6 + 1.5 u USD (first-order
## coupon); seller-side shifter: fee = 30 + 7 (-0.7 u + sqrt(0.51) w) USD per month (onboarding
## fee); market size log_pop ~ N(0, 1). The city plan pairs deeper coupons with fee waivers, so
## corr(coupon, fee) = -0.7; ad_budget is independent of both.
##   buyers    = 14.4 + 1.2 n_sellers + 0.10 ad_budget + 0.40 coupon + 2.0 log_pop + eps_B
##   n_sellers = 9.0 + 0.25 buyers - 0.10 fee + 0.6 log_pop + eps_S
## 100 more sellers bring 1,200 more buyers; 1,000 more buyers bring 25 more sellers; the feedback
## multiplier is 1 / (1 - 1.2 x 0.25) = 1 / 0.7. USD 1,000 of ads brings 100 buyers (USD 10 per buyer),
## USD 1 of coupon 400 buyers, USD 1 lower fee 10 sellers. Average market: 40,000 buyers, 1,600 sellers.
## (eps_B, eps_S): local demand shocks hit both sides, sd 2.5 (thousand buyers) and 0.8 (hundred
## sellers), correlation 0.8, homoskedastic (SIV.5 holds, so 3SLS is the efficient system estimator).
P1 <- list(n = 1500, ad0 = 40, ad_sd = 10, cp0 = 6, cp_sd = 1.5, fee0 = 30, fee_sd = 7, r = -0.7,
           EB = 40, ES = 16, bBS = 1.2, gad = 0.10, gcp = 0.40, lB = 2.0, sB = 2.5,
           bSB = 0.25, dfee = -0.10, lS = 0.6, sS = 0.8, rho = 0.8)
P1$aB <- P1$EB - P1$bBS * P1$ES - P1$gad * P1$ad0 - P1$gcp * P1$cp0           # 14.4
P1$aS <- P1$ES - P1$bSB * P1$EB - P1$dfee * P1$fee0                            # 9.0
set.seed(seed$S1)
S1 <- gen_S1(P1)

zn <- c("const", "ad_budget", "coupon", "fee", "log_pop"); yn <- c("buyers", "n_sellers")
struct <- function(par) {
  Gamma <- matrix(c(-1, par$bBS, par$bSB, -1), 2, dimnames = list(yn, c("buyer", "seller")))
  Delta <- cbind(buyer = c(par$aB, par$gad, par$gcp, 0, par$lB), seller = c(par$aS, 0, 0, par$dfee, par$lS))
  rownames(Delta) <- zn
  list(Gamma = Gamma, Delta = Delta, Pi = reduced_form(Gamma, Delta))
}
st <- struct(P1)
Omega <- matrix(c(P1$sB^2, P1$rho * P1$sB * P1$sS, P1$rho * P1$sB * P1$sS, P1$sS^2), 2, dimnames = list(yn, yn))
Gi <- solve(st$Gamma)
Omega_v <- t(Gi) %*% Omega %*% Gi                                             # Var(v_i), v_i' = -eps_i' Gamma^-1
cov_eps_v <- -Omega %*% Gi                                                    # E[eps_i v_i']
corr_eB_v <- cov_eps_v[1, 2] / sqrt(Omega[1, 1] * Omega_v[2, 2])              # endogeneity of n_sellers in the buyer equation

ident_of <- function(Pi_sub, n_endog) {
  rc <- rank_condition(Pi_sub)
  list(n_excluded = nrow(Pi_sub), n_endog = n_endog, order = order_condition(nrow(Pi_sub), n_endog),
       rank = rc$rank, rank_ok = rc$ok, overid = nrow(Pi_sub) - n_endog)
}
Pi_sub <- list(buyer = st$Pi["fee", "n_sellers", drop = FALSE],
               seller = st$Pi[c("ad_budget", "coupon"), "buyers", drop = FALSE])
Bmat <- rbind(st$Gamma, st$Delta)                                             # (M + r) x M structural parameters
rank_general <- list(buyer = rank_condition(Bmat["fee", , drop = FALSE])$rank,    # Theorem 9.2: rank R_m B = M - 1
                     seller = rank_condition(Bmat[c("ad_budget", "coupon"), , drop = FALSE])$rank)
st0 <- struct(modifyList(P1, list(dfee = 0)))                                 # counterfactual: the fee does not move sellers
nofee <- rank_condition(st0$Pi["fee", "n_sellers", drop = FALSE])

## Population asymptotic variances of 2SLS and 3SLS (homoskedastic) and the SE ratio for n_sellers.
## Q = E[z z'] from the instrument means and covariance; x_hat_im = Pi_fs' z_i with Pi_fs built from Pi.
pop_ratio <- function(par) {
  s <- struct(par)
  mu <- c(1, par$ad0, par$cp0, par$fee0, 0)
  Cz <- diag(c(0, par$ad_sd^2, par$cp_sd^2, par$fee_sd^2, 1))
  Cz[3, 4] <- Cz[4, 3] <- par$r * par$cp_sd * par$fee_sd
  Q <- Cz + tcrossprod(mu)
  e <- diag(5)
  PiB <- cbind(e[, 1], s$Pi[, "n_sellers"], e[, 2], e[, 3], e[, 5])
  PiS <- cbind(e[, 1], s$Pi[, "buyers"], e[, 4], e[, 5])
  pa <- pop_avar_sys(list(PiB, PiS), Q, Omega)
  list(se_ratio = sqrt(pa$V_3sls[2, 2] / pa$V_2sls[[1]][2, 2]), se_2sls_n = sqrt(pa$V_2sls[[1]][2, 2] / par$n),
       se_3sls_n = sqrt(pa$V_3sls[2, 2] / par$n), Q = Q)
}
pr <- pop_ratio(P1)
r_grid <- round(seq(-0.9, 0.9, by = 0.1), 10)
pop_avar <- list(se_ratio = pr$se_ratio, se_2sls_n = pr$se_2sls_n, se_3sls_n = pr$se_3sls_n,
                 se_ratio_independent = pop_ratio(modifyList(P1, list(r = 0)))$se_ratio,
                 by_r = data.frame(r = r_grid, se_ratio = vapply(r_grid, function(r) pop_ratio(modifyList(P1, list(r = r)))$se_ratio, numeric(1))))

## Large-sample check (one million markets, seed S1_large): probability limits of OLS, 2SLS and 3SLS.
plim_system <- function(d) {
  sy <- s1_system(d)
  tB <- tsls(sy$y$buyer, sy$X$buyer, sy$Z, vcov = "homo"); tS <- tsls(sy$y$seller, sy$X$seller, sy$Z, vcov = "homo")
  th <- three_sls(sy$y, sy$X, sy$Z)
  J <- hansen_j(sy$y$seller, sy$X$seller, sy$Z)
  list(N = nrow(d), ols_buyer = ols_fit(sy$y$buyer, sy$X$buyer)$b[["n_sellers"]],
       tsls_buyer = tB$b[["n_sellers"]], three_sls_buyer = th$b[["buyer:n_sellers"]],
       tsls_seller = tS$b[["buyers"]], three_sls_seller = th$b[["seller:buyers"]],
       se_ratio = sqrt(th$V["buyer:n_sellers", "buyer:n_sellers"] / tB$V["n_sellers", "n_sellers"]),
       rho_hat = th$Sigma[1, 2] / sqrt(th$Sigma[1, 1] * th$Sigma[2, 2]), J_seller_per_n = J$J / nrow(d))
}
set.seed(seed$S1_large)
plim_S1 <- plim_system(gen_S1(P1, n = N_large))

truth_S1 <- list(buyer = c(const = P1$aB, n_sellers = P1$bBS, ad_budget = P1$gad, coupon = P1$gcp, log_pop = P1$lB),
                 seller = c(const = P1$aS, buyers = P1$bSB, fee = P1$dfee, log_pop = P1$lS),
                 means = list(buyers = P1$EB, n_sellers = P1$ES),
                 Gamma = st$Gamma, Delta = st$Delta, Pi = st$Pi, Pi_sub = Pi_sub,
                 ident = list(buyer = ident_of(Pi_sub$buyer, 1), seller = ident_of(Pi_sub$seller, 1)),
                 rank_general = rank_general,
                 ident_no_fee = list(order = order_condition(1, 1), rank = nofee$rank, ok = nofee$ok),
                 Omega = Omega, Omega_v = Omega_v, corr_eB_v = corr_eB_v, multiplier = 1 / (1 - P1$bBS * P1$bSB),
                 pop_avar = pop_avar, plim = plim_S1)

## =================================================================
## S1m: the seller equation omits a coupon x buyers term
## =================================================================
## Same parameters and design as S1; the true seller equation is
##   n_sellers = 9.0 + (0.25 - 0.0045 (coupon - 6)) buyers - 0.10 fee + 0.6 log_pop + eps_S:
## coupon-driven buyers order less, so where the coupon is USD 1 deeper than the plan's USD 6, each 1,000
## extra buyers draw 0.45 fewer sellers (a 1.8% smaller cross-side response; at a coupon one sd above the
## mean, 2.7% smaller). The analyst fits the linear S1 seller equation. Because the omitted term contains
## coupon x E[buyers], the coupon is no longer a valid excluded instrument for that equation, and the
## seller equation's over-identifying restriction fails (Wooldridge 2010 section 9.3.1: one misspecified
## equation makes 3SLS inconsistent for all equations; 2SLS of the correct buyer equation stays consistent).
## WHICH buyer coefficient is contaminated depends on the same coupon-fee link (corr = -0.7) that produces the
## 3SLS efficiency gain in S1; 2SLS of the buyer equation is consistent for every coefficient in both designs.
## The violated moment is the seller equation's coupon moment. Through the error correlation (rho = 0.8) 3SLS
## carries it into the buyer equation, into the coefficient(s) that move with it: when the coupon is correlated
## with the fee, the instrument that identifies n_sellers, the n_sellers coefficient is contaminated (large-sample
## limits about 1.20 under 2SLS and 0.80 under 3SLS, truth 1.2); with independent shifters (corr(coupon, fee) = 0)
## the n_sellers limits coincide (about 1.20 both: no contagion on n_sellers, the same reason 3SLS gains nothing on
## n_sellers there) and the contamination sits on the buyer equation's own coupon coefficient (2SLS about 0.40,
## 3SLS about 0.73, truth 0.40). The coupon coefficient is contaminated under 3SLS in both designs. The
## corr(coupon, fee) = 0 limits are computed by 02b_contagion_r0.R (n = 1e6, its own seeds 9113 and 9114, stored
## in res_contagion_r0.rds, keys dgp_S1m_r0_* and dgp_S1m_r_minus0p7_* in numbers.json); the plim_S1m above
## is the corr = -0.7 limit of n_sellers only.
P1m <- c(P1, list(kappa = -0.0045))
set.seed(seed$S1m)
S1m <- gen_S1m(P1m)
set.seed(seed$S1m_large)
plim_S1m <- plim_system(gen_S1m(P1m, n = N_large))
truth_S1m <- list(buyer = truth_S1$buyer, seller = truth_S1$seller, kappa = P1m$kappa, plim = plim_S1m)

## =================================================================
## S2: SUR, two outcome equations for the same markets
## =================================================================
## 1,500 markets. free_ship_min ~ N(35, 5^2) USD (minimum basket for free delivery); premium_share =
## 0.30 + 0.008 (free_ship_min - 35) + 0.06 N(0, 1) (share of sellers with a premium badge); loyalty = 1 when
## 0.08 (free_ship_min - 35) + N(0, 1) > 0.25 (a loyalty program, on in about 41% of markets).
##   basket_usd = 18 + 0.50 free_ship_min + 25 premium_share + eps_1   (average order value, mean about USD 43)
##   orders_pb  = 3.6 - 0.04 free_ship_min + 0.5 loyalty + eps_2       (orders per active buyer per month, mean about 2.4)
## (eps_1, eps_2) normal, sd USD 4 and 0.45 orders, correlation 0.6 (local income shocks raise both).
## "Different regressors": each equation with its own variables (the truth); "same regressors": both
## equations on (1, free_ship_min, premium_share, loyalty), the true zero coefficients included.
P2 <- list(n = 1500, ship0 = 35, ship_sd = 5, prem0 = 0.30, prem_ship = 0.008, prem_sd = 0.06, loy_ship = 0.08,
           loy_cut = 0.25, b1 = c(18, 0.50, 25), b2 = c(3.6, -0.04, 0.5), s1 = 4, s2 = 0.45, rho = 0.6)
set.seed(seed$S2)
S2 <- gen_S2(P2)
p_loy <- 1 - pnorm(P2$loy_cut / sqrt(1 + P2$loy_ship^2 * P2$ship_sd^2))
E_basket <- P2$b1[1] + P2$b1[2] * P2$ship0 + P2$b1[3] * P2$prem0
E_orders <- P2$b2[1] + P2$b2[2] * P2$ship0 + P2$b2[3] * p_loy
truth_S2 <- list(beta = c(basket_const = P2$b1[1], basket_free_ship_min = P2$b1[2], basket_premium_share = P2$b1[3],
                          orders_const = P2$b2[1], orders_free_ship_min = P2$b2[2], orders_loyalty = P2$b2[3]),
                 beta_same_X = c(basket_const = P2$b1[1], basket_free_ship_min = P2$b1[2], basket_premium_share = P2$b1[3],
                                 basket_loyalty = 0, orders_const = P2$b2[1], orders_free_ship_min = P2$b2[2],
                                 orders_premium_share = 0, orders_loyalty = P2$b2[3]),
                 Omega = matrix(c(P2$s1^2, P2$rho * P2$s1 * P2$s2, P2$rho * P2$s1 * P2$s2, P2$s2^2), 2),
                 p_loyalty = p_loy, E_basket = E_basket, E_orders = E_orders,
                 # first-order effect of USD 1 on the free-shipping minimum on revenue per buyer (basket x orders)
                 revenue_effect = E_orders * P2$b1[2] + E_basket * P2$b2[2])

## =================================================================
## S3: the weak fee pilot
## =================================================================
## Same structural equations and errors as S1; coupon and fee independent (r = 0); the pilot nudges the
## fee by fee_sd ~ USD 0.35 around USD 30. The reduced-form fee coefficient for n_sellers is
## pi_f = -0.10 / 0.7 = -0.143 (hundred sellers per USD), the reduced-form error sd sigma_v = sqrt(Omega_v[2, 2]).
## Concentration parameter mu2 = n pi_f^2 Var(fee) / sigma_v^2 (E[F_fs] about 1 + mu2): fee_sd is set so that
## mu2 = 1 at n = 1,500. Monte Carlo (04_weak.R): R_mc = 2000 replications at mu2 = 1, R_grid = 1000 at each
## other mu2 on grid_mu2 (from 0.25, a nearly irrelevant instrument, to 128, a strong one).
pi_f <- st$Pi["fee", "n_sellers"]; sigma_v <- sqrt(Omega_v[2, 2])
fee_sd_for <- function(mu2, n) sqrt(mu2 * sigma_v^2 / (n * pi_f^2))
P3 <- modifyList(P1, list(r = 0, mu2 = 1, R_mc = 2000, R_grid = 1000,
                          grid_mu2 = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)))
P3$fee_sd <- fee_sd_for(P3$mu2, P3$n)
set.seed(seed$S3)
S3 <- gen_S3(P3)
grid3 <- data.frame(mu2 = P3$grid_mu2, fee_sd = fee_sd_for(P3$grid_mu2, P3$n))
# OLS probability limit in the S3 design (fee, coupon, ad_budget independent): beta + Cov(eps_B, v_S) / Var(n_sellers | W)
grid3$plim_ols <- P1$bBS + cov_eps_v[1, 2] / (pi_f^2 * grid3$fee_sd^2 + sigma_v^2)
truth_S3 <- list(beta = c(n_sellers = P1$bBS), mu2 = P3$mu2, fee_sd = P3$fee_sd, pi_f = pi_f, sigma_v = sigma_v,
                 corr_eB_v = corr_eB_v, plim_ols = P1$bBS + cov_eps_v[1, 2] / (pi_f^2 * P3$fee_sd^2 + sigma_v^2),
                 grid = grid3)

## =================================================================
## S4: an over-identified buyer equation, one invalid instrument
## =================================================================
## 1,500 markets, one quarter. Buyer-side: ad_budget ~ N(40, 10^2), coupon ~ 6 + 1.5 N(0, 1) (independent).
## Three seller-side programs randomized independently: fee ~ N(30, 7^2) USD, payout_days ~ U(1, 7)
## (days until a seller is paid), referral_bonus ~ U(20, 80) USD (paid per referred seller).
##   n_sellers = 9.2 + 0.25 buyers - 0.10 fee - 0.30 payout_days + 0.02 referral_bonus + 0.6 log_pop + eps_S
##   buyers    = 14.4 + 1.2 n_sellers + 0.10 ad_budget + 0.40 coupon + 2.0 log_pop + theta (referral_bonus - 50) + eps_B
## S4: theta = 0.05 (a referred seller brings regular customers along: USD 10 more bonus, 500 more buyers
## directly), so the referral bonus violates the exclusion restriction in the buyer equation; S4v: theta = 0.
## Errors as in S1 but scaled by h = exp(0.5 zr) / sqrt(E[exp(zr)]), zr the standardized referral bonus:
## referral waves make both sides more volatile (Var(eps | z) depends on an instrument, so the homoskedastic
## Sargan statistic is not valid; the robust Hansen J is). Monte Carlo (05_overid.R): R_mc = 2000 valid samples
## (size of the tests) and R_power = 1000 invalid ones (power).
P4 <- list(n = 1500, ad0 = 40, ad_sd = 10, cp0 = 6, cp_sd = 1.5, fee0 = 30, fee_sd = 7, pay_lo = 1, pay_hi = 7,
           ref_lo = 20, ref_hi = 80, ref_mid = 50, EB = 40, ES = 16, bBS = 1.2, gad = 0.10, gcp = 0.40, lB = 2.0, sB = 2.5,
           bSB = 0.25, dfee = -0.10, dpay = -0.30, dref = 0.02, lS = 0.6, sS = 0.8, rho = 0.8, eta = 0.5,
           theta = 0.05, R_mc = 2000, R_power = 1000)
P4$aB <- P4$EB - P4$bBS * P4$ES - P4$gad * P4$ad0 - P4$gcp * P4$cp0                                    # 14.4
P4$aS <- P4$ES - P4$bSB * P4$EB - P4$dfee * P4$fee0 - P4$dpay * (P4$pay_lo + P4$pay_hi) / 2 - P4$dref * P4$ref_mid   # 9.2
P4v <- modifyList(P4, list(theta = 0))
set.seed(seed$S4)
S4 <- gen_S4(P4)
set.seed(seed$S4v)
S4v <- gen_S4v(P4v)
truth_S4 <- list(beta = c(n_sellers = P4$bBS), theta = P4$theta, overid_df = 2,
                 buyer = c(const = P4$aB, n_sellers = P4$bBS, ad_budget = P4$gad, coupon = P4$gcp, log_pop = P4$lB))

## =================================================================
## Save
## =================================================================
saveRDS(S1, "bramble_S1.rds"); saveRDS(S1m, "bramble_S1m.rds"); saveRDS(S2, "bramble_S2.rds")
saveRDS(S3, "bramble_S3.rds"); saveRDS(S4, "bramble_S4.rds"); saveRDS(S4v, "bramble_S4v.rds")
tr <- list(S1 = truth_S1, S1m = truth_S1m, S2 = truth_S2, S3 = truth_S3, S4 = truth_S4,
           params = list(seed = seed, N_large = N_large, S1 = P1, S1m = P1m, S2 = P2, S3 = P3, S4 = P4, S4v = P4v))
saveRDS(tr, "truth09.rds")

cat(sprintf("S1: %d markets; mean buyers %.2f thousand, sellers %.2f hundred; corr(coupon, fee) %.3f; corr(eps_B, v_S) %.3f\n",
            nrow(S1), mean(S1$buyers), mean(S1$n_sellers), cor(S1$coupon, S1$fee), corr_eB_v))
cat(sprintf("S1 identification: buyer Pi_(m) = %.4f (rank %d), seller Pi_(m) = (%.4f, %.4f) (rank %d); no-fee rank %d\n",
            Pi_sub$buyer[1, 1], truth_S1$ident$buyer$rank, Pi_sub$seller[1, 1], Pi_sub$seller[2, 1], truth_S1$ident$seller$rank, nofee$rank))
cat(sprintf("S1 population SE ratio 3SLS / 2SLS for n_sellers: %.4f (independent shifters: %.6f)\n",
            pop_avar$se_ratio, pop_avar$se_ratio_independent))
cat(sprintf("S1 plim (N = %d): OLS %.4f, 2SLS %.4f, 3SLS %.4f; seller 2SLS %.4f, 3SLS %.4f; SE ratio %.4f\n",
            plim_S1$N, plim_S1$ols_buyer, plim_S1$tsls_buyer, plim_S1$three_sls_buyer, plim_S1$tsls_seller,
            plim_S1$three_sls_seller, plim_S1$se_ratio))
cat(sprintf("S1m plim: 2SLS %.4f, 3SLS %.4f; seller 2SLS %.4f; rho_hat %.3f; J_seller / n %.4f\n",
            plim_S1m$tsls_buyer, plim_S1m$three_sls_buyer, plim_S1m$tsls_seller, plim_S1m$rho_hat, plim_S1m$J_seller_per_n))
cat(sprintf("S3: fee_sd %.4f USD, mu2 %.2f, pi_f %.4f, sigma_v %.4f, OLS plim %.4f\n",
            P3$fee_sd, P3$mu2, pi_f, sigma_v, truth_S3$plim_ols))
cat("Saved bramble_S1.rds ... bramble_S4v.rds, truth09.rds\n")
