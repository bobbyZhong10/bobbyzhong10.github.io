# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Bramble data and truth09.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("iv_fns.R")
tr <- readRDS("truth09.rds")
S1 <- readRDS("bramble_S1.rds"); S1m <- readRDS("bramble_S1m.rds"); S2 <- readRDS("bramble_S2.rds")
S3 <- readRDS("bramble_S3.rds"); S4 <- readRDS("bramble_S4.rds"); S4v <- readRDS("bramble_S4v.rds")
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))
P <- tr$params

## ---- truth09 fields named on the chapter card -------------------------------------------
stopifnot(fin(tr$S1$buyer, 5L), fin(tr$S1$seller, 4L), fin(tr$S2$beta, 6L), fin(tr$S3$beta), fin(tr$S4$beta),
          identical(names(tr$S1$buyer), c("const", "n_sellers", "ad_budget", "coupon", "log_pop")),
          identical(names(tr$S1$seller), c("const", "buyers", "fee", "log_pop")),
          !is.null(names(tr$S2$beta)), identical(names(tr$S3$beta), "n_sellers"), identical(names(tr$S4$beta), "n_sellers"))
b <- tr$S1$buyer[["n_sellers"]]
stopifnot(tr$S3$beta[["n_sellers"]] == b, tr$S4$beta[["n_sellers"]] == b, identical(tr$S1m$buyer, tr$S1$buyer))
# every block has its own seed, all different integers
sd <- unlist(P$seed)
stopifnot(length(sd) >= 12, !anyDuplicated(sd), all(sd == round(sd)))

## ---- the generators reproduce the saved data from the stored seeds -----------------------
for (blk in c("S1", "S1m", "S2", "S3", "S4", "S4v")) {
  set.seed(P$seed[[blk]])
  regen <- get(paste0("gen_", blk))(P[[blk]])
  stopifnot(isTRUE(all.equal(regen, get(blk), tolerance = 0, check.attributes = TRUE)))
}

## ---- S1: two-sided equilibrium in 1,500 markets -------------------------------------------
P1 <- P$S1
cols1 <- c("market", "buyers", "n_sellers", "ad_budget", "coupon", "fee", "log_pop", "eps_buyer_true", "eps_seller_true")
for (dd in list(S1, S1m, S3)) stopifnot(is.data.frame(dd), all(cols1 %in% names(dd)),
                                        all(dd$buyers > 0), all(dd$n_sellers > 0), all(dd$ad_budget > 0),
                                        all(dd$coupon > 0), all(dd$fee > 0))
stopifnot(nrow(S1) == P1$n, nrow(S1m) == P$S1m$n, nrow(S3) == P$S3$n)
# both structural equations hold exactly with the true errors
sb <- tr$S1$buyer; ss <- tr$S1$seller
eqB <- function(d) d$buyers - (sb[["const"]] + sb[["n_sellers"]] * d$n_sellers + sb[["ad_budget"]] * d$ad_budget +
                                 sb[["coupon"]] * d$coupon + sb[["log_pop"]] * d$log_pop)
eqS <- function(d, kappa = 0) d$n_sellers - (ss[["const"]] + (ss[["buyers"]] + kappa * (d$coupon - P1$cp0)) * d$buyers +
                                               ss[["fee"]] * d$fee + ss[["log_pop"]] * d$log_pop)
stopifnot(max(abs(eqB(S1) - S1$eps_buyer_true)) < 1e-10, max(abs(eqS(S1) - S1$eps_seller_true)) < 1e-10,
          max(abs(eqB(S1m) - S1m$eps_buyer_true)) < 1e-10,
          max(abs(eqS(S1m, P$S1m$kappa) - S1m$eps_seller_true)) < 1e-10,
          max(abs(eqB(S3) - S3$eps_buyer_true)) < 1e-10, max(abs(eqS(S3) - S3$eps_seller_true)) < 1e-10)
# market scale: about 40,000 active buyers and 1,600 active sellers in the average market
stopifnot(abs(mean(S1$buyers) - 40) < 1, abs(mean(S1$n_sellers) - 16) < 0.5,
          abs(tr$S1$means$buyers - 40) < 1e-12, abs(tr$S1$means$n_sellers - 16) < 1e-12)
# the city plan links coupons and fee waivers; the two side-specific errors are correlated
stopifnot(abs(cor(S1$coupon, S1$fee) - P1$r) < 0.05, abs(cor(S1$eps_buyer_true, S1$eps_seller_true) - P1$rho) < 0.05,
          abs(cor(S1$ad_budget, S1$fee)) < 0.1, abs(cor(S3$coupon, S3$fee)) < 0.1, P$S3$r == 0)
# cross-side feedback is stable and the reduced form reproduces the equilibrium
stopifnot(sb[["n_sellers"]] * ss[["buyers"]] < 1)
Zs <- cbind(1, S1$ad_budget, S1$coupon, S1$fee, S1$log_pop)
Vs <- -cbind(S1$eps_buyer_true, S1$eps_seller_true) %*% solve(tr$S1$Gamma)
stopifnot(max(abs(Zs %*% tr$S1$Pi + Vs - cbind(S1$buyers, S1$n_sellers))) < 1e-10,
          max(abs(tr$S1$Pi - reduced_form(tr$S1$Gamma, tr$S1$Delta))) < 1e-14)
# identification of each equation (exclusion restrictions): order (Theorem 9.1) and rank (eq. 9.10)
idb <- tr$S1$ident$buyer; ids <- tr$S1$ident$seller
stopifnot(idb$n_excluded == 1, idb$n_endog == 1, isTRUE(idb$order), idb$rank == 1, isTRUE(idb$rank_ok), idb$overid == 0,
          ids$n_excluded == 2, ids$n_endog == 1, isTRUE(ids$order), ids$rank == 1, isTRUE(ids$rank_ok), ids$overid == 1,
          all(dim(tr$S1$Pi_sub$buyer) == c(1, 1)), all(dim(tr$S1$Pi_sub$seller) == c(2, 1)),
          tr$S1$rank_general$buyer == 1, tr$S1$rank_general$seller == 1)                  # Theorem 9.2: rank R_m B = M - 1
# without a fee shifter the buyer equation fails the rank condition although the order condition still holds
stopifnot(isTRUE(tr$S1$ident_no_fee$order), tr$S1$ident_no_fee$rank == 0, !tr$S1$ident_no_fee$ok)
# population asymptotic standard errors: 3SLS gains on n_sellers only because coupon and fee are correlated
pa <- tr$S1$pop_avar
stopifnot(pa$se_ratio < 0.95, abs(pa$se_ratio_independent - 1) < 1e-8, pa$se_3sls_n < pa$se_2sls_n,
          is.data.frame(pa$by_r), all(diff(pa$by_r$se_ratio[pa$by_r$r <= 0]) >= -1e-12))
# large-sample check (seed S1_large): 2SLS consistent, OLS not; 3SLS and 2SLS agree when both equations are right
pl <- tr$S1$plim
stopifnot(pl$N >= 1e6, abs(pl$tsls_buyer - b) < 0.01, abs(pl$three_sls_buyer - b) < 0.01, pl$ols_buyer - b > 0.5,
          abs(pl$tsls_seller - ss[["buyers"]]) < 0.005, abs(pl$three_sls_seller - pl$tsls_seller) < 1e-8,
          abs(pl$se_ratio - pa$se_ratio) < 0.01)

## ---- S1m: the seller equation has a coupon x buyers term the analyst omits -----------------
P1m <- P$S1m
stopifnot(P1m$kappa < 0, identical(P1m[names(P1)], P1), setequal(names(P1m), c(names(P1), "kappa")))   # only kappa added
slope_m <- ss[["buyers"]] + P1m$kappa * (S1m$coupon - P1$cp0)
stopifnot(all(sb[["n_sellers"]] * slope_m < 1), max(abs(slope_m - ss[["buyers"]])) < 0.1)    # a mild nonlinearity
plm <- tr$S1m$plim
stopifnot(plm$N >= 1e6, abs(plm$tsls_buyer - b) < 0.01, abs(plm$three_sls_buyer - b) > 0.2, plm$J_seller_per_n > 0.01)

## ---- S2: two outcome equations, correlated errors, different regressors --------------------------
P2 <- P$S2
stopifnot(is.data.frame(S2), nrow(S2) == P2$n,
          all(c("free_ship_min", "premium_share", "loyalty", "basket_usd", "orders_pb", "eps_basket_true", "eps_orders_true") %in% names(S2)),
          all(S2$premium_share > 0), all(S2$premium_share < 1), setequal(unique(S2$loyalty), 0:1),
          all(S2$basket_usd > 0), all(S2$orders_pb > 0))
bt <- tr$S2$beta
stopifnot(max(abs(S2$basket_usd - (bt[["basket_const"]] + bt[["basket_free_ship_min"]] * S2$free_ship_min +
                                    bt[["basket_premium_share"]] * S2$premium_share) - S2$eps_basket_true)) < 1e-10,
          max(abs(S2$orders_pb - (bt[["orders_const"]] + bt[["orders_free_ship_min"]] * S2$free_ship_min +
                                   bt[["orders_loyalty"]] * S2$loyalty) - S2$eps_orders_true)) < 1e-10)
bs <- tr$S2$beta_same_X
stopifnot(length(bs) == 8, bs[["basket_loyalty"]] == 0, bs[["orders_premium_share"]] == 0,
          abs(cor(S2$eps_basket_true, S2$eps_orders_true) - P2$rho) < 0.05,
          cor(S2$free_ship_min, S2$premium_share) > 0.2, cor(S2$free_ship_min, S2$loyalty) > 0.1)

## ---- S3: the weak fee pilot ------------------------------------------------------------------
P3 <- P$S3
stopifnot(abs(tr$S3$mu2 - P3$mu2) < 1e-12, P3$mu2 <= 2,
          abs(tr$S3$mu2 - P3$n * tr$S3$pi_f^2 * P3$fee_sd^2 / tr$S3$sigma_v^2) < 1e-12,
          abs(tr$S3$pi_f - tr$S1$Pi["fee", "n_sellers"]) < 1e-14, tr$S3$corr_eB_v > 0.9,
          abs(tr$S3$corr_eB_v - tr$S1$corr_eB_v) < 1e-14,
          P3$R_mc >= 1000, all(P3$R_grid >= 1000), P3$mu2 %in% P3$grid_mu2,
          abs(sd(S3$fee) / P3$fee_sd - 1) < 0.1, P3$fee_sd < 1)
# the grid's fee spreads give the stated concentration parameters
stopifnot(max(abs(P3$n * tr$S3$pi_f^2 * tr$S3$grid$fee_sd^2 / tr$S3$sigma_v^2 - tr$S3$grid$mu2)) < 1e-10)

## ---- S4: three seller-side programs; one also moves buyers directly in the invalid variant ----------
P4 <- P$S4; P4v <- P$S4v
cols4 <- c("buyers", "n_sellers", "ad_budget", "coupon", "fee", "payout_days", "referral_bonus", "log_pop",
           "eps_buyer_true", "eps_seller_true", "h_true")
stopifnot(all(cols4 %in% names(S4)), all(cols4 %in% names(S4v)), nrow(S4) == P4$n, nrow(S4v) == P4v$n,
          all(S4$buyers > 0), all(S4$n_sellers > 0), all(S4v$buyers > 0), all(S4v$n_sellers > 0),
          P4v$theta == 0, P4$theta > 0, identical(P4[setdiff(names(P4), "theta")], P4v[setdiff(names(P4v), "theta")]))
eqB4 <- function(d, par) d$buyers - (sb[["const"]] + sb[["n_sellers"]] * d$n_sellers + sb[["ad_budget"]] * d$ad_budget +
                                      sb[["coupon"]] * d$coupon + sb[["log_pop"]] * d$log_pop + par$theta * (d$referral_bonus - par$ref_mid))
stopifnot(max(abs(eqB4(S4, P4) - S4$eps_buyer_true)) < 1e-10, max(abs(eqB4(S4v, P4v) - S4v$eps_buyer_true)) < 1e-10)
# errors are heteroskedastic in the referral bonus (scale h, mean square one)
hi <- S4v$referral_bonus > median(S4v$referral_bonus)
stopifnot(var(S4v$eps_buyer_true[hi]) > 1.5 * var(S4v$eps_buyer_true[!hi]), abs(mean(S4v$h_true^2) - 1) < 0.1,
          P4$R_mc >= 1000)
