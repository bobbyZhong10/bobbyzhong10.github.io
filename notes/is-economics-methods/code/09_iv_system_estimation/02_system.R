# ------------------------------------------------------------------
# 02_system.R -- blocks S1 and S1m: the simultaneous buyer-seller system.
#
# S1 (both equations right):
#   OLS of the buyer equation (simultaneity bias), equation-by-equation 2SLS (tsls, checked against
#   ivreg and systemfit 2SLS), 3SLS (three_sls, checked against systemfit 3SLS with
#   methodResidCov = "noDfCor"), first-stage F_fs and F_eff, the control-function endogeneity test
#   (Wooldridge 2010 eq. 6.15), the over-identification tests of the seller equation (robust Hansen J and
#   Sargan), ivreg's diagnostics, the sample reduced form Pi_hat and its identification submatrices.
#   Ruling 14: the seller equation is over-identified, so 3SLS differs from 2SLS on the buyer equation
#   (the 3SLS - 2SLS differences are stored in `diff`) and coincides with 2SLS on the seller equation.
# S1m (the seller equation omits coupon x buyers): the same estimators; 3SLS carries the seller
#   equation's error into the buyer coefficient n_sellers, 2SLS of the buyer equation does not.
# Monte Carlo: R = 1000 replications of S1 (seed S1_mc) and of S1m (seed S1m_mc): sampling
#   distribution of the buyer coefficient n_sellers under 2SLS (HC0 interval) and 3SLS (eq. 8.45 interval).
#
# Standard errors: 2SLS "robust" = HC0 (eq. 5.34), "homo" = eq. 5.27 (n - k); 3SLS = eq. 8.45 with
# Omega_hat from 2SLS residuals divided by n (eq. 8.44), "robust" = the GMM sandwich of the 3SLS weight.
# The efficiency comparison se$tsls vs se$three_sls uses the homoskedastic formulas of both
# (Theorem 8.4 needs SIV.5, which holds in S1); the robust versions are stored too. The two differ in the
# degrees-of-freedom convention (2SLS n - k, 3SLS Omega_hat / n), so se$tsls_n_denominator is the 2SLS
# homoskedastic SE on the same n denominator as the 3SLS one (sigma_hat^2 = sum eps_hat^2 / n) and
# se$ratio_n_denominator the 3SLS SE over it (a like-for-like ratio). se$tsls_n_denominator = se$tsls
# sqrt((n - k) / n) with k = 5 regressors, so the like-for-like ratio is se$ratio sqrt(n / (n - k)).
#
# Output: res_system.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("iv_fns.R")
suppressPackageStartupMessages({ library(ivreg); library(systemfit) })

tr <- readRDS("truth09.rds")
S1 <- readRDS("bramble_S1.rds")
S1m <- readRDS("bramble_S1m.rds")
b <- tr$S1$buyer[["n_sellers"]]
bS <- tr$S1$seller[["buyers"]]
se_of <- function(V, j) sqrt(V[j, j])
inst <- ~ ad_budget + coupon + fee + log_pop
eqs <- list(buyer = buyers ~ n_sellers + ad_budget + coupon + log_pop, seller = n_sellers ~ buyers + fee + log_pop)

## ---- everything for one sample -------------------------------------------------------------
fit_sample <- function(d) {
  sy <- s1_system(d)
  tB <- tsls(sy$y$buyer, sy$X$buyer, sy$Z); tBh <- tsls(sy$y$buyer, sy$X$buyer, sy$Z, vcov = "homo")
  tS <- tsls(sy$y$seller, sy$X$seller, sy$Z); tSh <- tsls(sy$y$seller, sy$X$seller, sy$Z, vcov = "homo")
  th <- three_sls(sy$y, sy$X, sy$Z); thr <- three_sls(sy$y, sy$X, sy$Z, vcov = "robust")
  oB <- ols_fit(sy$y$buyer, sy$X$buyer); oBh <- ols_fit(sy$y$buyer, sy$X$buyer, vcov = "homo")
  oS <- ols_fit(sy$y$seller, sy$X$seller)
  jB <- "buyer:n_sellers"; jS <- "seller:buyers"
  # 2SLS homoskedastic variance with sigma_hat^2 = sum eps_hat^2 / n instead of / (n - k): V_homo (n - k) / n
  se_n_den <- sqrt(tBh$V["n_sellers", "n_sellers"] * (length(sy$y$buyer) - ncol(sy$X$buyer)) / length(sy$y$buyer))
  iB <- grep("^buyer:", names(th$b)); iS <- grep("^seller:", names(th$b))
  est <- list(
    ols = list(buyer_est = oB$b[["n_sellers"]], buyer_se = se_of(oB$V, "n_sellers"), buyer_se_homo = se_of(oBh$V, "n_sellers"),
               buyer_ci = ci95(oB$b[["n_sellers"]], se_of(oB$V, "n_sellers")),
               seller_est = oS$b[["buyers"]], seller_se = se_of(oS$V, "buyers"), seller_ci = ci95(oS$b[["buyers"]], se_of(oS$V, "buyers"))),
    tsls = list(buyer_est = tB$b[["n_sellers"]], buyer_se = se_of(tB$V, "n_sellers"), buyer_se_homo = se_of(tBh$V, "n_sellers"),
                buyer_ci = ci95(tB$b[["n_sellers"]], se_of(tB$V, "n_sellers")),
                seller_est = tS$b[["buyers"]], seller_se = se_of(tS$V, "buyers"), seller_se_homo = se_of(tSh$V, "buyers"),
                seller_ci = ci95(tS$b[["buyers"]], se_of(tS$V, "buyers")), buyer_b = tB$b, seller_b = tS$b),
    three_sls = list(buyer_est = th$b[[jB]], buyer_se = se_of(th$V, jB), buyer_se_robust = se_of(thr$V, jB),
                     buyer_ci = ci95(th$b[[jB]], se_of(th$V, jB)),
                     seller_est = th$b[[jS]], seller_se = se_of(th$V, jS), seller_ci = ci95(th$b[[jS]], se_of(th$V, jS)),
                     b = th$b, V = th$V, V_robust = thr$V, Sigma = th$Sigma,
                     rho_hat = th$Sigma[1, 2] / sqrt(th$Sigma[1, 1] * th$Sigma[2, 2])),
    diff = list(buyer = unname(th$b[iB] - tB$b), seller = unname(th$b[iS] - tS$b),
                buyer_n_sellers = th$b[[jB]] - tB$b[["n_sellers"]], seller_buyers = th$b[[jS]] - tS$b[["buyers"]]),
    se = list(tsls = se_of(tBh$V, "n_sellers"), three_sls = se_of(th$V, jB), tsls_robust = se_of(tB$V, "n_sellers"),
              three_sls_robust = se_of(thr$V, jB), ratio = se_of(th$V, jB) / se_of(tBh$V, "n_sellers"),
              tsls_n_denominator = se_n_den, ratio_n_denominator = se_of(th$V, jB) / se_n_den,
              ratio_robust = se_of(thr$V, jB) / se_of(tB$V, "n_sellers"),
              seller_tsls = se_of(tSh$V, "buyers"), seller_three_sls = se_of(th$V, jS)),
    J_seller = hansen_j(sy$y$seller, sy$X$seller, sy$Z), sargan_seller = sargan(sy$y$seller, sy$X$seller, sy$Z))
  est
}

## =================================================================
## S1
## =================================================================
sy <- s1_system(S1)
f1 <- fit_sample(S1)
ivB <- ivreg(buyers ~ n_sellers + ad_budget + coupon + log_pop | ad_budget + coupon + fee + log_pop, data = S1)
ivS <- ivreg(n_sellers ~ buyers + fee + log_pop | ad_budget + coupon + fee + log_pop, data = S1)
sf2 <- systemfit(eqs, method = "2SLS", inst = inst, data = S1, methodResidCov = "noDfCor")
sf3 <- systemfit(eqs, method = "3SLS", inst = inst, data = S1, methodResidCov = "noDfCor")
tsls_hand <- list(b = c(f1$tsls$buyer_b, f1$tsls$seller_b))
tsls_pkg <- list(b = c(coef(ivB), coef(ivS)), b_systemfit = coef(sf2))
three_sls_hand <- list(b = f1$three_sls$b, V = f1$three_sls$V)
three_sls_pkg <- list(b = coef(sf3), V = vcov(sf3))
n1 <- nrow(S1); k1 <- ncol(sy$X$buyer)
stopifnot(abs(f1$se$tsls_n_denominator^2 - vcov(ivB)["n_sellers", "n_sellers"] * (n1 - k1) / n1) < 1e-12,       # ivreg uses n - k
          abs(f1$se$ratio_n_denominator - f1$se$ratio * sqrt(n1 / (n1 - k1))) < 1e-12)
stopifnot(max(abs(tsls_hand$b - tsls_pkg$b)) < 1e-10, max(abs(tsls_hand$b - tsls_pkg$b_systemfit)) < 1e-10,
          max(abs(three_sls_hand$b - three_sls_pkg$b)) < 1e-6, max(abs(three_sls_hand$V - three_sls_pkg$V)) < 1e-8)

# first stages: n_sellers on fee (buyer equation), buyers on ad_budget and coupon (seller equation)
WB <- cbind(1, S1$ad_budget, S1$coupon, S1$log_pop); WS <- cbind(1, S1$fee, S1$log_pop)
fsB <- first_stage_F(S1$n_sellers, cbind(S1$fee), WB); fsS <- first_stage_F(S1$buyers, cbind(S1$ad_budget, S1$coupon), WS)
first_stage <- list(buyer_F = fsB$F, buyer_F_eff = eff_F(S1$n_sellers, cbind(S1$fee), WB), buyer_pi_fee = fsB$pi,
                    seller_F = fsS$F, seller_F_eff = eff_F(S1$buyers, cbind(S1$ad_budget, S1$coupon), WS),
                    seller_pi_ad = fsS$pi[1], seller_pi_coupon = fsS$pi[2])
# control-function endogeneity test of n_sellers in the buyer equation
cfB <- cf_test(sy$y$buyer, sy$X$buyer, sy$Z, endog = "n_sellers")
stopifnot(abs(cfB$b[["n_sellers"]] - f1$tsls$buyer_est) < 1e-10)
endog <- list(rho = cfB$rho, se = cfB$se, t = cfB$t, p = cfB$p)
# ivreg's diagnostics (homoskedastic weak-instrument F, Wu-Hausman, Sargan)
dgB <- summary(ivB, diagnostics = TRUE)$diagnostics; dgS <- summary(ivS, diagnostics = TRUE)$diagnostics
stopifnot(abs(dgS["Sargan", "statistic"] - f1$sargan_seller$stat) < 1e-8, abs(dgB["Weak instruments", "statistic"] - fsB$F) < 1e-8)
diagnostics <- list(buyer = dgB, seller = dgS)
overid <- list(J_seller = f1$J_seller$J, J_seller_p = f1$J_seller$p, J_seller_df = f1$J_seller$df,
               sargan_seller = f1$sargan_seller$stat, sargan_seller_p = f1$sargan_seller$p)
# sample reduced form and identification submatrices (the rank condition is about the population Pi;
# in a sample it is checked with the first-stage F)
Pi_hat <- solve(crossprod(sy$Z), crossprod(sy$Z, cbind(buyers = S1$buyers, n_sellers = S1$n_sellers)))
identification <- list(Pi_hat = Pi_hat, Pi_sub_buyer = Pi_hat["fee", "n_sellers", drop = FALSE],
                       Pi_sub_seller = Pi_hat[c("ad_budget", "coupon"), "buyers", drop = FALSE],
                       rank_buyer = rank_condition(Pi_hat["fee", "n_sellers", drop = FALSE]),
                       rank_seller = rank_condition(Pi_hat[c("ad_budget", "coupon"), "buyers", drop = FALSE]))

## =================================================================
## S1m
## =================================================================
fm <- fit_sample(S1m)
misspec <- list(tsls_est = fm$tsls$buyer_est, tsls_se = fm$tsls$buyer_se, tsls_ci = fm$tsls$buyer_ci,
                tsls_bias = fm$tsls$buyer_est - b, tsls_covers = covers_num(fm$tsls$buyer_ci, b),
                three_sls_est = fm$three_sls$buyer_est, three_sls_se = fm$three_sls$buyer_se, three_sls_ci = fm$three_sls$buyer_ci,
                three_sls_bias = fm$three_sls$buyer_est - b, three_sls_covers = covers_num(fm$three_sls$buyer_ci, b),
                three_sls_minus_tsls = fm$three_sls$buyer_est - fm$tsls$buyer_est,
                seller_tsls_est = fm$tsls$seller_est, seller_tsls_se = fm$tsls$seller_se, seller_tsls_ci = fm$tsls$seller_ci,
                seller_three_sls_est = fm$three_sls$seller_est, seller_diff = fm$diff$seller,
                ols_est = fm$ols$buyer_est, rho_hat = fm$three_sls$rho_hat,
                J_seller = fm$J_seller$J, J_seller_p = fm$J_seller$p, sargan_seller = fm$sargan_seller$stat,
                sargan_seller_p = fm$sargan_seller$p, se_three_sls = fm$se$three_sls, se_tsls_homo = fm$se$tsls)

## =================================================================
## Monte Carlo: S1 and S1m
## =================================================================
z <- qnorm(0.975)
one_rep <- function(par, kappa) {
  d <- gen_S1(par, kappa = kappa)
  s <- s1_system(d)
  tB <- tsls(s$y$buyer, s$X$buyer, s$Z)
  th <- three_sls(s$y, s$X, s$Z)
  c(tsls = tB$b[["n_sellers"]], se_tsls = sqrt(tB$V["n_sellers", "n_sellers"]),
    three_sls = th$b[["buyer:n_sellers"]], se_three_sls = sqrt(th$V["buyer:n_sellers", "buyer:n_sellers"]))
}
run_mc <- function(par, kappa, R, seed) {
  set.seed(seed)
  as.data.frame(t(vapply(seq_len(R), function(r) one_rep(par, kappa), numeric(4))))
}
summ <- function(sim) {
  list(R = nrow(sim), mean_tsls = mean(sim$tsls), sd_tsls = sd(sim$tsls), mean_three_sls = mean(sim$three_sls),
       sd_three_sls = sd(sim$three_sls), median_tsls = median(sim$tsls), median_three_sls = median(sim$three_sls),
       mean_se_tsls = mean(sim$se_tsls), mean_se_three_sls = mean(sim$se_three_sls),
       cover_tsls = mean(abs(sim$tsls - b) <= z * sim$se_tsls), cover_three_sls = mean(abs(sim$three_sls - b) <= z * sim$se_three_sls),
       sd_ratio = sd(sim$three_sls) / sd(sim$tsls))
}
R_mc <- 1000
t0 <- proc.time()[["elapsed"]]
sim1 <- run_mc(tr$params$S1, 0, R_mc, tr$params$seed$S1_mc)
sim1m <- run_mc(tr$params$S1m, tr$params$S1m$kappa, R_mc, tr$params$seed$S1m_mc)
mc_seconds <- proc.time()[["elapsed"]] - t0
mc <- list(S1 = summ(sim1), S1m = summ(sim1m), draws = list(S1 = sim1, S1m = sim1m))

## =================================================================
## Save
## =================================================================
res <- list(tsls_hand = tsls_hand, tsls_pkg = tsls_pkg, three_sls_hand = three_sls_hand, three_sls_pkg = three_sls_pkg,
            ols = f1$ols, tsls = f1$tsls[setdiff(names(f1$tsls), c("buyer_b", "seller_b"))],
            three_sls = f1$three_sls[setdiff(names(f1$three_sls), c("b", "V", "V_robust"))],
            se = f1$se, diff = f1$diff, first_stage = first_stage, endog = endog, overid = overid,
            diagnostics = diagnostics, identification = identification, misspec = misspec, mc = mc,
            n = nrow(S1))
saveRDS(res, "res_system.rds")

cat(sprintf("S1 buyer n_sellers: OLS %.4f [%.4f, %.4f]; 2SLS %.4f [%.4f, %.4f] (HC0 se %.4f, homo %.4f); 3SLS %.4f (se %.4f, robust %.4f)\n",
            f1$ols$buyer_est, f1$ols$buyer_ci[1], f1$ols$buyer_ci[2], f1$tsls$buyer_est, f1$tsls$buyer_ci[1], f1$tsls$buyer_ci[2],
            f1$tsls$buyer_se, f1$tsls$buyer_se_homo, f1$three_sls$buyer_est, f1$three_sls$buyer_se, f1$three_sls$buyer_se_robust))
cat(sprintf("S1 3SLS - 2SLS: buyer n_sellers %.5f, max |seller diff| %.2e; SE ratio 3SLS / 2SLS %.4f (robust %.4f)\n",
            f1$diff$buyer_n_sellers, max(abs(f1$diff$seller)), f1$se$ratio, f1$se$ratio_robust))
cat(sprintf("S1 seller buyers: 2SLS %.4f (se %.4f); J %.3f (p %.3f), Sargan %.3f (p %.3f); F_fs buyer %.1f (eff %.1f), seller %.1f (eff %.1f); CF t %.2f\n",
            f1$tsls$seller_est, f1$tsls$seller_se, overid$J_seller, overid$J_seller_p, overid$sargan_seller, overid$sargan_seller_p,
            first_stage$buyer_F, first_stage$buyer_F_eff, first_stage$seller_F, first_stage$seller_F_eff, endog$t))
cat(sprintf("S1m buyer n_sellers: 2SLS %.4f [%.4f, %.4f] (bias %.4f); 3SLS %.4f [%.4f, %.4f] (bias %.4f); seller J %.2f (p %.2e)\n",
            misspec$tsls_est, misspec$tsls_ci[1], misspec$tsls_ci[2], misspec$tsls_bias, misspec$three_sls_est,
            misspec$three_sls_ci[1], misspec$three_sls_ci[2], misspec$three_sls_bias, misspec$J_seller, misspec$J_seller_p))
for (k in c("S1", "S1m")) with(mc[[k]], cat(sprintf(
  "MC %s (R = %d): 2SLS mean %.4f sd %.4f cover %.3f | 3SLS mean %.4f sd %.4f cover %.3f | sd ratio %.3f\n",
  k, R, mean_tsls, sd_tsls, cover_tsls, mean_three_sls, sd_three_sls, cover_three_sls, sd_ratio)))
cat(sprintf("Monte Carlo time %.1f s\nSaved res_system.rds\n", mc_seconds))
