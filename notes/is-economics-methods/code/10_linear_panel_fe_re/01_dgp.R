# ============================================================================
# 01_dgp.R  --  Nimbus marketplace merchant panel, STATIC facet (chapter 10)
# ----------------------------------------------------------------------------
# "Nimbus" is a fictional marketplace. We observe merchants i over months t.
#
#  STATIC facet -- FE vs RE and the Hausman test.
#      y_it = beta * x_it + alpha_i + eps_it,   with Corr(alpha_i, x_it) > 0
#      (better merchants, high alpha_i, also run more promotions x_it). Pooled
#      OLS and random effects are biased for beta because they treat alpha_i as
#      uncorrelated with x; fixed effects (within) is consistent; Hausman
#      rejects RE. This teaches FE, strict exogeneity, and the Hausman test.
#
# All TRUE parameters known. Seed fixed. Scripts use relative paths; run with
# the working directory set to this folder. (The dynamic facet, Nickell bias
# and Arellano-Bond / Blundell-Bond, lives in code/11_dynamic_panel/.)
# ============================================================================

set.seed(737)

## ---- static panel: FE vs RE ------------------------------------------------
gen_static <- function(N = 600, T = 6, beta = 0.5, ca = 0.7) {
  alpha <- rnorm(N, 0, 1)                              # merchant quality (unobs.)
  id <- rep(1:N, each = T); tt <- rep(1:T, times = N)
  a  <- alpha[id]
  x  <- ca * a + rnorm(N*T, 0, 1)                      # promotions, corr with alpha
  y  <- beta * x + a + rnorm(N*T, 0, 1)
  data.frame(id, t = tt, y, x)
}
static <- gen_static()

saveRDS(list(static = static,
             truth = list(beta = 0.5, ca = 0.7)),
        "nimbus.rds")

cat(sprintf("static: N=600 T=6, beta=0.5, Corr(alpha,x) built via ca=0.7\n"))
cat(sprintf("  implied pooled-OLS beta plim = beta + Cov(x,a)/Var(x) = %.3f\n",
            0.5 + 0.7/(0.7^2*1 + 1)))

# ============================================================================
# Lines 1 to 37 above are the inherited legacy block P1, byte-identical (Phase 3: inherited blocks
# first and unchanged; tests/test_dgp.R checks the md5 of those lines). Everything below is new in v3.
# ----------------------------------------------------------------------------
# 01_dgp.R (v3) -- the Nimbus simulation world for chapter 10 (linear panels: FE, RE, FD and CRE).
#
# Nimbus is a (fictional) online marketplace; unit i is a merchant, period t a month. Money is in US
# dollars: Y_it is monthly GMV in thousands of USD, X_it is monthly promotion spend in hundreds of USD
# (the inherited block P1 keeps the legacy standardised units). One world, eight blocks:
#   P1   inherited static panel (above): FE vs RE vs pooled OLS, Hausman.            nimbus.rds
#   P2   FE vs FD: the same promotions and monthly shocks, iid errors vs random-walk errors;
#        Monte Carlo (04_fe_fd_feedback.R) for the sampling SDs.                      nimbus_P2.rds
#   P3   feedback: promotions respond to last month's shock (strict exogeneity fails). nimbus_P3.rds
#   P4   Hausman-Taylor: time-invariant Nimbus Verified status correlated with alpha_i,
#        fulfilment share an exogenous time-varying regressor (its mean is the instrument). nimbus_P4.rds
#   P5   random slopes correlated with the within-merchant promotion variability.    nimbus_P5.rds
#   P5b  random slopes correlated with the level Xbar_i only.                         nimbus_P5b.rds
#   P6   unbalanced panel: attrition on alpha_i vs attrition on the shock.            nimbus_P6.rds
#   P7   heteroskedastic, serially correlated errors (RE.1 holds, RE.3 fails); one sample here,
#        the size Monte Carlo of the classical and robust Hausman tests in 07_mc_hausman.R. nimbus_P7.rds
#   truth10.rds   list(P1, P2, P3, P4, P5, P5b, P6, P7, params): truth, closed-form probability limits,
#                 large-sample checks, parameters and seeds.
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 1.8): Y_it, X_it, alpha_i, eps_it; T_i observed periods; varphi the RE
# quasi-demeaning parameter; b_i random slopes, beta = E[b_i], beta_* the FE probability limit
# (P5$fe_plim); W_i and gamma for Hausman-Taylor; s_it the selection indicator of P6.
#
# Seeds (every block its own set.seed; P1 keeps the legacy 737): P2 1002, P3 1003, P4 1004, P5 1005,
# P5b 1015, P6 1006, P7 1007; Monte Carlo P2 1012, P7 1017; large-sample checks P1 1021, P3 1023,
# P5 1025, P5b 1035, P6 1026.
# Run with the working directory set to this folder (00_run_all.R does this), like the legacy block.
# ============================================================================
source("panel_fns.R")

seed <- list(P1 = 737, P2 = 1002, P3 = 1003, P4 = 1004, P5 = 1005, P5b = 1015, P6 = 1006, P7 = 1007,
             P2_mc = 1012, P7_mc = 1017, P1_large = 1021, P3_large = 1023, P5_large = 1025, P5b_large = 1035,
             P6_large = 1026)

## =================================================================
## P1: truth of the inherited block (closed form) and a large-sample check
## =================================================================
## Legacy DGP: alpha_i ~ N(0, 1), X_it = 0.7 alpha_i + N(0, 1), Y_it = 0.5 X_it + alpha_i + N(0, 1), N = 600,
## T = 6 (unitless, as in the legacy text). Probability limits (N -> infinity, T = 6):
##   pooled OLS: beta + ca / (ca^2 + 1)                       (the legacy 01_dgp.R prints 0.970)
##   between:    beta + ca / (ca^2 + 1 / T)                   (Cov(Xbar_i, alpha_i) / Var(Xbar_i))
##   Swamy-Arora: sigma2_eps -> 1; T sigma2_alpha + sigma2_eps -> T (sigma2_alpha + 1/T - ca^2 / (ca^2 + 1/T)),
##   the between residual variance, because alpha_i is partly absorbed by Xbar_i; varphi_plim follows, and
##   RE = w BE + (1 - w) FE with w = (1 - varphi)^2 T Var(Xbar) / ((T - 1) + (1 - varphi)^2 T Var(Xbar)).
##   varphi_components = 1 - sqrt(1 / (1 + T)) is the value RE would use with the true variance components.
P1 <- list(N = 600, T = 6, beta = 0.5, ca = 0.7, sigma_alpha = 1, sigma_x = 1, sigma_eps = 1)
v_xbar <- P1$ca^2 + 1 / P1$T
s21_plim <- P1$T * (1 + 1 / P1$T - P1$ca^2 / v_xbar)
varphi_plim <- 1 - sqrt(1 / s21_plim)
w_plim <- (1 - varphi_plim)^2 * P1$T * v_xbar / ((P1$T - 1) + (1 - varphi_plim)^2 * P1$T * v_xbar)
truth_P1 <- list(beta = P1$beta, plim_pols = P1$beta + P1$ca / (P1$ca^2 + 1), plim_be = P1$beta + P1$ca / v_xbar,
                 varphi_plim = varphi_plim, w_plim = w_plim,
                 plim_re = P1$beta + w_plim * P1$ca / v_xbar,
                 varphi_components = 1 - sqrt(1 / (1 + P1$T)))
set.seed(seed$P1_large)
big <- gen_static(N = 1e5)
reL <- re_fit(big$y, big$x, big$id)
truth_P1$large <- list(N = 1e5, pols = unname(coef(lm.fit(cbind(1, big$x), big$y))[2]),
                       be = unname(be_fit(big$y, big$x, big$id)$b["x"]), fe = unname(fe_fit(big$y, big$x, big$id)$b),
                       re = unname(reL$b["x"]), varphi = reL$varphi)
rm(big, reL)

## =================================================================
## P2: FE vs FD, iid errors vs random-walk errors
## =================================================================
## 500 merchants, 6 months. alpha_i ~ N(0, 1.5^2) (thousands of USD of monthly GMV). Promotion spend
## X_it = 6 + 0.4 alpha_i + u_it (600 USD a month on average; better merchants promote more), u_it a
## stationary AR(1) with rho 0.5 and unit variance (budgets are persistent). GMV
## Y_it = 10 + 0.5 X_it + alpha_i + eps_it: each extra 100 USD of promotion adds 500 USD of GMV.
## The same monthly shocks e_it ~ N(0, 1) enter in two ways: eps_it = e_it (shocks wash out, FE.3) or
## eps_it = e_i1 + ... + e_it (a random walk, shocks accumulate: a review or a stock-out changes the
## merchant's standing for good; FD.3). Monte Carlo: R_mc = 1000 samples of the same design.
P2 <- list(N = 500, T = 6, alpha_sd = 1.5, x_mean = 6, x_alpha = 0.4, rho_x = 0.5, beta = 0.5, mu0 = 10,
           sigma_e = 1, R_mc = 1000)
set.seed(seed$P2)
dP2 <- gen_P2(P2)
truth_P2 <- list(beta = P2$beta)

## =================================================================
## P3: feedback from last month's shock to this month's promotions
## =================================================================
## 1,000 merchants, 6 months. X_it = 6 + 0.4 alpha_i - 0.6 eps_{i,t-1} + u_it, u_it ~ N(0, 1): after a bad
## month (eps < 0) a merchant spends more on promotion. Y_it = 10 + 0.5 X_it + alpha_i + eps_it, eps ~ N(0, 1).
## X is contemporaneously exogenous but not strictly exogenous: Cov(X_{i,t+1}, eps_it) = -0.6.
## Probability limits (N -> infinity), with v = 0.6^2 + 1 the within variance of X:
##   FE: beta - feedback sigma2 / (T v) = 0.5 + 0.6 / (6 * 1.36)    (bias of order 1 / T)
##   FD: beta - feedback sigma2 / (2 v) = 0.5 + 0.6 / 2.72          (bias does not shrink with T)
P3 <- list(N = 1000, T = 6, alpha_sd = 1.5, x_mean = 6, x_alpha = 0.4, feedback = -0.6, x_sd = 1, beta = 0.5,
           mu0 = 10, sigma_e = 1)
set.seed(seed$P3)
dP3 <- gen_P3(P3)
vw <- P3$feedback^2 * P3$sigma_e^2 + P3$x_sd^2
truth_P3 <- list(beta = P3$beta, fe_plim = P3$beta - P3$feedback * P3$sigma_e^2 / (P3$T * vw),
                 fd_plim = P3$beta - P3$feedback * P3$sigma_e^2 / (2 * vw))
set.seed(seed$P3_large)
big <- gen_P3(modifyList(P3, list(N = 1e5)))
truth_P3$large <- list(N = 1e5, fe = unname(fe_fit(big$y, big$x, big$id)$b),
                       fd = unname(fd_fit(big$y, big$x, big$id, big$t)$b["x"]))
rm(big)

## =================================================================
## P4: Hausman-Taylor, certified merchants
## =================================================================
## 2,000 merchants, 6 months. alpha_i ~ N(0, 2^2). mu1_i ~ N(0, 1) is how close the merchant is to a Nimbus
## warehouse (geography, unrelated to alpha_i). x1_it = 5 + mu1_i + N(0, 0.6^2): the share of orders shipped
## by Nimbus Fulfilment, in tens of percentage points (50% on average). Nimbus Verified status needs a
## fulfilment record, and better merchants apply more: W_i = 1{1.2 mu1_i + 0.35 alpha_i + v_i > 0.3}, v ~ N(0, 1)
## (about 43% certified). x2_it = 6 + 0.3 alpha_i + N(0, 1) is promotion spend. Y_it = 8 + 0.3 x1_it + 0.5 x2_it
## + 1.5 W_i + alpha_i + eps_it, eps ~ N(0, 1.5^2): certification adds gamma = 1,500 USD of monthly GMV.
## Hausman-Taylor groups: X_1 = x1 (exogenous), X_2 = x2 (correlated with alpha_i), W_2 = W (correlated);
## Xbar_1i is the instrument for W_i (order condition dim X_1 = 1 >= dim W_2 = 1).
P4 <- list(N = 2000, T = 6, alpha_sd = 2, x1_mean = 5, x1_mu = 1.0, x1_sd = 0.6, x2_mean = 6, x2_alpha = 0.3,
           x2_sd = 1, w_mu = 1.2, w_alpha = 0.35, w_cut = 0.3, beta1 = 0.3, beta2 = 0.5, gamma = 1.5, mu0 = 8,
           sigma_e = 1.5)
set.seed(seed$P4)
dP4 <- gen_P4(P4)
truth_P4 <- list(gamma = P4$gamma, beta1 = P4$beta1, beta2 = P4$beta2)

## =================================================================
## P5 and P5b: random slopes
## =================================================================
## 800 merchants, 8 months. alpha_i ~ N(0, 1.5^2); promotion level mu_i = 10 + 0.4 alpha_i + N(0, 0.5^2)
## (1,000 USD a month); four merchant types by the within-merchant SD of their promotion spend,
## s_i in {0.5, 1, 1.5, 2} (50 to 200 USD), one quarter each: steady promoters to flash-sale merchants.
## X_it = mu_i + s_i u_it, u ~ N(0, 1). Y_it = 8 + b_i X_it + alpha_i + eps_it, eps ~ N(0, 1).
## P5:  b_i = 0.5 + 0.15 (s_i^2 - E[s^2]) + eta_i, eta ~ N(0, 0.1^2): the merchants who vary promotions most
##      are those for whom promotions work best. Average slope beta = E[b_i] = 0.5; FE converges to the
##      within-variance-weighted average beta_* = E[b_i s_i^2] / E[s_i^2] = 0.5 + 0.15 Var(s^2) / E[s^2]
##      (Wooldridge 2010 eqs. 11.84 to 11.85: E[Sum_t Xdd_it^2 | s_i] = (T - 1) s_i^2).
## P5b: b_i = 0.5 + 0.3 (mu_i - 10) + eta_i: slopes correlate with the level Xbar_i only, s_i independent of
##      b_i, so E[b_i | Xdd_it] = E[b_i] and FE converges to the average slope 0.5 (Wooldridge 2005).
P5 <- list(N = 800, T = 8, alpha_sd = 1.5, mu_mean = 10, mu_alpha = 0.4, mu_sd = 0.5, s_types = c(0.5, 1, 1.5, 2),
           s_probs = rep(0.25, 4), avg_slope = 0.5, kappa = 0.15, eta_sd = 0.1, design = "within_sd", mu0 = 8,
           sigma_e = 1)
P5b <- modifyList(P5, list(kappa = 0.3, design = "level"))
set.seed(seed$P5)
dP5 <- gen_P5(P5)
set.seed(seed$P5b)
dP5b <- gen_P5(P5b)
Es2 <- sum(P5$s_probs * P5$s_types^2); Vs2 <- sum(P5$s_probs * P5$s_types^4) - Es2^2
w_type <- P5$s_probs * P5$s_types^2 / Es2                          # FE weight of each type
truth_P5 <- list(avg_slope = P5$avg_slope, fe_plim = P5$avg_slope + P5$kappa * Vs2 / Es2, Es2 = Es2, Vs2 = Vs2,
                 type_slope = P5$avg_slope + P5$kappa * (P5$s_types^2 - Es2), type_weight_fe = w_type,
                 type_weight_avg = P5$s_probs)
truth_P5b <- list(avg_slope = P5b$avg_slope, fe_plim = P5b$avg_slope)
set.seed(seed$P5_large)
big <- gen_P5(modifyList(P5, list(N = 1e5)))
truth_P5$large <- list(N = 1e5, fe = unname(fe_fit(big$y, big$x, big$id)$b), mg = mg_fit(big$y, big$x, big$id)$est)
set.seed(seed$P5b_large)
big <- gen_P5(modifyList(P5b, list(N = 1e5)))
truth_P5b$large <- list(N = 1e5, fe = unname(fe_fit(big$y, big$x, big$id)$b), mg = mg_fit(big$y, big$x, big$id)$est)
rm(big)

## =================================================================
## P6: unbalanced panel, two attrition rules
## =================================================================
## 4,000 merchants, 8 months. alpha_i ~ N(0, 1.5^2); X_it = 6 + 0.5 alpha_i + N(0, 1);
## Y_it = 9 + 0.5 X_it + alpha_i + eps_it, eps ~ N(0, 1.5^2) (monthly GMV about 12,000 USD).
## s_alpha: from month 2 an active merchant leaves with probability pnorm(-1.5 - 0.4 alpha_i) (about 18% a
##          month at alpha_i = -1.5, 2% at alpha_i = 1.5): weaker merchants leave, FE stays consistent.
## s_eps:   a merchant leaves in the first month t >= 2 in which its GMV falls below 9,500 USD and does not
##          report that month: selection on eps_it through Y_it, FE is inconsistent.
P6 <- list(N = 4000, T = 8, alpha_sd = 1.5, x_mean = 6, x_alpha = 0.5, beta = 0.5, mu0 = 9, sigma_e = 1.5,
           haz0 = -1.5, haz_alpha = -0.4, floor = 9.5)
set.seed(seed$P6)
dP6 <- gen_P6(P6)
truth_P6 <- list(beta = P6$beta)
set.seed(seed$P6_large)
big <- gen_P6(modifyList(P6, list(N = 1e5)))
ka <- big$s_alpha == 1; ke <- big$s_eps == 1
truth_P6$large <- list(N = 1e5, fe_alpha = unname(fe_fit(big$y[ka], big$x[ka], big$id[ka])$b),
                       fe_eps = unname(fe_fit(big$y[ke], big$x[ke], big$id[ke])$b),
                       keep_alpha = mean(ka), keep_eps = mean(ke))
rm(big, ka, ke)

## =================================================================
## P7: RE.1 holds, RE.3 fails (one sample; the Monte Carlo is in 07_mc_hausman.R)
## =================================================================
## 1,000 merchants, 5 months. X_it = 6 + mu_i + u_it, mu_i ~ N(0, 1), u_it a stationary AR(1) (rho 0.6, unit
## variance). alpha_i = (0.3 + 1.0 |Xbar_i - 6|) z_i, z ~ N(0, 1): mean zero given X_i (RE.1(b) holds), but
## merchants whose average promotion is far from typical are more heterogeneous. eps_it = (0.5 + 0.5 |X_it - 6|)
## e_it, e_it a stationary AR(1) (rho 0.6): heteroskedastic and serially correlated (RE.3 fails).
## Y_it = 22 + 0.5 X_it + alpha_i + eps_it. RE is consistent but not efficient, so the classical Hausman
## statistic is not chi^2 (Wooldridge 2010 §10.7.3, p. 333). Monte Carlo: R_mc = 2000 samples.
P7 <- list(N = 1000, T = 5, x_mean = 6, mu_sd = 1, rho_x = 0.6, alpha_base = 0.3, alpha_slope = 1.0,
           eps_base = 0.5, eps_slope = 0.5, rho_e = 0.6, beta = 0.5, mu0 = 22, R_mc = 2000)
set.seed(seed$P7)
dP7 <- gen_P7(P7)
truth_P7 <- list(beta = P7$beta)

## =================================================================
## Save
## =================================================================
saveRDS(dP2, "nimbus_P2.rds"); saveRDS(dP3, "nimbus_P3.rds"); saveRDS(dP4, "nimbus_P4.rds")
saveRDS(dP5, "nimbus_P5.rds"); saveRDS(dP5b, "nimbus_P5b.rds"); saveRDS(dP6, "nimbus_P6.rds")
saveRDS(dP7, "nimbus_P7.rds")
tr <- list(P1 = truth_P1, P2 = truth_P2, P3 = truth_P3, P4 = truth_P4, P5 = truth_P5, P5b = truth_P5b,
           P6 = truth_P6, P7 = truth_P7,
           params = list(seed = seed, P1 = P1, P2 = P2, P3 = P3, P4 = P4, P5 = P5, P5b = P5b, P6 = P6, P7 = P7))
saveRDS(tr, "truth10.rds")

cat(sprintf("P1 plims: pooled %.4f, between %.4f, RE %.4f (varphi %.4f); large sample: %.4f %.4f %.4f %.4f\n",
            tr$P1$plim_pols, tr$P1$plim_be, tr$P1$plim_re, tr$P1$varphi_plim,
            tr$P1$large$pols, tr$P1$large$be, tr$P1$large$re, tr$P1$large$varphi))
cat(sprintf("P3 plims: FE %.4f FD %.4f; large sample FE %.4f FD %.4f\n",
            tr$P3$fe_plim, tr$P3$fd_plim, tr$P3$large$fe, tr$P3$large$fd))
cat(sprintf("P4: certified share %.3f; P5: avg slope %.3f, FE plim %.4f (large %.4f, MG %.4f); P5b large FE %.4f\n",
            mean(dP4$w[dP4$t == 1]), tr$P5$avg_slope, tr$P5$fe_plim, tr$P5$large$fe, tr$P5$large$mg, tr$P5b$large$fe))
cat(sprintf("P6: kept %.3f (alpha) %.3f (shock); large-sample FE %.4f (alpha) %.4f (shock)\n",
            mean(dP6$s_alpha), mean(dP6$s_eps), tr$P6$large$fe_alpha, tr$P6$large$fe_eps))
cat("Saved nimbus_P2.rds ... nimbus_P7.rds, truth10.rds\n")
