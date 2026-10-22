# ------------------------------------------------------------------
# 01_dgp.R -- the Marigold simulation world for chapter 13 (quantile regression)
#
# Marigold is a (fictional) e-commerce platform that hosts online storefronts for independent
# merchants. Money is in US dollars. One world, seven data sets (all new in v3; nothing is
# inherited from a legacy chapter). Quantile levels taus = (0.1, 0.25, 0.5, 0.75, 0.9).
#   Q1   AutoPrice experiment: an algorithmic pricing tool is switched on for a random half of
#        5,002 merchants. Monthly revenue (USD thousands). The tool leaves mean revenue unchanged
#        (ATE = 0 exactly) but spreads the distribution: small merchants lose, large ones gain.
#        Covariate X: standardised log revenue of the previous year.
#   Q1r  the same merchants, the same marginal distributions of Y(0) and Y(1) as Q1 (so the same
#        QTE at every tau), with Y(1) re-paired to Y(0) (rank invariance fails): the quantiles of
#        the individual effects tau_i differ from the QTE.
#   Q2   Fulfilled by Marigold (FBM), a covariate whose share could change: log monthly revenue
#        rises by 0.5 at every conditional quantile (a pure location shift given log catalogue
#        size), yet the unconditional quantile effect of a larger FBM share (UQPE) varies with tau.
#   Q3   AdBoost (sponsored listings) adoption is endogenous (merchants with better unobserved
#        prospects adopt); a randomised ad-credit offer Z is the instrument. Log monthly revenue.
#   Q4   monthly merchant panel with merchant fixed effects: log revenue on the share of the
#        catalogue on promotion, panel lengths \mathcal{T} in {5, 20} (single samples here; the Monte
#        Carlo is in 05_panel_mc.R).
#   Q5   One-Click Checkout launched in one region: two regions (S_i), two quarters (T_i), repeated
#        cross sections; the time effect is nonlinear in revenue (holiday growth is larger for
#        larger merchants) and the treated region's merchants are larger, so QDID is biased and CIC
#        is not. Monthly revenue (USD thousands).
#   Q5f  the fuzzy version: adoption is voluntary (always-takers, switchers, never-takers; a stable
#        adoption rate in the control region), for the switchers' LQTE (dCDH 2018).
# Every number in chapter 13 comes from these data and from truth13.rds.
#
# Outputs (in this directory):
#   marigold_Q1.rds    data.frame: merchant x D y y0_true y1_true
#   marigold_Q1r.rds   data.frame: merchant x D y y0_true y1_true
#   marigold_Q2.rds    data.frame: merchant fbm x logy y
#   marigold_Q3.rds    data.frame: merchant x Z D logy y D0_true D1_true type_true logy0_true logy1_true zu_true
#   marigold_Q4.rds    list(T5, T20): data.frames id t x y alpha_true
#   marigold_Q5.rds    data.frame: merchant S T D y y0_true y1_true z_true
#   marigold_Q5f.rds   data.frame: merchant S T D y type_true z_true
#   truth13.rds        list(Q1, Q1r, Q2, Q3, Q4, Q5, params): truth and parameters
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 1.11, 2.1, 2.2): tau is the quantile level; tau_i = Y_i(1) - Y_i(0) the
# individual effect (in formulas with tau_i the level is written u); QTE(tau), QTT(tau), LQTE(tau),
# UQPE(tau); beta(tau) the QR coefficient, beta_D(tau) the IV-QR structural effect; S_i the group
# (region), T_i the period indicator (as in chapter 20), \mathcal{T} the panel length (written Tn or
# T_short / T_long here); Q_d the CIC map (Q_{d=0} for the untreated). Field names follow the card.
#
# Seeds (every block its own set.seed, nothing inherited): Q1 1301, Q1r 1302, Q2 1303, Q3 1304,
# Q4 single samples 1305 (\mathcal{T} = 5) and 1306 (\mathcal{T} = 20), Q5 1307, Q5f 1308; population
# checks Q1 1311, Q3 1314, Q5 1317, Q5f 1318. Seeds used by the estimation scripts are stored here too:
# bootstraps Q1 13101, Q1r 13102, Q2 13201, Q3 13301, Q5 13501, Q5f 13502, toy data for the chapter-20 cross-check
# 13601, and the Q4 Monte Carlo 134000 + r for replication r.
# Seeds were fixed before any result was seen and are not changed afterwards (ruling 25).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("qr_fns.R")
if (file.exists("res_extra.rds")) invisible(file.remove("res_extra.rds"))      # results from older data are stale

taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
taus_fine <- round(seq(0.05, 0.95, by = 0.05), 2)                     # grid for the figures (contains taus)
seed <- list(Q1 = 1301, Q1r = 1302, Q2 = 1303, Q3 = 1304, Q4_T5 = 1305, Q4_T20 = 1306, Q5 = 1307, Q5f = 1308,
             Q1_pop = 1311, Q3_pop = 1314, Q5_pop = 1317, Q5f_pop = 1318,
             Q1_boot = 13101, Q1r_boot = 13102, Q2_boot = 13201, Q3_boot = 13301, Q5_boot = 13501, Q5f_boot = 13502,
             cross_toy = 13601, Q4_mc = 134000)
stopifnot(all(taus %in% taus_fine))
qtile_root <- function(cdf, t, lo, hi) uniroot(function(q) cdf(q) - t, c(lo, hi), tol = 1e-13)$root

## =================================================================
## Q1: the AutoPrice experiment (revenue in USD thousands)
## =================================================================
## 5,002 merchants, complete randomisation: AutoPrice (D = 1) for exactly 2,501 (both arms odd, so
## every within-arm sample quantile is unique). X ~ N(0, 1) is the standardised log of last year's
## revenue. Without the tool log Y(0) = m + b X + s e: median 10 (USD 10,000 a month), log sd
## sqrt(v) = 0.5 (v = b^2 + s^2). With the tool Y(1) = c Y(0)^k, k = 1.4: AutoPrice raises prices
## where demand is strong and cuts them where it is thin, so revenue moves away from the middle;
## c is chosen so that E[Y(1)] = E[Y(0)] (ATE = 0 exactly):
##   E[Y(0)] = exp(m + v / 2),  E[Y(1)] = c exp(k m + k^2 v / 2),  log c = (1 - k) m + (1 - k^2) v / 2.
## Marginal quantiles: F^{-1}_{Y(0)}(tau) = exp(m + sqrt(v) z_tau), F^{-1}_{Y(1)}(tau) = c F^{-1}_{Y(0)}(tau)^k,
##   QTE(tau) = c q0^k - q0 (USD thousands; about -1.65 at tau = 0.1 and +2.77 at tau = 0.9).
## Rank invariance holds (Y(1) is increasing in Y(0)), so QTE(tau) is the effect on the merchant at
## rank tau; the individual effect tau_i = c Y0^k - Y0 is U-shaped in Y0 (smallest near Y0 = 5.8), so
## even here the quantiles of tau_i are not the QTE (they coincide only when QTE(u) is monotone in u).
## Log scale: QTE in log points = log c + (k - 1)(m + sqrt(v) z_tau); conditional on X = x the
## tau-quantile of log Y(d) is linear in x, and the conditional QTE at X = 0 is
## log c + (k - 1)(m + s z_tau): smaller in absolute value, because the within-X spread s < sqrt(v).
P1 <- list(n = 5002, n_treat = 2501, m = log(10), b = 0.4, s = 0.3, k = 1.4, rho_r = 0.5, n_pop = 1e6)
P1$v <- P1$b^2 + P1$s^2
P1$logc <- (1 - P1$k) * P1$m + (1 - P1$k^2) * P1$v / 2
set.seed(seed$Q1)
Q1 <- gen_Q1(P1)

q0_of <- function(t) exp(P1$m + sqrt(P1$v) * qnorm(t))
q1_of <- function(t) exp(P1$logc + P1$k * log(q0_of(t)))                 # c q0^k
qte_Q1 <- function(t) q1_of(t) - q0_of(t)
ate_check <- exp(P1$logc + P1$k * P1$m + P1$k^2 * P1$v / 2) - exp(P1$m + P1$v / 2)
stopifnot(abs(ate_check) < 1e-12)
# quantiles of tau_i = g(Y0) = c Y0^k - Y0 under rank invariance (Q1)
g1 <- function(y) exp(P1$logc) * y^P1$k - y
ystar <- (1 / (exp(P1$logc) * P1$k))^(1 / (P1$k - 1)); gmin <- g1(ystar)
F0 <- function(y) pnorm((log(y) - P1$m) / sqrt(P1$v))
G1 <- function(t) {                                                      # P(g(Y0) <= t)
  if (t <= gmin) return(0)
  hi <- uniroot(function(y) g1(y) - t, c(ystar, 1e4), tol = 1e-14)$root
  lo <- if (t >= 0) 0 else uniroot(function(y) g1(y) - t, c(1e-12, ystar), tol = 1e-14)$root
  F0(hi) - (if (lo > 0) F0(lo) else 0)
}
ind_Q1 <- function(t) qtile_root(Vectorize(G1), t, gmin + 1e-12, 200)
truth_Q1 <- list(taus = taus, qte = qte_Q1(taus), ate = 0, ate_check = ate_check,
                 q0 = q0_of(taus), q1 = q1_of(taus), mean_y = exp(P1$m + P1$v / 2),
                 sd_y0 = exp(P1$m + P1$v / 2) * sqrt(exp(P1$v) - 1), sd_y1 = exp(P1$m + P1$v / 2) * sqrt(exp(P1$k^2 * P1$v) - 1),
                 qte_log = P1$logc + (P1$k - 1) * (P1$m + sqrt(P1$v) * qnorm(taus)),
                 cqte_log = P1$logc + (P1$k - 1) * (P1$m + P1$s * qnorm(taus)),
                 ind_effect_quantiles = vapply(taus, ind_Q1, numeric(1)),
                 ind_min = gmin, y_at_ind_min = ystar, share_lose = G1(0),
                 fine = data.frame(tau = taus_fine, qte = qte_Q1(taus_fine), ind_Q1 = vapply(taus_fine, ind_Q1, numeric(1))))

## =================================================================
## Q1r: the same marginals, effects re-paired
## =================================================================
## The same 5,002 merchants, X, D and Y(0); Y(1) has exactly Q1's marginal distribution but is tied to
## Y(0) through a Gaussian copula with correlation rho_r = 0.5: the tool reshuffles merchants' ranks.
## The data of an experiment (one potential outcome per merchant) have the same distribution in Q1 and
## Q1r, so no estimator can tell them apart; the QTE is identical, the distribution of tau_i is not.
## Joint law: (log Y0, log Y1r) bivariate normal; given Y0 = exp(m + sqrt(v) z0),
##   log Y1r ~ N(log c + k (m + sqrt(v) rho_r z0), k^2 v (1 - rho_r^2)), and
##   P(tau_i <= t) = E[ Phi((log(t + Y0) - mu_{1|0}) / sigma_{1|0}) 1{t + Y0 > 0} ] (one-dimensional integral).
set.seed(seed$Q1r)
Q1r <- gen_Q1r(P1, Q1)
q1r_of <- function(t) exp(P1$logc + P1$k * (P1$m + sqrt(P1$v) * qnorm(t)))  # the marginal of Y(1) under Q1r
sig10 <- P1$k * sqrt(P1$v) * sqrt(1 - P1$rho_r^2)
G1r <- function(t) {
  lo <- if (t < 0) (log(-t) - P1$m) / sqrt(P1$v) else -12
  integrate(function(z0) {
    y0 <- exp(P1$m + sqrt(P1$v) * z0)
    mu10 <- P1$logc + P1$k * (P1$m + sqrt(P1$v) * P1$rho_r * z0)
    dnorm(z0) * pnorm((log(pmax(t + y0, 1e-300)) - mu10) / sig10)
  }, max(lo, -12), 12, rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L)$value
}
ind_Q1r <- function(t) qtile_root(G1r, t, -200, 400)
truth_Q1r <- list(taus = taus, qte = q1r_of(taus) - q0_of(taus), rho_r = P1$rho_r,
                  ind_effect_quantiles = vapply(taus, ind_Q1r, numeric(1)), share_lose = G1r(0),
                  fine = data.frame(tau = taus_fine, ind_Q1r = vapply(taus_fine, ind_Q1r, numeric(1))))
stopifnot(max(abs(truth_Q1r$qte - truth_Q1$qte)) < 1e-12)
tau_curve <- round(seq(0.02, 0.98, by = 0.01), 2)                      # curves for figure (a)
truth_Q1$curve <- data.frame(tau = tau_curve, qte = qte_Q1(tau_curve), ind_Q1 = vapply(tau_curve, ind_Q1, numeric(1)),
                             ind_Q1r = vapply(tau_curve, ind_Q1r, numeric(1)))

## Population check: one million merchants with both potential outcomes, under both couplings.
set.seed(seed$Q1_pop)
pop <- gen_Q1(modifyList(P1, list(n = P1$n_pop, n_treat = P1$n_pop / 2)))
popr <- gen_Q1r(P1, pop)
truth_Q1$pop <- list(n = nrow(pop), qte = emp_quantile(pop$y1_true, taus) - emp_quantile(pop$y0_true, taus),
                     ate = mean(pop$y1_true) - mean(pop$y0_true),
                     ind_effect_quantiles = emp_quantile(pop$y1_true - pop$y0_true, taus))
truth_Q1r$pop <- list(n = nrow(popr), qte = emp_quantile(popr$y1_true, taus) - emp_quantile(popr$y0_true, taus),
                      ind_effect_quantiles = emp_quantile(popr$y1_true - popr$y0_true, taus),
                      share_lose = mean(popr$y1_true < popr$y0_true))
rm(pop, popr)

## =================================================================
## Q2: Fulfilled by Marigold (log monthly revenue in USD)
## =================================================================
## 6,002 merchants; exactly 1,201 (p = 0.2) use FBM (W = 1; an observational covariate, independent
## of the catalogue size; the fixed count keeps both groups odd so every QR solution is unique).
## X ~ N(0, 1): standardised log catalogue size.
##   log Y = mu + g W + bx X + s_W e,  mu = log(8000), g = 0.5, bx = 0.5, s0 = s1 = 0.35.
## Conditional quantiles: Quant_tau(log Y | W, X) = mu + g W + bx X + s_W z_tau, linear, so the QR
## coefficient on W is g + (s1 - s0) z_tau = 0.5 at every tau (a pure location shift).
## Unconditional: log Y | W = w ~ N(mu + g w, bx^2 + s_w^2); with the share p of FBM merchants,
##   F(q; p) = p F1(q) + (1 - p) F0(q). Raising p (holding F(Y | W, X) and F(X) fixed) moves the
##   tau-quantile by UQPE(tau) = dq/dp = (F0(q) - F1(q)) / f(q) at q = F^{-1}(tau; p) (FFL 2009,
##   Corollary 3), the target of the RIF regression coefficient on W (X independent of W, so the
##   population OLS coefficient on W equals E[RIF | W = 1] - E[RIF | W = 0]).
## truth13$Q2$uqe is UQPE at tau_uqe = 0.9 (about 0.665, against the conditional 0.5).
## The discrete counterfactual: the share rises by dp = 0.1 (to 0.3); cf = quantiles at p + dp.
P2 <- list(n = 6002, n_fbm = 1201, p = 0.2, mu = log(8000), g = 0.5, bx = 0.5, s0 = 0.35, s1 = 0.35,
           tau_uqe = 0.9, dp = 0.1)
set.seed(seed$Q2)
Q2 <- gen_Q2(P2)
mixq <- function(t, p) qtile_root(function(q) mix_cdf_Q2(q, P2, p), t, P2$mu - 10, P2$mu + 10)
uqpe_Q2 <- function(t) {
  q <- mixq(t, P2$p)
  (pnorm(q, P2$mu, sqrt(P2$bx^2 + P2$s0^2)) - pnorm(q, P2$mu + P2$g, sqrt(P2$bx^2 + P2$s1^2))) / mix_pdf_Q2(q, P2, P2$p)
}
truth_Q2 <- list(taus = taus, q = vapply(taus, mixq, numeric(1), p = P2$p),
                 uqpe = vapply(taus, uqpe_Q2, numeric(1)), cqr = P2$g + (P2$s1 - P2$s0) * qnorm(taus),
                 marg_diff = P2$g + (sqrt(P2$bx^2 + P2$s1^2) - sqrt(P2$bx^2 + P2$s0^2)) * qnorm(taus),
                 cf = vapply(taus, mixq, numeric(1), p = P2$p + P2$dp),
                 fine = data.frame(tau = taus_fine, uqpe = vapply(taus_fine, uqpe_Q2, numeric(1)),
                                   cqr = P2$g + (P2$s1 - P2$s0) * qnorm(taus_fine)))
truth_Q2$uqe <- truth_Q2$uqpe[taus == P2$tau_uqe]
truth_Q2$uqe_discrete <- (truth_Q2$cf - truth_Q2$q) / P2$dp
truth_Q2$pct_cf <- exp(truth_Q2$cf - truth_Q2$q) - 1          # proportional change of the revenue quantiles

## =================================================================
## Q3: AdBoost adoption with a randomised ad-credit offer (log monthly revenue in USD)
## =================================================================
## 6,002 merchants; exactly 3,001 get a USD 200 ad-credit offer (Z = 1; both arms odd). Adoption
## D(z) = 1{Phi(Z_V) > 1 - pi_z}, pi0 = 0.3, pi1 = 0.7: 30% always-takers, 40% compliers (adopt only
## with the credit), 30% never-takers; monotone (no defiers). The revenue rank Z_U and the adoption
## index Z_V are standard normal with correlation rho = 0.5: merchants with better unobserved
## prospects adopt, so D is endogenous. X ~ N(0, 1) (standardised log followers), independent.
##   log Y(d) = mu0 + d dmu + s_d Z_U + bx X,  mu0 = log(8000), dmu = 0.10, s0 = 0.5, s1 = 0.6, bx = 0.3.
## Chernozhukov-Hansen structure (rank invariance: the same Z_U in both potential outcomes, U = Phi(Z_U)
## independent of Z): the structural quantile effect is beta_D(tau) = dmu + (s1 - s0) z_tau.
## The complier band (1 - pi1, 1 - pi0] = (0.3, 0.7] is symmetric around 0.5, so Z_U | complier is
## symmetric around 0; with X symmetric and independent, every complier distribution of log Y(d) is
## symmetric around mu0 + d dmu and LQTE(0.5) = dmu = beta_D(0.5) exactly (lqte_median). Away from the
## median the complier QTE and beta_D(tau) differ (compliers' ranks are concentrated in the middle).
## Complier cdf: P(log Y(d) <= y | C) = int_a^b Phi((y - mu_d - s_d rho v) / sigma_d) phi(v) dv / (Phi(b) - Phi(a)),
##   a = qnorm(1 - pi1), b = qnorm(1 - pi0), sigma_d^2 = s_d^2 (1 - rho^2) + bx^2.
P3 <- list(n = 6002, n_z1 = 3001, mu0 = log(8000), dmu = 0.10, s0 = 0.5, s1 = 0.6, bx = 0.3, rho = 0.5,
           pi0 = 0.3, pi1 = 0.7, tau = 0.5, n_pop = 1e6, grid_step = 0.002)
set.seed(seed$Q3)
Q3 <- gen_Q3(P3)
aC <- qnorm(1 - P3$pi1); bC <- qnorm(1 - P3$pi0)
FC <- function(y, d) {
  mu_d <- P3$mu0 + d * P3$dmu; s_d <- if (d == 1) P3$s1 else P3$s0
  sig_d <- sqrt(s_d^2 * (1 - P3$rho^2) + P3$bx^2)
  integrate(function(v) pnorm((y - mu_d - s_d * P3$rho * v) / sig_d) * dnorm(v), aC, bC, rel.tol = 1e-13)$value /
    (pnorm(bC) - pnorm(aC))
}
lqte_Q3 <- function(t) qtile_root(function(y) FC(y, 1), t, P3$mu0 - 5, P3$mu0 + 5) -
  qtile_root(function(y) FC(y, 0), t, P3$mu0 - 5, P3$mu0 + 5)
truth_Q3 <- list(taus = taus, lqte = vapply(taus, lqte_Q3, numeric(1)), struct = P3$dmu + (P3$s1 - P3$s0) * qnorm(taus),
                 share_always = P3$pi0, share_complier = P3$pi1 - P3$pi0, share_never = 1 - P3$pi1,
                 late_log = P3$dmu,                                    # E[log Y(1) - log Y(0) | C] = dmu (symmetry)
                 fine = data.frame(tau = taus_fine, lqte = vapply(taus_fine, lqte_Q3, numeric(1)),
                                   struct = P3$dmu + (P3$s1 - P3$s0) * qnorm(taus_fine)))
truth_Q3$lqte_median <- truth_Q3$lqte[taus == 0.5]
stopifnot(abs(truth_Q3$lqte_median - P3$dmu) < 1e-8)
set.seed(seed$Q3_pop)
pop <- gen_Q3(modifyList(P3, list(n = P3$n_pop, n_z1 = P3$n_pop / 2)))
cp <- pop$type_true == "complier"
truth_Q3$pop <- list(n = nrow(pop), n_complier = sum(cp),
                     lqte = emp_quantile(pop$logy1_true[cp], taus) - emp_quantile(pop$logy0_true[cp], taus),
                     late_log = mean(pop$logy1_true[cp] - pop$logy0_true[cp]))
rm(pop, cp)

## =================================================================
## Q4: merchant panel with fixed effects (log monthly revenue in USD)
## =================================================================
## N = 200 merchants observed for \mathcal{T} = 5 or 20 months. alpha_i ~ N(0, 0.6^2) (merchant size);
## X_it = plogis(0.5 alpha_i / 0.6 + nu_it) in (0, 1), the share of the catalogue on promotion (larger
## merchants promote more: X is correlated with alpha_i, so pooled QR without the effects is biased).
##   Y_it = mu + alpha_i + beta X_it + sig (1 + gam X_it) e_it,  mu = log(8000), beta = 0.4, sig = 0.3, gam = 1.
## beta(tau) = beta + sig gam z_tau: promotions move the top of the conditional distribution more
## (0.016 at tau = 0.1, 0.784 at tau = 0.9). truth13$Q4$beta_tau is beta(0.9).
## The Monte Carlo (05_panel_mc.R, R = 1000 replications, seed 134000 + r) compares QR with merchant
## dummies and Canay's two-step estimator at \mathcal{T} = 5 and 20. Ruling 17: this is a finite-sample
## illustration, not a check of an asymptotic theorem: Kato, Galvao and Montes-Rojas (2012) need
## n^2 (log n)^3 / \mathcal{T} -> 0 for asymptotic normality of the dummy-variable estimator and Canay's (2011)
## rate condition is disputed (Besstremyannaya and Golovan 2019).
P4 <- list(N = 200, T_short = 5, T_long = 20, mu = log(8000), sa = 0.6, cx = 0.5, beta = 0.4, sig = 0.3, gam = 1,
           tau = 0.9, R = 1000)
set.seed(seed$Q4_T5)
Q4_T5 <- gen_Q4(P4, P4$T_short)
set.seed(seed$Q4_T20)
Q4_T20 <- gen_Q4(P4, P4$T_long)
truth_Q4 <- list(taus = taus, beta_tau = P4$beta + P4$sig * P4$gam * qnorm(P4$tau),
                 beta_grid = P4$beta + P4$sig * P4$gam * qnorm(taus), beta_mean = P4$beta)

## =================================================================
## Q5: sharp changes-in-changes, One-Click Checkout (revenue in USD thousands)
## =================================================================
## Two regions (S = 1: launch region), two quarters (T = 0: Q3, T = 1: the holiday quarter), 2,500
## merchants per cell, repeated cross sections. Rank z: N(0, 1) in the control region, N(0.6, 0.8^2) in
## the launch region (larger, more homogeneous merchants), the same distribution in both quarters.
##   Y(0) = h(z, t) = A exp(g t + s (1 + kap t) z),  A = 10, g = 0.05, s = 0.5, kap = 0.4;
##   Y(1) = (1 + delta) Y(0), delta = 0.12 (One-Click Checkout raises revenue by 12%); D = S T.
## Athey and Imbens (2006) Assumptions 3.1 to 3.4 hold (h strictly increasing in z; z | S time-invariant;
## full support), so CIC identifies QTT(tau) = delta h(mu1 + sd1 z_tau, 1). The time effect
## y -> h(h^{-1}(y, 0), 1) is a power function, and the launch region's rank distribution differs, so the
## quantile DiD (QDID) confounds the region difference with the nonlinear time effect.
P5 <- list(nc = 2500, A = 10, g = 0.05, s = 0.5, kap = 0.4, mu1 = 0.6, sd1 = 0.8, delta = 0.12, nc_pop = 250000)
set.seed(seed$Q5)
Q5 <- gen_Q5(P5)
z1_of <- function(t) P5$mu1 + P5$sd1 * qnorm(t)
qtt_Q5 <- function(t) P5$delta * h_Q5(z1_of(t), 1, P5)
qdid_Q5 <- function(t) (1 + P5$delta) * h_Q5(z1_of(t), 1, P5) - h_Q5(z1_of(t), 0, P5) - (h_Q5(qnorm(t), 1, P5) - h_Q5(qnorm(t), 0, P5))
truth_Q5 <- list(taus = taus, qtt = qtt_Q5(taus), qdid_pop = qdid_Q5(taus),
                 q11 = (1 + P5$delta) * h_Q5(z1_of(taus), 1, P5), cf = h_Q5(z1_of(taus), 1, P5),
                 att = P5$delta * P5$A * exp(P5$g + P5$s * (1 + P5$kap) * P5$mu1 + (P5$s * (1 + P5$kap) * P5$sd1)^2 / 2),
                 fine = data.frame(tau = taus_fine, qtt = qtt_Q5(taus_fine), qdid = qdid_Q5(taus_fine)))
set.seed(seed$Q5_pop)
pop <- gen_Q5(modifyList(P5, list(nc = P5$nc_pop)))
c11 <- pop$S == 1 & pop$T == 1
truth_Q5$pop <- list(n_cell = P5$nc_pop, qtt = emp_quantile(pop$y1_true[c11], taus) - emp_quantile(pop$y0_true[c11], taus))
rm(pop, c11)

## =================================================================
## Q5f: fuzzy changes-in-changes, voluntary One-Click Checkout (revenue in USD thousands)
## =================================================================
## The same regions, quarters, h and delta as Q5, 10,000 merchants per cell. Adoption is voluntary:
## in the launch region 10% adopted in a beta before the launch (always-takers), 60% adopt in the
## holiday quarter (switchers), 30% never adopt; in the control region a stable 15% use an older
## version (always-takers). Adopters have higher ranks: z = mu_S + sd_S (rho Z_V + sqrt(1 - rho^2) eta),
## rho = 0.4, type from Phi(Z_V). dCDH (2018) Assumptions 1, 2, 6 (stable control-group rate), 7
## (Y(d) = h_d(U_d, T) with h_0 = h, h_1 = (1 + delta) h, ranks time-invariant within group and type)
## and 8 (support) hold, so Theorem 2.3 identifies the switchers' period-1 distributions of Y(0) and
## Y(1) and LQTE_sw(tau) = delta h(mu1 + sd1 w_tau, 1), with w_tau the tau-quantile of
## rho Z_V + sqrt(1 - rho^2) eta given Z_V in (qnorm(0.3), qnorm(0.9)]. truth13$Q5$lqte is LQTE_sw(0.5);
## late_sw = delta E[h(z, 1) | switcher] is the switchers' LATE (chapter 20's W_CIC target).
P5f <- modifyList(P5, list(nc = 10000, rho = 0.4, p_at1 = 0.10, p_sw1 = 0.60, p_at0 = 0.15, tau = 0.5, nc_pop = 250000))
set.seed(seed$Q5f)
Q5f <- gen_Q5f(P5f)
aS <- qnorm(1 - P5f$p_at1 - P5f$p_sw1); bS <- qnorm(1 - P5f$p_at1)
FzS <- function(w) integrate(function(v) pnorm((w - P5f$rho * v) / sqrt(1 - P5f$rho^2)) * dnorm(v), aS, bS,
                             rel.tol = 1e-13)$value / (pnorm(bS) - pnorm(aS))
lqte_sw <- function(t) P5f$delta * h_Q5(P5f$mu1 + P5f$sd1 * qtile_root(FzS, t, -8, 8), 1, P5f)
lam <- P5f$s * (1 + P5f$kap) * P5f$sd1                                 # E[exp(lam w) | band], w = rho Z_V + sqrt(1 - rho^2) eta
Eexp <- exp(lam^2 / 2) * (pnorm(bS - lam * P5f$rho) - pnorm(aS - lam * P5f$rho)) / (pnorm(bS) - pnorm(aS))
truth_Q5$lqte_grid <- vapply(taus, lqte_sw, numeric(1))
truth_Q5$lqte <- truth_Q5$lqte_grid[taus == P5f$tau]
truth_Q5$late_sw <- P5f$delta * P5f$A * exp(P5f$g + P5f$s * (1 + P5f$kap) * P5f$mu1) * Eexp
truth_Q5$shares_f <- list(always1 = P5f$p_at1, switcher1 = P5f$p_sw1, never1 = 1 - P5f$p_at1 - P5f$p_sw1, always0 = P5f$p_at0)
set.seed(seed$Q5f_pop)
pop <- gen_Q5f(modifyList(P5f, list(nc = P5f$nc_pop)))
sw <- pop$type_true == "switcher" & pop$T == 1
y0s <- h_Q5(pop$z_true[sw], 1, P5f); y1s <- (1 + P5f$delta) * y0s
truth_Q5$pop$lqte_grid <- emp_quantile(y1s, taus) - emp_quantile(y0s, taus)
truth_Q5$pop$late_sw <- mean(y1s - y0s)
truth_Q5$pop$n_switcher <- sum(sw)
rm(pop, sw, y0s, y1s)

## =================================================================
## Save
## =================================================================
saveRDS(Q1, "marigold_Q1.rds"); saveRDS(Q1r, "marigold_Q1r.rds"); saveRDS(Q2, "marigold_Q2.rds")
saveRDS(Q3, "marigold_Q3.rds"); saveRDS(list(T5 = Q4_T5, T20 = Q4_T20), "marigold_Q4.rds")
saveRDS(Q5, "marigold_Q5.rds"); saveRDS(Q5f, "marigold_Q5f.rds")
truth13 <- list(Q1 = truth_Q1, Q1r = truth_Q1r, Q2 = truth_Q2, Q3 = truth_Q3, Q4 = truth_Q4, Q5 = truth_Q5,
                params = list(taus = taus, taus_fine = taus_fine, seed = seed,
                              Q1 = P1, Q2 = P2, Q3 = P3, Q4 = P4, Q5 = P5, Q5f = P5f))
saveRDS(truth13, "truth13.rds")
cat("Saved marigold_Q1 ... Q5f and truth13.rds\n")
cat(sprintf("Q1 QTE: %s\nQ1r individual-effect quantiles: %s\nQ2 UQPE: %s (uqe %.4f)\nQ3 LQTE: %s\nQ4 beta(0.9) = %.4f\nQ5 QTT: %s; LQTE_sw(0.5) = %.4f\n",
            paste(round(truth_Q1$qte, 3), collapse = " "), paste(round(truth_Q1r$ind_effect_quantiles, 3), collapse = " "),
            paste(round(truth_Q2$uqpe, 3), collapse = " "), truth_Q2$uqe, paste(round(truth_Q3$lqte, 4), collapse = " "),
            truth_Q4$beta_tau, paste(round(truth_Q5$qtt, 3), collapse = " "), truth_Q5$lqte))
