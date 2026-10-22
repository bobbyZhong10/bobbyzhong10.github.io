# ============================================================================
# 01_dgp.R  --  Helix developer-platform adoption + usage DGP
# ----------------------------------------------------------------------------
# "Helix" is a fictional B2B developer platform. Each client firm i has an
# observed standardized "engagement" covariate x_i (past API activity, z-scored).
# Two outcomes are generated, both driven by x through a nonlinear link:
#
#   (A) ADOPTION  D_i in {0,1}: whether firm i adopts a paid "Pro" tier.
#       True model is LOGIT:  P(D=1|x) = Lambda(a0 + a1 x),  Lambda = plogis.
#       => the object of MLE. A naive linear-probability OLS gives nonsense
#          (predicted probabilities outside [0,1]); logit MLE fixes it and
#          lets us exhibit Jensen consistency + Information-Matrix equality.
#
#   (B) USAGE  y_i in {0,1,2,...}: monthly count of Pro API calls (in 1000s).
#       True conditional MEAN is  E[y|x] = exp(b0 + b1 x)  (log link), BUT the
#       true distribution is NEGATIVE BINOMIAL (overdispersed), NOT Poisson.
#       => a data scientist who fits POISSON MLE is running QMLE: the mean is
#          correctly specified so beta is CONSISTENT, but the Information-Matrix
#          equality FAILS (Var != mean), so naive Fisher-information SEs are
#          wrong and the Huber-White SANDWICH is needed. The same conditional-
#          mean restriction  E[(y - exp(b0+b1 x)) h(x)] = 0  becomes a GMM
#          system: h=(1,x) just-identifies (= Poisson score), h=(1,x,x^2)
#          over-identifies -> efficient weighting + Hansen J test of the mean.
#
# All TRUE parameters are known; oracle targets are computed exactly / by a
# giant reference sample. Seed fixed. Scripts read/write RELATIVE paths, so
# run with the working directory set to this folder.
# ============================================================================

set.seed(512)

## ---- true parameters -------------------------------------------------------
a0 <- -0.40; a1 <- 1.10          # logit adoption:  P(D=1|x) = plogis(a0 + a1 x)
b0 <-  1.00; b1 <- 0.50          # log-mean usage:  E[y|x]  = exp(b0 + b1 x)
theta_nb <- 2.0                  # NB size (dispersion): Var(y|x) = mu + mu^2/theta

n_main <- 4000                   # main analysis sample

gen_helix <- function(n) {
  x  <- rnorm(n)                                   # standardized engagement
  ## (A) adoption -- logit
  pD <- plogis(a0 + a1 * x)
  D  <- rbinom(n, 1, pD)
  ## (B) usage -- negative binomial with log mean (overdispersed)
  mu <- exp(b0 + b1 * x)
  y  <- rnbinom(n, size = theta_nb, mu = mu)       # Var = mu + mu^2/theta_nb
  data.frame(x = x, D = D, pD = pD, y = y, mu = mu)
}

dat <- gen_helix(n_main)

## ---- oracle truths (from a 5,000,000-row reference draw) -------------------
set.seed(9001)
big <- gen_helix(5e6)

# average marginal effect (AME) of x on adoption probability, logit truth:
#   d/dx P = a1 * lambda(a0 + a1 x),  lambda = dLogistic
ame_logit <- mean(a1 * dlogis(a0 + a1 * big$x))
# marginal effect at the mean (MEM), x = 0:
mem_logit <- a1 * dlogis(a0 + a1 * 0)
base_rate <- mean(big$D)

oracle <- list(
  a0 = a0, a1 = a1, b0 = b0, b1 = b1, theta_nb = theta_nb,
  ame_logit = ame_logit, mem_logit = mem_logit, base_rate = base_rate
)

saveRDS(dat,    "helix_data.rds")
saveRDS(oracle, "oracle.rds")

cat(sprintf("n_main            = %d\n", n_main))
cat(sprintf("adoption rate     = %.4f (oracle base rate %.4f)\n", mean(dat$D), base_rate))
cat(sprintf("mean(y)           = %.4f\n", mean(dat$y)))
cat(sprintf("var(y)/mean(y)    = %.4f  (>1 => overdispersion vs Poisson)\n",
            var(dat$y) / mean(dat$y)))
cat(sprintf("oracle AME (logit)= %.4f\n", ame_logit))
cat(sprintf("oracle MEM (logit)= %.4f\n", mem_logit))

# ============================================================================
# Lines 1 to 77 above are the inherited legacy block B1, byte-identical (Phase 3: inherited blocks
# first and unchanged; tests/test_dgp.R checks the md5 of those lines). Chapter 12 uses the adoption
# outcome D (Y_i in the text) and x (X_i in the text); the usage count y of the legacy block is left in
# helix_data.rds and is not used in this chapter. Everything below is new in v3.
# ----------------------------------------------------------------------------
# 01_dgp.R (v3) -- the Helix simulation world for chapter 12 (binary and multinomial response).
#
# Helix is a (fictional) B2B developer platform; unit i is a client firm. Money is in US dollars.
# One world, seven blocks:
#   B1   inherited (above): adoption of the paid Pro tier, logit, n = 4,000 (LPM 121%, logit MLE, AME,
#        consistency simulation in 03_mle_logit.R). Truth extended below by quadrature.    helix_data.rds
#   B2   two features (API engagement x1, Helix Deploy usage x2) and their interaction, probit:
#        beta_jj' > 0 but the average cross-partial effect is negative (Ai and Norton 2003).  helix_B2.rds
#   B3   endogenous trial engagement D_i, randomised onboarding reminder Z_i as instrument,
#        cloud spend (USD) as exogenous regressor; bivariate normal errors (Rivers and Vuong 1988). helix_B3.rds
#   B4   an omitted covariate (budget flexibility a_i) independent of X_i: the probit coefficient
#        shrinks by 1 / sqrt(1 + sigma_a^2), the APE does not (Wooldridge 2010 eqs. 15.35 to 15.38). helix_B4.rds
#   B5   monthly adoption panel, N = 1,000 firms, T = 8 months, firm effects alpha_i correlated with
#        X_it (Mundlak form), FE logit correct.                                              helix_B5.rds
#        Two designs at T = 4 (controller ruling 42), kept apart everywhere:
#          T4        a 4-month DGP: alpha_i built on the 4-month mean of X_it (the Monte Carlo of 06_panel.R
#                    with T = 4, truths dgp_B5_T4_*, dgp_B5_ape_T4);
#          T8first4  months 1 to 4 of the T = 8 panel: alpha_i built on the 8-month mean (the one-sample
#                    "T = 4" panel, panel_T4_* keys; truths dgp_B5_T8first4_*).
#   B6   1 to 5 satisfaction rating, ordered probit with a heteroskedastic latent error for
#        managed accounts (Bond and Lang 2019 warning).                                       helix_B6.rds
#   B7   plan choice Free / Pro / Enterprise: MNL (IIA holds) and a nested-logit variant with
#        the paid nest {Pro, Enterprise} (IIA fails).                                         helix_B7.rds, helix_B7n.rds
#   truth12.rds  list(B1, ..., B7, params): truth (closed form or quadrature), large-sample checks,
#                parameters and seeds.
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 1.10): Y_i (code: D in B1, y or Y elsewhere), X_i, link F with density F';
# beta_jj' the interaction coefficient (B2); D_i, beta_D, Z_i, pi = (pi_1, pi_2), v_i, rho (B3);
# a_i, sigma_a (B4); alpha_i = xi_0 + xi Xbar_i + a_i, T (B5); mu_1 < ... < mu_4 the cut points (B6);
# beta_(m), labels beta_Pro and beta_Ent with Free the base, lambda the nest parameter (B7).
#
# Seeds (every block its own set.seed; B1 keeps the legacy 512, 9001 and 778): B2 1202, B3 1203,
# B4 1204, B5 1205, B6 1206, B7 1207, B7 nested variant 1217; large-sample checks of the truth
# B2 1222, B3 1223, B4 1224, B5 1225 (T = 8) and 1245 (T = 4), B6 1226, B7 1227; bootstrap of B3 1233;
# the large-sample check of the T8first4 truth of B5 (months 1 to 4 of a T = 8 panel) has its own seed 1235,
# added in fix round 1 (no existing seed changed);
# Monte Carlo seeds used by 06 and 07 are stored here too (B5: 1254000 + r at T = 4 and 1258000 + r at
# T = 8; B7: 1270000 + r). They were fixed before any result was seen. B6_mc (1260000 + r) was added
# afterwards under ruling 29: the single-sample B6 interval misses its truth (an ordinary sampling miss,
# kept as it is), and the chapter card's assertion is replaced by its Monte Carlo counterpart.
# Run with the working directory set to this folder (00_run_all.R does this), like the legacy block.
# ============================================================================
rm(big)
source("bin_fns.R")

seed <- list(B1 = 512, B1_oracle = 9001, B1_consistency = 778, B2 = 1202, B3 = 1203, B4 = 1204, B5 = 1205,
             B6 = 1206, B7 = 1207, B7_nested = 1217, B2_large = 1222, B3_large = 1223, B4_large = 1224,
             B5_large = 1225, B5_large_T4 = 1245, B5_large_first4 = 1235, B6_large = 1226, B7_large = 1227, B3_boot = 1233,
             B5_mc_T4 = 1254000, B5_mc_T8 = 1258000, B7_mc = 1270000, B6_mc = 1260000)
ng <- normal_grid()                                       # quadrature nodes for N(0, 1)

## =================================================================
## B1: truth of the inherited block by quadrature (no random numbers)
## =================================================================
## P(Y = 1 | X) = Lambda(beta0 + beta1 X), beta = (a0, a1) = (-0.40, 1.10), X ~ N(0, 1). The legacy
## oracle (5,000,000 draws, seed 9001) is kept as B1$ame and B1$base_rate; the exact values are added.
## Because X is normal, the population linear projection slope of Y on X is Cov(X, Y) = E[X Lambda(.)]
## = E[Lambda'(.)] beta1 (Stein's lemma): the LPM slope estimates the APE exactly in this design.
## The probit pseudo-true value maximises the expected probit log likelihood under the logit truth.
P1 <- list(n = n_main, beta = c(a0, a1))
EX <- function(f) sum(ng$w * f(ng$z))
ame_exact <- EX(function(x) a1 * dlogis(a0 + a1 * x))
base_exact <- EX(function(x) plogis(a0 + a1 * x))
proj_slope <- EX(function(x) x * plogis(a0 + a1 * x))
qprob <- function(b) -EX(function(x) { pl <- plogis(a0 + a1 * x); u <- b[1] + b[2] * x
  pl * pnorm(u, log.p = TRUE) + (1 - pl) * pnorm(-u, log.p = TRUE) })
pq <- optim(c(-0.25, 0.65), qprob, method = "BFGS", control = list(reltol = 1e-15, maxit = 1000))$par
truth_B1 <- list(beta = c(a0, a1), ame = ame_logit, mem = mem_logit, base_rate = base_rate,
                 ame_exact = ame_exact, base_rate_exact = base_exact,
                 lpm_slope_pop = proj_slope, lpm_intercept_pop = base_exact,
                 share_above_pop = pnorm(-(1 - base_exact) / proj_slope), share_below_pop = pnorm(-base_exact / proj_slope),
                 probit_pseudo = pq, probit_ape_pseudo = EX(function(x) pq[2] * dnorm(pq[1] + pq[2] * x)),
                 logit_over_probit_pseudo = a1 / pq[2])

## =================================================================
## B2: two features and their interaction (probit)
## =================================================================
## 4,000 firms in the Helix Deploy beta. x1 = API engagement, x2 = Deploy pipeline usage (both z-scores,
## Corr 0.3). P(Y = 1 | x) = Phi(1.0 + 1.0 x1 + 0.7 x2 + 0.2 x1 x2): the interaction coefficient
## beta_jj' = 0.2 > 0 (the features are complements in the latent index), but about 72% of these firms
## adopt and most sit where Phi is concave, so the cross-partial
##   d^2 P / dx1 dx2 = beta_jj' phi(u) + (beta1 + beta_jj' x2)(beta2 + beta_jj' x1) phi'(u),  phi'(u) = -u phi(u),
## is negative on average (truth -0.019). Truth by two-dimensional quadrature; large-sample check below.
P2 <- list(n = 4000, rho_x = 0.3, beta = c(1.0, 1.0, 0.7, 0.2))
set.seed(seed$B2)
B2 <- gen_B2(P2)
q2 <- expand.grid(i = seq_along(ng$z), j = seq_along(ng$z))
x1q <- ng$z[q2$i]; x2q <- P2$rho_x * x1q + sqrt(1 - P2$rho_x^2) * ng$z[q2$j]; wq <- ng$w[q2$i] * ng$w[q2$j]
Xq2 <- cbind(1, x1q, x2q, x1q * x2q)
uq <- drop(Xq2 %*% P2$beta); cq <- cross_partial_i(P2$beta, Xq2, 2, 3, "probit")
truth_B2 <- list(beta = setNames(P2$beta, c("(Intercept)", "x1", "x2", "x1:x2")),
                 avg_cross = sum(wq * cq), share_pos = sum(wq * (cq > 0)), rate = sum(wq * pnorm(uq)),
                 ape1 = sum(wq * dnorm(uq) * (P2$beta[2] + P2$beta[4] * x2q)),
                 ape2 = sum(wq * dnorm(uq) * (P2$beta[3] + P2$beta[4] * x1q)),
                 naive_int_me = P2$beta[4] * sum(wq * dnorm(uq)))
rm(q2, x1q, x2q, wq, Xq2, uq, cq)
set.seed(seed$B2_large)
L2 <- gen_B2(P2, n = 2e6)
XL2 <- cbind(1, L2$x1, L2$x2, L2$x1 * L2$x2)
cL2 <- cross_partial_i(P2$beta, XL2, 2, 3, "probit")
truth_B2$large <- list(N = nrow(L2), avg_cross = mean(cL2), share_pos = mean(cL2 > 0), rate = mean(L2$y))
rm(L2, XL2, cL2)

## =================================================================
## B3: endogenous trial engagement, randomised onboarding reminder (control function)
## =================================================================
## 4,000 firms in a 30-day trial. X = log(monthly cloud spend / USD 5,000), sd 0.8 (median spend
## USD 5,000 a month). Helix sends a sequence of onboarding reminders to a random half (Z = 1).
## Trial engagement (z-score) D = -0.4 + 0.3 X + 0.8 Z + v, v ~ N(0, 1): the reminder raises engagement by
## 0.8 standard deviations. Adoption Y = 1{-0.6 + 0.3 X + 0.5 D + eps > 0}, eps = 0.5 v + sqrt(0.75) e:
## an urgent project raises both engagement and the value of Pro (Corr(eps, v) = 0.5, Var(eps) = 1), so a
## probit of Y on (X, D) overstates the engagement effect. The reminder moves adoption only through D.
## Truth: APE of D = E[beta_D phi(beta0 + betaX X + beta_D D)] over the population of (X, D), in closed form
## (given Z the index is normal); the control-function probit estimates (beta, beta_D, rho) / sigma_e with
## rho = 0.5 / sigma_v and sigma_e = sqrt(1 - 0.25).
P3 <- list(n = 4000, spend_median_kusd = 5, sd_logspend = 0.8, p_reminder = 0.5, pi = c(-0.4, 0.3, 0.8),
           sigma_v = 1, beta = c(-0.6, 0.3), beta_D = 0.5, r = 0.5)
set.seed(seed$B3)
B3 <- gen_B3(P3)
ape_B3 <- function(p) {
  c0 <- p$beta[1] + p$beta_D * p$pi[1]; c1 <- p$beta[2] + p$beta_D * p$pi[2]; c2 <- p$beta_D * p$pi[3]
  s2 <- c1^2 * p$sd_logspend^2 + p$beta_D^2 * p$sigma_v^2
  m <- c0 + c2 * c(0, 1); pz <- c(1 - p$p_reminder, p$p_reminder)
  # the ASF integrates eps over its marginal N(0, 1); the adoption rate uses Cov(eps, beta_D v) = beta_D r sigma_v
  list(ape = p$beta_D * sum(pz * dnorm(m / sqrt(1 + s2)) / sqrt(1 + s2)),
       rate = sum(pz * pnorm(m / sqrt(1 + s2 + 2 * p$beta_D * p$r * p$sigma_v))))
}
a3 <- ape_B3(P3)
sigma_e <- sqrt(1 - P3$r^2); rho <- P3$r / P3$sigma_v
truth_B3 <- list(beta = P3$beta, beta_D = P3$beta_D, pi = P3$pi, rho = rho, sigma_e = sigma_e,
                 scaled = c(P3$beta, P3$beta_D, rho) / sigma_e, ape = a3$ape, rate = a3$rate,
                 partial_r2_Z = P3$pi[3]^2 * P3$p_reminder * (1 - P3$p_reminder) /
                   (P3$pi[3]^2 * P3$p_reminder * (1 - P3$p_reminder) + P3$sigma_v^2))
set.seed(seed$B3_large)
L3 <- gen_B3(P3, n = 2e6)
gL3 <- glm.fit(cbind(1, L3$X, L3$D), L3$Y, family = binomial("probit"))
bL3 <- gL3$coefficients
truth_B3$large <- list(N = nrow(L3), ape = mean(P3$beta_D * dnorm(P3$beta[1] + P3$beta[2] * L3$X + P3$beta_D * L3$D)),
                       rate = mean(L3$Y), naive_coef = unname(bL3),
                       naive_ape = mean(dnorm(drop(cbind(1, L3$X, L3$D) %*% bL3))) * unname(bL3[3]))
rm(L3, gL3)

## =================================================================
## B4: omitted heterogeneity independent of X (attenuation; the APE survives)
## =================================================================
## 4,000 firms. Y = 1{-0.2 + 0.9 X + a + eps > 0}: a = budget flexibility, N(0, 1), unobserved and independent
## of engagement X; eps ~ N(0, 1). A probit of Y on X estimates beta / sqrt(1 + sigma_a^2) = 0.9 / sqrt(2)
## = 0.636 (eq. 15.35), and its APE is the APE of the average structural function (eqs. 15.37, 15.38):
##   APE = E[(beta1 / s) phi((beta0 + beta1 X) / s)],  s = sqrt(1 + sigma_a^2),
## the same number as E[beta1 phi(beta0 + beta1 X + a)] (the derivative averaged over X and a).
P4 <- list(n = 4000, beta = c(-0.2, 0.9), sigma_a = 1)
set.seed(seed$B4)
B4 <- gen_B4(P4)
s4 <- sqrt(1 + P4$sigma_a^2); k4 <- sqrt(1 + (P4$beta[2] / s4)^2)
truth_B4 <- list(beta = P4$beta[2], beta_vec = P4$beta, sigma_a = P4$sigma_a, scale = s4,
                 beta_att = P4$beta[2] / s4, beta_att_vec = P4$beta / s4,
                 ape = (P4$beta[2] / s4) * dnorm(P4$beta[1] / s4 / k4) / k4,
                 rate = pnorm(P4$beta[1] / sqrt(1 + P4$sigma_a^2 + P4$beta[2]^2)))
truth_B4$lpm_slope_pop <- truth_B4$ape                     # X normal: Stein's lemma, as in B1
set.seed(seed$B4_large)
L4 <- gen_B4(P4, n = 2e6)
truth_B4$large <- list(N = nrow(L4), ape = mean(P4$beta[2] * dnorm(P4$beta[1] + P4$beta[2] * L4$X + L4$a_true)),
                       rate = mean(L4$Y),
                       short_coef = unname(glm.fit(cbind(1, L4$X), L4$Y, family = binomial("probit"))$coefficients))
rm(L4)

## =================================================================
## B5: monthly adoption panel with firm effects correlated with X
## =================================================================
## 1,000 firms over T = 8 months. Monthly engagement X_it = c_i + w_it, c_i ~ N(0, 0.7^2),
## w_it ~ N(0, 1). Firm effect alpha_i = -0.3 + 1.0 Xbar_i + a_i, a_i ~ N(0, 1): firms that engage more on
## average also like Pro more, whatever their engagement this month. Y_it = 1 (on Pro in month t) with
## probability Lambda(alpha_i + 0.8 X_it): FE logit (conditional MLE) is correctly specified for beta = 0.8.
## Pooled logit ignores alpha_i; its probability limit (below) is far from 0.8.
## Truth APE = E[beta Lambda'(alpha_i + beta X_it)] over all firm-months.
## EXACT truths (controller ruling 44): (X_it, alpha_i) is bivariate normal with the same law in every month t
## (Var X = sigma_c^2 + 1, Cov(X, alpha) = xi (sigma_c^2 + 1/T), Var alpha = xi^2 (sigma_c^2 + 1/T) + sigma_a^2), so the
## true APE, the rate, the pooled-logit limit and the pooled-APE limit are two-dimensional Gauss-Hermite integrals
## (b5_quadrature in bin_fns.R, 120 nodes per dimension, stable at 60). The large draws (N_large firms, own seeds)
## are kept as CHECKS only, under $draw, and tests/test_dgp.R asserts them within their Monte Carlo error (sd about
## 0.003 for the pooled limit, 0.00013 for the APE). The CRE-probit limits and the share of movers are not
## quadrature values: they stay large-draw values (they are not needed as exact truths).
## Three designs, because "T = 4" means two different things (controller ruling 42):
##   T8        the T = 8 panel (the one-sample panel, seed B5): pooled limit 0.97690, APE 0.124527, pooled-APE limit 0.18877;
##   T4        a 4-month DGP, gen_B5(T = 4): alpha_i = xi0 + xi Xbar_i^(4) + a_i with the 4-month mean (the Monte
##             Carlo at T = 4 draws from this DGP; seed B5_large_T4): pooled limit 1.04498, APE 0.121431;
##   T8first4  months 1 to 4 of the T = 8 panel: alpha_i keeps the 8-month mean Xbar_i^(8). Because the law of
##             (X_it, alpha_i) does not depend on t, its pooled-logit limit, true APE, pooled-APE limit and rate EQUAL
##             the T = 8 values exactly (0.97690, 0.124527, 0.18877); they differ from the T4 truths. The one-sample
##             "T = 4" panel (months 1 to 4 of helix_B5.rds) is compared with these. Its large draw (T = 8, keep t <= 4)
##             has its own seed B5_large_first4 and is a check, like the other two.
P5 <- list(N = 1000, T = 8, T_short = 4, sigma_c = 0.7, sigma_a = 1, xi0 = -0.3, xi = 1.0, beta = 0.8,
           R_mc = 1000, N_large = 2e5)
set.seed(seed$B5)
B5 <- gen_B5(P5)
truth_large_B5 <- function(T, first = NULL) {
  L <- gen_B5(P5, N = P5$N_large, T = T)
  if (!is.null(first)) L <- L[L$t <= first, ]                  # months 1 to `first` of the T-month panel (alpha_i built on all T months)
  Tn <- if (is.null(first)) T else first
  pooled <- glm.fit(cbind(1, L$X), L$Y, family = binomial("logit"))$coefficients
  xb <- ave(L$X, L$id)
  cre <- glm.fit(cbind(1, L$X, xb), L$Y, family = binomial("probit"))$coefficients
  list(N = P5$N_large, T = Tn, ape = mean(P5$beta * dlogis(L$alpha_true + P5$beta * L$X)), rate = mean(L$Y),
       pooled_plim = unname(pooled[2]),
       pooled_ape_plim = mean(dlogis(pooled[1] + pooled[2] * L$X)) * unname(pooled[2]),
       cre_probit_plim = unname(cre), cre_probit_ape_plim = mean(dnorm(cre[1] + cre[2] * L$X + cre[3] * xb)) * unname(cre[2]),
       share_movers = mean(tapply(L$Y, L$id, function(z) sum(z) > 0 && sum(z) < Tn)))
}
set.seed(seed$B5_large)
T8 <- truth_large_B5(P5$T)
set.seed(seed$B5_large_T4)
T4 <- truth_large_B5(P5$T_short)
set.seed(seed$B5_large_first4)
T8first4 <- truth_large_B5(P5$T, first = P5$T_short)
T8first4$T_dgp <- P5$T                                       # the DGP behind the four months
# the truths are the quadrature values; the large draws move to $draw (checks)
with_exact <- function(L, q) {
  L$draw <- L[c("ape", "rate", "pooled_plim", "pooled_ape_plim")]
  L[c("ape", "rate", "pooled_plim", "pooled_ape_plim")] <- q[c("ape", "rate", "pooled_plim", "pooled_ape_plim")]
  L$pooled_int <- q$pooled_int; L$quad_nodes <- 120L
  L
}
q8 <- b5_quadrature(P5, P5$T, 120L); q4 <- b5_quadrature(P5, P5$T_short, 120L)
T8 <- with_exact(T8, q8); T4 <- with_exact(T4, q4)
T8first4 <- with_exact(T8first4, q8)                         # months 1 to 4 of the T = 8 panel: the T = 8 values (same law in every month)
truth_B5 <- list(beta = P5$beta, ape = T8$ape, ape_T4 = T4$ape, T8 = T8, T4 = T4, T8first4 = T8first4)

## =================================================================
## B6: 1 to 5 satisfaction rating (ordered probit, heteroskedastic latent error)
## =================================================================
## 3,000 firms answer the quarterly survey. Y* = 0.5 x1 - 0.3 bill + 0.2 G + sigma_G eps: x1 = engagement,
## bill = monthly Helix bill in thousands of USD (lognormal, median USD 1,200, log sd 0.5),
## G = 1 for managed accounts (30%, an account manager), sigma_G = 1 (self-serve) or 1.5 (managed).
## Cut points mu = (-1.9, -1.1, -0.2, 0.8). The chapter's target is the APE of engagement on the share of
## 5-star ratings, E[(beta1 / sigma_G) phi((mu_4 - x'beta - delta G) / sigma_G)] (by quadrature: x1 integrates
## out in closed form, bill on a grid). Bond and Lang (2019): managed accounts have the higher latent mean
## (+0.2) and the larger spread; under the increasing transformation h(Y*) = -exp(-Y*) the ranking flips.
## The ratio is CONDITIONAL on the covariates: E[h | x1, bill, G = 1] / E[h | x1, bill, G = 0] =
## exp(-0.2 + (1.5^2 - 1) / 2) = exp(0.425) > 1 (h is negative, so managed accounts are lower on the h scale
## at every (x1, bill)). The unconditional E[h | G] is minus infinity in both groups: bill is lognormal and
## enters Y* with coefficient -0.3, so E[exp(0.3 bill)] diverges, and the unconditional ratio does not exist.
## The ratings cannot tell the two latent scales apart.
P6 <- list(n = 3000, beta = c(0.5, -0.3), delta_G = 0.2, sigma = c(1, 1.5), cuts = c(-1.9, -1.1, -0.2, 0.8),
           bill_median_kusd = 1.2, sd_logbill = 0.5, p_managed = 0.3, R_mc = 1000)
set.seed(seed$B6)
B6 <- gen_B6(P6)
bill_q <- exp(log(P6$bill_median_kusd) + P6$sd_logbill * ng$z)
cdf6 <- function(mu, g) {                                  # P(Y* <= mu | G = g), x1 and bill integrated out
  sd <- sqrt(P6$sigma[g + 1]^2 + P6$beta[1]^2)
  sum(ng$w * pnorm((mu - P6$beta[2] * bill_q - P6$delta_G * g) / sd))
}
shares_G <- t(sapply(0:1, function(g) diff(c(0, sapply(P6$cuts, cdf6, g = g), 1))))
dimnames(shares_G) <- list(c("self_serve", "managed"), paste0("r", 1:5))
pG <- c(1 - P6$p_managed, P6$p_managed)
ape_top6 <- sum(pG * sapply(0:1, function(g) {
  sd <- sqrt(P6$sigma[g + 1]^2 + P6$beta[1]^2)
  P6$beta[1] * sum(ng$w * dnorm((P6$cuts[4] - P6$beta[2] * bill_q - P6$delta_G * g) / sd)) / sd
}))
truth_B6 <- list(beta = P6$beta, delta_G = P6$delta_G, sigma = P6$sigma, cuts = P6$cuts, ape_top = ape_top6,
                 shares = drop(pG %*% shares_G), shares_G = shares_G, mean_rating_G = drop(shares_G %*% 1:5),
                 mean_bill_kusd = P6$bill_median_kusd * exp(P6$sd_logbill^2 / 2),
                 bl_latent_diff = P6$delta_G, bl_log_ratio_h = -P6$delta_G + (P6$sigma[2]^2 - P6$sigma[1]^2) / 2)
set.seed(seed$B6_large)
L6 <- gen_B6(P6, n = 2e6)
sg <- P6$sigma[L6$G + 1]
truth_B6$large <- list(N = nrow(L6), ape_top = mean(P6$beta[1] / sg * dnorm((P6$cuts[4] - P6$beta[1] * L6$x1 -
                         P6$beta[2] * L6$bill_kusd - P6$delta_G * L6$G) / sg)),
                       shares = as.vector(table(L6$rating)) / nrow(L6))
rm(L6, sg)

## =================================================================
## B7: plan choice Free / Pro / Enterprise (MNL and nested variant)
## =================================================================
## 5,000 firms choose a plan. x = engagement z-score, s = standardised log seats. Free is the base:
## V_Pro = 0.7 + 1.0 x + 0.3 s, V_Ent = 0.0 + 0.4 x + 1.1 s (beta_Pro, beta_Ent). Shares about 27% / 47% / 26%.
## IIA variant: iid extreme-value errors, MNL is correct. Nested variant (seed B7_nested): the same V
## with the paid nest {Pro, Enterprise}, lambda = 0.4 (Pro and Enterprise are close substitutes), drawn from
## the nested-logit probabilities. The Hausman and McFadden test drops Free: among paying firms, IIA says
## the Pro-versus-Enterprise odds keep the full-model coefficients beta_Ent - beta_Pro; in the nested
## variant they are (beta_Ent - beta_Pro) / lambda. MNL truth (shares, APE of x and s on each plan's
## probability) by two-dimensional quadrature.
P7 <- list(n = 5000, levels = c("Free", "Pro", "Enterprise"),
           coef = rbind(Pro = c(0.7, 1.0, 0.3), Enterprise = c(0.0, 0.4, 1.1)),
           lambda_nested = 0.4, drop = "Free", R_mc = 4000)
colnames(P7$coef) <- c("(Intercept)", "x", "s")
set.seed(seed$B7)
B7 <- gen_B7(P7)
set.seed(seed$B7_nested)
B7n <- gen_B7(P7, lambda = P7$lambda_nested)
q7 <- expand.grid(i = seq_along(ng$z), j = seq_along(ng$z))
Wq <- cbind(1, ng$z[q7$i], ng$z[q7$j]); wq7 <- ng$w[q7$i] * ng$w[q7$j]
Pq <- mnl_probs(P7$coef, Wq); Pqn <- nl_probs(P7$coef, Wq, P7$lambda_nested)
ape7 <- function(j) { bj <- c(0, P7$coef[, j]); colSums(wq7 * Pq * (matrix(bj, nrow(Pq), 3, byrow = TRUE) - drop(Pq %*% bj))) }
truth_B7 <- list(coef = P7$coef, lambda_nested = P7$lambda_nested,
                 shares = setNames(colSums(wq7 * Pq), P7$levels), shares_nested = setNames(colSums(wq7 * Pqn), P7$levels),
                 ape_x = setNames(ape7(2), P7$levels), ape_s = setNames(ape7(3), P7$levels),
                 within_paid_nested = (P7$coef["Enterprise", ] - P7$coef["Pro", ]) / P7$lambda_nested)
rm(q7, Wq, wq7, Pq, Pqn)
set.seed(seed$B7_large)
L7 <- gen_B7(P7, n = 1e6)
truth_B7$large <- list(N = nrow(L7), shares = setNames(as.vector(table(L7$choice)) / nrow(L7), P7$levels),
                       ape_x = setNames(mnl_ape(P7$coef, cbind(1, L7$x, L7$s), 2), P7$levels))
rm(L7)

## =================================================================
## save
## =================================================================
saveRDS(B2, "helix_B2.rds"); saveRDS(B3, "helix_B3.rds"); saveRDS(B4, "helix_B4.rds"); saveRDS(B5, "helix_B5.rds")
saveRDS(B6, "helix_B6.rds"); saveRDS(B7, "helix_B7.rds"); saveRDS(B7n, "helix_B7n.rds")
truth12 <- list(B1 = truth_B1, B2 = truth_B2, B3 = truth_B3, B4 = truth_B4, B5 = truth_B5, B6 = truth_B6, B7 = truth_B7,
                params = list(seed = seed, B1 = P1, B2 = P2, B3 = P3, B4 = P4, B5 = P5, B6 = P6, B7 = P7))
saveRDS(truth12, "truth12.rds")

cat("\n=== v3 blocks (truth12.rds) ===\n")
cat(sprintf("B1 exact AME %.5f (legacy oracle %.5f), LPM projection slope %.5f, probit pseudo-true (%.4f, %.4f)\n",
            ame_exact, ame_logit, proj_slope, pq[1], pq[2]))
cat(sprintf("B2 adoption %.3f, avg cross-partial %.5f (large %.5f), share positive %.3f\n",
            mean(B2$y), truth_B2$avg_cross, truth_B2$large$avg_cross, truth_B2$share_pos))
cat(sprintf("B3 adoption %.3f, APE of D %.5f (large %.5f), naive probit APE plim %.5f\n",
            mean(B3$Y), truth_B3$ape, truth_B3$large$ape, truth_B3$large$naive_ape))
cat(sprintf("B4 adoption %.3f, beta %.2f -> attenuated %.4f, APE %.5f (large %.5f)\n",
            mean(B4$Y), truth_B4$beta, truth_B4$beta_att, truth_B4$ape, truth_B4$large$ape))
cat(sprintf("B5 Pro share %.3f, exact APE %.6f (4-month DGP: %.6f; months 1-4 of the T = 8 panel: %.6f), pooled logit plim %.5f (4-month DGP: %.5f; months 1-4 of T = 8: %.5f)\n",
            mean(B5$Y), truth_B5$ape, truth_B5$ape_T4, T8first4$ape, T8$pooled_plim, T4$pooled_plim, T8first4$pooled_plim))
cat(sprintf("   large-draw checks: APE %.5f / %.5f / %.5f, pooled plim %.4f / %.4f / %.4f (T8 / 4-month DGP / months 1-4)\n",
            T8$draw$ape, T4$draw$ape, T8first4$draw$ape, T8$draw$pooled_plim, T4$draw$pooled_plim, T8first4$draw$pooled_plim))
cat(sprintf("B6 rating shares %s, APE on P(5) %.5f (large %.5f)\n",
            paste(sprintf("%.3f", truth_B6$shares), collapse = " "), truth_B6$ape_top, truth_B6$large$ape_top))
cat(sprintf("B7 shares (MNL) %s; nested %s\n", paste(sprintf("%.3f", truth_B7$shares), collapse = " "),
            paste(sprintf("%.3f", truth_B7$shares_nested), collapse = " ")))
