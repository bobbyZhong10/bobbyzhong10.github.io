# ------------------------------------------------------------------
# 01_dgp.R -- the Larkspur simulation world for chapter 5
#             (M-estimation and nonlinear regression)
#
# Larkspur is a (fictional) online marketplace for small independent sellers of
# home goods, apparel, jewellery and stationery. Money is in US dollars.
# One world, four blocks (all new in v3; nothing is inherited from the legacy chapter):
#   M1   Spotlight experiment: a storefront-promotion feature is switched on for a
#        random half of 20,000 merchants; 30-day sales, many merchants sell nothing.
#        Spotlight changes both whether a merchant sells at all and how much.
#        Sales are recorded in cents and in USD (the same numbers / 100).
#   M2   exponential conditional mean, estimated by NLS (Gauss-Newton and stats::nls);
#        M2m the same kind of data with a misspecified mean (pseudo-true value).
#   M3   two-step estimation: a traffic index estimated in a first step (OLS of a noisy
#        week of log storefront visits on listing signals) is the regressor of a
#        second-step exponential NLS for sales.
#   M4   overdispersed order counts for Poisson QMLE; the Monte Carlo (05_mc_poisson.R)
#        replaces the retired legacy Helix Poisson demonstration (ruling 6).
# Every number in chapter 5 comes from these data and from truth05.rds.
#
# Outputs (in this directory):
#   larkspur_M1.rds   data.frame: merchant category tenure_std D cents usd y0_cents_true y1_cents_true
#   larkspur_M2.rds   data.frame: listings x1 free_ship y m_true
#   larkspur_M2m.rds  data.frame: x y m_true
#   larkspur_M3.rds   data.frame: z1 z2 v y t_true
#   larkspur_M4.rds   data.frame: x y mu_true
#   truth05.rds       list(M1, M2, M2m, M3, M4, params): truth and parameters
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 1.4): D_i, y_i(1), y_i(0) in M1; the unit multiplier a
# (a = 100 turns USD into cents); theta_* is a pseudo-true value; gamma is the
# first-step parameter of M3. Field names follow the chapter card.
#
# Seeds (every block its own set.seed, nothing inherited): M1 5101, M1 population
# check 5111, M2 5102, M2m 5103, M2m large-sample check 5113, M3 5104, M4 5105;
# the Monte Carlo seeds used by 04 and 05 are stored here too (M3 5114 and
# 52400 + j for grid point j, M4 5115 and 52500 + j).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("m_fns.R")

# Seeds were fixed before any result was seen and were not tuned afterwards. With these draws the
# M2 intercept interval and the M4 sample slope interval (2.07 sandwich SEs from beta1) miss their
# truth: expected 5% events. Coverage claims rest on the Monte Carlo (04, 05), not on one sample.
seed <- list(M1 = 5101, M1_pop = 5111, M2 = 5102, M2m = 5103, M2m_large = 5113, M3 = 5104, M4 = 5105,
             M3_mc = 5114, M3_grid = 52400, M4_mc = 5115, M4_grid = 52500)

## =================================================================
## M1: the Spotlight experiment (merchant cross-section, many zeros)
## =================================================================
## 20,000 merchants, complete randomisation: Spotlight (D = 1) for exactly half.
## Outcome: sales in the 30 days after assignment, recorded in whole cents.
## Under control a merchant in category c makes at least one sale with probability
## p0[c] (40%, 48%, 36%, 30%: about 60% of merchants sell nothing in a month);
## Spotlight multiplies that probability by k = 1.25 (the same latent draw decides,
## so nobody who sells without Spotlight stops selling with it). Given a sale, log
## sales are log(median0[c]) + 0.25 tenure_std + 1.1 eps (median 120 to 320 USD at
## average tenure, mean about 450 USD among sellers), and Spotlight adds
## delta = 0.10 (+10.5% for sellers). Then
##   E[y(d) | c, tenure] = p0[c] k^d exp(log median0[c] + delta d + 0.25 tenure + 1.1^2 / 2),
## an exponential mean in (D, category, tenure) with coefficient on D
##   tau = log(E[y(1)] / E[y(0)]) = log(1.25) + 0.10 = 0.3231 (a 38.1% increase),
## the same in every category: Poisson PML with category and tenure regressors is
## correctly specified, and so is PPML on D alone (randomisation).
P1 <- list(n = 20000, categories = c("home_decor", "apparel", "jewelry", "stationery"),
           pi = c(0.35, 0.25, 0.25, 0.15), p0 = c(0.40, 0.48, 0.36, 0.30), k = 1.25,
           median0 = c(180, 260, 320, 120), sigma = 1.1, beta_tenure = 0.25, delta = 0.10)
set.seed(seed$M1)
M1 <- gen_M1(P1)

share_sell0 <- sum(P1$pi * P1$p0)
mean_y0 <- sum(P1$pi * P1$p0 * P1$median0) * exp(P1$sigma^2 / 2 + P1$beta_tenure^2 / 2)
truth_M1 <- list(log_ratio = log(P1$k) + P1$delta, ratio = P1$k * exp(P1$delta),
                 share_sell0 = share_sell0, share_sell1 = P1$k * share_sell0, ext_effect = (P1$k - 1) * share_sell0,
                 mean_y0_usd = mean_y0, mean_y1_usd = P1$k * exp(P1$delta) * mean_y0)

## Population check and the log(1 + a y) estimands: one million merchants with both
## potential outcomes. a is the unit multiplier: y_a = a * y_usd (a = 0.001: thousands
## of USD, a = 1: USD, a = 100: cents). The population ATE on log(1 + a y) is what the
## OLS coefficient of log(1 + a y) on D estimates, and it moves with a.
set.seed(seed$M1_pop)
pop <- gen_M1(P1, n = 1e6)
y0p <- pop$y0_cents_true / 100
y1p <- pop$y1_cents_true / 100
a_grid <- 10^seq(-3, 4, by = 0.25)
ate_a <- vapply(a_grid, function(a) mean(log1p(a * y1p)) - mean(log1p(a * y0p)), numeric(1))
truth_M1$pop <- list(n = nrow(pop), log_ratio = log(mean(y1p) / mean(y0p)), mean_y0_usd = mean(y0p), mean_y1_usd = mean(y1p),
                     share_sell0 = mean(y0p > 0), share_sell1 = mean(y1p > 0),
                     ate_log1p_kusd = mean(log1p(y1p / 1000)) - mean(log1p(y0p / 1000)),
                     ate_log1p_usd = mean(log1p(y1p)) - mean(log1p(y0p)),
                     ate_log1p_cents = mean(log1p(100 * y1p)) - mean(log1p(100 * y0p)),
                     ate_asinh_usd = mean(asinh(y1p)) - mean(asinh(y0p)),
                     ate_asinh_cents = mean(asinh(100 * y1p)) - mean(asinh(100 * y0p)),
                     sweep = data.frame(a = a_grid, ate_log1p = ate_a))
rm(pop, y0p, y1p)

## =================================================================
## M2: exponential mean, NLS
## =================================================================
## 2,000 active merchants. listings: catalogue size, median 40, log sd 0.9;
## x1 = log(listings / 40); free_ship: 40% offer free shipping.
## Monthly sales (USD thousands): E[y | x] = exp(log 1.2 + 0.45 x1 + 0.20 free_ship):
## 1,200 USD at the median catalogue without free shipping, elasticity 0.45 with
## respect to catalogue size, +22% with free shipping. The error is multiplicative
## (mean-one lognormal, log sd 0.7), so Var(y | x) is proportional to E[y | x]^2:
## NLS is consistent, the homoskedastic (model-based) NLS variance is not.
P2 <- list(n = 2000, listings_median = 40, sd_log_listings = 0.9, p_free_ship = 0.4,
           theta = c(log(1.2), 0.45, 0.20), sigma_eta = 0.7)
set.seed(seed$M2)
M2 <- gen_M2(P2)
truth_M2 <- list(theta = P2$theta)

## =================================================================
## M2m: misspecified mean, pseudo-true value
## =================================================================
## 2,000 active merchants, x = log(catalogue size / 40) ~ N(0, 0.9^2). True mean with
## diminishing returns: E[y | x] = exp(0.2 + 0.55 x - 0.10 x^2) (USD thousands); the
## analyst fits exp(b0 + b1 x). NLS converges to the minimiser of
## E[(E[y | x] - exp(b0 + b1 x))^2] (White 1981), the pseudo-true value theta_*;
## Poisson PML converges to a different pseudo-true value. Both by quadrature
## (pseudo_true_exp), checked against NLS on one million noise-free draws.
## truth05$M2m$pseudo is the slope b1 of the NLS pseudo-true value (the coefficient the
## chapter's interval is about); pseudo_vec is (b0, b1).
P2m <- list(n = 2000, sd_x = 0.9, a0 = 0.2, a1 = 0.55, a2 = -0.10, sigma_eta = 0.5)
set.seed(seed$M2m)
M2m <- gen_M2m(P2m)
m0 <- function(x) exp(P2m$a0 + P2m$a1 * x + P2m$a2 * x^2)
pv <- pseudo_true_exp(m0, P2m$sd_x, "nls")
pv_pp <- pseudo_true_exp(m0, P2m$sd_x, "ppml")
set.seed(seed$M2m_large)
xl <- rnorm(1e6, 0, P2m$sd_x)
gl <- gauss_newton(m0(xl), cbind(1, xl), exp_mean, exp_jac, pv)
stopifnot(gl$converged)
truth_M2m <- list(pseudo = pv[2], pseudo_vec = pv, pseudo_ppml = pv_pp[2], pseudo_ppml_vec = pv_pp,
                  large = list(N = length(xl), pseudo_vec = gl$theta))
rm(xl)

## =================================================================
## M3: two-step M-estimation with a generated regressor
## =================================================================
## 2,000 merchants. z1: search-placement score, z2: log followers (both standardised).
## Traffic index t = 0.30 z1 + 0.20 z2 (log storefront visits relative to the median);
## the platform sees one week of log visits v = t + 0.6 e, a noisy measure (first-stage
## R^2 about 0.27). Monthly sales (USD thousands): E[y | z] = exp(log 1.5 + 0.8 t),
## an elasticity of 0.8 with respect to traffic, mean-one lognormal error (log sd 0.35).
## Step 1: OLS of v on (1, z1, z2) gives gamma_hat and t_hat = z'gamma_hat.
## Step 2: NLS of y on exp(theta0 + theta1 t_hat). Because theta1 != 0, the second-step
## score depends on gamma (Wooldridge 2010 eq. 12.36: F_o != 0, condition 12.37 fails),
## so the naive variance that treats t_hat as data is too small.
## Monte Carlo (04_twostep.R): R_mc = 2000 replications at theta1 = 0.8, and R_grid = 1000
## at each other theta1 on grid_theta1 (theta1 = 0 is the case where 12.37 holds).
P3 <- list(n = 2000, gamma = c(0, 0.30, 0.20), theta = c(log(1.5), 0.8), sigma_e = 0.6, sigma_eta = 0.35,
           R_mc = 2000, grid_theta1 = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), R_grid = 1000)
set.seed(seed$M3)
M3 <- gen_M3(P3)
var_t <- sum(P3$gamma[2:3]^2)
truth_M3 <- list(theta = P3$theta, gamma = P3$gamma, r2_first_stage = var_t / (var_t + P3$sigma_e^2))

## =================================================================
## M4: overdispersed counts for Poisson QMLE
## =================================================================
## 2,000 active merchants. x: standardised seller rating. Monthly orders y with
## E[y | x] = exp(1.4 + 0.45 x) (about 4 orders at the average rating) and negative
## binomial dispersion, size 3: Var(y | x) = mu + mu^2 / 3. Averaged over x the
## variance-to-mean ratio E[Var(y | x)] / E[y] = 1 + E[mu^2] / (3 E[mu]) is about 2.8.
## The mean is right, so Poisson QMLE is consistent (Gourieroux, Monfort and Trognon
## 1984), but the information-matrix equality fails: the model-based variance A_0^-1 / n
## is too small and the sandwich A_0^-1 B_0 A_0^-1 / n is right.
## Monte Carlo (05_mc_poisson.R): R_mc = 2000 replications at size 3, R_grid = 1000 at each
## other size on grid_size (Inf: Poisson, where the sandwich collapses to the model-based variance).
P4 <- list(n = 2000, beta = c(1.4, 0.45), size = 3, R_mc = 2000,
           grid_size = c(Inf, 12, 6, 3, 1.5, 0.75), R_grid = 1000)
set.seed(seed$M4)
M4 <- gen_M4(P4)
## Three different overdispersion quantities (population values; they are not comparable with each other):
##   avg_cond_var_to_mean  E[Var(y | x)] / E[y] = 1 + E[mu^2] / (size E[mu])        (about 2.83; Monte Carlo x-axis)
##   pearson_dispersion    E[Var(y | x) / E[y | x]] = 1 + E[mu] / size, the target of the Pearson
##                         statistic sum (y - mu_hat)^2 / mu_hat / (n - 2)           (about 2.50)
##   marginal_var_to_mean  Var(y) / E[y], which also contains Var(mu) from x         (about 3.84)
## with E[mu] = exp(beta0 + beta1^2 / 2) and E[mu^2] = exp(2 beta0 + 2 beta1^2) for x ~ N(0, 1).
Emu <- exp(P4$beta[1] + P4$beta[2]^2 / 2)
Emu2 <- exp(2 * P4$beta[1] + 2 * P4$beta[2]^2)
avg_cond <- function(size) 1 + Emu2 / (size * Emu)
truth_M4 <- list(beta = P4$beta, size = P4$size, avg_cond_var_to_mean = avg_cond(P4$size),
                 pearson_dispersion = 1 + Emu / P4$size,
                 marginal_var_to_mean = (Emu + Emu2 / P4$size + Emu2 - Emu^2) / Emu,
                 grid = data.frame(size = P4$grid_size, avg_cond_var_to_mean = avg_cond(P4$grid_size)))

## =================================================================
## Save
## =================================================================
saveRDS(M1, "larkspur_M1.rds"); saveRDS(M2, "larkspur_M2.rds"); saveRDS(M2m, "larkspur_M2m.rds")
saveRDS(M3, "larkspur_M3.rds"); saveRDS(M4, "larkspur_M4.rds")
tr <- list(M1 = truth_M1, M2 = truth_M2, M2m = truth_M2m, M3 = truth_M3, M4 = truth_M4,
           params = list(seed = seed, M1 = P1, M2 = P2, M2m = P2m, M3 = P3, M4 = P4))
saveRDS(tr, "truth05.rds")

cat(sprintf("M1: %d merchants; zero sales %.3f (control) %.3f (Spotlight); mean USD %.2f / %.2f; true log ratio %.4f\n",
            nrow(M1), mean(M1$usd[M1$D == 0] == 0), mean(M1$usd[M1$D == 1] == 0),
            mean(M1$usd[M1$D == 0]), mean(M1$usd[M1$D == 1]), tr$M1$log_ratio))
cat(sprintf("M1 population: log ratio %.4f; ATE on log(1 + y): %.4f (USD thousands) %.4f (USD) %.4f (cents)\n",
            tr$M1$pop$log_ratio, tr$M1$pop$ate_log1p_kusd, tr$M1$pop$ate_log1p_usd, tr$M1$pop$ate_log1p_cents))
cat(sprintf("M2m: NLS pseudo-true (%.4f, %.4f); PPML pseudo-true (%.4f, %.4f); large-sample NLS (%.4f, %.4f)\n",
            pv[1], pv[2], pv_pp[1], pv_pp[2], gl$theta[1], gl$theta[2]))
cat(sprintf("M3: first-stage R^2 (population) %.3f; M4: E[Var(y | x)] / E[y] %.3f\n",
            tr$M3$r2_first_stage, tr$M4$avg_cond_var_to_mean))
cat("Saved larkspur_M1.rds ... larkspur_M4.rds, truth05.rds\n")
