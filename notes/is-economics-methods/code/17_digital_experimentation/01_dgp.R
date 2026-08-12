# ------------------------------------------------------------------
# 01_dgp.R -- Lumen A/B platform DGP for Module 14 (Digital experiments).
#
# Setting: Lumen is a delivery marketplace running online controlled
# experiments (A/B tests). The flagship metric is 14-day revenue per
# user, Y. It is noisy: a mean around 5 currency units with an SD around
# 6, because a minority of heavy users dominate spending (the Lewis-Reiley
# problem: the signal you want is tiny next to user-level noise).
#
# Each user carries a PRE-experiment covariate X (their revenue in the
# 14 days before the test), strongly correlated with Y. This is the raw
# material for CUPED variance reduction.
#
# The feature under test lifts revenue by a true 2% (delta = 0.10 on a
# base of 5). Everything downstream (power/MDE, CUPED, peeking) is
# calibrated to these primitives. Writes lumen_params.rds.
# ------------------------------------------------------------------
set.seed(1414)

## ---- population primitives ---------------------------------------
mu0      <- 5.00     # mean control revenue per user
beta_v   <- 5.00     # loading on persistent user value v ~ N(0,1)
sd_idio  <- 3.317    # idiosyncratic SD so that Var(Y0) = 25 + 11 = 36
sd_Y0    <- sqrt(beta_v^2 + sd_idio^2)                 # = 6
delta    <- 0.10     # true additive treatment effect (2% of mu0)
rel_lift <- delta / mu0

# implied pre/post correlation (used by CUPED)
rho <- beta_v^2 / (beta_v^2 + sd_idio^2)               # = 0.6944 as corr^2? see below
rho_corr <- beta_v^2 / sqrt((beta_v^2 + sd_idio^2) * (beta_v^2 + sd_idio^2))

## a generator we reuse everywhere -----------------------------------
# returns data.frame with X (pre), Y0 (control potential), Y1 (treated
# potential). Treatment is a simple additive lift; heavy-user structure
# comes from the persistent value v plus a lognormal spike.
gen_users <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mu0 <- 5.00; beta_v <- 5.00; sd_idio <- 3.317; delta <- 0.10
  v   <- rnorm(n)                                  # persistent user value
  X   <- mu0 + beta_v * v + rnorm(n, 0, sd_idio)   # pre-period revenue
  Y0  <- mu0 + beta_v * v + rnorm(n, 0, sd_idio)   # control potential
  Y1  <- Y0 + delta                                # treated potential
  data.frame(X = X, Y0 = Y0, Y1 = Y1)
}

## verify the moments on a large sample ------------------------------
chk <- gen_users(2e5, seed = 1)
emp_sd   <- sd(chk$Y0)
emp_rho  <- cor(chk$X, chk$Y0)
emp_mean <- mean(chk$Y0)

params <- list(mu0 = mu0, beta_v = beta_v, sd_idio = sd_idio,
               sd_Y0 = sd_Y0, delta = delta, rel_lift = rel_lift,
               rho_corr = emp_rho, gen = gen_users,
               emp_sd = emp_sd, emp_mean = emp_mean)
saveRDS(params, "lumen_params.rds")

cat(sprintf("Control revenue: mean %.3f, SD %.3f (target 5, 6)\n",
            emp_mean, emp_sd))
cat(sprintf("Pre/post correlation rho = %.3f (CUPED reduction 1-rho^2 = %.3f)\n",
            emp_rho, 1 - emp_rho^2))
cat(sprintf("True effect delta = %.2f  (relative lift %.1f%%)\n",
            delta, 100 * rel_lift))
