# ============================================================================
# 04_mte.R  --  Marginal Treatment Effects with the continuous instrument
# ----------------------------------------------------------------------------
# The binary encouragement only identified LATE on one complier window.  With a
# CONTINUOUS instrument (randomized onboarding INTENSITY) whose propensity spans
# (0,1), we trace the whole MTE curve over unobserved resistance U_D, and show
# that ATE / ATT / LATE are all weighted integrals of the SAME MTE curve.
#
# Three estimators:
#   (A) Local IV (Heckman-Vytlacil): E[Y|P=p] = E[Y0] + integral_0^p MTE(u)du,
#       so MTE(u) = d/dp E[Y|P=p].  With a linear MTE this makes E[Y|P] a
#       QUADRATIC in P; fit Y ~ P + P^2 and differentiate.
#   (B) Semiparametric derivative: loess of Y on P_hat, numerical derivative.
#   (C) Control-function two-step: plug E[U_D | D, P] into the outcome equation
#       and recover the structural MTE (t0, t1) directly.
#
# Truth:  MTE(u) = 6 - 8u.
# ============================================================================

suppressMessages(library(fixest))

d <- readRDS("meridian_data.rds")
dat <- d$dat; truth <- d$truth; par <- d$params
t0 <- par$t0; t1 <- par$t1

## ---- first stage for the CONTINUOUS design: propensity score ---------------
# True propensity is probit in Z_cont, so a probit first stage is consistent.
ps_fit <- glm(D_cont ~ Z_cont, data = dat, family = binomial(link = "probit"))
P_hat  <- predict(ps_fit, type = "response")

## ---- (A) Local IV: quadratic E[Y|P], differentiate -------------------------
loc <- lm(Y_cont ~ P_hat + I(P_hat^2), data = dat)
b1 <- coef(loc)["P_hat"]            # -> estimates t0
b2 <- coef(loc)["I(P_hat^2)"]       # -> estimates -t1/2
mte_a_int   <- unname(b1)           # MTE intercept
mte_a_slope <- unname(2 * b2)       # MTE slope
mte_a_slope_se <- unname(2 * summary(loc)$coef["I(P_hat^2)", "Std. Error"])
MTE_A <- function(u) mte_a_int + mte_a_slope * u

## ---- (B) Semiparametric derivative of E[Y|P] -------------------------------
lo <- loess(Y_cont ~ P_hat, data = dat, span = 0.5, degree = 2)
ug <- seq(0.02, 0.98, by = 0.01)
Eyp <- predict(lo, newdata = data.frame(P_hat = ug))
mte_semi <- c(NA, diff(Eyp) / diff(ug))     # numerical derivative on grid

## ---- (C) Control-function two-step -----------------------------------------
# Y = b_y0 + b_w W - g0 U_D + D(t0 - t1 U_D) + e.  With U_D ~ U(0,1) and
# D = 1{U_D < P}: E[U_D|D=1,P]=P/2, E[U_D|D=0,P]=(1+P)/2.  Substitute:
#   Y = b_y0 + b_w W + t0 D - g0*cf_main - t1*cf_int + resid
# where cf_main = D*P/2 + (1-D)*(1+P)/2 , cf_int = D*P/2.
cf_main <- dat$D_cont * (P_hat / 2) + (1 - dat$D_cont) * (1 + P_hat) / 2
cf_int  <- dat$D_cont * (P_hat / 2)
cf <- lm(Y_cont ~ W_size + D_cont + cf_main + cf_int, data = dat)
cf_t0 <- unname(coef(cf)["D_cont"])       # -> t0
cf_t1 <- unname(-coef(cf)["cf_int"])      # -> t1
cf_ate <- cf_t0 - cf_t1 / 2
MTE_C <- function(u) cf_t0 - cf_t1 * u

## bootstrap SE for the control-function ATE (accounts for two-step)
set.seed(99)
Bboot <- 300; n <- nrow(dat); ate_bs <- numeric(Bboot)
for (b in seq_len(Bboot)) {
  ix <- sample.int(n, n, replace = TRUE)
  db <- dat[ix, ]
  ph <- predict(glm(D_cont ~ Z_cont, data = db, family = binomial("probit")),
                type = "response")
  m1 <- db$D_cont * (ph / 2) + (1 - db$D_cont) * (1 + ph) / 2
  m2 <- db$D_cont * (ph / 2)
  cb <- lm(Y_cont ~ W_size + D_cont + m1 + m2, data = db)
  ate_bs[b] <- coef(cb)["D_cont"] + coef(cb)["m2"] / 2   # t0 - t1/2
}
cf_ate_se <- sd(ate_bs)

## ---- ATE / ATT / LATE as WEIGHTED INTEGRALS of the MTE ---------------------
# We integrate the CONTROL-FUNCTION MTE line: with the propensity piling up near
# 0 and 1, the fully nonparametric local-IV slope is consistent but imprecise at
# this n (slope SE ~ 1.3), whereas the control-function estimate, which imposes
# linearity in U_D, pins the curve down tightly.  The two agree on shape.
gg <- seq(0.0025, 0.9975, by = 0.0025); du <- gg[2] - gg[1]
mte_grid <- MTE_C(gg)

# ATE: uniform weight
ate_int <- sum(mte_grid) * du

# ATT weight w_TT(u) = P(P > u) / E[P]   (continuous-design propensity dist.)
EP <- mean(P_hat)
S_P <- sapply(gg, function(u) mean(P_hat > u))
w_tt <- S_P / EP
att_int <- sum(mte_grid * w_tt) * du
att_oracle_cont <- mean(dat$tau[dat$D_cont == 1])

# LATE on the binary-instrument complier window [P0_main, P1_main]: uniform there
win <- gg >= par$P0_main & gg < par$P1_main
late_int <- mean(mte_grid[win])

cat("=== MTE via the continuous instrument ===\n")
cat(sprintf("First-stage probit propensity support: [%.3f, %.3f]\n",
            min(P_hat), max(P_hat)))
cat("\n--- (A) Local IV (quadratic E[Y|P], nonparametric cross-check) ---\n")
cat(sprintf("MTE(u) = %.3f + (%.3f)*u   slope SE %.3f   [true: 6.000 - 8.000*u]\n",
            mte_a_int, mte_a_slope, mte_a_slope_se))
cat("   (consistent but imprecise: propensity concentrates near 0 and 1)\n")
cat("\n--- (C) Control-function two-step ---\n")
cat(sprintf("recovered t0 = %.3f (true 6), t1 = %.3f (true 8)\n", cf_t0, cf_t1))
cat(sprintf("CF ATE = t0 - t1/2 = %.3f  (bootstrap SE %.3f)  [true %.3f]\n",
            cf_ate, cf_ate_se, truth$ATE))
cat("\n--- ATE/ATT/LATE as weighted integrals of the (control-function) MTE ---\n")
cat(sprintf("ATE  = uniform-weighted integral   = %.3f   [true ATE  %.3f]\n",
            ate_int, truth$ATE))
cat(sprintf("ATT  = P(P>u)/E[P]-weighted        = %.3f   [oracle ATT(cont) %.3f]\n",
            att_int, att_oracle_cont))
cat(sprintf("LATE = uniform over [%.3f,%.3f]     = %.3f   [binary-IV 2SLS 0.792, true LATE %.3f]\n",
            par$P0_main, par$P1_main, late_int, truth$LATE))

## ---- save curve + markers for the figure -----------------------------------
saveRDS(list(
  mte_a_int = mte_a_int, mte_a_slope = mte_a_slope, mte_a_slope_se = mte_a_slope_se,
  cf_t0 = cf_t0, cf_t1 = cf_t1, cf_ate = cf_ate, cf_ate_se = cf_ate_se,
  ate_int = ate_int, att_int = att_int, late_int = late_int,
  att_oracle_cont = att_oracle_cont,
  P0 = par$P0_main, P1 = par$P1_main,
  grid = gg, mte_true = 6 - 8 * gg,
  mte_cf = mte_grid,               # control-function line (used for integrals)
  mte_A_line = MTE_A(gg),          # local-IV line (nonparametric cross-check)
  semi_u = ug, semi_mte = mte_semi,
  w_tt = w_tt),
  "res_mte.rds")
