# ------------------------------------------------------------------
# 01_dgp.R -- "Cadence" payments-platform switching panel.
#
# Each of N small merchants routes card payments through ONE of two
# competing platforms, A or B, in each of T months. Two forces make
# a merchant's choices persistent over time:
#   (1) persistent unobserved heterogeneity theta_i : an intrinsic,
#       time-invariant taste for A vs B (fit with the merchant's
#       business). High-theta merchants sit on A almost regardless
#       of fees.
#   (2) TRUE state dependence delta : having used a platform last
#       month raises its utility this month (integration, staff
#       habit, stored data) -- a genuine switching cost / lock-in.
#
# The identification problem the chapter is built around: a naive
# dynamic logit that regresses this month's choice on last month's
# choice, without controlling for theta_i, loads theta_i's
# persistence onto the lagged-choice coefficient and BADLY overstates
# delta (Heckman's spurious state dependence). The initial choice is
# itself correlated with theta_i (the initial-conditions problem),
# which makes it worse.
#
# Latent index for choosing A over B at month t (t >= 2):
#   Ustar_it = theta_i + beta*dFee_t + delta*L_{i,t-1} + nu_it
#     dFee_t      = fee_At - fee_Bt        (dollars/month, exogenous)
#     L_{i,t-1}   = 1{last=A} - 1{last=B}  in {-1,+1}
#     nu_it       ~ Logistic(0,1)
#   choose A iff Ustar_it > 0.
# theta_i ~ N(0, sigma_theta^2). The month-1 choice (initial
# condition) is drawn from theta_i plus noise, so it is correlated
# with theta_i by construction.
# ------------------------------------------------------------------
library(data.table)

set.seed(20)

## ---- dimensions --------------------------------------------------
N_merch  <- 5000
T_months <- 10              # month 1 = initial condition; 2..10 estimated

## ---- structural parameters --------------------------------------
beta_true    <- -0.55       # fee sensitivity (per dollar/month)
delta_true   <-  0.90       # TRUE state dependence (switching cost)
sigma_theta  <-  1.40       # sd of persistent taste heterogeneity

## implied structural switching cost in dollars/month = delta/|beta|
sc_dollars_true <- delta_true / abs(beta_true)

## ---- persistent heterogeneity -----------------------------------
theta <- rnorm(N_merch, 0, sigma_theta)

## ---- exogenous monthly fee gap dFee_t (A minus B) ---------------
## Platforms move their merchant fees around a bit month to month;
## treat the gap as exogenous AR(1) time variation common to all.
dfee <- numeric(T_months)
dfee[1] <- rnorm(1, 0, 1.2)
for (t in 2:T_months) dfee[t] <- 0.5 * dfee[t - 1] + rnorm(1, 0, 1.2)

## ---- simulate choices -------------------------------------------
choiceA <- matrix(NA_integer_, N_merch, T_months)

## initial condition: month-1 choice correlated with theta_i
##   (merchants tend to START on the platform they intrinsically
##    prefer, plus onboarding noise) -> initial-conditions problem
u0 <- theta + rlogis(N_merch, 0, 1) + beta_true * dfee[1]
choiceA[, 1] <- as.integer(u0 > 0)

for (t in 2:T_months) {
  L_prev <- 2L * choiceA[, t - 1] - 1L               # {-1,+1}
  ustar  <- theta + beta_true * dfee[t] + delta_true * L_prev +
            rlogis(N_merch, 0, 1)
  choiceA[, t] <- as.integer(ustar > 0)
}

## ---- long panel --------------------------------------------------
dt <- data.table(
  merchant = rep(1:N_merch, each = T_months),
  month    = rep(1:T_months, times = N_merch),
  choiceA  = as.vector(t(choiceA)),
  theta    = rep(theta, each = T_months),
  dfee     = rep(dfee, times = N_merch)
)
setkey(dt, merchant, month)
dt[, lagA  := shift(choiceA), by = merchant]          # 0/1 last choice
dt[, Lprev := 2L * lagA - 1L]                          # {-1,+1}
dt[, initA := choiceA[1], by = merchant]               # month-1 choice
dt[, dfee_bar := mean(dfee), by = merchant]            # time-mean (constant here)

## ---- observed inertia (raw persistence in choices) --------------
est <- dt[month >= 2]
stay_rate <- est[, mean(choiceA == lagA)]              # share who repeat
switch_rate <- 1 - stay_rate

## ---- decomposition of persistence via counterfactual re-sim -----
## Re-simulate two counterfactual worlds holding the SAME shocks
## structure but shutting down one channel at a time, to quantify
## how much observed inertia is TRUE state dependence vs heterogeneity.
sim_stay <- function(delta_c, sigma_c, seed) {
  set.seed(seed)
  th <- rnorm(N_merch, 0, sigma_c)
  cA <- matrix(NA_integer_, N_merch, T_months)
  u0 <- th + rlogis(N_merch, 0, 1) + beta_true * dfee[1]
  cA[, 1] <- as.integer(u0 > 0)
  for (t in 2:T_months) {
    Lp <- 2L * cA[, t - 1] - 1L
    us <- th + beta_true * dfee[t] + delta_c * Lp + rlogis(N_merch, 0, 1)
    cA[, t] <- as.integer(us > 0)
  }
  mean(cA[, 2:T_months] == cA[, 1:(T_months - 1)])
}
stay_full      <- sim_stay(delta_true, sigma_theta, 20)   # both channels
stay_no_sd     <- sim_stay(0,          sigma_theta, 20)   # heterogeneity only
stay_no_hetero <- sim_stay(delta_true, 0,           20)   # state dependence only
stay_neither   <- sim_stay(0,          0,           20)   # neither (pure fee/noise)

cat("== Cadence switching-panel DGP ==\n")
cat(sprintf("Merchants: %d, months: %d (1 = initial condition)\n", N_merch, T_months))
cat(sprintf("True beta = %.3f, delta = %.3f, sigma_theta = %.3f\n",
            beta_true, delta_true, sigma_theta))
cat(sprintf("Implied structural switching cost = delta/|beta| = $%.2f / month\n",
            sc_dollars_true))
cat(sprintf("Overall share on A (months 2-10): %.3f\n", est[, mean(choiceA)]))
cat(sprintf("Raw stay rate (choice == last choice): %.3f  (switch %.3f)\n",
            stay_rate, switch_rate))
cat("Persistence decomposition (stay rate under counterfactual worlds):\n")
cat(sprintf("  both channels on            : %.3f\n", stay_full))
cat(sprintf("  heterogeneity only (delta=0): %.3f\n", stay_no_sd))
cat(sprintf("  state dep only (sigma=0)    : %.3f\n", stay_no_hetero))
cat(sprintf("  neither (fees+noise only)   : %.3f\n", stay_neither))

## ---- save --------------------------------------------------------
out_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
saveRDS(dt, file.path(out_dir, "cadence_panel.rds"))
fwrite(dt, file.path(out_dir, "cadence_panel.csv"))
saveRDS(list(beta = beta_true, delta = delta_true, sigma_theta = sigma_theta,
             sc_dollars = sc_dollars_true, dfee = dfee,
             stay_rate = stay_rate,
             stay_full = stay_full, stay_no_sd = stay_no_sd,
             stay_no_hetero = stay_no_hetero, stay_neither = stay_neither),
        file.path(out_dir, "truth24.rds"))
cat("Saved cadence_panel.rds / .csv / truth24.rds\n")
