# ============================================================================
# 01_dgp.R  --  Meridian analytics-dashboard randomized-encouragement experiment
# ----------------------------------------------------------------------------
# A B2B SaaS platform ("Meridian") wants to know whether client firms that ADOPT
# a newly launched analytics dashboard (D) end up with higher engagement /
# retention (Y).  Adoption is voluntary and ENDOGENOUS: an unobserved firm
# "sophistication" (low latent resistance U_D) drives BOTH adoption and the
# outcome, so OLS -- even adjusted for observed covariates -- is upward biased.
#
# Meridian ran a randomized ENCOURAGEMENT: a random subset of clients got an
# extra customer-success push (onboarding invite + tutorial), Z in {0,1}.
# Compliance is imperfect -> compliers / always-takers / never-takers, and
# monotonicity holds BY CONSTRUCTION (Z enters the adoption index with a
# positive coefficient, so raising Z never pushes anyone OUT of adoption:
# no defiers).
#
# Framework: Heckman-Vytlacil threshold-crossing / MTE model.
#   U_D ~ Uniform(0,1)  = unobserved resistance to adoption (0 = sophisticated,
#                         adopts at the slightest push; 1 = stubborn hold-out).
#   D(z) = 1{ P(z) > U_D }, propensity P(z) = Phi(a0 + a1*z).
#   Types with binary Z:
#       always-taker : U_D <  P(0)          (adopts regardless)
#       complier     : P(0) <= U_D < P(1)   (adopts iff encouraged)
#       never-taker  : U_D >= P(1)          (never adopts)
#   Outcome (correlated random coefficient, so LATE != ATE):
#       Y0 = b_y0 + b_w*W - g0*U_D + e0
#       tau_i = t0 - t1*U_D + e_tau   (gain FALLS with resistance)
#       Y1 = Y0 + tau ;  Y = Y0 + D*tau
#   => MTE(u) = E[tau | U_D = u] = t0 - t1*u   (a clean straight line over U_D).
#
# We ALSO build a CONTINUOUS instrument variant (encouragement INTENSITY, a
# randomized onboarding dosage Z_cont ~ U(0,1)) whose propensity spans (0,1),
# giving full support of U_D so the whole MTE curve is traced; and a WEAK-IV
# variant (tiny a1) to exhibit the weak-instrument pathology.
#
# All TRUE estimands (ATE, ATT, complier-LATE) are computed exactly from the
# generated potential outcomes.  Seed is fixed.
# ============================================================================

set.seed(404)

## ---- sample size -----------------------------------------------------------
n <- 30000L

## ---- structural parameters -------------------------------------------------
# adoption index (main binary design): P(z) = Phi(a0 + a1*z)
a0     <- -0.20   # -> P(0) = Phi(-0.20) = 0.4207  (always-taker share)
a1     <-  1.40   # -> P(1) = Phi( 1.20) = 0.8849  (strong first stage)

# adoption index (WEAK-IV variant): tiny encouragement effect. Chosen so the
# first-stage F stays low (~4) even at this large n -> weak-instrument pathology.
a1_weak <- 0.028  # -> P(1) = Phi(-0.172) = 0.4317, complier share ~ 0.011

# adoption index (CONTINUOUS instrument variant): full-support propensity
ac0    <- -2.40
ac1    <-  4.80   # index in [-2.4, 2.4] -> P in [0.008, 0.992]

# outcome parameters
b_y0   <- 60.0    # baseline engagement level
b_w    <-  5.0    # observed covariate (firm size) effect on outcome
g0     <- 12.0    # baseline FALLS with resistance U_D (sophistication -> high Y0)
t0     <-  6.0    # MTE intercept
t1     <-  8.0    # MTE slope: gain falls with resistance  => MTE(u) = 6 - 8u
sd_e0  <-  3.0
sd_tau <-  1.5

# covariate <-> unobserved resistance link (so covariate-adjusted OLS still
# fails: W proxies sophistication only partially, residual selection remains)
rho_wu <- -0.70   # corr(W, latent resistance) : bigger firms adopt more easily

## ---- draw the unobservables and covariate ----------------------------------
U_D    <- runif(n)                       # unobserved resistance, ~ U(0,1)
latent <- qnorm(U_D)                     # standard-normal version of resistance
W      <- rho_wu * latent + sqrt(1 - rho_wu^2) * rnorm(n)   # firm size (std)

## ---- potential outcomes (same for every instrument variant) ---------------
e0   <- rnorm(n, 0, sd_e0)
etau <- rnorm(n, 0, sd_tau)
Y0   <- b_y0 + b_w * W - g0 * U_D + e0
tau  <- t0 - t1 * U_D + etau
Y1   <- Y0 + tau

## ---- instruments -----------------------------------------------------------
Z       <- rbinom(n, 1, 0.5)             # main binary encouragement
Z_cont  <- runif(n)                      # continuous encouragement intensity

## propensities
P0_main <- pnorm(a0 + a1 * 0)
P1_main <- pnorm(a0 + a1 * 1)
Pz      <- pnorm(a0 + a1 * Z)            # per-firm propensity, main design
Pz_weak <- pnorm(a0 + a1_weak * Z)       # weak-IV design
Pz_cont <- pnorm(ac0 + ac1 * Z_cont)     # continuous design

## adoption decisions (monotone in each instrument)
D      <- as.integer(U_D < Pz)           # main design
D_weak <- as.integer(U_D < Pz_weak)      # weak-IV design
D_cont <- as.integer(U_D < Pz_cont)      # continuous design

## observed outcomes under each design
Y      <- Y0 + D      * tau
Y_weak <- Y0 + D_weak * tau
Y_cont <- Y0 + D_cont * tau

## ---- compliance types (main binary design) --------------------------------
# always-taker: adopts under both z ; complier: adopts only if encouraged ;
# never-taker : adopts under neither.  Monotonicity => no defiers.
type <- ifelse(U_D <  P0_main, "always",
        ifelse(U_D <  P1_main, "complier", "never"))
type <- factor(type, levels = c("always", "complier", "never"))
is_complier <- U_D >= P0_main & U_D < P1_main

## ---- TRUE estimands (exact, from generated potential outcomes) -------------
TRUE_ATE  <- mean(tau)                    # population average
TRUE_ATT  <- mean(tau[D == 1])            # among the treated (main design)
TRUE_LATE <- mean(tau[is_complier])       # complier average = what IV targets

## ---- assemble and persist --------------------------------------------------
dat <- data.frame(
  Y, D, Z, W,
  Y_weak, D_weak,
  Y_cont, D_cont, Z_cont,
  U_D, W_size = W, Y0, Y1, tau,
  Pz, Pz_weak, Pz_cont,
  type, is_complier
)

params <- list(
  n = n, seed = 404,
  a0 = a0, a1 = a1, a1_weak = a1_weak, ac0 = ac0, ac1 = ac1,
  b_y0 = b_y0, b_w = b_w, g0 = g0, t0 = t0, t1 = t1,
  sd_e0 = sd_e0, sd_tau = sd_tau, rho_wu = rho_wu,
  P0_main = P0_main, P1_main = P1_main,
  P0_weak = pnorm(a0), P1_weak = pnorm(a0 + a1_weak),
  complier_share = P1_main - P0_main,
  always_share = P0_main, never_share = 1 - P1_main,
  complier_share_weak = pnorm(a0 + a1_weak) - pnorm(a0)
)

truth <- list(ATE = TRUE_ATE, ATT = TRUE_ATT, LATE = TRUE_LATE,
              mte_intercept = t0, mte_slope = -t1)

out_dir <- "/Users/bobbyzhong/Library/Mobile Documents/com~apple~CloudDocs/Research/Method/Micro Theory_Microeconometrics/notes_v2/code/12_iv_late_mte"
saveRDS(list(dat = dat, params = params, truth = truth),
        file.path(out_dir, "meridian_data.rds"))

## ---- console summary -------------------------------------------------------
cat("=== Meridian DGP (seed", params$seed, ", n =", n, ") ===\n")
cat(sprintf("Propensities (main): P(0)=%.4f  P(1)=%.4f\n", P0_main, P1_main))
cat(sprintf("Type shares (main): always=%.4f  complier=%.4f  never=%.4f\n",
            mean(type == "always"), mean(type == "complier"), mean(type == "never")))
cat(sprintf("First-stage E[D|Z=1]-E[D|Z=0] (main) = %.4f\n",
            mean(D[Z == 1]) - mean(D[Z == 0])))
cat(sprintf("First-stage (weak)                    = %.4f  (complier share %.4f)\n",
            mean(D_weak[Z == 1]) - mean(D_weak[Z == 0]),
            params$complier_share_weak))
cat(sprintf("D_cont support of P: [%.3f, %.3f]\n", min(Pz_cont), max(Pz_cont)))
cat("\n--- TRUE estimands (oracle) ---\n")
cat(sprintf("TRUE ATE           = %.4f\n", TRUE_ATE))
cat(sprintf("TRUE ATT (treated) = %.4f\n", TRUE_ATT))
cat(sprintf("TRUE complier-LATE = %.4f\n", TRUE_LATE))
cat(sprintf("MTE(u) = %.2f - %.2f*u   (ATE = intercept - slope/2 = %.3f)\n",
            t0, t1, t0 - t1/2))
cat("\nSaved: meridian_data.rds\n")
