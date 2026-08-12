# ============================================================================
# 01_dgp.R  —  Northwind AI assistant adoption: data-generating process
# ----------------------------------------------------------------------------
# Outcome  : Y = log(tickets resolved per agent) at a client firm.
# Covariates X (pre-treatment, all observed):
#   size  = firm size (standardized log employees)
#   dig   = existing digitization level (standardized)
#   eff   = baseline ticketing efficiency (standardized)
#   ind   = industry, 4 categories (factor)
# Structural (potential-outcome) model, SHARED across all datasets:
#   Y_i(0) = mu(X_i) + phi * U_i + eps_i
#   tau_i  = tau_base + g_obs(X_i) + delta * U_i          (individual effect)
#   Y_i(1) = Y_i(0) + tau_i
# where U_i is an UNOBSERVED firm-quality / growth factor that raises the
# baseline outcome (phi) and the gain (delta); g_obs is the observable gain.
#
# Two datasets from the SAME structural effects:
#   (A) RCT subsample     : D assigned by a coin flip  -> clean ATE.
#   (B) Observational     : D ~ propensity in X AND the (observable) gain
#                           -> selection on gains, naive comparison biased up.
# CIA toggle: the observational sample carries TWO adoption vectors built from
#   the same firms and same U:
#     D_cia  (lambda = 0)  selection loads only on X and g_obs(X); U is left
#                          out of selection so U _|_ D | X  -> CIA holds,
#                          adjustment recovers the ATT.
#     D_fail (lambda > 0)  selection also loads on the unobserved U, which
#                          itself moves Y(0); adopters then differ in an
#                          unobserved way  -> CIA fails, residual bias remains
#                          (this is the endogeneity that motivates IV, ch.10).
# ============================================================================

set.seed(20260711)

# ---- sample sizes ----------------------------------------------------------
n_obs <- 6000L   # observational study (Study B)
n_rct <- 2000L   # randomized pilot   (Study A)

# ---- structural parameters -------------------------------------------------
# baseline outcome  mu(X) = a0 + b*size + b*dig + b*eff + industry shift
a0     <- 1.50
b_size <- 0.30
b_dig  <- 0.25
b_eff  <- 0.20
ind_shift <- c(A = 0.00, B = 0.10, C = 0.20, D = 0.30)
sd_eps <- 0.40   # idiosyncratic outcome noise
phi    <- 0.30   # loading of the UNOBSERVED factor U on the baseline Y(0)

# individual treatment effect  tau_i = tau_base + g_obs(X) + delta*U
tau_base <- 0.20                       # => true ATE ~ 0.20 (X, U mean 0)
gg_dig   <- 0.12                       # gain rises with digitization
gg_size  <- 0.10                       # gain rises with size
gg_eff   <- -0.05                      # low-efficiency firms gain a bit more
delta    <- 0.15                       # loading of UNOBSERVED gain U on tau

# selection index (Study B): P(adopt) = plogis(theta0 + X terms + gain + lam*U)
theta0    <- -0.20
s_size    <- 0.40
s_dig     <- 0.50
s_eff     <- 0.40
s_gain    <- 1.50                      # selection ON the observable gain
lambda    <- 1.60                      # extra loading on U in the CIA-FAILS arm

# ---- helper: draw the shared structural pieces for n firms ------------------
draw_firms <- function(n) {
  size <- rnorm(n)
  dig  <- rnorm(n)
  eff  <- rnorm(n)
  ind  <- factor(sample(c("A","B","C","D"), n, replace = TRUE),
                 levels = c("A","B","C","D"))
  U    <- rnorm(n)                      # unobserved gain component
  eps  <- rnorm(n, 0, sd_eps)

  mu    <- a0 + b_size*size + b_dig*dig + b_eff*eff + ind_shift[as.character(ind)]
  g_obs <- gg_dig*dig + gg_size*size + gg_eff*eff
  tau   <- tau_base + g_obs + delta*U   # individual treatment effect

  Y0 <- as.numeric(mu + phi*U + eps)
  Y1 <- as.numeric(Y0 + tau)

  data.frame(size, dig, eff, ind, U, eps,
             mu = as.numeric(mu), g_obs = g_obs, tau = tau,
             Y0 = Y0, Y1 = Y1)
}

observed_Y <- function(df, D) df$Y0 + D * (df$Y1 - df$Y0)

# ============================================================================
# Study A: randomized pilot  (D ~ Bernoulli(0.5), independent of everything)
# ============================================================================
rct <- draw_firms(n_rct)
rct$D <- rbinom(n_rct, 1, 0.5)
rct$Y <- observed_Y(rct, rct$D)

# ============================================================================
# Study B: observational sample  (selection on X and on gains)
# ============================================================================
obs <- draw_firms(n_obs)

lin_common <- with(obs, theta0 + s_size*size + s_dig*dig + s_eff*eff + s_gain*g_obs)
p_cia  <- plogis(lin_common)                       # CIA holds  (lambda = 0)
p_fail <- plogis(lin_common + lambda * obs$U)      # CIA fails  (loads on U)

obs$e_true_cia  <- p_cia                            # true propensity, CIA arm
obs$e_true_fail <- p_fail
obs$D_cia  <- rbinom(n_obs, 1, p_cia)
obs$D_fail <- rbinom(n_obs, 1, p_fail)

obs$Y_cia  <- observed_Y(obs, obs$D_cia)
obs$Y_fail <- observed_Y(obs, obs$D_fail)

# ============================================================================
# Truths computable from the individual effects
# ============================================================================
true_ATE_pop  <- mean(obs$tau)                             # population ATE (Study B frame)
true_ATE_rct  <- mean(rct$tau)                             # ATE in the RCT sample
true_ATT_cia  <- mean(obs$tau[obs$D_cia  == 1])            # ATT, CIA arm
true_ATT_fail <- mean(obs$tau[obs$D_fail == 1])            # ATT, CIA-fails arm
adopt_rate_cia  <- mean(obs$D_cia)
adopt_rate_fail <- mean(obs$D_fail)

truths <- list(
  true_ATE_pop = true_ATE_pop,
  true_ATE_rct = true_ATE_rct,
  true_ATT_cia = true_ATT_cia,
  true_ATT_fail = true_ATT_fail,
  adopt_rate_cia = adopt_rate_cia,
  adopt_rate_fail = adopt_rate_fail
)

# quick diagnostic: naive gaps
naive_cia  <- mean(obs$Y_cia[obs$D_cia==1])  - mean(obs$Y_cia[obs$D_cia==0])
naive_fail <- mean(obs$Y_fail[obs$D_fail==1]) - mean(obs$Y_fail[obs$D_fail==0])
rct_dim    <- mean(rct$Y[rct$D==1]) - mean(rct$Y[rct$D==0])

cat("---- DGP truths ----\n")
print(round(unlist(truths), 4))
cat(sprintf("RCT diff-in-means      : %.4f\n", rct_dim))
cat(sprintf("naive obs (CIA arm)    : %.4f  (bias vs ATT = %.4f)\n",
            naive_cia,  naive_cia  - true_ATT_cia))
cat(sprintf("naive obs (CIA-fails)  : %.4f  (bias vs ATT = %.4f)\n",
            naive_fail, naive_fail - true_ATT_fail))
cat(sprintf("adopt rate CIA / fail  : %.3f / %.3f\n",
            adopt_rate_cia, adopt_rate_fail))

# ---- save for downstream scripts ------------------------------------------
outdir <- dirname(sub("^--file=", "",
                      grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(outdir) || length(outdir) == 0 || outdir == "") outdir <- "."
saveRDS(list(rct = rct, obs = obs, truths = truths,
             params = list(a0=a0,b_size=b_size,b_dig=b_dig,b_eff=b_eff,
                           ind_shift=ind_shift,sd_eps=sd_eps,phi=phi,
                           tau_base=tau_base,gg_dig=gg_dig,gg_size=gg_size,
                           gg_eff=gg_eff,delta=delta,theta0=theta0,
                           s_size=s_size,s_dig=s_dig,s_eff=s_eff,
                           s_gain=s_gain,lambda=lambda,
                           n_obs=n_obs,n_rct=n_rct)),
        file.path(outdir, "northwind.rds"))
cat("\nsaved:", file.path(outdir, "northwind.rds"), "\n")
