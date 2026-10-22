# ------------------------------------------------------------------
# 01_dgp.R -- Halcyon "Growth Playbook" DGP for Module 13 (Causal ML).
#
# Setting: Halcyon is a marketplace platform. It rolled out an AI
# merchant-advisory tool ("Growth Playbook") to a subset of merchants.
# Rollout was OBSERVATIONAL, not randomized: merchants who were already
# on a growth trajectory were more likely to be offered / to adopt the
# tool. Unconfoundedness holds given a high-dimensional covariate vector
# X (p = 20 merchant features), but propensity e(X) and baseline g(X)
# are nonlinear, so credible ATE/CATE estimation needs ML for the
# nuisance functions plus Neyman-orthogonal, cross-fit moments.
#
# Outcome Y : merchant 90-day log GMV growth (in log points).
# Treatment D: adopted Growth Playbook (1) or not (0).
# CATE tau(X): heterogeneous; largest for SMALL, NEW merchants, can be
#   negative for large established merchants (tool distracts them).
#
# Writes: halcyon.rds (full data.frame), truth.rds (population targets).
# Run from this directory.
# ------------------------------------------------------------------
set.seed(1313)

n <- 10000L
p <- 20L

## ---- covariates X: correlated merchant features -------------------
# X1  = merchant size (standardized log catalog size); small = negative
# X2  = tenure on platform (standardized); new = negative
# X3  = pre-period demand momentum (standardized)
# X4..X20 = other operational / category features (many irrelevant)
Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
L <- chol(Sigma)
X <- matrix(rnorm(n * p), n, p) %*% L
colnames(X) <- paste0("X", 1:p)
X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]

## ---- shared nonlinear confounding structure -----------------------
# The SAME nonlinear terms (an interaction and a curvature) drive both
# adoption and baseline outcome, so linear-in-X controls cannot remove
# the confounding but a flexible ML learner can. This is exactly what
# makes ML nuisances necessary rather than optional.
h_conf <- 0.7 * X1 * X3 + 0.6 * (X2^2 - 1)

## ---- propensity e(X): nonlinear, confounded -----------------------
# Bigger, higher-momentum merchants are more likely to adopt.
ps_index <- 0.40 * X1 + 0.32 * X3 + 0.50 * h_conf
e_true   <- plogis(ps_index)
e_true   <- pmin(pmax(e_true, 0.05), 0.95)   # enforce overlap
D <- rbinom(n, 1, e_true)

## ---- baseline g0(X) = E[Y(0)|X]: nonlinear ------------------------
g0 <- 0.20 * X3 + 0.15 * X1 + 0.50 * h_conf + 0.10 * X[, 6]

## ---- CATE tau(X): heterogeneous treatment effect ------------------
# Small (X1 low) and new (X2 low) merchants gain most; large + old lose.
tau <- function(x1, x2) 0.12 - 0.10 * x1 - 0.06 * x2 + 0.03 * x1 * x2
tau_i <- tau(X1, X2)

## ---- assemble outcome --------------------------------------------
sigma_y <- 0.30
Y <- g0 + tau_i * D + rnorm(n, 0, sigma_y)

dat <- data.frame(Y = Y, D = D, X)
dat$e_true <- e_true
dat$tau_true <- tau_i

## ---- population / sample truths ----------------------------------
ATE_true <- mean(tau_i)                 # sample-average treatment effect
ATT_true <- mean(tau_i[D == 1])
# oracle best linear projection of tau on (X1, X2)
blp_oracle <- coef(lm(tau_i ~ X1 + X2))

truth <- list(
  n = n, p = p,
  ATE = ATE_true, ATT = ATT_true,
  blp_oracle = blp_oracle,
  tau_fun = tau,
  share_treated = mean(D),
  mean_e = mean(e_true),
  overlap_min = min(e_true), overlap_max = max(e_true)
)

saveRDS(dat, "halcyon.rds")
saveRDS(truth, "truth.rds")

cat(sprintf("n = %d, p = %d, share treated = %.3f\n", n, p, mean(D)))
cat(sprintf("True SATE = %.4f   ATT = %.4f\n", ATE_true, ATT_true))
cat(sprintf("Overlap: e in [%.3f, %.3f]\n", min(e_true), max(e_true)))
cat(sprintf("Oracle BLP tau ~ 1 + X1 + X2: b0=%.4f b1=%.4f b2=%.4f\n",
            blp_oracle[1], blp_oracle[2], blp_oracle[3]))
cat(sprintf("tau range: [%.3f, %.3f]; share tau<0 = %.3f\n",
            min(tau_i), max(tau_i), mean(tau_i < 0)))
