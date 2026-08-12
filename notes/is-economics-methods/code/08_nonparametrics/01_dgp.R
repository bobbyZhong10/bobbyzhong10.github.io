# ============================================================================
# 01_dgp.R  --  Solstice marketplace price-demand DGP (two facets, one seed)
# ----------------------------------------------------------------------------
# "Solstice" is a fictional marketplace. For each listing we see a price x and
# the resulting log-demand y. The true demand curve is NONLINEAR: demand is
# flat (inelastic) at low and high prices and drops steeply through a middle
# "sensitive zone", so the local price sensitivity varies enormously with the
# price. A single linear / constant-elasticity fit reports ONE slope and hides
# all of this.
#
#   g(x) = 3 - 1.8 * plogis(1.5 * (x - 5.5))          (the true E[log-demand|x])
#
#  (A) DEMAND-CURVE facet:  y = g(x) + eps.  Nonparametric regression (kernel,
#      local linear, splines) recovers g; CV picks the bandwidth; bias-variance
#      is visible as under- vs over-smoothing measured against the known g.
#
#  (B) PARTIALLY-LINEAR facet:  y = beta*z + g(x) + eps, with ad-exposure z
#      CORRELATED with price. A linear-in-price model biases beta; Robinson's
#      double-residual estimator partials x out of y and z nonparametrically
#      and recovers beta at root-n despite g being nonparametric.
#
# All truths known. Seed fixed. Scripts use relative paths; run from this dir.
# ============================================================================

set.seed(848)

# true log-demand: a gentle overall decline plus a steep "sensitive zone" near 5
g  <- function(x) 4 - 0.18 * x - 1.2 * plogis(2 * (x - 5))
gp <- function(x) -0.18 - 1.2 * 2 * dlogis(2 * (x - 5))   # its derivative (slope)

sigma <- 0.35

## ---- (A) demand curve ------------------------------------------------------
n <- 600
xA <- runif(n, 1, 9)
yA <- g(xA) + rnorm(n, 0, sigma)
demand <- data.frame(x = xA, y = yA)

## ---- (B) partially linear: z NONLINEARLY related to price ------------------
## ad exposure rises with price but in an S-shaped, nonlinear way, so a linear-
## in-price control cannot remove its confounding; Robinson's flexible control can.
beta <- 0.8
mz <- function(x) 0.15 * x + 1.1 * plogis(1.8 * (x - 5))   # nonlinear E[z|x]
z <- mz(xA) + rnorm(n, 0, 0.7)
yB <- beta * z + g(xA) + rnorm(n, 0, sigma)
partial <- data.frame(x = xA, z = z, y = yB)

## reference grid for evaluating estimators against the known truth
grid <- seq(1, 9, length = 200)

## true local slope and elasticity at a few prices (elasticity = x * g'(x))
prices <- c(2, 5.5, 8)
truth_slope <- gp(prices)
truth_elas  <- prices * gp(prices)

saveRDS(list(demand = demand, partial = partial, grid = grid,
             g = g, gp = gp, sigma = sigma,
             truth = list(beta = beta, g_grid = g(grid),
                          prices = prices, slope = truth_slope, elas = truth_elas)),
        "solstice.rds")

cat("=== Solstice price-demand DGP ===\n")
cat(sprintf("n = %d listings, price in [1,9], sigma = %.2f\n", n, sigma))
cat("true local slope g'(x) at prices 2 / 5.5 / 8:\n  ")
cat(sprintf("%.3f  %.3f  %.3f\n", truth_slope[1], truth_slope[2], truth_slope[3]))
cat("true point elasticity x*g'(x) at 2 / 5.5 / 8:\n  ")
cat(sprintf("%.3f  %.3f  %.3f\n", truth_elas[1], truth_elas[2], truth_elas[3]))
cat(sprintf("partially-linear truth: beta = %.2f (ad exposure), z corr with price\n", beta))
cat(sprintf("  cor(z, x) = %.3f\n", cor(z, xA)))
