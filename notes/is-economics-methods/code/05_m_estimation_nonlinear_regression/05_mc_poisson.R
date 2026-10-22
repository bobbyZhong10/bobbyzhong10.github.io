# ------------------------------------------------------------------
# 05_mc_poisson.R -- block M4: Poisson QMLE on overdispersed order counts.
#
# Replaces the retired legacy Helix Poisson demonstration (ruling 6) on Larkspur's own data.
# The conditional mean exp(beta0 + beta1 x) is right, the variance is negative binomial
# (Var = mu + mu^2 / size), so Poisson PML is consistent but the information-matrix equality
# A_0 = B_0 fails. Three standard errors for beta1_hat:
#   model     A^-1 / n, the Poisson (Fisher information) variance, vcov(glm(family = poisson));
#   sandwich  A^-1 B A^-1 / n (HC0, m_vcov on the Poisson scores; = sandwich::sandwich(glm));
#   quasi     model times the Pearson dispersion (glm(family = quasipoisson)): right only when
#             Var(y | x) is proportional to the mean, which the negative binomial is not.
# Monte Carlo: R_mc = 2000 replications at size 3 (seed M4_mc) and R_grid = 1000 at each other size
# on the grid (seed M4_grid + j for grid point j); size = Inf is Poisson, where B_0 = A_0 and the
# sandwich collapses to the model-based variance. Coverage of est -/+ qnorm(0.975) se.
# Each replication fits with glm.fit (the same IRLS as glm) and hand-coded scores, for speed.
#
# Output: res_mc.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("m_fns.R")
suppressPackageStartupMessages(library(sandwich))

tr <- readRDS("truth05.rds")
P4 <- tr$params$M4
M4 <- readRDS("larkspur_M4.rds")
z <- qnorm(0.975)
b1 <- P4$beta[2]

## ---- the sample ------------------------------------------------------------------------
g <- glm(y ~ x, family = poisson(), data = M4, control = glm.control(epsilon = 1e-12, maxit = 100))
X <- cbind(1, M4$x)
pp <- poisson_parts(M4$y, X, coef(g))
n <- nrow(M4)
se_model <- sqrt(vcov(g)[2, 2])
se_sand <- sqrt(m_vcov(pp$S, pp$A)[2, 2])
stopifnot(abs(se_sand / sqrt(sandwich::sandwich(g)[2, 2]) - 1) < 1e-6, abs(se_model / sqrt(solve(pp$A)[2, 2] / n) - 1) < 1e-6)
phi <- sum((M4$y - pp$mu)^2 / pp$mu) / (n - 2)                     # Pearson dispersion (1 under Poisson; target truth05$M4$pearson_dispersion)
gq <- glm(y ~ x, family = quasipoisson(), data = M4, control = glm.control(epsilon = 1e-12, maxit = 100))
se_quasi <- sqrt(vcov(gq)[2, 2])
est <- unname(coef(g)[2])
# with seed 5105 the sandwich interval misses beta1 (2.07 SEs): an expected 5% event, seed not tuned (see 01_dgp.R seeds)
# marginal_var_to_mean is Var(y) / E[y] over the whole sample (contains Var(mu) from x): not comparable with
# the Monte Carlo x-axis E[Var(y | x)] / E[y]; its population value is truth05$M4$marginal_var_to_mean
sample <- list(est = est, beta = unname(coef(g)), se_model = se_model, se_sandwich = se_sand, se_quasi = se_quasi,
               ci_model = ci95(est, se_model), ci_sandwich = ci95(est, se_sand), pearson_dispersion = phi,
               se_ratio = se_sand / se_model, marginal_var_to_mean = var(M4$y) / mean(M4$y), mean_y = mean(M4$y),
               share_zero = mean(M4$y == 0))

## ---- Monte Carlo --------------------------------------------------------------------------
one_rep <- function(size) {
  d <- gen_M4(P4, size = size)
  Xr <- cbind(1, d$x)
  fit <- glm.fit(Xr, d$y, family = poisson(), control = list(epsilon = 1e-10, maxit = 50, trace = FALSE))
  p <- poisson_parts(d$y, Xr, fit$coefficients)
  Ai <- solve(p$A)
  v_model <- Ai[2, 2] / P4$n
  dispersion <- sum((d$y - p$mu)^2 / p$mu) / (P4$n - 2)
  c(est = fit$coefficients[[2]], se_model = sqrt(v_model), se_sandwich = sqrt(m_vcov(p$S, p$A)[2, 2]),
    se_quasi = sqrt(v_model * dispersion), converged = fit$converged)
}
summarise_mc <- function(sim, size) {
  ok <- sim[, "converged"] == 1
  s <- sim[ok, , drop = FALSE]
  cov <- function(se) mean(abs(s[, "est"] - b1) <= z * s[, se])
  list(R = nrow(sim), n_converged = sum(ok), size = size,
       avg_cond_var_to_mean = 1 + exp(P4$beta[1] + 1.5 * P4$beta[2]^2) / size,     # E[Var(y | x)] / E[y]
       cover_model = cov("se_model"), cover_sandwich = cov("se_sandwich"), cover_quasi = cov("se_quasi"),
       mean_est = mean(s[, "est"]), sd_est = sd(s[, "est"]),
       mean_se_model = mean(s[, "se_model"]), mean_se_sandwich = mean(s[, "se_sandwich"]), mean_se_quasi = mean(s[, "se_quasi"]))
}
run_mc <- function(size, R, seed) {
  set.seed(seed)
  t(vapply(seq_len(R), function(r) one_rep(size), numeric(5)))
}
t0 <- proc.time()[["elapsed"]]
mc <- summarise_mc(run_mc(P4$size, P4$R_mc, tr$params$seed$M4_mc), P4$size)
grid <- do.call(rbind, lapply(seq_along(P4$grid_size), function(j) {
  sz <- P4$grid_size[j]
  s <- if (identical(sz, P4$size)) mc else summarise_mc(run_mc(sz, P4$R_grid, tr$params$seed$M4_grid + j), sz)
  as.data.frame(s)
}))
mc_seconds <- proc.time()[["elapsed"]] - t0
stopifnot(mc$n_converged == mc$R, all(grid$n_converged == grid$R))

res <- c(mc, list(sample = sample, grid = grid))   # run time is printed, not saved (byte-reproducible .rds)
saveRDS(res, "res_mc.rds")

cat(sprintf("M4 sample: beta1_hat %.4f; SE model %.4f, sandwich %.4f, quasi %.4f; Pearson dispersion %.3f\n",
            est, se_model, se_sand, se_quasi, phi))
cat(sprintf("M4 Monte Carlo (%d reps, size %g): coverage model %.4f, sandwich %.4f, quasi %.4f; sd %.4f; mean SE %.4f / %.4f (%.1f s incl. grid)\n",
            mc$R, mc$size, mc$cover_model, mc$cover_sandwich, mc$cover_quasi, mc$sd_est, mc$mean_se_model, mc$mean_se_sandwich, mc_seconds))
print(grid[, c("size", "avg_cond_var_to_mean", "R", "cover_model", "cover_sandwich", "cover_quasi", "sd_est", "mean_se_model", "mean_se_sandwich")],
      row.names = FALSE)
cat("Saved res_mc.rds\n")
