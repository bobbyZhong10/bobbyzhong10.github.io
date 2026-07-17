# ============================================================================
# 04_qmle_poisson.R  --  Poisson QMLE on overdispersed (NB) usage counts.
# The conditional MEAN E[y|x]=exp(b0+b1 x) is correct, so Poisson MLE is a
# QMLE and is CONSISTENT for (b0,b1). But the true Var(y|x) != E[y|x], so the
# Information-Matrix equality FAILS: the naive Fisher-information ("model") SE
# is wrong (too small), and the Huber-White SANDWICH SE is the right one.
# We verify by Monte-Carlo: the empirical SD of b1_hat across many draws is the
# truth the SEs are trying to estimate.
# ============================================================================

suppressMessages({ library(sandwich); library(lmtest) })
dat    <- readRDS("helix_data.rds")
oracle <- readRDS("oracle.rds")

## Poisson MLE (QMLE here). glm gives the model/Fisher SE by default.
pm <- glm(y ~ x, data = dat, family = poisson())

se_model <- sqrt(diag(vcov(pm)))                 # Fisher-info SE (assumes Var=mean)
se_sand  <- sqrt(diag(sandwich(pm)))             # Huber-White sandwich (robust)

cat("=== Poisson QMLE on overdispersed counts ===\n")
cat(sprintf("b0_hat = %.4f  b1_hat = %.4f   (truth %.2f, %.2f)\n",
            coef(pm)[1], coef(pm)[2], oracle$b0, oracle$b1))
cat(sprintf("b1  model/Fisher SE = %.4f\n", se_model[2]))
cat(sprintf("b1  sandwich   SE   = %.4f\n", se_sand[2]))
cat(sprintf("sandwich / model ratio (b1) = %.3f\n", se_sand[2]/se_model[2]))
cat(sprintf("b0  model SE = %.4f   sandwich SE = %.4f  (ratio %.3f)\n",
            se_model[1], se_sand[1], se_sand[1]/se_model[1]))

## overdispersion diagnostic: Pearson dispersion phi = sum r^2 / (n-k)
mu_hat <- fitted(pm); pearson <- (dat$y - mu_hat)/sqrt(mu_hat)
phi <- sum(pearson^2)/(nrow(dat) - 2)
cat(sprintf("Pearson dispersion phi = %.4f  (=1 under Poisson; >1 => overdispersed)\n", phi))

## ---- Monte-Carlo truth: SD of b1_hat across 2000 fresh NB draws ------------
set.seed(313)
genNB <- function(n){
  x  <- rnorm(n); mu <- exp(oracle$b0 + oracle$b1*x)
  y  <- rnbinom(n, size = oracle$theta_nb, mu = mu); data.frame(x, y)
}
R <- 2000; n <- nrow(dat)
mc <- t(replicate(R, {
  d  <- genNB(n); m <- glm(y ~ x, d, family = poisson())
  c(b1 = coef(m)[2], se_model = sqrt(vcov(m)[2,2]), se_sand = sqrt(sandwich(m)[2,2]))
}))
sd_true      <- sd(mc[,"b1.x"])
mean_se_mod  <- mean(mc[,"se_model"])
mean_se_sand <- mean(mc[,"se_sand"])
cov_model <- mean(abs(mc[,"b1.x"] - oracle$b1) <= 1.96*mc[,"se_model"])
cov_sand  <- mean(abs(mc[,"b1.x"] - oracle$b1) <= 1.96*mc[,"se_sand"])

cat("\n=== Monte-Carlo (2000 reps): which SE tells the truth? ===\n")
cat(sprintf("empirical SD of b1_hat (the truth) = %.4f\n", sd_true))
cat(sprintf("avg model/Fisher SE   = %.4f   (understates truth)\n", mean_se_mod))
cat(sprintf("avg sandwich   SE     = %.4f   (matches truth)\n", mean_se_sand))
cat(sprintf("95%% CI coverage: model = %.3f   sandwich = %.3f\n", cov_model, cov_sand))
cat(sprintf("mean b1_hat across reps = %.4f  (consistent for truth %.2f)\n",
            mean(mc[,"b1.x"]), oracle$b1))

saveRDS(list(coef=coef(pm), se_model=se_model, se_sand=se_sand, phi=phi,
             sd_true=sd_true, mean_se_mod=mean_se_mod, mean_se_sand=mean_se_sand,
             cov_model=cov_model, cov_sand=cov_sand), "res_qmle.rds")
