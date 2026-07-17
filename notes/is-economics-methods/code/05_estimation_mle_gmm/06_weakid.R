# ============================================================================
# 06_weakid.R  --  weak identification as a FLAT objective. Identification asks
# the population objective to have a well-separated peak at theta_0, i.e. the
# moment Jacobian D = E[dg/dtheta] to be far from zero (a steep criterion).
# When the identifying variation is thin, D -> 0, the criterion goes nearly
# flat, and even though theta is "technically" identified the estimator has a
# wide, non-normal sampling distribution and standard errors explode.
# We contrast a STRONG design (regressor varies) with a WEAK one (regressor
# barely varies), both estimating b1 in the Poisson conditional mean.
# ============================================================================

oracle <- readRDS("oracle.rds")
b0 <- oracle$b0; b1 <- oracle$b1; th <- oracle$theta_nb

## one estimate of (b0,b1) by just-identified GMM = Poisson score, given data
fitb <- function(x, y) coef(glm(y ~ x, family = poisson()))

## curvature of the criterion in b1: |D_{b1}| ~ E[x^2 mu] (steepness proxy)
curv <- function(x) mean(x^2 * exp(b0 + b1*x))

sim <- function(sd_x, n = 4000, R = 1000, seed = 1) {
  set.seed(seed)
  est <- t(replicate(R, {
    x  <- rnorm(n, 0, sd_x)
    mu <- exp(b0 + b1*x)
    y  <- rnbinom(n, size = th, mu = mu)
    m  <- glm(y ~ x, family = poisson())
    c(b1 = coef(m)[2], se = sqrt(sandwich::sandwich(m)[2,2]))
  }))
  list(b1 = est[,1], se = est[,2], curv = curv(rnorm(2e5, 0, sd_x)))
}

strong <- sim(sd_x = 1.00, seed = 41)
weak   <- sim(sd_x = 0.05, seed = 42)

report <- function(tag, s) {
  cov <- mean(abs(s$b1 - b1) <= 1.96 * s$se)
  cat(sprintf("%-8s: median b1 = %+.4f  mean = %+.4f  SD = %.4f  IQR = %.3f\n",
              tag, median(s$b1), mean(s$b1), sd(s$b1),
              IQR(s$b1)))
  cat(sprintf("          avg reported SE = %.4f   curvature |D_b1| = %.4f   95%% cover = %.3f\n",
              mean(s$se), s$curv, cov))
}
cat("=== weak identification: flat criterion, exploding sampling variance ===\n")
cat(sprintf("truth b1 = %.2f\n", b1))
report("STRONG", strong)
report("WEAK",   weak)
cat(sprintf("\nSD(weak)/SD(strong) = %.1fx ;  curvature strong/weak = %.1fx\n",
            sd(weak$b1)/sd(strong$b1), strong$curv/weak$curv))

## objective-function slice Q(b1) at b0 fixed, for one strong vs one weak draw
slice <- function(sd_x, seed) {
  set.seed(seed); n <- 4000
  x  <- rnorm(n, 0, sd_x); mu <- exp(b0 + b1*x); y <- rnbinom(n, size=th, mu=mu)
  grid <- seq(b1 - 0.8, b1 + 0.8, length = 121)
  Q <- sapply(grid, function(bb){
    g1 <- mean(y - exp(b0 + bb*x)); g2 <- mean(x*(y - exp(b0 + bb*x)))
    g1^2 + g2^2
  })
  data.frame(b1 = grid, Q = Q / max(Q))    # normalized
}
sl_strong <- slice(1.00, 41); sl_weak <- slice(0.05, 42)

saveRDS(list(strong_sd=sd(strong$b1), weak_sd=sd(weak$b1),
             strong_curv=strong$curv, weak_curv=weak$curv,
             strong_cov=mean(abs(strong$b1-b1)<=1.96*strong$se),
             weak_cov=mean(abs(weak$b1-b1)<=1.96*weak$se),
             weak_median=median(weak$b1), weak_meanse=mean(weak$se),
             sl_strong=sl_strong, sl_weak=sl_weak), "res_weak.rds")
