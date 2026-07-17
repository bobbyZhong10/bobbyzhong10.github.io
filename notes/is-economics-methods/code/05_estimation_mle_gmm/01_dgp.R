# ============================================================================
# 01_dgp.R  --  Helix developer-platform adoption + usage DGP
# ----------------------------------------------------------------------------
# "Helix" is a fictional B2B developer platform. Each client firm i has an
# observed standardized "engagement" covariate x_i (past API activity, z-scored).
# Two outcomes are generated, both driven by x through a nonlinear link:
#
#   (A) ADOPTION  D_i in {0,1}: whether firm i adopts a paid "Pro" tier.
#       True model is LOGIT:  P(D=1|x) = Lambda(a0 + a1 x),  Lambda = plogis.
#       => the object of MLE. A naive linear-probability OLS gives nonsense
#          (predicted probabilities outside [0,1]); logit MLE fixes it and
#          lets us exhibit Jensen consistency + Information-Matrix equality.
#
#   (B) USAGE  y_i in {0,1,2,...}: monthly count of Pro API calls (in 1000s).
#       True conditional MEAN is  E[y|x] = exp(b0 + b1 x)  (log link), BUT the
#       true distribution is NEGATIVE BINOMIAL (overdispersed), NOT Poisson.
#       => a data scientist who fits POISSON MLE is running QMLE: the mean is
#          correctly specified so beta is CONSISTENT, but the Information-Matrix
#          equality FAILS (Var != mean), so naive Fisher-information SEs are
#          wrong and the Huber-White SANDWICH is needed. The same conditional-
#          mean restriction  E[(y - exp(b0+b1 x)) h(x)] = 0  becomes a GMM
#          system: h=(1,x) just-identifies (= Poisson score), h=(1,x,x^2)
#          over-identifies -> efficient weighting + Hansen J test of the mean.
#
# All TRUE parameters are known; oracle targets are computed exactly / by a
# giant reference sample. Seed fixed. Scripts read/write RELATIVE paths, so
# run with the working directory set to this folder.
# ============================================================================

set.seed(512)

## ---- true parameters -------------------------------------------------------
a0 <- -0.40; a1 <- 1.10          # logit adoption:  P(D=1|x) = plogis(a0 + a1 x)
b0 <-  1.00; b1 <- 0.50          # log-mean usage:  E[y|x]  = exp(b0 + b1 x)
theta_nb <- 2.0                  # NB size (dispersion): Var(y|x) = mu + mu^2/theta

n_main <- 4000                   # main analysis sample

gen_helix <- function(n) {
  x  <- rnorm(n)                                   # standardized engagement
  ## (A) adoption -- logit
  pD <- plogis(a0 + a1 * x)
  D  <- rbinom(n, 1, pD)
  ## (B) usage -- negative binomial with log mean (overdispersed)
  mu <- exp(b0 + b1 * x)
  y  <- rnbinom(n, size = theta_nb, mu = mu)       # Var = mu + mu^2/theta_nb
  data.frame(x = x, D = D, pD = pD, y = y, mu = mu)
}

dat <- gen_helix(n_main)

## ---- oracle truths (from a 5,000,000-row reference draw) -------------------
set.seed(9001)
big <- gen_helix(5e6)

# average marginal effect (AME) of x on adoption probability, logit truth:
#   d/dx P = a1 * lambda(a0 + a1 x),  lambda = dLogistic
ame_logit <- mean(a1 * dlogis(a0 + a1 * big$x))
# marginal effect at the mean (MEM), x = 0:
mem_logit <- a1 * dlogis(a0 + a1 * 0)
base_rate <- mean(big$D)

oracle <- list(
  a0 = a0, a1 = a1, b0 = b0, b1 = b1, theta_nb = theta_nb,
  ame_logit = ame_logit, mem_logit = mem_logit, base_rate = base_rate
)

saveRDS(dat,    "helix_data.rds")
saveRDS(oracle, "oracle.rds")

cat(sprintf("n_main            = %d\n", n_main))
cat(sprintf("adoption rate     = %.4f (oracle base rate %.4f)\n", mean(dat$D), base_rate))
cat(sprintf("mean(y)           = %.4f\n", mean(dat$y)))
cat(sprintf("var(y)/mean(y)    = %.4f  (>1 => overdispersion vs Poisson)\n",
            var(dat$y) / mean(dat$y)))
cat(sprintf("oracle AME (logit)= %.4f\n", ame_logit))
cat(sprintf("oracle MEM (logit)= %.4f\n", mem_logit))
