# ============================================================================
# 05_testing.R  --  the testing trinity (Wald / LR / LM) on a Halcyon logit,
# their asymptotic agreement and finite-sample gap, plus the Wald statistic's
# non-invariance to reparameterization.
#   repeat_i ~ Bernoulli( Lambda(g0 + g1 * tenure_i) )
# H0: g1 = 0  (tenure does not predict repeat ordering).
# ============================================================================

set.seed(818)
n <- 1500
tenure <- rnorm(n)                                   # standardized customer tenure
g0 <- -0.3; g1 <- 0.35
p  <- plogis(g0 + g1 * tenure)
rep_order <- rbinom(n, 1, p)

## unrestricted and restricted logit fits
mu <- glm(rep_order ~ tenure, family = binomial())   # unrestricted
mr <- glm(rep_order ~ 1,      family = binomial())    # restricted: g1 = 0

## --- Wald: (g1_hat)^2 / Var(g1_hat) ---
g1_hat <- unname(coef(mu)["tenure"]); v1 <- vcov(mu)["tenure","tenure"]
W  <- g1_hat^2 / v1

## --- Likelihood ratio: 2 (ll_u - ll_r) ---
LR <- 2 * (as.numeric(logLik(mu)) - as.numeric(logLik(mr)))

## --- Lagrange multiplier / score: score at restricted fit ---
# score of tenure coefficient, evaluated at restricted model, / info
pr <- predict(mr, type = "response")                 # fitted prob under H0
X  <- cbind(1, tenure)
s_r  <- colSums(X * (rep_order - pr))                 # score vector at restricted
I_r  <- t(X) %*% (X * (pr*(1-pr)))                    # info at restricted
LM <- as.numeric(t(s_r) %*% solve(I_r) %*% s_r)

cat("=== testing trinity for H0: g1 = 0 (chi^2_1; truth g1 = 0.35) ===\n")
cat(sprintf("g1_hat = %.4f\n", g1_hat))
cat(sprintf("Wald W  = %.4f   p = %.5f\n", W,  pchisq(W,  1, lower.tail = FALSE)))
cat(sprintf("LR      = %.4f   p = %.5f\n", LR, pchisq(LR, 1, lower.tail = FALSE)))
cat(sprintf("LM      = %.4f   p = %.5f\n", LM, pchisq(LM, 1, lower.tail = FALSE)))
cat("(all three ~ chi^2_1, agree closely; small finite-sample gaps)\n")

## --- Wald non-invariance: test H0: g1 = 0.5 two equivalent ways -------------
## (a) directly on g1 ;  (b) on the reparameterization phi = log(g1)
g1v <- g1_hat; se1 <- sqrt(v1)
W_direct <- (g1v - 0.5)^2 / v1
# delta-method SE of log(g1): d/dg1 log(g1) = 1/g1
se_log <- se1 / g1v
W_repar <- (log(g1v) - log(0.5))^2 / se_log^2
cat("\n=== Wald non-invariance: same null g1 = 0.5, two parameterizations ===\n")
cat(sprintf("test on g1        : W = %.4f  p = %.4f\n", W_direct,
            pchisq(W_direct, 1, lower.tail = FALSE)))
cat(sprintf("test on log(g1)   : W = %.4f  p = %.4f\n", W_repar,
            pchisq(W_repar, 1, lower.tail = FALSE)))
cat("(different answers for an identical hypothesis -> Wald is not invariant)\n")

saveRDS(list(g1_hat=g1_hat, W=W, LR=LR, LM=LM,
             pW=pchisq(W,1,lower.tail=FALSE), pLR=pchisq(LR,1,lower.tail=FALSE),
             pLM=pchisq(LM,1,lower.tail=FALSE),
             W_direct=W_direct, W_repar=W_repar,
             pWd=pchisq(W_direct,1,lower.tail=FALSE),
             pWr=pchisq(W_repar,1,lower.tail=FALSE)), "res_test.rds")
