# ============================================================================
# 03_mle_logit.R  --  logit MLE: hand-rolled Newton-Raphson vs glm; consistency
# across sample sizes; and the sandwich = Fisher-information SE agreement under
# correct specification (Information-Matrix equality holds for a true logit).
# ============================================================================

suppressMessages({ library(sandwich); library(lmtest) })
dat    <- readRDS("helix_data.rds")
oracle <- readRDS("oracle.rds")

## ---- hand-rolled logit MLE via Newton-Raphson ------------------------------
# log-lik: sum y*log(F) + (1-y)*log(1-F), F = plogis(Xb)
# score:   sum (y - F(Xb)) X        (Conlon: nabla ell = sum (y - F(Z)) X)
# Hessian: -sum F(1-F) X X'
X <- cbind(1, dat$x); y <- dat$D
newton_logit <- function(X, y, tol = 1e-10, maxit = 100) {
  b <- rep(0, ncol(X))
  for (it in 1:maxit) {
    Xb <- as.vector(X %*% b); F <- plogis(Xb)
    score <- crossprod(X, y - F)                 # k x 1
    W <- F * (1 - F)
    H <- -crossprod(X, X * W)                     # -sum F(1-F) X X'
    step <- solve(H, score)
    b <- b - step
    if (max(abs(step)) < tol) break
  }
  list(coef = b, iter = it, F = plogis(as.vector(X %*% b)), W = W)
}
fit <- newton_logit(X, y)

## cross-check against glm
g <- glm(D ~ x, data = dat, family = binomial())

## Fisher-information SE (model-based) vs Huber-White sandwich SE (robust)
Fh   <- fit$F; Wh <- Fh * (1 - Fh)
I_hat <- crossprod(X, X * Wh)                      # observed = expected info (logit)
V_fisher <- solve(I_hat)
se_fisher <- sqrt(diag(V_fisher))
# sandwich: (sum F(1-F)XX')^{-1} (sum (y-F)^2 XX') (.)^{-1}
B <- crossprod(X, X * (y - Fh)^2)
V_sand <- V_fisher %*% B %*% V_fisher
se_sand <- sqrt(diag(V_sand))

cat("=== logit MLE (Newton-Raphson, hand-rolled) ===\n")
cat(sprintf("iterations to converge = %d\n", fit$iter))
cat(sprintf("beta0 = %.4f   beta1 = %.4f   (truth: %.2f, %.2f)\n",
            fit$coef[1], fit$coef[2], oracle$a0, oracle$a1))
cat(sprintf("glm    beta0 = %.4f   beta1 = %.4f  (max abs diff %.2e)\n",
            coef(g)[1], coef(g)[2], max(abs(fit$coef - coef(g)))))
cat("=== SEs: Fisher-information vs sandwich (should agree under correct spec) ===\n")
cat(sprintf("beta1 Fisher SE   = %.4f\n", se_fisher[2]))
cat(sprintf("beta1 sandwich SE = %.4f\n", se_sand[2]))
cat(sprintf("beta1 glm/vcovHC  = %.4f / %.4f\n",
            sqrt(vcov(g)[2,2]), sqrt(vcovHC(g, type="HC0")[2,2])))

## ---- average marginal effect (AME) recovered from logit --------------------
ame_hat <- mean(dlogis(as.vector(X %*% fit$coef)) * fit$coef[2])
cat(sprintf("AME (logit MLE)   = %.4f   (oracle %.4f)\n", ame_hat, oracle$ame_logit))

## ---- consistency across sample sizes (fresh draws, same DGP) ---------------
set.seed(778)
gen <- function(n){ x<-rnorm(n); D<-rbinom(n,1,plogis(oracle$a0+oracle$a1*x)); data.frame(x,D) }
ns  <- c(100, 400, 1600, 6400, 25600, 102400)
R   <- 400
cons <- sapply(ns, function(n){
  est <- replicate(R, coef(glm(D~x, gen(n), family=binomial()))[2])
  c(mean = mean(est), sd = sd(est), rmse = sqrt(mean((est-oracle$a1)^2)))
})
colnames(cons) <- ns
cat("\n=== consistency of beta1_hat (truth 1.10): mean / sd / rmse across n ===\n")
print(round(cons, 4))

saveRDS(list(coef = fit$coef, se_fisher = se_fisher, se_sand = se_sand,
             ame = ame_hat, consistency = cons), "res_mle.rds")
