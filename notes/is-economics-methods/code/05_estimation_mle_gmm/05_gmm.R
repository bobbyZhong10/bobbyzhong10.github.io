# ============================================================================
# 05_gmm.R  --  the same conditional-mean restriction, as GMM.
#   E[(y - exp(b0+b1 x)) h(x)] = 0.
#   h=(1,x)      -> JUST-identified: these ARE the Poisson score FOCs.
#   h=(1,x,x^2)  -> OVER-identified (q=3, k=2): identity vs efficient weighting;
#                   Hansen J test of the exp-mean specification (df = q-k = 1).
# Hand-rolled so the weight matrix and two-step update are fully visible.
# ============================================================================

dat    <- readRDS("helix_data.rds")
oracle <- readRDS("oracle.rds")
y <- dat$y; x <- dat$x; n <- length(y)

## moment matrix G_i(theta) = h(x_i) * (y_i - exp(b0+b1 x_i)),  n x q
gmat <- function(b, H) {
  r <- y - exp(b[1] + b[2]*x)        # residual
  H * r                              # n x q
}
gbar   <- function(b, H) colMeans(gmat(b, H))
Qobj   <- function(b, H, W) { g <- gbar(b, H); as.numeric(t(g) %*% W %*% g) }
Jac    <- function(b, H) {           # D = d gbar / d theta  (q x k)
  mu <- exp(b[1] + b[2]*x)
  cbind(colMeans(H * (-mu)), colMeans(H * (-mu * x)))
}
gmm_fit <- function(H, W, b0 = c(0,0)) {
  o <- optim(b0, Qobj, H = H, W = W, method = "BFGS",
             control = list(reltol = 1e-12, maxit = 500))
  o$par
}
Smat <- function(b, H) {             # S = (1/n) sum g_i g_i'  (demeaned)
  g <- gmat(b, H); gc <- sweep(g, 2, colMeans(g)); crossprod(gc)/n
}
vcov_gmm <- function(b, H, W) {      # sandwich (DWD')^{-1}(DWSW'D')(DWD')^{-1}/n
  D <- t(Jac(b, H)); S <- Smat(b, H) # D here is k x q (t of q x k Jacobian)
  bread <- solve(D %*% W %*% t(D))
  bread %*% (D %*% W %*% S %*% W %*% t(D)) %*% bread / n
}

## ---- just-identified: h=(1,x) == Poisson score ----------------------------
H2 <- cbind(1, x)
b_just <- gmm_fit(H2, diag(2))
cat("=== just-identified GMM h=(1,x) (should equal Poisson QMLE) ===\n")
cat(sprintf("b0 = %.4f  b1 = %.4f\n", b_just[1], b_just[2]))

## ---- over-identified: h=(1,x,x^2) -----------------------------------------
H3 <- cbind(1, x, x^2)
# step 1: identity weight
b_id <- gmm_fit(H3, diag(3))
se_id <- sqrt(diag(vcov_gmm(b_id, H3, diag(3))))
# step 2: efficient weight W = S(b_id)^{-1}
Whes <- solve(Smat(b_id, H3))
b_eff <- gmm_fit(H3, Whes)
V_eff <- vcov_gmm(b_eff, H3, Whes)
se_eff <- sqrt(diag(V_eff))

cat("\n=== over-identified GMM h=(1,x,x^2): identity vs efficient weight ===\n")
cat(sprintf("identity   : b1 = %.4f  (SE %.4f)\n", b_id[2], se_id[2]))
cat(sprintf("efficient  : b1 = %.4f  (SE %.4f)\n", b_eff[2], se_eff[2]))
cat(sprintf("efficiency gain: SE ratio eff/identity (b1) = %.3f\n", se_eff[2]/se_id[2]))

## Hansen J at the efficient estimate: N * g' W g  ->  chi^2_{q-k}, df = 1
Jstat <- n * Qobj(b_eff, H3, solve(Smat(b_eff, H3)))
Jp <- pchisq(Jstat, df = 1, lower.tail = FALSE)
cat(sprintf("Hansen J = %.4f  (df=1)  p = %.4f  (mean correctly specified => no reject)\n",
            Jstat, Jp))

## ---- misspecified variant: TRUE mean has curvature exp(b0+b1 x + b2 x^2) ---
## but we fit the same 2-param model. The x^2 moment should now detect it.
set.seed(20)
b2_true <- 0.25
mu_mis  <- exp(oracle$b0 + oracle$b1*x + b2_true*x^2)
y_mis   <- rnbinom(n, size = oracle$theta_nb, mu = mu_mis)
ymem <- y; y <- y_mis                       # temporarily point moment fns at y_mis
b_mis  <- gmm_fit(H3, diag(3)); b_mis <- gmm_fit(H3, solve(Smat(b_mis, H3)))
Jmis <- n * Qobj(b_mis, H3, solve(Smat(b_mis, H3)))
Jmisp <- pchisq(Jmis, df = 1, lower.tail = FALSE)
y <- ymem
cat("\n=== misspecified mean (omitted x^2 curvature, b2=0.25): J should reject ===\n")
cat(sprintf("Hansen J = %.4f  (df=1)  p = %.6f\n", Jmis, Jmisp))

saveRDS(list(b_just=b_just, b_id=b_id, se_id=se_id, b_eff=b_eff, se_eff=se_eff,
             Jstat=Jstat, Jp=Jp, Jmis=Jmis, Jmisp=Jmisp), "res_gmm.rds")
