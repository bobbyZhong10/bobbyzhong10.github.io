# ============================================================================
# 05_robinson.R  --  semiparametrics: the partially linear model
#   y = beta*z + g(x) + eps,   z (ad exposure) correlated with price x.
# A linear-in-price model biases beta because g is nonlinear. Robinson's (1988)
# double-residual estimator partials x out of y and z nonparametrically, then
# regresses the residuals, recovering beta at root-n despite g being unknown.
# ============================================================================

S <- readRDS("solstice.rds"); truth <- S$truth
d <- S$partial

## naive 1: regress y on z only (omits the price effect entirely)
b_naive1 <- coef(lm(y ~ z, data = d))["z"]
## naive 2: y on z and a LINEAR price term (wrong functional form for g)
b_naive2 <- coef(lm(y ~ z + x, data = d))["z"]
## naive 3: y on z and a cubic-in-price (better, but still parametric)
b_naive3 <- coef(lm(y ~ z + poly(x, 3), data = d))["z"]

## Robinson double residual: partial x out of BOTH y and z nonparametrically
##   e_y = y - E[y|x],  e_z = z - E[z|x];  beta = regression of e_y on e_z
m_yx <- smooth.spline(d$x, d$y)                    # E[y|x] nonparametric
m_zx <- smooth.spline(d$x, d$z)                    # E[z|x] nonparametric
ey <- d$y - predict(m_yx, d$x)$y
ez <- d$z - predict(m_zx, d$x)$y
b_rob <- coef(lm(ey ~ ez - 1))["ez"]
se_rob <- sqrt(vcov(lm(ey ~ ez - 1))["ez","ez"])

cat("=== partially linear: beta on ad exposure (truth 0.80) ===\n")
cat(sprintf("naive OLS  y ~ z           = %.4f  (omits price; biased)\n", b_naive1))
cat(sprintf("naive OLS  y ~ z + x       = %.4f  (linear-in-price; still biased)\n", b_naive2))
cat(sprintf("naive OLS  y ~ z + poly3   = %.4f  (cubic price; better)\n", b_naive3))
cat(sprintf("Robinson double-residual   = %.4f (SE %.4f)  (consistent)\n", b_rob, se_rob))

## Monte-Carlo: Robinson is root-n consistent and ~unbiased for beta
set.seed(55)
g <- S$g; sigma <- S$sigma
mz <- function(x) 0.15*x + 1.1*plogis(1.8*(x-5))
mc <- replicate(300, {
  X <- runif(600,1,9); Z <- mz(X)+rnorm(600,0,0.7); Y <- 0.8*Z + g(X) + rnorm(600,0,sigma)
  dd <- data.frame(x=X, z=Z, y=Y)
  eyr <- dd$y - predict(smooth.spline(dd$x, dd$y), dd$x)$y
  ezr <- dd$z - predict(smooth.spline(dd$x, dd$z), dd$x)$y
  c(lin = coef(lm(y ~ z + x, dd))["z"], rob = coef(lm(eyr ~ ezr - 1))[1])
})
cat(sprintf("\nMC (300 draws): mean beta  linear-in-x = %.4f  Robinson = %.4f  (truth 0.80)\n",
            mean(mc["lin.z",]), mean(mc[2,])))
cat(sprintf("MC SD of Robinson beta = %.4f  (root-n precise)\n", sd(mc[2,])))

saveRDS(list(b_naive1=b_naive1, b_naive2=b_naive2, b_naive3=b_naive3,
             b_rob=b_rob, se_rob=se_rob,
             mc_lin=mean(mc["lin.z",]), mc_rob=mean(mc[2,]), mc_sd=sd(mc[2,])),
        "res_robinson.rds")
