# ============================================================================
# 04_series_cv.R  --  the bias-variance tradeoff and how cross-validation picks
# the smoothing. Undersmoothing (small h) -> low bias, high variance; over-
# smoothing (large h) -> high bias, low variance; MSE is U-shaped in h. CV finds
# the bottom. Series/splines tell the same story with the number of knots.
# ============================================================================

suppressMessages(library(np))
S <- readRDS("solstice.rds"); truth <- S$truth
d <- S$demand; grid <- S$grid; g <- S$g; sigma <- S$sigma

Kg <- function(u) dnorm(u)
ll_fit <- function(x0, X, Y, h){ w<-Kg((x0-X)/h); Xm<-cbind(1,X-x0)
  solve(t(Xm)%*%(Xm*w), t(Xm)%*%(w*Y))[1] }

## ---- bias-variance decomposition over a grid of bandwidths -----------------
## at a set of interior evaluation points, Monte-Carlo the estimator and split
## integrated MSE into (bias^2) + (variance).
set.seed(31)
epts <- seq(2, 8, length = 25)
hs   <- c(0.15, 0.25, 0.4, 0.6, 0.9, 1.3, 2.0)
R    <- 250
bv <- sapply(hs, function(h){
  M <- t(replicate(R, {
    Xd <- runif(600,1,9); Yd <- g(Xd)+rnorm(600,0,sigma)
    sapply(epts, function(x0) ll_fit(x0, Xd, Yd, h))
  }))                                              # R x |epts|
  mean_est <- colMeans(M)
  bias2 <- mean((mean_est - g(epts))^2)
  varc  <- mean(apply(M, 2, var))
  c(bias2 = bias2, variance = varc, mse = bias2 + varc)
})
colnames(bv) <- hs
cat("=== bias-variance vs bandwidth h (integrated over interior prices) ===\n")
print(round(bv, 5))
h_opt <- hs[which.min(bv["mse",])]
cat(sprintf("MSE-minimizing h on this grid = %.2f\n", h_opt))

## ---- cross-validation picks the bandwidth automatically --------------------
bw <- npregbw(y ~ x, data = d, regtype = "ll", bwmethod = "cv.ls",
              ckertype = "gaussian")
cat(sprintf("least-squares CV-optimal bandwidth = %.3f\n", bw$bw))

## ---- series / splines: number of knots is the same dial -------------------
## natural cubic spline with df chosen by GCV vs deliberately over/under-smoothed
sp_cv  <- smooth.spline(d$x, d$y)                   # GCV picks df
fit_df <- function(df) predict(smooth.spline(d$x, d$y, df = df), grid)$y
mse_grid <- function(f) mean((f - truth$g_grid)^2)
cat(sprintf("\nsmoothing spline GCV effective df = %.1f\n", sp_cv$df))
cat(sprintf("MSE vs truth: df=3 (oversmooth)=%.5f  GCV df=%.1f=%.5f  df=30 (overfit)=%.5f\n",
            mse_grid(fit_df(3)), sp_cv$df,
            mse_grid(predict(sp_cv, grid)$y), mse_grid(fit_df(30))))

saveRDS(list(bv=bv, hs=hs, h_opt=h_opt, cv_bw=bw$bw, epts=epts,
             gcv_df=sp_cv$df, sp_grid=predict(sp_cv, grid)$y,
             mse_over=mse_grid(fit_df(3)), mse_gcv=mse_grid(predict(sp_cv,grid)$y),
             mse_under=mse_grid(fit_df(30))), "res_cv.rds")
