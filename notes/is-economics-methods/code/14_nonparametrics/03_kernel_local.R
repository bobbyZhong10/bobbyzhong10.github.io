# ============================================================================
# 03_kernel_local.R  --  Nadaraya-Watson (local constant) vs local linear.
# Hand-rolled so the weighted-average and locally-weighted-least-squares forms
# are visible. The single-dataset fits feed the figure; a Monte-Carlo isolates
# the BIAS (a repeated-sampling property), exposing NW's boundary bias which
# local linear removes. Local linear also delivers the slope (price sensitivity).
# ============================================================================

S <- readRDS("solstice.rds"); truth <- S$truth
d <- S$demand; grid <- S$grid; g <- S$g; sigma <- S$sigma

Kg <- function(u) dnorm(u)                          # Gaussian kernel
nw <- function(x0, X, Y, h) { w <- Kg((x0-X)/h); sum(w*Y)/sum(w) }
ll <- function(x0, X, Y, h) {                       # returns c(m_hat, m'_hat)
  w <- Kg((x0-X)/h); Xm <- cbind(1, X - x0)
  as.numeric(solve(t(Xm) %*% (Xm*w), t(Xm) %*% (w*Y)))
}
h <- 0.6

## ---- single-dataset fitted curves (for the figure) ------------------------
mhat_nw <- sapply(grid, function(x0) nw(x0, d$x, d$y, h))
mhat_ll <- sapply(grid, function(x0) ll(x0, d$x, d$y, h)[1])
mse <- function(fit) mean((fit - truth$g_grid)^2)
cat("=== single fit (h = 0.6): overall MSE vs the known curve ===\n")
cat(sprintf("Nadaraya-Watson MSE = %.5f   local-linear MSE = %.5f\n",
            mse(mhat_nw), mse(mhat_ll)))

## local-linear slope recovers the varying price sensitivity
sl <- sapply(truth$prices, function(x0) ll(x0, d$x, d$y, h)[2])
cat("local-linear slope at prices 2 / 5.5 / 8 = ")
cat(sprintf("%.3f / %.3f / %.3f  (truth %.3f / %.3f / %.3f)\n",
            sl[1], sl[2], sl[3], truth$slope[1], truth$slope[2], truth$slope[3]))

## ---- Monte-Carlo bias: NW boundary bias vs local-linear fix ---------------
set.seed(9)
bias_at <- function(x0, R = 400) {
  est <- t(replicate(R, {
    X <- runif(600, 1, 9); Y <- g(X) + rnorm(600, 0, sigma)
    c(nw(x0, X, Y, h), ll(x0, X, Y, h)[1])
  }))
  c(nw = mean(est[,1]) - g(x0), ll = mean(est[,2]) - g(x0))
}
cat("\n=== Monte-Carlo bias (400 draws): NW vs local linear ===\n")
bt <- sapply(c(1.2, 5.0, 8.8), bias_at)
colnames(bt) <- c("x=1.2 (left edge)", "x=5.0 (interior)", "x=8.8 (right edge)")
print(round(bt, 4))
cat("=> NW is biased at the boundaries; local linear is not.\n")

saveRDS(list(mhat_nw=mhat_nw, mhat_ll=mhat_ll, h=h,
             mse_nw=mse(mhat_nw), mse_ll=mse(mhat_ll),
             ll_slopes=sl, bias_tbl=bt), "res_kernel.rds")
