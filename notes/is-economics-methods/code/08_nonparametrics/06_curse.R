# ============================================================================
# 06_curse.R  --  the curse of dimensionality: local methods degrade fast as
# the number of regressors grows. (a) the edge-length a neighborhood must span
# to hold a fixed fraction of the data; (b) a kernel-regression MSE that blows
# up with dimension at fixed n.
# ============================================================================

## (a) edge length e_d(r) = r^(1/d) to capture fraction r of a unit hypercube
r <- 0.10
ds <- c(1, 2, 5, 10)
edge <- r^(1/ds)
cat("=== to capture 10% of the data, a cube must span this fraction per side ===\n")
for (i in seq_along(ds))
  cat(sprintf("  d = %2d :  edge = %.3f  (of the [0,1] range)\n", ds[i], edge[i]))
cat("  => in d=10 you must span 80% of every axis to grab 10% of the points:\n")
cat("     'local' is no longer local.\n")

## (b) kernel regression MSE at the center vs dimension, fixed n
## truth m(x) = sum_j sin(x_j); estimate m(0) by Gaussian-product-kernel NW.
set.seed(72)
n <- 2000; h <- 0.8; R <- 200
mse_dim <- sapply(c(1, 2, 3, 5), function(d) {
  err <- replicate(R, {
    X <- matrix(runif(n*d, -1, 1), n, d)
    m <- rowSums(sin(X)); Y <- m + rnorm(n, 0, 0.3)
    x0 <- rep(0, d)
    w <- exp(-0.5 * rowSums(sweep(X, 2, x0)^2) / h^2)   # product Gaussian kernel
    mhat <- sum(w * Y) / sum(w)                          # NW at 0; truth m(0)=0
    mhat^2
  })
  mean(err)
})
cat("\n=== NW regression MSE at x=0 vs dimension d (n = 2000 fixed) ===\n")
for (i in seq_along(c(1,2,3,5)))
  cat(sprintf("  d = %d :  MSE = %.4f\n", c(1,2,3,5)[i], mse_dim[i]))
cat(sprintf("  => MSE grows %.0fx from d=1 to d=5 at the same sample size.\n",
            mse_dim[4]/mse_dim[1]))

saveRDS(list(ds=ds, edge=edge, dims=c(1,2,3,5), mse_dim=mse_dim), "res_curse.rds")
