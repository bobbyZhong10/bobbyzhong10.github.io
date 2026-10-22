# ------------------------------------------------------------------
# 03_sur.R -- block S2: seemingly unrelated regressions.
#
# Two outcome equations for the same 1,500 markets: the average basket (USD) and orders per active
# buyer, with correlated errors.
#   same_X: both equations on (1, free_ship_min, premium_share, loyalty). FGLS equals equation-by-equation
#           OLS exactly (Wooldridge 2010 Theorem 7.6), but the joint covariance of the two equations'
#           estimates still matters: the first-order effect of the free-shipping minimum on revenue per
#           buyer, E[orders] beta_ship,basket + E[basket] beta_ship,orders, needs Cov(beta_hat_1, beta_hat_2).
#   diff_X: each equation on its own regressors (the true exclusions). FGLS uses the other equation's
#           residuals and is more efficient; with the same Omega_hat, its standard errors are below the
#           OLS ones for every coefficient (Gauss-Markov with Omega_hat known; the OLS variance of
#           equation m is Omega_hat_mm (X_m'X_m)^-1, the same divisor n as in Omega_hat).
# sur_fgls is checked against systemfit(method = "SUR", methodResidCov = "noDfCor").
#
# Output: res_sur.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("iv_fns.R")
suppressPackageStartupMessages(library(systemfit))

tr <- readRDS("truth09.rds")
S2 <- readRDS("bramble_S2.rds")
y <- list(basket = S2$basket_usd, orders = S2$orders_pb)
X_same <- cbind(const = 1, free_ship_min = S2$free_ship_min, premium_share = S2$premium_share, loyalty = S2$loyalty)
X_diff <- list(basket = X_same[, c("const", "free_ship_min", "premium_share")],
               orders = X_same[, c("const", "free_ship_min", "loyalty")])

fit_spec <- function(Xlist) {
  s <- sur_fgls(y, Xlist)
  o <- mapply(function(yy, X) ols_fit(yy, X, vcov = "homo"), y, Xlist, SIMPLIFY = FALSE)
  b_ols <- unlist(lapply(o, `[[`, "b")); names(b_ols) <- names(s$b)
  # OLS variances with the same Omega_hat (divisor n): Omega_hat_mm (X_m'X_m)^-1
  se_ols <- unlist(lapply(seq_along(Xlist), function(m) sqrt(s$Sigma[m, m] * diag(solve(crossprod(Xlist[[m]]))))))
  se_lm <- unlist(lapply(o, function(f) sqrt(diag(f$V))))           # conventional lm standard errors (n - k_m)
  names(se_ols) <- names(se_lm) <- names(s$b)
  list(sur = s$b, ols = b_ols, se_sur = sqrt(diag(s$V)), se_ols = se_ols, se_lm = se_lm, V = s$V, Sigma = s$Sigma,
       rho_hat = s$Sigma[1, 2] / sqrt(s$Sigma[1, 1] * s$Sigma[2, 2]))
}
same_X <- fit_spec(list(basket = X_same, orders = X_same))
diff_X <- fit_spec(X_diff)
same_X$max_abs_sur_minus_ols <- max(abs(same_X$sur - same_X$ols))
diff_X$se_ratio <- diff_X$se_sur / diff_X$se_ols
diff_X$sur_minus_ols <- diff_X$sur - diff_X$ols

## package checks
fs_same <- systemfit(list(basket = basket_usd ~ free_ship_min + premium_share + loyalty,
                          orders = orders_pb ~ free_ship_min + premium_share + loyalty),
                     method = "SUR", data = S2, methodResidCov = "noDfCor")
fs_diff <- systemfit(list(basket = basket_usd ~ free_ship_min + premium_share, orders = orders_pb ~ free_ship_min + loyalty),
                     method = "SUR", data = S2, methodResidCov = "noDfCor")
pkg <- list(same_max_abs = max(abs(same_X$sur - coef(fs_same))), diff_max_abs = max(abs(diff_X$sur - coef(fs_diff))),
            diff_V_rel = max(abs(diff_X$V - vcov(fs_diff))) / max(abs(vcov(fs_diff))))
stopifnot(pkg$same_max_abs < 1e-10, pkg$diff_max_abs < 1e-10, pkg$diff_V_rel < 1e-8)

## coverage of the true coefficients by the diff_X SUR intervals (single sample)
bt <- tr$S2$beta
covers_diff <- vapply(seq_along(bt), function(j) covers_num(ci95(diff_X$sur[j], diff_X$se_sur[j]), bt[[j]]), numeric(1))
names(covers_diff) <- names(bt)

## revenue per buyer: first-order effect of USD 1 on the free-shipping minimum,
## w' beta with w = mean(orders) on beta_ship,basket and mean(basket) on beta_ship,orders
lin_comb <- function(fit) {
  nm <- names(fit$sur)
  w <- setNames(numeric(length(nm)), nm)
  w["basket:free_ship_min"] <- mean(S2$orders_pb); w["orders:free_ship_min"] <- mean(S2$basket_usd)
  V_blk <- fit$V; i1 <- grep("^basket:", nm); i2 <- grep("^orders:", nm)
  V_blk[i1, i2] <- 0; V_blk[i2, i1] <- 0
  est <- sum(w * fit$sur); se_full <- sqrt(drop(t(w) %*% fit$V %*% w)); se_blk <- sqrt(drop(t(w) %*% V_blk %*% w))
  list(est = est, se_full = se_full, se_blockdiag = se_blk, se_ratio = se_full / se_blk, ci = ci95(est, se_full),
       covers_truth = covers_num(ci95(est, se_full), tr$S2$revenue_effect),
       cov_ship = fit$V["basket:free_ship_min", "orders:free_ship_min"],
       corr_ship = fit$V["basket:free_ship_min", "orders:free_ship_min"] /
         sqrt(fit$V["basket:free_ship_min", "basket:free_ship_min"] * fit$V["orders:free_ship_min", "orders:free_ship_min"]),
       w_orders = mean(S2$orders_pb), w_basket = mean(S2$basket_usd))
}
revenue <- list(same_X = lin_comb(same_X), diff_X = lin_comb(diff_X))

res <- list(same_X = same_X, diff_X = diff_X, pkg = pkg, covers_diff = covers_diff, revenue = revenue,
            sample = list(mean_basket = mean(S2$basket_usd), mean_orders = mean(S2$orders_pb), share_loyalty = mean(S2$loyalty),
                          mean_premium = mean(S2$premium_share), n = nrow(S2)))
saveRDS(res, "res_sur.rds")

cat(sprintf("S2 same regressors: max |SUR - OLS| = %.2e (Theorem 7.6); rho_hat %.3f\n", same_X$max_abs_sur_minus_ols, same_X$rho_hat))
cat("S2 different regressors: SE ratio SUR / OLS (same Omega_hat):\n"); print(round(diff_X$se_ratio, 4))
cat(sprintf("Revenue effect (same X): %.4f, se with cross-equation covariance %.4f, ignoring it %.4f (ratio %.3f); truth %.4f\n",
            revenue$same_X$est, revenue$same_X$se_full, revenue$same_X$se_blockdiag, revenue$same_X$se_ratio, tr$S2$revenue_effect))
cat("Saved res_sur.rds\n")
