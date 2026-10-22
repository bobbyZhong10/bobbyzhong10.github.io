# ------------------------------------------------------------------
# 02_naive.R -- The estimators that fail, for Module 13.
#   (a) naive difference in means           (confounding)
#   (b) OLS with linear controls            (misspecified nuisance)
#   (c) single-ML plug-in (no orthogonality, no cross-fitting)
#       -> regularization bias
# Writes res_naive.rds
# ------------------------------------------------------------------
set.seed(1313)
suppressMessages({library(ranger)})

dat   <- readRDS("halcyon.rds")
truth <- readRDS("truth.rds")
Xnames <- paste0("X", 1:truth$p)
X <- as.matrix(dat[, Xnames]); Y <- dat$Y; D <- dat$D

## ---- (a) naive difference in means -------------------------------
dim_est <- mean(Y[D == 1]) - mean(Y[D == 0])
dim_se  <- sqrt(var(Y[D == 1]) / sum(D == 1) + var(Y[D == 0]) / sum(D == 0))

## ---- (b) OLS with linear controls --------------------------------
ols <- lm(Y ~ D + ., data = data.frame(Y, D, X))
ols_est <- coef(ols)["D"]
ols_se  <- summary(ols)$coefficients["D", "Std. Error"]

## ---- (c) single-ML plug-in ---------------------------------------
# Fit ONE random forest for the outcome on (D, X) on the full sample,
# then read off the treatment "effect" by prediction differencing.
# This is the tempting-but-wrong plug-in: no orthogonalization, no
# sample splitting. Regularization in the forest biases the effect.
rf_full <- ranger(Y ~ ., data = data.frame(Y, D, X), num.trees = 1000,
                  mtry = 8, min.node.size = 5, seed = 1313)
d1 <- data.frame(D = 1, X); d0 <- data.frame(D = 0, X)
plugin_est <- mean(predict(rf_full, d1)$predictions -
                   predict(rf_full, d0)$predictions)

## ---- a second plug-in flavor: residualize with ML but no crossfit
# g-hat from RF of Y on X (ignoring D), then OLS of (Y - ghat) on D.
# This is the "naive FWL with ML" the DML slides warn against.
rf_g <- ranger(Y ~ ., data = data.frame(Y, X), num.trees = 1000,
               mtry = 8, min.node.size = 5, seed = 1313)
ghat <- predict(rf_g, data.frame(X))$predictions
naive_fwl <- coef(lm(I(Y - ghat) ~ D))["D"]

res_naive <- list(
  dim = c(est = dim_est, se = dim_se),
  ols = c(est = unname(ols_est), se = unname(ols_se)),
  plugin = c(est = plugin_est),
  naive_fwl = c(est = unname(naive_fwl)),
  ATE_true = truth$ATE, ATT_true = truth$ATT
)
saveRDS(res_naive, "res_naive.rds")

cat(sprintf("True SATE                 = %.4f\n", truth$ATE))
cat(sprintf("True ATT                  = %.4f\n", truth$ATT))
cat(sprintf("(a) diff in means         = %.4f (SE %.4f)\n", dim_est, dim_se))
cat(sprintf("(b) OLS linear controls   = %.4f (SE %.4f)\n", ols_est, ols_se))
cat(sprintf("(c) single-RF plug-in     = %.4f\n", plugin_est))
cat(sprintf("(c') naive ML-FWL no CF   = %.4f\n", naive_fwl))
