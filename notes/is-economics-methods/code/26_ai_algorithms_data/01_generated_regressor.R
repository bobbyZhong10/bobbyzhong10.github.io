# ------------------------------------------------------------------
# 01_generated_regressor.R -- the LLM-generated-regressor problem.
#
# A researcher wants the effect of a product's true review "sentiment"
# x_star on its sales y. x_star is latent. They run a large language
# model over each product's reviews to score sentiment, producing a
# NOISY measure x_hat = x_star + u (the LLM is imperfect). They then
# regress y on x_hat and read off the coefficient.
#
# Because x_hat is a mismeasured regressor, OLS is attenuated:
#   plim beta_hat = beta * reliability,  reliability = var(x*)/var(x_hat).
# And the naive standard error ignores the generation error, so it is
# too small (Pagan 1984 generated-regressors problem).
#
# Fix: a small GOLD-STANDARD validation subsample of size n_L where a
# human hand-labels the true x_star. From it we (a) estimate the
# reliability and correct the attenuation (errors-in-variables), and
# (b) form a prediction-powered / measurement-error-corrected estimate
# that combines the large LLM-scored sample with the small labeled set.
# ------------------------------------------------------------------
set.seed(22)

## ---- dimensions & truth -----------------------------------------
N     <- 4000        # products with LLM-scored sentiment
n_L   <- 400         # gold-standard hand-labeled subsample
b0    <- 1.0
b1    <- 2.0         # TRUE effect of sentiment on sales
sig_e <- 1.5         # outcome noise
sig_u <- 1.0         # LLM measurement-error sd (x_hat = x* + u)

## ---- generate latent sentiment, outcome, LLM measure ------------
x_star <- rnorm(N, 0, 1)
y      <- b0 + b1 * x_star + rnorm(N, 0, sig_e)
x_hat  <- x_star + rnorm(N, 0, sig_u)          # LLM score = truth + noise

reliab_true <- var(x_star) / (var(x_star) + sig_u^2)   # ~ 1/(1+1) = 0.5

## ---- (1) naive OLS on the LLM measure ---------------------------
m_naive <- lm(y ~ x_hat)
b_naive <- coef(m_naive)["x_hat"]
se_naive <- summary(m_naive)$coefficients["x_hat", "Std. Error"]

## ---- (2) gold-standard-only OLS (small labeled set) -------------
lab <- sample(N, n_L)
m_gold <- lm(y[lab] ~ x_star[lab])
b_gold <- coef(m_gold)[2]
se_gold <- summary(m_gold)$coefficients[2, "Std. Error"]

## ---- (3) errors-in-variables correction using validation --------
## Estimate reliability from the labeled set (where both x* and x_hat
## are observed): reliability = var(x*)/var(x_hat), estimated by the
## slope of x* on x_hat (regression-calibration).
lab_df <- data.frame(x_star = x_star[lab], x_hat = x_hat[lab])
reliab_hat <- coef(lm(x_star ~ x_hat, data = lab_df))[2]   # ~ var(x*)/var(x_hat)
b_eiv <- b_naive / reliab_hat

## ---- (4) measurement-error-corrected (regression calibration) ---
## Predict x* from x_hat using the labeled set, impute on the full
## sample, then regress y on the calibrated regressor. Equivalent to
## the EIV correction but framed as combining labeled + unlabeled data
## (the prediction-powered / post-prediction idea).
cal <- lm(x_star ~ x_hat, data = lab_df)
x_cal <- predict(cal, newdata = data.frame(x_hat = x_hat))
m_rc <- lm(y ~ x_cal)
b_rc <- coef(m_rc)["x_cal"]

## bootstrap SE for the corrected estimate (accounts for the
## generation/calibration step that the naive SE ignores)
set.seed(2201)
B <- 500; bs <- numeric(B)
for (b in 1:B) {
  ib <- sample(N, N, replace = TRUE)
  il <- sample(lab, n_L, replace = TRUE)
  rl <- coef(lm(x_star ~ x_hat, data = data.frame(x_star = x_star[il], x_hat = x_hat[il])))[2]
  bs[b] <- coef(lm(y[ib] ~ x_hat[ib]))[2] / rl
}
se_eiv <- sd(bs)

cat("== LLM-generated-regressor problem (true beta1 = 2.0) ==\n")
cat(sprintf("Reliability (true) = %.3f, estimated from labels = %.3f\n",
            reliab_true, reliab_hat))
cat(sprintf("(1) naive OLS on x_hat      beta = %.3f (SE %.3f)  [attenuated ~%.0f%%]\n",
            b_naive, se_naive, 100 * (1 - b_naive / b1)))
cat(sprintf("(2) gold-only OLS (n_L=%d)  beta = %.3f (SE %.3f)  [unbiased, imprecise]\n",
            n_L, b_gold, se_gold))
cat(sprintf("(3) EIV / reliability-corrected beta = %.3f (bootstrap SE %.3f)\n", b_eiv, se_eiv))
cat(sprintf("(4) regression-calibration      beta = %.3f\n", b_rc))
cat(sprintf("Naive SE (%.3f) vs corrected SE (%.3f): naive understates uncertainty\n",
            se_naive, se_eiv))

saveRDS(list(b1 = b1, reliab_true = reliab_true, reliab_hat = reliab_hat,
             b_naive = b_naive, se_naive = se_naive,
             b_gold = b_gold, se_gold = se_gold,
             b_eiv = b_eiv, se_eiv = se_eiv, b_rc = b_rc,
             N = N, n_L = n_L),
        file.path(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))),
          "res26_genreg.rds"))
cat("Saved res26_genreg.rds\n")
