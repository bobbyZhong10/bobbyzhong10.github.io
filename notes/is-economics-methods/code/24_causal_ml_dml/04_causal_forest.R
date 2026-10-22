# ------------------------------------------------------------------
# 04_causal_forest.R -- CATE with a causal forest (Wager-Athey), Mod 13.
#   - grf causal_forest: honest, out-of-bag CATE
#   - AIPW ATE (doubly robust)
#   - calibration test (test_calibration)
#   - best_linear_projection: recover the oracle BLP of tau on (X1,X2)
#   - GATES: AIPW group ATE by CATE quintile
#   - CLAN: covariate means in top vs bottom CATE group
#   - RATE / TOC (rank_average_treatment_effect)
# Writes res_cf.rds
# ------------------------------------------------------------------
set.seed(1313)
suppressMessages(library(grf))

dat   <- readRDS("halcyon.rds")
truth <- readRDS("truth.rds")
Xnames <- paste0("X", 1:truth$p)
X <- as.matrix(dat[, Xnames]); Y <- dat$Y; D <- dat$D

## ---- honest causal forest ----------------------------------------
cf <- causal_forest(X, Y, D, num.trees = 4000, honesty = TRUE,
                    tune.parameters = "all", seed = 1313)
tau_hat <- predict(cf)$predictions          # out-of-bag CATE

## ---- doubly robust ATE (AIPW) ------------------------------------
ate <- average_treatment_effect(cf, target.sample = "all")

## ---- calibration test --------------------------------------------
cal <- test_calibration(cf)

## ---- best linear projection of CATE on (X1, X2) ------------------
blp <- best_linear_projection(cf, X[, c(1, 2)])

## ---- quality of CATE recovery vs oracle --------------------------
tau_true_i <- dat$tau_true
cate_cor  <- cor(tau_hat, tau_true_i)
cate_rmse <- sqrt(mean((tau_hat - tau_true_i)^2))

## ---- GATES: AIPW scores, groups by CATE quintile -----------------
# doubly robust scores Gamma_i (AIPW), then average within CATE quintile
DR <- get_scores(cf)                          # AIPW/doubly-robust scores
q  <- cut(tau_hat, quantile(tau_hat, seq(0, 1, .2)),
          include.lowest = TRUE, labels = 1:5)
gates <- sapply(1:5, function(g) {
  gi <- q == g
  c(est = mean(DR[gi]), se = sd(DR[gi]) / sqrt(sum(gi)),
    tau_true = mean(tau_true_i[gi]))
})
gates <- t(gates)

## ---- CLAN: covariate means, top vs bottom CATE quintile ----------
clan <- rbind(
  bottom = colMeans(X[q == 1, c(1, 2, 3)]),
  top    = colMeans(X[q == 5, c(1, 2, 3)])
)
colnames(clan) <- c("X1_size", "X2_tenure", "X3_momentum")

## ---- RATE (TOC-based) : is targeting by CATE valuable? -----------
rate <- rank_average_treatment_effect(cf, tau_hat)

res_cf <- list(
  ate = ate, cal = cal, blp = blp,
  blp_oracle = truth$blp_oracle,
  cate_cor = cate_cor, cate_rmse = cate_rmse,
  gates = gates, clan = clan,
  rate = c(estimate = rate$estimate, se = rate$std.err),
  ATE_true = truth$ATE
)
## ---- doubly robust score matrix for policy learning -------------
gamma_dr <- policytree::double_robust_scores(cf)   # n x 2: control, treat

saveRDS(res_cf, "res_cf.rds")
saveRDS(data.frame(tau_hat = tau_hat, tau_true = tau_true_i,
                   X1 = X[, 1], X2 = X[, 2], quint = as.integer(q)),
        "cate_frame.rds")
saveRDS(list(gamma = gamma_dr, X = X, tau_true = tau_true_i,
             tau_hat = tau_hat), "dr_scores.rds")

cat(sprintf("AIPW ATE = %.4f (SE %.4f)  [true %.4f]\n",
            ate[1], ate[2], truth$ATE))
cat(sprintf("CATE cor with truth = %.3f ; RMSE = %.4f\n",
            cate_cor, cate_rmse))
cat("Calibration test:\n"); print(round(cal, 4))
cat("\nBest linear projection of CATE on (X1,X2):\n"); print(blp)
cat(sprintf("Oracle BLP: b0=%.4f b1=%.4f b2=%.4f\n",
            truth$blp_oracle[1], truth$blp_oracle[2], truth$blp_oracle[3]))
cat("\nGATES (AIPW group ATE by CATE quintile):\n")
print(round(gates, 4))
cat("\nCLAN (covariate means bottom vs top CATE quintile):\n")
print(round(clan, 3))
cat(sprintf("\nRATE (AUTOC) = %.4f (SE %.4f)\n", rate$estimate, rate$std.err))
