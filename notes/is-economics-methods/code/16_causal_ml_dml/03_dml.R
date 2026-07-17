# ------------------------------------------------------------------
# 03_dml.R -- Double/Debiased ML for the ATE, Module 13.
#   (1) hand-rolled PLR: cross-fit Neyman-orthogonal residual-on-residual
#       - contrast NO cross-fitting vs 5-fold cross-fitting
#   (2) DoubleML package: PLR and IRM (AIPW) with honest SEs
# Writes res_dml.rds
# ------------------------------------------------------------------
set.seed(1313)
suppressMessages({
  library(ranger); library(data.table)
})

dat   <- readRDS("halcyon.rds")
truth <- readRDS("truth.rds")
Xnames <- paste0("X", 1:truth$p)
X <- as.matrix(dat[, Xnames]); Y <- dat$Y; D <- dat$D
n <- nrow(dat)

## overlap-weighted true effect (the PLR estimand under heterogeneity):
w  <- dat$e_true * (1 - dat$e_true)
theta_plr_true <- sum(w * dat$tau_true) / sum(w)

rf_fit <- function(y, Xtr, Xnew, seed = 1313) {
  df <- data.frame(y = y, Xtr)
  m  <- ranger(y ~ ., data = df, num.trees = 1000, mtry = 8,
               min.node.size = 5, seed = seed)
  predict(m, data.frame(Xnew))$predictions
}

## ---- (1a) orthogonal moment, NO cross-fitting ---------------------
# Fit l(X)=E[Y|X], m(X)=E[D|X] on the SAME data used to form residuals.
lhat_nocf <- rf_fit(Y, X, X)
mhat_nocf <- rf_fit(D, X, X)
Yt <- Y - lhat_nocf; Dt <- D - mhat_nocf
theta_nocf <- sum(Dt * Yt) / sum(Dt * Dt)

## ---- (1b) orthogonal moment WITH 5-fold cross-fitting ------------
K <- 5
fold <- sample(rep(1:K, length.out = n))
lhat <- mhat <- numeric(n)
for (k in 1:K) {
  tr <- fold != k; te <- fold == k
  lhat[te] <- rf_fit(Y[tr], X[tr, ], X[te, ], seed = 1000 + k)
  mhat[te] <- rf_fit(D[tr], X[tr, ], X[te, ], seed = 2000 + k)
}
Yt_cf <- Y - lhat; Dt_cf <- D - mhat
theta_cf <- sum(Dt_cf * Yt_cf) / sum(Dt_cf * Dt_cf)
# influence-function SE for the PLR estimator
psi_a <- -(Dt_cf^2)
psi_b <- Dt_cf * (Yt_cf - theta_cf * Dt_cf)
J <- mean(psi_a)
se_cf <- sqrt(mean(psi_b^2) / J^2 / n)

## ---- (2) DoubleML package: PLR and IRM (AIPW) --------------------
have_dml <- requireNamespace("DoubleML", quietly = TRUE) &&
            requireNamespace("mlr3", quietly = TRUE) &&
            requireNamespace("mlr3learners", quietly = TRUE)
dml_plr <- dml_irm <- NULL
if (have_dml) {
  suppressMessages({
    library(DoubleML); library(mlr3); library(mlr3learners)
    lgr::get_logger("mlr3")$set_threshold("error")
  })
  df <- data.frame(Y = Y, D = D, X)
  dml_data <- DoubleMLData$new(df, y_col = "Y", d_cols = "D",
                               x_cols = Xnames)
  ml_g <- lrn("regr.ranger", num.trees = 1000, mtry = 8, min.node.size = 5)
  ml_m <- lrn("regr.ranger", num.trees = 1000, mtry = 8, min.node.size = 5)
  ml_m_cl <- lrn("classif.ranger", num.trees = 1000, mtry = 8,
                 min.node.size = 5, predict_type = "prob")

  set.seed(1313)
  plr <- DoubleMLPLR$new(dml_data, ml_l = ml_g, ml_m = ml_m,
                         n_folds = 5, n_rep = 1)
  plr$fit()
  dml_plr <- c(est = plr$coef[["D"]], se = plr$se[["D"]])

  set.seed(1313)
  irm <- DoubleMLIRM$new(dml_data, ml_g = ml_g, ml_m = ml_m_cl,
                         n_folds = 5, n_rep = 1,
                         trimming_threshold = 0.02)
  irm$fit()
  dml_irm <- c(est = irm$coef[["D"]], se = irm$se[["D"]])
}

res_dml <- list(
  theta_plr_true = theta_plr_true,
  ATE_true = truth$ATE,
  nocf = theta_nocf,
  cf   = c(est = theta_cf, se = se_cf),
  dml_plr = dml_plr,
  dml_irm = dml_irm
)
saveRDS(res_dml, "res_dml.rds")

cat(sprintf("True SATE                        = %.4f\n", truth$ATE))
cat(sprintf("True overlap-weighted (PLR) eff. = %.4f\n", theta_plr_true))
cat(sprintf("orthogonal, NO cross-fit         = %.4f\n", theta_nocf))
cat(sprintf("orthogonal, 5-fold cross-fit     = %.4f (SE %.4f)\n",
            theta_cf, se_cf))
if (have_dml) {
  cat(sprintf("DoubleML PLR                     = %.4f (SE %.4f)\n",
              dml_plr["est"], dml_plr["se"]))
  cat(sprintf("DoubleML IRM (AIPW, ATE)         = %.4f (SE %.4f)\n",
              dml_irm["est"], dml_irm["se"]))
} else cat("DoubleML not available\n")
