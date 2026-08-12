# ------------------------------------------------------------------
# 03b_montecarlo.R -- sampling distribution of naive-ML vs DML, Mod 13.
# Repeats the Halcyon DGP R times and, on each draw, computes:
#   (A) naive plug-in : one forest of Y on (D,X), difference predictions
#                       (the "just use machine learning" estimator)
#   (B) DML           : 5-fold cross-fit Neyman-orthogonal PLR, with
#                       an influence-function 95% CI
# Reports bias, SD, RMSE, and CI coverage for the SATE.
# Lightweight learners (n and trees reduced) keep it fast.
# Writes mc_frame.rds (per-rep estimates) and res_mc.rds (summary).
# ------------------------------------------------------------------
suppressMessages({library(ranger)})

R  <- 300L
n  <- 4000L
p  <- 20L
Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
Lc <- chol(Sigma)

draw <- function(seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p) %*% Lc
  colnames(X) <- paste0("X", 1:p)
  X1 <- X[,1]; X2 <- X[,2]; X3 <- X[,3]
  h  <- 0.7 * X1 * X3 + 0.6 * (X2^2 - 1)
  e  <- pmin(pmax(plogis(0.40*X1 + 0.32*X3 + 0.50*h), 0.05), 0.95)
  D  <- rbinom(n, 1, e)
  g0 <- 0.20*X3 + 0.15*X1 + 0.50*h + 0.10*X[,6]
  tau <- 0.12 - 0.10*X1 - 0.06*X2 + 0.03*X1*X2
  Y  <- g0 + tau*D + rnorm(n, 0, 0.30)
  list(X = X, D = D, Y = Y, SATE = mean(tau))
}

rf <- function(y, Xtr, Xnew, seed)
  predict(ranger(y ~ ., data = data.frame(y = y, Xtr), num.trees = 300,
                 mtry = 8, min.node.size = 5, seed = seed),
          data.frame(Xnew))$predictions

est_naive <- function(d) {                       # plug-in, no orthogonality
  df <- data.frame(Y = d$Y, D = d$D, d$X)
  m  <- ranger(Y ~ ., data = df, num.trees = 300, mtry = 8,
               min.node.size = 5, seed = 7)
  mean(predict(m, data.frame(D = 1, d$X))$predictions -
       predict(m, data.frame(D = 0, d$X))$predictions)
}

est_dml <- function(d, K = 5) {                  # cross-fit orthogonal PLR
  nn <- length(d$Y); fold <- sample(rep(1:K, length.out = nn))
  lh <- mh <- numeric(nn)
  for (k in 1:K) {
    tr <- fold != k; te <- fold == k
    lh[te] <- rf(d$Y[tr], d$X[tr, ], d$X[te, ], 100 + k)
    mh[te] <- rf(d$D[tr], d$X[tr, ], d$X[te, ], 200 + k)
  }
  Yt <- d$Y - lh; Dt <- d$D - mh
  th <- sum(Dt * Yt) / sum(Dt * Dt)
  psi_b <- Dt * (Yt - th * Dt); J <- mean(-(Dt^2))
  se <- sqrt(mean(psi_b^2) / J^2 / nn)
  c(est = th, se = se)
}

naive <- dmlv <- dmlse <- sate <- numeric(R)
for (r in 1:R) {
  d <- draw(5000 + r)
  sate[r]  <- d$SATE
  naive[r] <- est_naive(d)
  b <- est_dml(d); dmlv[r] <- b["est"]; dmlse[r] <- b["se"]
  if (r %% 50 == 0) cat("rep", r, "\n")
}

cover <- mean(abs(dmlv - sate) <= 1.96 * dmlse)
summ <- data.frame(
  method = c("naive plug-in", "DML cross-fit"),
  bias   = c(mean(naive - sate), mean(dmlv - sate)),
  sd     = c(sd(naive), sd(dmlv)),
  rmse   = c(sqrt(mean((naive - sate)^2)), sqrt(mean((dmlv - sate)^2)))
)
res_mc <- list(summ = summ, coverage = cover,
               mean_sate = mean(sate), R = R, n = n)
saveRDS(data.frame(naive = naive, dml = dmlv, dml_se = dmlse, sate = sate),
        "mc_frame.rds")
saveRDS(res_mc, "res_mc.rds")

cat(sprintf("\nMonte Carlo: R = %d draws, n = %d each, mean SATE = %.4f\n",
            R, n, mean(sate)))
summ_num <- summ; summ_num[-1] <- round(summ_num[-1], 4)
print(summ_num)
cat(sprintf("\nDML 95%% CI coverage of SATE = %.3f\n", cover))
