# ------------------------------------------------------------------
# 04_peeking.R -- continuous monitoring, peeking, and always-valid
#                 inference (mSPRT), Module 14.
#   - naive fixed-0.05 peeking: type-I error inflates with # looks
#   - mSPRT (Johari-Pekelis-Walsh-Koomen): an always-valid test whose
#     type-I error stays <= alpha under continuous monitoring
#   - power / expected stopping time of mSPRT under the true effect
# Writes res_peek.rds
# ------------------------------------------------------------------
set.seed(1414)
p     <- readRDS("lumen_params.rds")
gen   <- p$gen; sigma <- p$sd_Y0; delta <- p$delta
alpha <- 0.05

n_max   <- 40000L          # per-arm sample at the final look
n_looks <- 14L             # one look per day for two weeks
looks   <- round(seq(n_max / n_looks, n_max, length.out = n_looks))

## mSPRT for a two-sample normal mean difference, mixing N(0, tau2). ---
tau2 <- 0.01
msprt_lambda <- function(dbar, n_arm, s2 = sigma^2, tau2 = 0.01) {
  vD <- 2 * s2 / n_arm                        # variance of the mean diff
  sqrt(vD / (vD + tau2)) * exp(tau2 * dbar^2 / (2 * vD * (vD + tau2)))
}

## ---- naive peeking type-I error vs number of looks (A/A) ---------
naive_fpr <- function(K, R = 4000) {
  lk <- round(seq(n_max / K, n_max, length.out = K))
  set.seed(42); crossed <- logical(R)
  for (r in 1:R) {
    d <- gen(2 * n_max)                       # A/A: both arms ~ Y0
    yT <- d$Y0[1:n_max]; yC <- d$Y0[(n_max + 1):(2 * n_max)]
    hit <- FALSE
    for (nk in lk) {
      db <- mean(yT[1:nk]) - mean(yC[1:nk])
      se <- sqrt(var(yT[1:nk]) / nk + var(yC[1:nk]) / nk)
      if (abs(db / se) > qnorm(1 - alpha / 2)) { hit <- TRUE; break }
    }
    crossed[r] <- hit
  }
  mean(crossed)
}
naive_curve <- sapply(c(1, 2, 5, 10, 14), naive_fpr)

## ---- mSPRT type-I error under continuous monitoring (A/A) --------
msprt_run <- function(R = 4000, effect = 0, seed = 99) {
  set.seed(seed)
  reject <- logical(R); stop_n <- rep(n_max, R)
  for (r in 1:R) {
    d <- gen(2 * n_max)
    yT <- d$Y0[1:n_max] + effect; yC <- d$Y0[(n_max + 1):(2 * n_max)]
    hit <- FALSE
    for (nk in looks) {
      db <- mean(yT[1:nk]) - mean(yC[1:nk])
      if (msprt_lambda(db, nk, tau2 = tau2) >= 1 / alpha) {
        hit <- TRUE; stop_n[r] <- nk; break
      }
    }
    reject[r] <- hit
  }
  list(rej = mean(reject), stop_med = median(stop_n))
}
msprt_aa <- msprt_run(effect = 0)               # type-I error
msprt_ab <- msprt_run(effect = delta, seed = 100)  # power + stopping

## for reference: naive fixed-0.05 peeking rejection under the effect
naive_ab_pow <- {
  set.seed(101); R <- 2000; rej <- logical(R)
  for (r in 1:R) {
    d <- gen(2 * n_max)
    yT <- d$Y0[1:n_max] + delta; yC <- d$Y0[(n_max + 1):(2 * n_max)]
    hit <- FALSE
    for (nk in looks) {
      db <- mean(yT[1:nk]) - mean(yC[1:nk])
      se <- sqrt(var(yT[1:nk]) / nk + var(yC[1:nk]) / nk)
      if (abs(db / se) > qnorm(1 - alpha / 2)) { hit <- TRUE; break }
    }
    rej[r] <- hit
  }
  mean(rej)
}

res_peek <- list(
  n_max = n_max, n_looks = n_looks, tau2 = tau2, alpha = alpha,
  naive_curve = data.frame(looks = c(1, 2, 5, 10, 14), fpr = naive_curve),
  msprt_typeI = msprt_aa$rej,
  msprt_power = msprt_ab$rej, msprt_stop_median = msprt_ab$stop_med,
  naive_power_peek = naive_ab_pow
)
saveRDS(res_peek, "res_peek.rds")

cat("Naive fixed-0.05 peeking, type-I error vs number of looks (A/A):\n")
print(res_peek$naive_curve)
cat(sprintf("\nmSPRT type-I error under continuous monitoring (A/A) = %.3f\n",
            msprt_aa$rej))
cat(sprintf("mSPRT power under the true effect = %.3f (median stop n = %d/arm)\n",
            msprt_ab$rej, msprt_ab$stop_med))
cat(sprintf("(naive peeking 'power' under the effect = %.3f, but its A/A error is %.3f)\n",
            naive_ab_pow, tail(naive_curve, 1)))
