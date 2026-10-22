# ------------------------------------------------------------------
# 04_contract.R -- block V3: the LEN contract for Vireo's ad-review moderators, with and without an extra signal.
#
# Closed forms (Holmstrom and Milgrom 1987; len_contract in mech_fns.R): output x = e + eps, Var(eps) = 4
# (thousand USD a month), CARA A = 0.5, c(e) = e^2 / 2. Four contracts: x alone; x and the user-report rate
# (y = -0.5 e + eta, Var(eta) = 1, independent of eps); x and a peer team's output (no effort content, noise
# correlated 0.6 with eps, variance 4); x and a pure-noise metric (variance 1, independent). The principal's
# expected profit is e - c(e) - risk premium (the fixed salary keeps the moderator's certainty equivalent at 0).
# Small simulation (seed 3103, Ruling P23): 500 replications of 1,000 moderator-months under the contract on x
# alone and the contract on x and the report rate (common draws); the mean and the 2.5% to 97.5% range of the
# realised profit gain, and the moderator's realised certainty equivalent.
# Output: res_len.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("mech_fns.R")
tr <- readRDS("truth03.rds"); P3 <- tr$params$V3; seed <- tr$params$seed
rng <- function(x) unname(stats::quantile(x, c(0.025, 0.975), type = 7))

## ---- closed forms ----------------------------------------------------------------------------
sig_matrix <- function(s) {
  cv <- s$rho * sqrt(P3$sigma2 * s$var)
  matrix(c(P3$sigma2, cv, cv, s$var), 2)
}
K <- list(without = len_contract(P3$A, P3$c2, matrix(P3$sigma2), 1))
for (nm in names(P3$signals)) K[[nm]] <- len_contract(P3$A, P3$c2, sig_matrix(P3$signals[[nm]]), c(1, P3$signals[[nm]]$k))
cases <- rbindlist(lapply(names(K), function(nm) {
  k <- K[[nm]]; s <- P3$signals[[nm]]
  data.table(case = nm, k = if (is.null(s)) NA_real_ else s$k, var = if (is.null(s)) NA_real_ else s$var,
             rho = if (is.null(s)) NA_real_ else s$rho, sigma_s2 = k$sigma_s2, beta = k$beta, b_x = k$b[1],
             b_signal = if (length(k$b) > 1) k$b[2] else 0, effort = k$effort, effort_cost = k$effort_cost,
             risk_premium = k$risk_premium, profit = k$profit, pay_sd = k$pay_sd, w0 = k$w0)
}))
stopifnot(abs(cases$beta[1] - len_bstar(P3$A, P3$sigma2, P3$c2)) < 1e-15)

## ---- small simulation (seed 3103) ---------------------------------------------------------------
set.seed(seed$V3_sim)
k0 <- K$without; k1 <- K$report; s <- P3$signals$report; L <- t(chol(sig_matrix(s)))
R <- P3$R_sim; N <- P3$N_sim
sim <- matrix(NA_real_, R, 4, dimnames = list(NULL, c("profit_without", "profit_with", "ce_without", "ce_with")))
for (r in seq_len(R)) {
  z <- L %*% matrix(stats::rnorm(2 * N), 2, N)                       # (eps, eta), correlated as Sigma says
  eps <- z[1, ]; eta <- z[2, ]
  # contract on x alone
  x0 <- k0$effort + eps; w_0 <- k0$w0 + k0$b[1] * x0
  # contract on x and the report rate
  x1 <- k1$effort + eps; y1 <- s$k * k1$effort + eta; w_1 <- k1$w0 + k1$b[1] * x1 + k1$b[2] * y1
  ce <- function(w, e) -log(mean(exp(-P3$A * w))) / P3$A - P3$c2 * e^2 / 2
  sim[r, ] <- c(mean(x0 - w_0), mean(x1 - w_1), ce(w_0, k0$effort), ce(w_1, k1$effort))
}
gain <- sim[, "profit_with"] - sim[, "profit_without"]
simres <- list(R = R, N = N, gain_mean = mean(gain), gain_sd = stats::sd(gain), gain_range = rng(gain),
               share_gain_positive = mean(gain > 0),
               profit_without_mean = mean(sim[, "profit_without"]), profit_without_range = rng(sim[, "profit_without"]),
               profit_with_mean = mean(sim[, "profit_with"]), profit_with_range = rng(sim[, "profit_with"]),
               ce_without_mean = mean(sim[, "ce_without"]), ce_with_mean = mean(sim[, "ce_with"]),
               ce_without_range = rng(sim[, "ce_without"]), ce_with_range = rng(sim[, "ce_with"]))

## ---- save ----------------------------------------------------------------------------------------
k3 <- list(cases = cases, b_without = K$without$beta, b_with_signal = K$report$beta,
           profit_without = K$without$profit, profit_with_signal = K$report$profit,
           b_peer = K$peer$beta, profit_peer = K$peer$profit, b_noise = K$noise$beta, profit_noise = K$noise$profit,
           sim = simres)
saveRDS(k3, "res_len.rds")
for (i in seq_len(nrow(cases))) with(cases[i], cat(sprintf(
  "%-8s sigma_s^2 %.4f  beta %.4f  b_x %.4f  b_signal %+.4f  effort %.4f  risk premium %.4f  profit %.4f (thousand USD a month)\n",
  case, sigma_s2, beta, b_x, b_signal, effort, risk_premium, profit)))
cat(sprintf("simulation, %d x %d moderator-months: profit gain mean %.4f range [%.4f, %.4f] (closed form %.4f); CE %.4f / %.4f\n",
            R, N, simres$gain_mean, simres$gain_range[1], simres$gain_range[2], k3$profit_with_signal - k3$profit_without,
            simres$ce_without_mean, simres$ce_with_mean))
cat("Saved res_len.rds\n")
