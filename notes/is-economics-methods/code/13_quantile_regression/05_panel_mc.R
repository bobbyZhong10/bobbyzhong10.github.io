# ------------------------------------------------------------------
# 05_panel_mc.R -- block Q4: panel QR with merchant fixed effects (chapter 13)
# Ruling 17: a finite-sample Monte Carlo illustration, NOT a demonstration of an asymptotic theorem.
# Kato, Galvao and Montes-Rojas (2012) need n^2 (log n)^3 / \mathcal{T} -> 0 for asymptotic normality of the
# dummy-variable QR estimator, and Canay's (2011) rate condition (n / \mathcal{T}^s -> 0) is disputed
# (Besstremyannaya and Golovan 2019); with \mathcal{T} in {5, 20} neither theorem is being checked.
# Estimators of beta(tau) at every tau, for \mathcal{T} = 5 and 20 (N = 200 merchants):
#   pooled   QR of Y_it on (1, X_it), ignoring alpha_i (biased: X is correlated with alpha_i)
#   dummies  QR with one dummy per merchant (fe_dummy_qr, Koenker 2004 without penalty)
#   canay    Canay's two-step estimator (canay_qr)
# single: the two saved samples (01_dgp.R); mc: R = 1000 replications, replication r drawn after
# set.seed(134000 + r) (first the \mathcal{T} = 5 panel, then the \mathcal{T} = 20 panel), run in parallel with
# parallel::mclapply (results do not depend on the number of cores).
# x$Q4$canay_T20 and x$Q4$dummies_T5 (the card's test) are Monte Carlo means at tau = 0.9.
# Output: res_extra.rds element Q4 (save_extra).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("qr_fns.R")
suppressPackageStartupMessages(library(quantreg))
tr <- readRDS("truth13.rds"); P <- tr$params; P4 <- P$Q4; taus <- P$taus
Q4 <- readRDS("marigold_Q4.rds")

est_panel <- function(d) {
  X <- cbind(x = d$x)
  rbind(pooled = vapply(taus, function(t) .muffle_nonunique(rq.fit(cbind(1, d$x), d$y, tau = t, method = "br"))$coefficients[2], numeric(1)),
        dummies = vapply(taus, function(t) fe_dummy_qr(d$y, X, d$id, t)$beta[["x"]], numeric(1)),
        canay = vapply(taus, function(t) canay_qr(d$y, X, d$id, t)$beta[["x"]], numeric(1)))
}
single <- list(T5 = est_panel(Q4$T5), T20 = est_panel(Q4$T20))
for (s in names(single)) colnames(single[[s]]) <- taus

t0 <- proc.time()[["elapsed"]]
cores <- max(1L, min(8L, parallel::detectCores() - 1L))
reps <- parallel::mclapply(seq_len(P4$R), function(r) {
  set.seed(P$seed$Q4_mc + r)
  d5 <- gen_Q4(P4, P4$T_short); d20 <- gen_Q4(P4, P4$T_long)
  c(est_panel(d5), est_panel(d20))                                     # column-major: estimator fastest, then tau, then T
}, mc.cores = cores, mc.set.seed = FALSE)
bad <- vapply(reps, function(v) !is.numeric(v) || length(v) != 30 || anyNA(v), logical(1))
if (any(bad)) stop("05_panel_mc.R: ", sum(bad), " Monte Carlo replications failed")
A <- array(do.call(rbind, reps), dim = c(P4$R, 3, length(taus), 2),
           dimnames = list(NULL, c("pooled", "dummies", "canay"), taus, c("T5", "T20")))
truth_grid <- tr$Q4$beta_grid
mc_mean <- apply(A, c(2, 3, 4), mean)
mc_sd <- apply(A, c(2, 3, 4), sd)
mc_bias <- sweep(mc_mean, 2, truth_grid)
mc_rmse <- sqrt(mc_bias^2 + mc_sd^2 * (P4$R - 1) / P4$R)
k <- which(taus == P4$tau)
cat(sprintf("Monte Carlo: %d replications on %d cores in %.1f s\n", P4$R, cores, proc.time()[["elapsed"]] - t0))

Q4res <- list(tau = P4$tau, R = P4$R,
              canay_T20 = mc_mean["canay", k, "T20"], dummies_T5 = mc_mean["dummies", k, "T5"],
              canay_T5 = mc_mean["canay", k, "T5"], dummies_T20 = mc_mean["dummies", k, "T20"],
              pooled_T5 = mc_mean["pooled", k, "T5"], pooled_T20 = mc_mean["pooled", k, "T20"],
              mc_mean = mc_mean, mc_sd = mc_sd, mc_bias = mc_bias, mc_rmse = mc_rmse,
              mcse_mean = mc_sd / sqrt(P4$R), single = single, truth_grid = truth_grid)
save_extra("Q4", Q4res)
cat("Saved res_extra.rds: Q4\n")
cat("truth beta(tau):", round(truth_grid, 4), "\n")
for (Tn in c("T5", "T20")) { cat("Monte Carlo bias,", Tn, "\n"); print(round(mc_bias[, , Tn], 4)) }
cat("Monte Carlo sd, T5 and T20\n"); print(round(mc_sd[, , "T5"], 4)); print(round(mc_sd[, , "T20"], 4))
