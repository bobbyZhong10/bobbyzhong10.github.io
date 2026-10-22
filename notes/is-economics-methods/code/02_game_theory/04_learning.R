# ------------------------------------------------------------------
# 04_learning.R -- block G3: regret matching and fictitious play.
#
# Regret matching (Hart and Mas-Colell 2000): the empirical joint distribution of play converges to the set of
# correlated equilibria. For the two 50,000-period runs, ce_violation of the empirical distribution after t
# periods (the checkpoints of 01_dgp.R: about 200, log-spaced) and at T, the cell each run settles on and when; for the 100
# independent runs of each game, the mean and largest violation at each checkpoint and where the runs end.
# Fictitious play on matching pennies: the empirical marginal frequencies (Robinson 1951: they converge to
# the mixed equilibrium (1/2, 1/2) in two-player zero-sum games) along the run, and the share of periods on
# each cell.
# Output: res_learn.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("games_fns.R")
tr <- readRDS("truth02.rds"); P3 <- tr$params$G3
G3 <- readRDS("osprey_G3.rds")
Tn <- P3$T; cp <- P3$checkpoints

# ce_violation of the empirical joint distribution of the first t periods of one run, at each checkpoint t
viol_path <- function(a1, a2, A, B) {
  cell <- (a1 - 1L) * 2L + a2
  cum <- vapply(1:4, function(k) cumsum(cell == k), numeric(length(cell)))
  vapply(cp, function(t) ce_violation(matrix(cum[t, ] / t, 2, 2, byrow = TRUE), A, B), numeric(1))
}
games <- list(compat = list(A = P3$compat_A, B = P3$compat_B, run = G3$rm_compat, reps = G3$reps_compat),
              chicken = list(A = P3$chicken_A, B = P3$chicken_B, run = G3$rm_chicken, reps = G3$reps_chicken))
paths <- do.call(rbind, lapply(names(games), function(g) {
  z <- games[[g]]; rp <- attr(z$reps, "path")
  data.frame(game = g, t = cp, viol_run = viol_path(attr(z$run, "a1"), attr(z$run, "a2"), z$A, z$B),
             viol_reps_mean = colMeans(rp), viol_reps_max = apply(rp, 2, max),
             share_reps_zero = colMeans(rp == 0))
}))
mu_T <- lapply(games, function(z) { m <- z$run; attributes(m) <- list(dim = c(2L, 2L)); m })
# where a run ends: the cell holding at least 99% of the last 10,000 periods (0 if none)
last_cell <- function(a1, a2) {
  idx <- (Tn - 9999L):Tn; tab <- tabulate((a1[idx] - 1L) * 2L + a2[idx], 4) / 10000
  if (max(tab) >= 0.99) which.max(tab) else 0L
}
# first period from which the run never leaves its final cell
lock_in <- function(a1, a2) {
  cell <- (a1 - 1L) * 2L + a2; ch <- which(cell != cell[Tn])
  if (length(ch)) max(ch) + 1L else 1L
}
fp <- G3$fp_mp
fp_path <- data.frame(t = cp, p = cumsum(fp$a1 == 1L)[cp] / cp, q = cumsum(fp$a2 == 1L)[cp] / cp)
fp_joint <- matrix(tabulate((fp$a1 - 1L) * 2L + fp$a2, 4) / Tn, 2, 2, byrow = TRUE)

res <- list(
  T = Tn, R = P3$R, checkpoints = cp, paths = paths,
  ce_violation_T = vapply(games, function(z) ce_violation(z$run, z$A, z$B), numeric(1)),
  mu_T = mu_T,
  final_cell = vapply(games, function(z) last_cell(attr(z$run, "a1"), attr(z$run, "a2")), integer(1)),
  lock_in_period = vapply(games, function(z) lock_in(attr(z$run, "a1"), attr(z$run, "a2")), integer(1)),
  reps_max_violation = vapply(games, function(z) max(z$reps$viol), numeric(1)),
  reps_mean_violation = vapply(games, function(z) mean(z$reps$viol), numeric(1)),
  reps_share_zero_T = vapply(games, function(z) mean(z$reps$viol == 0), numeric(1)),
  # share of runs whose empirical distribution puts at least 0.99 on one cell, by cell (11, 12, 21, 22)
  reps_share_cell = lapply(games, function(z) colMeans(as.matrix(z$reps[, c("mu11", "mu12", "mu21", "mu22")]) >= 0.99)),
  fp_freq = fp$freq, fp_path = fp_path, fp_joint = fp_joint,
  fp_ce_violation = ce_violation(fp_joint, P3$mp_A, P3$mp_B),
  fp_max_dev_last_half = max(abs(c(fp_path$p[fp_path$t >= Tn / 2], fp_path$q[fp_path$t >= Tn / 2]) - 0.5)))
saveRDS(res, "res_learn.rds")

for (g in names(games))
  cat(sprintf("%-8s CE violation at T = %d: %.6f; final cell %d (from period %d); 100 runs: max violation at T %.6f, mean %.6f; shares on cells %s\n",
              g, Tn, res$ce_violation_T[[g]], res$final_cell[[g]], res$lock_in_period[[g]], res$reps_max_violation[[g]],
              res$reps_mean_violation[[g]], paste(sprintf("%.2f", res$reps_share_cell[[g]]), collapse = " ")))
cat(sprintf("fictitious play, matching pennies: p = %.4f, q = %.4f; joint %s; CE violation %.4f\n",
            res$fp_freq[1], res$fp_freq[2], paste(round(as.vector(t(fp_joint)), 4), collapse = " "), res$fp_ce_violation))
cat("Saved res_learn.rds\n")
