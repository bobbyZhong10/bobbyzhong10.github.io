# ------------------------------------------------------------------
# 05_interference.R -- SUTVA violation in a marketplace and how cluster
#                      / switchback designs fix it, Module 14.
#
# Mechanism: each market-slot has a fixed courier capacity C and N users.
# The feature raises a user's order propensity (p0 -> p1). Orders share
# capacity: if demand exceeds C, every wanter is served with prob C/demand.
# A user-level A/B splits users within a slot, so treated users' extra
# demand cannibalizes capacity from control users -> the naive
# treated-minus-control gap OVERSTATES the global rollout effect.
# Switchback (randomize whole slots) and cluster (randomize whole markets)
# assignments each let a unit experience one global world, removing the
# cross-user contamination.
# Writes res_interf.rds
# ------------------------------------------------------------------
set.seed(1414)

M       <- 80L        # markets
Tslots  <- 12L        # time slots per market
N       <- 600L       # users per market-slot
C       <- 150L       # courier capacity per market-slot
p0      <- 0.20       # control order propensity
p1      <- 0.35       # treated order propensity
r       <- 1.0        # value of a served order

## slot-level congestion multiplier (some slots busy, some quiet) -----
draw_mult <- function() matrix(runif(M * Tslots, 0.6, 1.6), M, Tslots)

## serve one slot: given wanter count, everyone served w.p. min(1,C/dem)
# returns per-slot average realized outcome given each user's want-prob
slot_outcome <- function(want_prob_vec) {
  wants  <- rbinom(length(want_prob_vec), 1, want_prob_vec)
  demand <- sum(wants)
  serve_p <- if (demand > C) C / demand else 1
  served  <- wants * rbinom(length(wants), 1, serve_p)
  served * r
}

## ---- TRUE global ATE: all-treat vs all-control (large sim) --------
truth_pass <- function(prob, mult) {
  ys <- numeric(0)
  for (j in 1:M) for (t in 1:Tslots) {
    wp <- pmin(prob * mult[j, t], 0.95)
    ys <- c(ys, slot_outcome(rep(wp, N)))
  }
  mean(ys)
}
set.seed(7)
mult_truth <- draw_mult()
avg1 <- mean(replicate(6, truth_pass(p1, draw_mult())))
avg0 <- mean(replicate(6, truth_pass(p0, draw_mult())))
true_ate <- avg1 - avg0

## ---- Design A: user-level A/B (interference-contaminated) --------
sim_userAB <- function(mult) {
  yT <- yC <- numeric(0)
  for (j in 1:M) for (t in 1:Tslots) {
    z  <- rbinom(N, 1, 0.5)                    # within-slot split
    wp <- pmin(ifelse(z == 1, p1, p0) * mult[j, t], 0.95)
    y  <- slot_outcome(wp)
    yT <- c(yT, y[z == 1]); yC <- c(yC, y[z == 0])
  }
  mean(yT) - mean(yC)
}

## ---- Design B: switchback (whole slot is T or C) -----------------
sim_switchback <- function(mult) {
  yT <- yC <- numeric(0)
  for (j in 1:M) for (t in 1:Tslots) {
    zt <- rbinom(1, 1, 0.5)
    wp <- pmin((if (zt == 1) p1 else p0) * mult[j, t], 0.95)
    y  <- slot_outcome(rep(wp, N))
    if (zt == 1) yT <- c(yT, y) else yC <- c(yC, y)
  }
  mean(yT) - mean(yC)
}

## ---- Design C: cluster (whole market is T or C) ------------------
sim_cluster <- function(mult) {
  yT <- yC <- numeric(0)
  zc <- rbinom(M, 1, 0.5)
  for (j in 1:M) for (t in 1:Tslots) {
    wp <- pmin((if (zc[j] == 1) p1 else p0) * mult[j, t], 0.95)
    y  <- slot_outcome(rep(wp, N))
    if (zc[j] == 1) yT <- c(yT, y) else yC <- c(yC, y)
  }
  mean(yT) - mean(yC)
}

## ---- replicate each design to get bias and SE --------------------
R <- 200L
set.seed(2024)
est_ab <- est_sb <- est_cl <- numeric(R)
for (r_i in 1:R) {
  est_ab[r_i] <- sim_userAB(draw_mult())
  est_sb[r_i] <- sim_switchback(draw_mult())
  est_cl[r_i] <- sim_cluster(draw_mult())
}

summ <- function(v) c(mean = mean(v), sd = sd(v), bias = mean(v) - true_ate)
res_interf <- list(
  true_ate = true_ate, avg1 = avg1, avg0 = avg0,
  userAB = summ(est_ab), switchback = summ(est_sb), cluster = summ(est_cl),
  params = c(M = M, Tslots = Tslots, N = N, C = C, p0 = p0, p1 = p1),
  est = data.frame(userAB = est_ab, switchback = est_sb, cluster = est_cl)
)
saveRDS(res_interf, "res_interf.rds")

cat(sprintf("Global worlds: avg outcome all-control %.4f, all-treat %.4f\n",
            avg0, avg1))
cat(sprintf("TRUE global ATE = %.4f\n\n", true_ate))
cat(sprintf("User-level A/B : mean %.4f  bias %+.4f  (SD %.4f)  <- overstates\n",
            summ(est_ab)["mean"], summ(est_ab)["bias"], summ(est_ab)["sd"]))
cat(sprintf("Switchback     : mean %.4f  bias %+.4f  (SD %.4f)\n",
            summ(est_sb)["mean"], summ(est_sb)["bias"], summ(est_sb)["sd"]))
cat(sprintf("Cluster        : mean %.4f  bias %+.4f  (SD %.4f)\n",
            summ(est_cl)["mean"], summ(est_cl)["bias"], summ(est_cl)["sd"]))
cat(sprintf("\nNaive A/B overstates the global effect by %.1fx\n",
            summ(est_ab)["mean"] / true_ate))
