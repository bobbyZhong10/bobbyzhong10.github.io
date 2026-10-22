# ------------------------------------------------------------------
# 02_algo_collusion.R -- two Q-learning pricing algorithms in a
# repeated Bertrand-logit game (Calvano-Calzolari-Denicolo-Pastorello
# 2020 style). Neither bot is told to collude and they never
# communicate; each just maximizes its own discounted profit by
# tabular Q-learning. The question is whether they nonetheless learn
# supra-competitive (collusive) prices.
#
# Demand (logit): q_i = exp((a - p_i)/mu) / (sum_j exp((a - p_j)/mu) + exp(a0/mu))
# Profit: pi_i = (p_i - c) q_i.
# State s_t = (p_{1,t-1}, p_{2,t-1}); action = own price on a grid.
# Q-update: Q(s,a) <- (1-lr) Q(s,a) + lr [pi + delta * max_a' Q(s',a')]
# Exploration: epsilon-greedy with epsilon decaying as exp(-beta * t).
#
# Collusion index Delta = (pbar - pNash) / (pMonopoly - pNash), the
# fraction of the way from the competitive to the collusive price.
# ------------------------------------------------------------------
set.seed(22)

## ---- demand & profit primitives ---------------------------------
a <- 2; a0 <- 0; mu <- 0.25; cost <- 1
dem <- function(p) {                      # p = c(p1, p2)
  e <- exp((a - p) / mu); e0 <- exp(a0 / mu)
  e / (sum(e) + e0)
}
prof <- function(p) (p - cost) * dem(p)

## ---- Bertrand-Nash price (symmetric fixed point) ----------------
br <- function(pj) optimize(function(pi) -prof(c(pi, pj))[1], c(cost, 3))$minimum
pN <- 1.5
for (it in 1:200) pN <- br(pN)
pNash <- pN

## ---- symmetric monopoly (joint-profit) price --------------------
pMon <- optimize(function(p) -sum(prof(c(p, p))), c(cost, 3))$minimum

## ---- price grid (Calvano: xi below Nash to above monopoly) ------
m  <- 11
xi <- 0.1
p_lo <- pNash - xi * (pMon - pNash)
p_hi <- pMon  + xi * (pMon - pNash)
grid <- seq(p_lo, p_hi, length.out = m)

## ---- precompute per-action profits for every state pair ---------
## profit matrix PI[i,j] = profit to firm choosing grid[i] vs rival grid[j]
PI <- matrix(0, m, m)
for (i in 1:m) for (j in 1:m) PI[i, j] <- prof(c(grid[i], grid[j]))[1]

## ---- Q-learning session ------------------------------------------
lr    <- 0.125     # learning rate
delta <- 0.95      # discount
beta  <- 6e-6      # exploration decay
Tmax  <- 1000000

run_session <- function(seed, keep_path = FALSE) {
  set.seed(seed)
  Q1 <- array(mean(PI) / (1 - delta), dim = c(m, m, m))  # Q1[p1prev,p2prev,a1]
  Q2 <- array(mean(PI) / (1 - delta), dim = c(m, m, m))
  s1 <- sample(m, 1); s2 <- sample(m, 1)
  path <- if (keep_path) numeric(Tmax) else NULL
  tail_sum <- 0; tail_n <- 0; t0 <- 0.95 * Tmax
  for (t in 1:Tmax) {
    eps <- exp(-beta * t)
    a1 <- if (runif(1) < eps) sample(m, 1) else which.max(Q1[s1, s2, ])
    a2 <- if (runif(1) < eps) sample(m, 1) else which.max(Q2[s1, s2, ])
    r1 <- PI[a1, a2]; r2 <- PI[a2, a1]
    Q1[s1, s2, a1] <- (1 - lr) * Q1[s1, s2, a1] + lr * (r1 + delta * max(Q1[a1, a2, ]))
    Q2[s1, s2, a2] <- (1 - lr) * Q2[s1, s2, a2] + lr * (r2 + delta * max(Q2[a1, a2, ]))
    s1 <- a1; s2 <- a2
    pm <- (grid[a1] + grid[a2]) / 2
    if (keep_path) path[t] <- pm
    if (t > t0) { tail_sum <- tail_sum + pm; tail_n <- tail_n + 1 }
  }
  list(pbar = tail_sum / tail_n, path = path)
}

## ---- run several independent sessions ---------------------------
seeds <- 1:6
pbars <- numeric(length(seeds)); path1 <- NULL
for (k in seq_along(seeds)) {
  r <- run_session(seeds[k], keep_path = (k == 1))
  pbars[k] <- r$pbar
  if (k == 1) path1 <- r$path
}
Delta_vec <- (pbars - pNash) / (pMon - pNash)
Delta_mean <- mean(Delta_vec)

cat("== Algorithmic pricing: two Q-learners, Bertrand-logit ==\n")
cat(sprintf("Bertrand-Nash price = %.3f, monopoly price = %.3f\n", pNash, pMon))
cat(sprintf("Grid: %d prices from %.3f to %.3f\n", m, p_lo, p_hi))
cat(sprintf("Sessions: %d, converged prices = %s\n", length(seeds),
            paste(sprintf("%.3f", pbars), collapse = ", ")))
cat(sprintf("Collusion index per session = %s\n",
            paste(sprintf("%.3f", Delta_vec), collapse = ", ")))
cat(sprintf("Mean collusion index Delta = %.3f (range %.3f to %.3f)\n",
            Delta_mean, min(Delta_vec), max(Delta_vec)))
cat(sprintf("Supra-competitive sessions: %d of %d\n",
            sum(Delta_vec > 0.1), length(seeds)))

saveRDS(list(pNash = pNash, pMon = pMon, grid = grid, path1 = path1,
             pbars = pbars, Delta_vec = Delta_vec, Delta_mean = Delta_mean,
             Tmax = Tmax, m = m,
             params = list(a=a,a0=a0,mu=mu,cost=cost,lr=lr,delta=delta,beta=beta)),
        file.path(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))),
          "res26_collusion.rds"))
cat("Saved res26_collusion.rds\n")
