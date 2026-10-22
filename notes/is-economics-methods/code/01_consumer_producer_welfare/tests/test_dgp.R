# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Plover data and truth01.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("welfare_fns.R")
tr <- readRDS("truth01.rds")
W1 <- readRDS("plover_W1.rds"); W1g <- readRDS("plover_W1g.rds")
W2 <- readRDS("plover_W2.rds"); W2b <- readRDS("plover_W2_bronars.rds")
P <- tr$params
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))

## ---- truth01 fields named on the chapter card ------------------------------------------
N <- P$W1$N
stopifnot(fin(tr$W1$cv_i, N), fin(tr$W1$ev_i, N), fin(tr$W1$dcs_i, N), fin(tr$W1$cv_mean), fin(tr$W1g$cv_mean),
          is.character(tr$W2$type), length(tr$W2$type) == P$W2$n_users, setequal(unique(tr$W2$type), c("rational", "noisy")))
stopifnot(abs(tr$W1$cv_mean - mean(tr$W1$cv_i)) < 1e-12, abs(tr$W1g$cv_mean - mean(tr$W1g$cv_i)) < 1e-12)
# seeds: integers, all different (W1g reuses W1's seed by design: same subscribers, only beta differs)
sd <- unlist(P$seed)
stopifnot(length(sd) == 3, !anyDuplicated(sd), all(sd == round(sd)))

## ---- the generators reproduce the saved data from the stored seeds ----------------------
set.seed(P$seed$W1); stopifnot(identical(gen_W1(P$W1), W1))
set.seed(P$seed$W1); stopifnot(identical(gen_W1(P$W1, common_beta = TRUE), W1g))
set.seed(P$seed$W2); stopifnot(identical(gen_W2(P$W2), W2))
set.seed(P$seed$W2_bronars); stopifnot(identical(gen_W2_bronars(W2, P$W2$n_bronars_per_user), W2b))

## ---- W0: hand block, no random numbers ---------------------------------------------------
stopifnot(is.list(tr$W0), fin(tr$W0$cv), fin(tr$W0$ev), fin(tr$W0$dcs), tr$W0$ev < tr$W0$dcs, tr$W0$dcs < tr$W0$cv)

## ---- W1: price experiment, Stone-Geary subscribers ------------------------------------------
stopifnot(N >= 5000, nrow(W1) == N,
          identical(names(W1), c("subscriber", "p1", "w", "x1", "beta_true", "g1_true", "g2_true")))
grid <- c(2, 2.5, 3, 3.5, 4)
stopifnot(identical(P$W1$p_grid, grid), setequal(unique(W1$p1), grid), all(table(W1$p1) == N / 5))  # complete randomisation, equal arms
stopifnot(max(abs(W1$x1 - sg_demand(W1$p1, W1$w, W1$beta_true, W1$g1_true, W1$g2_true))) < 1e-12)  # observed usage = demand
stopifnot(all(W1$beta_true > 0 & W1$beta_true < 1), all(W1$g1_true > 0), all(W1$g2_true >= 0),
          all(W1$x1 > W1$g1_true),                                                          # interior
          all(W1$w - P$W1$p_out * W1$g1_true - W1$g2_true > 0))                            # utility defined up to 12 USD
# income: lognormal, median about 100 USD; marginal budget share independent of income
stopifnot(abs(median(W1$w) / P$W1$w_median - 1) < 0.03, abs(sd(log(W1$w)) / P$W1$w_sdlog - 1) < 0.05,
          abs(cor(W1$beta_true, W1$w)) < 0.05, abs(mean(W1$beta_true) - P$W1$beta_mean) < 0.01,
          sd(W1$beta_true) > 0.05)                                                          # heterogeneous beta: Gorman fails
# randomisation balances income across the price arms
arm_w <- tapply(W1$w, W1$p1, mean)
stopifnot(max(abs(arm_w / mean(W1$w) - 1)) < 0.03)
# usage falls with price; monthly spending on premium is a plausible share of the budget
arm_x <- tapply(W1$x1, W1$p1, mean)
stopifnot(all(diff(arm_x) < 0), mean(W1$p1 * W1$x1 / W1$w) > 0.2, mean(W1$p1 * W1$x1 / W1$w) < 0.5)
# W1g: the same subscribers, incomes, subsistence quantities and prices; beta common (Gorman holds)
keep <- c("subscriber", "p1", "w", "g1_true", "g2_true")
stopifnot(identical(W1g[keep], W1[keep]), all(W1g$beta_true == P$W1$beta_mean),
          max(abs(W1g$x1 - sg_demand(W1g$p1, W1g$w, W1g$beta_true, W1g$g1_true, W1g$g2_true))) < 1e-12)

## ---- W1 truth: closed forms agree with cv_ev_exact, vartia_cv and dcs_numeric; ordering ------------
for (i in c(1L, 2L, 777L, N)) {
  bi <- W1$beta_true[i]; g1i <- W1$g1_true[i]; g2i <- W1$g2_true[i]; wi <- W1$w[i]
  e_i <- function(p1, u) sg_expenditure(p1, u, bi, g1i, g2i)
  x_i <- function(p1, m) sg_demand(p1, m, bi, g1i, g2i)
  ce <- cv_ev_exact(e_i, 2, 4, wi)
  stopifnot(abs(tr$W1$cv_i[i] - ce[["cv"]]) < 1e-9, abs(tr$W1$ev_i[i] - ce[["ev"]]) < 1e-9,
            abs(tr$W1$cv_i[i] - vartia_cv(x_i, 2, 4, wi)) < 1e-9,
            abs(tr$W1$dcs_i[i] - dcs_numeric(x_i, 2, 4, wi)) < 1e-9,
            abs(tr$W1$cv12_i[i] - vartia_cv(x_i, 2, P$W1$p_out, wi)) < 1e-8)
}
stopifnot(all(tr$W1$ev_i < tr$W1$dcs_i), all(tr$W1$dcs_i < tr$W1$cv_i),
          all(tr$W1g$ev_i < tr$W1g$dcs_i), all(tr$W1g$dcs_i < tr$W1g$cv_i))
# the representative consumer of the population (LES with b = mean beta, fitted to per-capita demand):
# exact under Gorman (W1g), off under heterogeneous beta (W1)
stopifnot(abs(tr$W1g$cv_rep_pop - tr$W1g$cv_mean) < 1e-9 * tr$W1g$cv_mean,
          abs(tr$W1$cv_rep_pop - tr$W1$cv_mean) > 0.005 * tr$W1$cv_mean)
# the per-capita true demand curve stored for the figure is the sample mean of individual demand
dc <- tr$W1$demand_curve
j <- which.min(abs(dc$p1 - 3))
stopifnot(abs(dc$x1_mean[j] - mean(sg_demand(dc$p1[j], W1$w, W1$beta_true, W1$g1_true, W1$g2_true))) < 1e-12)

## ---- W2: revealed-preference panel ----------------------------------------------------------
nu <- P$W2$n_users; Tn <- P$W2$n_budgets
stopifnot(nu >= 300, Tn == 8, nrow(W2) == nu * Tn,
          all(c("user", "t", "type_true", "alpha_true", "rho_true", "p1", "p2", "m", "x1", "x2", "x1_opt_true", "x2_opt_true") %in% names(W2)))
ty <- W2$type_true[W2$t == 1]
stopifnot(identical(ty, tr$W2$type), sum(ty == "rational") == round(0.7 * nu), sum(ty == "noisy") == round(0.3 * nu))
stopifnot(all(W2$x1 > 0), all(W2$x2 > 0), max(abs(W2$p1 * W2$x1 + W2$p2 * W2$x2 - W2$m)) < 1e-9)   # on the budget line
rat <- W2$type_true == "rational"
stopifnot(identical(W2$x1[rat], W2$x1_opt_true[rat]), identical(W2$x2[rat], W2$x2_opt_true[rat]),
          mean(abs(W2$x1[!rat] / W2$x1_opt_true[!rat] - 1)) > 0.05)                      # noisy users leave the optimum
# CES first-order condition at the optimum: MRS = p1 / p2
mrs <- with(W2, alpha_true * x1_opt_true^(rho_true - 1) / ((1 - alpha_true) * x2_opt_true^(rho_true - 1)))
stopifnot(max(abs(mrs / (W2$p1 / W2$p2) - 1)) < 1e-9)
# Bronars choosers: the users' own budgets, uniform on the budget line
nb <- P$W2$n_bronars_per_user
stopifnot(nrow(W2b) == nu * nb * Tn, length(unique(W2b$chooser)) == nu * nb)
key <- function(d) paste(d$user, d$t)
mt <- match(key(W2b), key(W2))
stopifnot(!anyNA(mt), identical(W2b$p1, W2$p1[mt]), identical(W2b$p2, W2$p2[mt]), identical(W2b$m, W2$m[mt]),
          max(abs(W2b$p1 * W2b$x1 + W2b$p2 * W2b$x2 - W2b$m)) < 1e-9)
sh <- W2b$p1 * W2b$x1 / W2b$m
stopifnot(abs(mean(sh) - 0.5) < 0.01, abs(var(sh) - 1 / 12) < 0.005, suppressWarnings(ks.test(sh, "punif")$p.value) > 0.001)

## ---- W3: deterministic risk block ---------------------------------------------------------------
P3 <- P$W3
stopifnot(P3$prob * P3$loss < min(P3$prices), max(P3$prices) < P3$loss, P3$buy_at < P3$decline_at,
          all(c(P3$buy_at, P3$decline_at) %in% P3$prices))
