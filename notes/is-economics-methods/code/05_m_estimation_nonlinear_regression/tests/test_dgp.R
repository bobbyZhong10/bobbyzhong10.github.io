# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Larkspur data and truth05.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("m_fns.R")
tr <- readRDS("truth05.rds")
M1 <- readRDS("larkspur_M1.rds"); M2 <- readRDS("larkspur_M2.rds"); M2m <- readRDS("larkspur_M2m.rds")
M3 <- readRDS("larkspur_M3.rds"); M4 <- readRDS("larkspur_M4.rds")
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))

## ---- truth05 fields named on the chapter card ------------------------------------
stopifnot(fin(tr$M1$log_ratio), fin(tr$M2$theta, 3L), fin(tr$M2m$pseudo), fin(tr$M3$theta, 2L), fin(tr$M4$beta, 2L))
# every block has its own seed, all different
sd <- unlist(tr$params$seed)
stopifnot(length(sd) >= 11, !anyDuplicated(sd), all(sd == round(sd)))

## ---- the generators reproduce the saved data from the stored seeds ------------------
for (b in c("M1", "M2", "M2m", "M3", "M4")) {
  set.seed(tr$params$seed[[b]])
  regen <- get(paste0("gen_", b))(tr$params[[b]])
  stopifnot(isTRUE(all.equal(regen, get(b), tolerance = 0, check.attributes = TRUE)))
}

## ---- M1: merchant cross-section, many zero-sales merchants, sales in USD and cents ----
P1 <- tr$params$M1
stopifnot(is.data.frame(M1), nrow(M1) == P1$n,
          all(c("merchant", "category", "tenure_std", "D", "cents", "usd", "y0_cents_true", "y1_cents_true") %in% names(M1)))
stopifnot(all(M1$cents >= 0), all(M1$cents == round(M1$cents)), max(abs(100 * M1$usd - M1$cents)) < 1e-6)
stopifnot(sum(M1$D) == P1$n / 2, setequal(unique(M1$D), 0:1))                         # complete randomisation, half treated
stopifnot(identical(M1$cents, ifelse(M1$D == 1, M1$y1_cents_true, M1$y0_cents_true)))  # observed = D y(1) + (1 - D) y(0)
stopifnot(all(M1$y1_cents_true >= M1$y0_cents_true),                                   # the feature never lowers sales
          any(M1$y0_cents_true == 0 & M1$y1_cents_true > 0))                           # and it moves the extensive margin
z0 <- mean(M1$usd[M1$D == 0] == 0); z1 <- mean(M1$usd[M1$D == 1] == 0)
stopifnot(z0 > 0.45, z0 < 0.75, z1 < z0 - 0.05)                                        # plausible share of zero-sales merchants
pos <- M1$usd[M1$usd > 0]
stopifnot(median(pos) > 50, median(pos) < 1000, mean(pos) < 2000, quantile(pos, 0.99) < 20000)   # monthly USD of small merchants
# truth: E[y(d) | category, tenure] = exp(alpha_c + tau d + beta tenure) with tau = log(k) + delta
stopifnot(all(P1$k * P1$p0 < 1), abs(sum(P1$pi) - 1) < 1e-12,
          abs(tr$M1$log_ratio - (log(P1$k) + P1$delta)) < 1e-12,
          abs(tr$M1$ratio - exp(tr$M1$log_ratio)) < 1e-12,
          abs(tr$M1$mean_y1_usd / tr$M1$mean_y0_usd - tr$M1$ratio) < 1e-12)
# large population sample (seed M1_pop) agrees with the closed form
stopifnot(tr$M1$pop$n >= 1e6, abs(tr$M1$pop$log_ratio - tr$M1$log_ratio) < 0.01,
          abs(tr$M1$pop$mean_y0_usd / tr$M1$mean_y0_usd - 1) < 0.02,
          abs(tr$M1$pop$share_sell1 - tr$M1$pop$share_sell0 - tr$M1$ext_effect) < 0.005)
# the log(1 + y) estimand itself depends on the unit (population ATE on log(1 + a y))
stopifnot(tr$M1$pop$ate_log1p_cents - tr$M1$pop$ate_log1p_usd > 0.2 * tr$M1$pop$ate_log1p_usd,
          tr$M1$pop$ate_log1p_kusd < tr$M1$pop$ate_log1p_usd)

## ---- M2: exponential mean, heteroskedastic multiplicative error --------------------
P2 <- tr$params$M2
stopifnot(is.data.frame(M2), nrow(M2) == P2$n, all(c("x1", "free_ship", "y", "m_true") %in% names(M2)),
          all(M2$y > 0), setequal(unique(M2$free_ship), 0:1), identical(tr$M2$theta, P2$theta))
stopifnot(max(abs(M2$m_true - exp(P2$theta[1] + P2$theta[2] * M2$x1 + P2$theta[3] * M2$free_ship))) < 1e-12,
          median(M2$y) > 0.3, median(M2$y) < 5)                                          # USD thousands per month
# Var(y | x) grows with the mean (constant coefficient of variation): model-based NLS variance is wrong
lo <- M2$m_true < median(M2$m_true)
stopifnot(var(M2$y[!lo] - M2$m_true[!lo]) > 2 * var(M2$y[lo] - M2$m_true[lo]))

## ---- M2m: misspecified mean, pseudo-true value --------------------------------------
P2m <- tr$params$M2m
stopifnot(is.data.frame(M2m), nrow(M2m) == P2m$n, all(c("x", "y", "m_true") %in% names(M2m)), all(M2m$y > 0))
stopifnot(abs(tr$M2m$pseudo - tr$M2m$pseudo_vec[2]) < 1e-12, fin(tr$M2m$pseudo_vec, 2L), fin(tr$M2m$pseudo_ppml_vec, 2L))
m0 <- function(x) exp(P2m$a0 + P2m$a1 * x + P2m$a2 * x^2)
# quadrature is accurate (half the grid gives the same answer) and agrees with a large noise-free sample
stopifnot(max(abs(pseudo_true_exp(m0, P2m$sd_x, "nls", n_grid = 20001) - tr$M2m$pseudo_vec)) < 1e-8,
          max(abs(pseudo_true_exp(m0, P2m$sd_x, "ppml", n_grid = 20001) - tr$M2m$pseudo_ppml_vec)) < 1e-8,
          tr$M2m$large$N >= 1e6, max(abs(tr$M2m$large$pseudo_vec - tr$M2m$pseudo_vec)) < 0.01)
# misspecification matters: the pseudo-true slope is not the true elasticity at x = 0, and
# NLS and PPML have different pseudo-true values (each estimator defines its own estimand)
stopifnot(abs(tr$M2m$pseudo - P2m$a1) > 0.05, abs(tr$M2m$pseudo_ppml - tr$M2m$pseudo) > 0.02)

## ---- M3: two-step, generated traffic index --------------------------------------------
P3 <- tr$params$M3
stopifnot(is.data.frame(M3), nrow(M3) == P3$n, all(c("z1", "z2", "v", "y", "t_true") %in% names(M3)), all(M3$y > 0),
          identical(tr$M3$theta, P3$theta), identical(tr$M3$gamma, P3$gamma), P3$theta[2] != 0)   # theta1 != 0: Gamma_0 != 0
r2 <- summary(lm(v ~ z1 + z2, data = M3))$r.squared
stopifnot(r2 > 0.1, r2 < 0.5, abs(tr$M3$r2_first_stage - r2) < 0.05)                    # weekly visits are a noisy measure
stopifnot(P3$R_mc >= 1000, all(P3$R_grid >= 1000), 0 %in% P3$grid_theta1, P3$theta[2] %in% P3$grid_theta1)

## ---- M4: overdispersed counts ------------------------------------------------------------
P4 <- tr$params$M4
stopifnot(is.data.frame(M4), nrow(M4) == P4$n, all(c("x", "y", "mu_true") %in% names(M4)),
          all(M4$y >= 0), all(M4$y == round(M4$y)), identical(tr$M4$beta, P4$beta))
stopifnot(abs(tr$M4$avg_cond_var_to_mean - (1 + exp(P4$beta[1] + 1.5 * P4$beta[2]^2) / P4$size)) < 1e-12,
          tr$M4$avg_cond_var_to_mean > 2,
          abs(tr$M4$pearson_dispersion - (1 + exp(P4$beta[1] + P4$beta[2]^2 / 2) / P4$size)) < 1e-12,
          tr$M4$marginal_var_to_mean > tr$M4$avg_cond_var_to_mean)            # the marginal ratio adds Var(mu)
# the three population overdispersion quantities agree with a large noise-free draw of x
set.seed(1); xb <- rnorm(1e6); mub <- exp(P4$beta[1] + P4$beta[2] * xb); vb <- mub + mub^2 / P4$size
stopifnot(abs(mean(vb) / mean(mub) / tr$M4$avg_cond_var_to_mean - 1) < 0.02,
          abs(mean(vb / mub) / tr$M4$pearson_dispersion - 1) < 0.02,
          abs((mean(vb) + var(mub)) / mean(mub) / tr$M4$marginal_var_to_mean - 1) < 0.02)
stopifnot(var(M4$y) / mean(M4$y) > 2, mean(M4$y == 0) < 0.3, mean(M4$y) > 2, mean(M4$y) < 10)   # monthly orders
stopifnot(P4$R_mc >= 1000, all(P4$R_grid >= 1000), Inf %in% P4$grid_size, P4$size %in% P4$grid_size)
