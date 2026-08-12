# ------------------------------------------------------------------
# 02_power_mde.R -- power, MDE, and the underpowered opening run, Mod 14.
#   - one seeded A/B run at n = 15000/arm: a non-significant p-value
#   - analytic power for the true 2% effect at that n
#   - empirical power (simulation) to confirm the formula
#   - sample size needed for 80% power
#   - power curve vs n
# Writes res_power.rds
# ------------------------------------------------------------------
set.seed(1414)
p <- readRDS("lumen_params.rds")
gen <- p$gen; delta <- p$delta; sigma <- p$sd_Y0

za <- qnorm(0.975); zb <- qnorm(0.80)

## ---- the opening run: n = 15000 per arm --------------------------
n_run <- 15000L
d <- gen(2 * n_run, seed = 24)
Z <- rep(0:1, each = n_run)                 # balanced A/B assignment
Yobs <- ifelse(Z == 1, d$Y1, d$Y0)
tt <- t.test(Yobs[Z == 1], Yobs[Z == 0], var.equal = FALSE)
diff_run <- unname(diff(rev(tt$estimate)))  # treat - control
se_run   <- tt$stderr; p_run <- tt$p.value
t_run    <- unname(tt$statistic)

## ---- analytic MDE and power at n = 15000 -------------------------
se_theory <- sigma * sqrt(2 / n_run)
mde_run   <- (za + zb) * se_theory          # 80%-power MDE at this n
ncp       <- delta / se_theory
power_run <- pnorm(ncp - za) + pnorm(-ncp - za)

## ---- empirical power at n = 15000 (2000 sims) --------------------
sim_power <- function(n, del, R = 2000, seed = 7) {
  set.seed(seed); rej <- logical(R)
  for (r in 1:R) {
    dd <- gen(2 * n)
    zz <- rep(0:1, each = n)
    yy <- ifelse(zz == 1, dd$Y0 + del, dd$Y0)
    rej[r] <- t.test(yy[zz == 1], yy[zz == 0])$p.value < 0.05
  }
  mean(rej)
}
emp_power_run <- sim_power(n_run, delta)

## ---- n for 80% power to detect the true effect -------------------
n_needed <- ceiling(2 * sigma^2 * (za + zb)^2 / delta^2)

## ---- power curve -------------------------------------------------
ns <- c(5e3, 1e4, 1.5e4, 2.5e4, 4e4, 5.65e4, 8e4, 1.2e5)
pw <- sapply(ns, function(n) {
  s <- sigma * sqrt(2 / n); nc <- delta / s
  pnorm(nc - za) + pnorm(-nc - za)
})

res_power <- list(
  n_run = n_run, diff_run = diff_run, se_run = se_run,
  t_run = t_run, p_run = p_run,
  mde_run = mde_run, mde_run_rel = mde_run / p$mu0,
  power_run = power_run, emp_power_run = emp_power_run,
  n_needed = n_needed, sigma = sigma, delta = delta,
  curve = data.frame(n = ns, power = pw)
)
saveRDS(res_power, "res_power.rds")

cat(sprintf("Opening run (n = %d/arm): diff = %+.4f, SE = %.4f, t = %.2f, p = %.3f\n",
            n_run, diff_run, se_run, t_run, p_run))
cat(sprintf("  80%%-power MDE at this n = %.4f (%.1f%% lift); true effect = %.2f (2%%)\n",
            mde_run, 100 * mde_run / p$mu0, delta))
cat(sprintf("  analytic power for the true effect = %.3f ; empirical = %.3f\n",
            power_run, emp_power_run))
cat(sprintf("  n for 80%% power to detect delta = %.2f : %d per arm\n",
            delta, n_needed))
cat("\nPower curve:\n"); print(round(res_power$curve, 3))
