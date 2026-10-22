# ------------------------------------------------------------------
# 03_menu.R -- block V2: Vireo's SaaS usage menus (screening, deterministic).
#
# Continuous types (Mussa and Rosen 1978): for the regular and the bimodal type density on the 2,001-point
# grid of theta in [2, 10] USD per thousand calls, the optimal menu q(theta) = max(0, Psi_bar) / cost2 with
# rent U = int q and payment t = theta q - U; also the menu without ironing (q = max(0, Psi) / cost2), which
# is not incentive compatible when Psi falls. Two types: theta_L = 4, theta_H = 7, mu = 0.4.
# Output: res_menu.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("mech_fns.R")
tr <- readRDS("truth03.rds"); P2 <- tr$params$V2
th <- tr$V2$grid; c2 <- tr$V2$cost2

at_theta <- function(m, x, col) stats::approx(m$theta, m[[col]], xout = x)$y
run <- function(par) {
  dens <- type_density(par, th)
  m <- optimal_menu(th, dens, c2); m0 <- optimal_menu(th, dens, c2, iron = FALSE)
  sm <- menu_summary(m, c2)
  k <- which(m$psi_bar > 0)[1]                                  # exclusion threshold: Psi_bar crosses 0
  theta_psi0 <- th[k - 1] + (0 - m$psi_bar[k - 1]) / (m$psi_bar[k] - m$psi_bar[k - 1]) * (th[k] - th[k - 1])
  dec <- which(diff(m$psi) < 0)
  c(sm, list(menu = m, menu_naive = m0, q = m$q, theta_psi0 = theta_psi0,
             mean_theta = sum(m$w * m$theta), psi_monotone = length(dec) == 0L,
             psi_fall_lo = if (length(dec)) th[min(dec)] else NA_real_,
             psi_fall_hi = if (length(dec)) th[max(dec) + 1L] else NA_real_,
             profit_naive_if_truthful = sum(m0$w * (m0$t - c2 * m0$q^2 / 2)), ic_gain_naive = menu_ic_gain(m0),
             rent_share = sm$rent / sm$surplus, profit_share_eff = sm$profit / sm$surplus_efficient,
             q_at = setNames(at_theta(m, c(4, 6, 8), "q"), c("t4", "t6", "t8")),
             q_eff_at = setNames(at_theta(m, c(4, 6, 8), "q_eff"), c("t4", "t6", "t8")),
             t_at = setNames(at_theta(m, c(4, 6, 8), "t"), c("t4", "t6", "t8")),
             U_at = setNames(at_theta(m, c(4, 6, 8), "U"), c("t4", "t6", "t8"))))
}
s2 <- list(regular = run(P2$regular), irregular = run(P2$irregular))
tt <- P2$two_type
s2$two_type <- c(tt, screening_two_type(tt$theta_L, tt$theta_H, tt$mu, tt$cost2))
saveRDS(s2, "res_menu.rds")

for (nm in c("regular", "irregular")) with(s2[[nm]], cat(sprintf(
  "%-9s: excluded below theta = %.4f (Psi_bar = 0 at %.4f); q_top %.4f = q_eff_top %.4f; profit %.4f, rent %.4f, efficient surplus %.4f; bunch [%s, %s] at q = %.4f (mass %.4f); IC gain %.2e, without ironing %.4f\n",
  nm, theta_exclusion, theta_psi0, q_top, q_eff_top, profit, rent, surplus_efficient,
  format(bunch_lo, digits = 5), format(bunch_hi, digits = 5), bunch_q, bunch_mass, ic_gain, ic_gain_naive)))
with(s2$two_type, cat(sprintf("two types: q_L %.4f (first best %.4f), q_H %.4f; t_L %.4f, t_H %.4f; rent %.4f; profit %.4f (first best %.4f)\n",
                              q_L, q_L_fb, q_H, t_L, t_H, rent_H, profit, profit_fb)))
cat("Saved res_menu.rds\n")
