# ------------------------------------------------------------------
# 05_global_game.R -- block G4: the global-game threshold as private noise shrinks.
#
# theta_star: the threshold of the unique equilibrium, solved numerically from the marginal player's
# indifference condition (gg_threshold) for sigma = 0.1, 0.05, 0.01 with the public prior N(0.25, 0.5^2).
# Checks from the simulated data of 01_dgp.R: theta_star_sim, the state at which the simulated share of
# 2,000 merchants adopting crosses 1/2 (linear interpolation on the theta grid); mc_indiff, the marginal
# player's expected payoff from adopting, m(x*) + E[l | x*] - 1, with E[l | x*] simulated (200,000 draws).
# Uniqueness: gg_min_slope > 0 at each sigma, and sigma_max, the largest private noise for which the
# threshold equilibrium is unique whatever the public signal (here prior sd 0.5). Counter-example: with the
# public signal as precise as the private one (prior sd = sigma = 0.1) and centred at 1/2, three thresholds.
# Output: res_gg.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("games_fns.R")
tr <- readRDS("truth02.rds"); P4 <- tr$params$G4
G4 <- readRDS("osprey_G4.rds")

theta_star <- gg_threshold(P4$sigma, prior_mean = P4$prior_mean, prior_sd = P4$prior_sd)
crossing <- function(th, share) {
  i <- which(share[-length(share)] < 0.5 & share[-1] >= 0.5)[1]
  th[i] + (0.5 - share[i]) * (th[i + 1] - th[i]) / (share[i + 1] - share[i])
}
pop <- G4$pop; pop$share <- pop$n_adopt / pop$n_pop
theta_star_sim <- vapply(P4$sigma, function(s) { d <- pop[pop$sigma == s, ]; crossing(d$theta, d$share) }, numeric(1))
mc <- G4$mc
a <- 1 / P4$prior_sd^2
m_post <- vapply(seq_along(P4$sigma), function(k) { b <- 1 / P4$sigma[k]^2; (a * P4$prior_mean + b * theta_star[k]) / (a + b) }, numeric(1))
El_exact <- 1 - m_post                                                             # from the indifference condition
multi <- list(sigma = 0.1, prior_mean = 0.5, prior_sd = 0.1)
multi$roots <- gg_roots(multi$sigma, multi$prior_mean, multi$prior_sd)
multi$sigma_max <- gg_sigma_max(multi$prior_sd)

res <- list(
  sigma = P4$sigma, theta_star = theta_star, theta_star_limit = tr$G4$theta_star_limit,
  dist_to_limit = abs(theta_star - tr$G4$theta_star_limit),
  theta_star_sim = theta_star_sim, m_post = m_post, El_exact = El_exact, El_sim = mc$El_sim,
  mc_indiff = mc$m_post + mc$El_sim - 1,
  min_slope = gg_min_slope(P4$sigma, P4$prior_sd), sigma_max = gg_sigma_max(P4$prior_sd),
  # share adopting at theta = 1/2 and at theta = 0.55: how sharp the tipping is at each sigma
  share_at_half = stats::pnorm((0.5 - theta_star) / P4$sigma),
  share_at_055 = stats::pnorm((0.55 - theta_star) / P4$sigma),
  flat_prior_threshold = gg_threshold(P4$sigma),
  prior_mean = P4$prior_mean, prior_sd = P4$prior_sd, n_pop = P4$n_pop, n_mc = P4$n_mc, pop = pop,
  multiple_example = multi)
saveRDS(res, "res_gg.rds")

for (k in seq_along(P4$sigma))
  cat(sprintf("sigma %.2f  theta* %.6f  |theta* - 1/2| %.6f  simulated crossing %.6f  E[l|x*] exact %.6f sim %.6f  min slope %.4f\n",
              P4$sigma[k], theta_star[k], res$dist_to_limit[k], theta_star_sim[k], El_exact[k], mc$El_sim[k], res$min_slope[k]))
cat(sprintf("unique for every public signal when sigma < %.6f (prior sd %.2f); prior sd 0.1, mean 0.5, sigma 0.1: thresholds %s\n",
            res$sigma_max, P4$prior_sd, paste(sprintf("%.6f", multi$roots), collapse = " ")))
cat("Saved res_gg.rds\n")
