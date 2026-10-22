# ------------------------------------------------------------------
# 06_initial_conditions.R -- block D4: an onboarding cohort that starts at 80% of its long-run GMV
# (2,000 merchants, 8 months, rho = 0.6), and its stationary twin (the same draws, start_share = 1).
#   Arellano-Bond (difference GMM) needs no assumption on the initial condition; Blundell-Bond (system
#   GMM) adds the level moments E[Delta Y_{i,t-1}(alpha_i + eps_it)] = 0, which need mean stationarity.
# Both with month effects (effect = "twoways"), every available lag, one-step with the robust SE (main),
# two-step with the Windmeijer SE alongside; FE for reference.
# Reads nimbus_D4.rds, truth11.rds; writes res_d4.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("dyn_fns.R")
suppressPackageStartupMessages(library(plm))
tr <- readRDS("truth11.rds"); P <- tr$params$D4
d <- readRDS("nimbus_D4.rds")
pd <- pdata.frame(d, index = c("id", "t"))

gm <- function(v, tf, model = "onestep") {
  f <- stats::as.formula(sprintf("%s ~ lag(%s, 1) | lag(%s, 2:99)", v, v, v))
  gmm_summary(pgmm(f, pd, effect = "twoways", model = model, transformation = tf), par = sprintf("lag(%s, 1)", v))
}
fe <- fe_ar1(d$y, d$id, d$t)
res <- list(ab = gm("y", "d"), bb = gm("y", "ld"), ab_twostep = gm("y", "d", "twosteps"), bb_twostep = gm("y", "ld", "twosteps"),
            ab_stat = gm("y_stat_true", "d"), bb_stat = gm("y_stat_true", "ld"),
            fe = list(est = unname(fe$est["rho"]), se = unname(fe$se["rho"]), ci = ci95(unname(fe$est["rho"]), unname(fe$se["rho"]))),
            N = P$N, T = P$T)
saveRDS(res, "res_d4.rds")

cat(sprintf("=== D4: %d merchants, %d months, start at %.0f%% of the long-run level; truth rho %.2f ===\n",
            P$N, P$T, 100 * P$start_share, P$rho))
for (k in c("ab", "bb", "ab_twostep", "bb_twostep", "ab_stat", "bb_stat")) {
  s <- res[[k]]
  cat(sprintf("%-11s %.4f (SE %.4f) [%.4f, %.4f]; r = %d, Hansen p %.4f, m2 p %.3f\n", k, s$est, s$se, s$ci[1], s$ci[2],
              s$r, s$J_p, s$m2_p))
}
cat(sprintf("FE          %.4f (SE %.4f)\n", res$fe$est, res$fe$se))
cat("Saved res_d4.rds\n")
