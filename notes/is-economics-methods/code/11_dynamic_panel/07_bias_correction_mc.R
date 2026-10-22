# ------------------------------------------------------------------
# 07_bias_correction_mc.R -- block D5: bias-corrected FE at T_obs = 8, 16, 32 months (T = 7, 15, 31
# regression periods), the D1 design: the legacy gen_dyn(N = 500, T = T_obs, rho = 0.70) from nimbus.rds.
# R_mc = 1,000 samples per length (seeds D5_T8, D5_T16, D5_T32), in each:
#   FE                  rho_hat_FE (fe_ar1), CR0 SE
#   Hahn-Kuersteiner    hk_correct(rho_hat_FE, T), SE ((T + 1)/T) SE_FE (the correction is linear)
#   half-panel jackknife spj_fe (Dhaene and Jochmans 2015, eqs. 2.2 and 2.4: the halves weighted by their shares
#                       |S|/T of the regression periods; fix round 2), with the FE SE (their Theorem 2.1: the same
#                       asymptotic variance as FE as N, T grow)
# Reports per length the mean, SD, RMSE and the coverage of the nominal 95% interval.
# Reads nimbus.rds, truth11.rds; writes res_d5.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("dyn_fns.R")
gen_dyn <- readRDS("nimbus.rds")$gen_dyn                       # the legacy generator of blocks D1 and D1p
tr <- readRDS("truth11.rds"); P <- tr$params$D5
z <- qnorm(0.975)

rows <- lapply(seq_along(P$T_obs), function(j) {
  Tobs <- P$T_obs[j]; TT <- Tobs - 1
  set.seed(tr$params$seed[[paste0("D5_T", Tobs)]])
  m <- t(vapply(seq_len(P$R_mc), function(r) {
    d <- gen_dyn(N = P$N, T = Tobs, rho = P$rho, burn = P$burn, sigma_a = P$sigma_a, sigma_e = P$sigma_e)
    f <- fe_ar1(d$y, d$id, d$t)
    c(fe = unname(f$est["rho"]), se = unname(f$se["rho"]), spj = spj_fe(d$y, d$id, d$t))
  }, numeric(3)))
  hk <- hk_correct(m[, "fe"], TT); se_hk <- (TT + 1) / TT * m[, "se"]
  cov <- function(est, se) mean(abs(est - P$rho) <= z * se)
  data.frame(T_obs = Tobs, T = TT, R = P$R_mc,
             fe = mean(m[, "fe"]), hk = mean(hk), spj = mean(m[, "spj"]),
             sd_fe = sd(m[, "fe"]), sd_hk = sd(hk), sd_spj = sd(m[, "spj"]),
             rmse_fe = sqrt(mean((m[, "fe"] - P$rho)^2)), rmse_hk = sqrt(mean((hk - P$rho)^2)),
             rmse_spj = sqrt(mean((m[, "spj"] - P$rho)^2)),
             mean_se_fe = mean(m[, "se"]), cover_fe = cov(m[, "fe"], m[, "se"]), cover_hk = cov(hk, se_hk),
             cover_spj = cov(m[, "spj"], m[, "se"]))
})
tab <- do.call(rbind, rows)
res <- c(as.list(tab[setdiff(names(tab), "R")]), list(R = P$R_mc, N = P$N, table = tab))
saveRDS(res, "res_d5.rds")

cat(sprintf("=== D5: %d samples of %d merchants per length; truth rho %.2f ===\n", P$R_mc, P$N, P$rho))
print(cbind(tab[, c("T_obs", "T", "fe", "hk", "spj", "sd_fe", "sd_spj", "rmse_fe", "rmse_hk", "rmse_spj",
                    "cover_fe", "cover_hk", "cover_spj")],
            fe_plim = tr$D5$fe_plim, hk_plim = tr$D5$hk_plim, spj_plim = tr$D5$spj_plim), digits = 4, row.names = FALSE)
cat("Saved res_d5.rds\n")
