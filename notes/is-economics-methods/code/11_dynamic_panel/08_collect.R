# ------------------------------------------------------------------
# 08_collect.R -- res_extra.rds: the results of 02 and 04 to 07 under the field names of the chapter card
# (tests/test_conclusions.R reads it together with res_dynamic.rds and truth11.rds).
#   ah, ah_diff, ab_onestep, ab_twostep, bb_twostep   block D1 (02): Anderson-Hsiao with the level and the
#                                                     differenced instrument, Arellano-Bond one-step (robust
#                                                     SE) and two-step, Blundell-Bond two-step (Windmeijer SE)
#   D1, D1p                                           02: pooled OLS, FE, HK, SPJ and the GMM fits per block;
#                                                     D1p$ah_diff_mc: the Anderson-Hsiao (differenced instrument)
#                                                     Monte Carlo at rho = 0.90 (cover, median, mean, draws)
#   D2                                                04: J_p_full, est_full, est_collapsed, ols_transformed,
#                                                     fe, n_inst_full, n_inst_collapsed, path, mc, ...
#   D3                                                05: fe and ab (the coefficient beta on X_it; fe_rho,
#                                                     ab_rho for rho), ab_twostep, ab_endog, ab_exog
#   D4                                                06: ab, bb (one-step), *_twostep, ab_stat, bb_stat, fe
#   D5                                                07: per T_obs = 8, 16, 32: fe, hk, spj (Monte Carlo
#                                                     means), sd_*, rmse_*, cover_*
# Reads res_d1.rds, res_d2.rds, res_d3.rds, res_d4.rds, res_d5.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
d1 <- readRDS("res_d1.rds")
x <- c(d1[c("ah", "ah_diff", "ab_onestep", "ab_twostep", "bb_twostep", "D1", "D1p")],
       list(D2 = readRDS("res_d2.rds"), D3 = readRDS("res_d3.rds"), D4 = readRDS("res_d4.rds"), D5 = readRDS("res_d5.rds")))
saveRDS(x, "res_extra.rds")
cat("Saved res_extra.rds\n")
