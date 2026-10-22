# ------------------------------------------------------------------
# 08_collect.R -- res_extra.rds: the results of 05 to 07 under the field names of the chapter card
# (tests/test_conclusions.R reads it together with truth12.rds).
#   B2   coef_interaction (beta_jj'), avg_cross and avg_cross_ci (average cross-partial effect), ...
#   B3   naive (probit of Y on X, D), cf (Rivers-Vuong), ivmle, tsls, ols_lpm
#   B4   beta (short-probit coefficient on X), ape, ape_ci, long (with a_i), ...
#   B5   clogit, pooled, fe, bc (one panel with T = 8; fe and bc carry the Monte Carlo biases bias_T4 and
#        bias_T8), cmle (hand-written CMLE with robust SE), T4 (the same estimators on months 1 to 4),
#        cre (pooled CRE probit), cre_re (CRE random-effects probit), mc (Monte Carlo summaries), sample
#   B6   ape_top, ape_top_ci (ordered probit, P(rating = 5)), mc (ruling 29), homo, Bond and Lang numbers
#   B7   iia (MNL, HM test, its Monte Carlo rejection rate mc_reject), nested (HM test, hm_p)
# Reads res_cs.rds, res_panel.rds, res_mult.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
cs <- readRDS("res_cs.rds"); pn <- readRDS("res_panel.rds"); mu <- readRDS("res_mult.rds")
mc <- lapply(pn$mc, function(m) m[setdiff(names(m), "draws")])
B5 <- list(clogit = pn$T8$clogit, pooled = pn$T8$pooled,
           fe = c(pn$T8$fe, list(bias_T4 = mc$T4$fe_bias, bias_T8 = mc$T8$fe_bias)),
           bc = c(pn$T8$bc, list(bias_T4 = mc$T4$bc_bias, bias_T8 = mc$T8$bc_bias)),
           cmle = pn$T8$cmle, T8 = pn$T8, T4 = pn$T4, cre = pn$cre, cre_re = pn$cre_re, mc = mc, sample = pn$sample)
x <- list(B2 = cs$B2, B3 = cs$B3, B4 = cs$B4, B5 = B5, B6 = mu$B6, B7 = mu$B7)
saveRDS(x, "res_extra.rds")
cat("Saved res_extra.rds\n")
