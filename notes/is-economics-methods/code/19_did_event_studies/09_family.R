# ------------------------------------------------------------------
# 09_family.R -- The imputation / regression-based family on the Kestrel panel:
#   (1) did2s (Gardner 2022; Butts & Gardner 2022): first stage fits
#       city and month FE on untreated observations only, second stage
#       regresses the residualised outcome on D; GMM-corrected SEs
#       clustered by city.
#   (2) etwfe (Wooldridge 2021/2023, Wong et al.): saturated cohort-by-time
#       interactions in one regression, aggregated with
#       marginaleffects via emfx(type = "simple") (equal weights on all
#       treated observations).
#
# Both target the same estimand as BJS imputation (04_modern.R): the
# equal-weight average of the effect over all 930 treated city-months,
# with all untreated observations (never-treated and not-yet-treated)
# as the comparison. With city and month FE and no covariates the point
# estimates coincide with BJS up to numerical precision (did2s: exactly;
# etwfe with cgroup = "notyet": ~2e-8, from marginaleffects); only the
# variance estimators differ.
#
# etwfe codes never-treated cities with gvar = 0 (our g = 0). cgroup =
# "notyet" (the package default) uses not-yet-treated cities as controls
# and matches BJS; cgroup = "never" restricts the controls to the 8
# never-treated cities and matches the Sun-Abraham ATT of 04_modern.R.
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(did2s)
  library(etwfe)
})

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dt <- readRDS(file.path(code_dir, "kestrel_panel.rds"))
m  <- readRDS(file.path(code_dir, "res_modern.rds"))
df <- as.data.frame(dt)

warn_log <- character()
with_log <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    warn_log <<- c(warn_log, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}
pack <- function(est, se) {
  list(est = unname(est), se = unname(se),
       lb = unname(est - qnorm(0.975) * se), ub = unname(est + qnorm(0.975) * se))
}

## ---- (1) did2s ---------------------------------------------------
fit_d2s <- with_log(did2s(df, yname = "lgmv",
                          first_stage  = ~ 0 | city_id + month,
                          second_stage = ~ i(D, ref = 0),
                          treatment = "D", cluster_var = "city_id"))
cat("== did2s ==\n"); print(fit_d2s)
d2s <- pack(coef(fit_d2s)[["D::1"]], se(fit_d2s)[["D::1"]])

## ---- (2) etwfe ---------------------------------------------------
fit_etwfe <- function(cgroup) {
  mod <- with_log(etwfe(fml = lgmv ~ 0, tvar = month, gvar = g, data = df,
                        ivar = city_id, cgroup = cgroup, vcov = ~city_id))
  ef <- with_log(emfx(mod, type = "simple"))
  cat(sprintf("\n== etwfe (cgroup = %s), emfx type = simple ==\n", cgroup)); print(ef)
  pack(ef$estimate, ef$std.error)
}
etw       <- fit_etwfe("notyet")
etw_never <- fit_etwfe("never")

## ---- Comparison with 04_modern.R ---------------------------------
cmp <- data.frame(
  estimator = c("did2s", "etwfe (notyet)", "etwfe (never)", "BJS imputation", "SA average ATT"),
  estimate  = c(d2s$est, etw$est, etw_never$est, m$imp$est, m$sa_att$est),
  se        = c(d2s$se, etw$se, etw_never$se, m$imp$se, m$sa_att$se))
cat("\n== Family vs BJS / SA ==\n"); print(cmp, digits = 6)
cat(sprintf("\n|did2s - BJS| = %.2e;  |etwfe(notyet) - BJS| = %.2e;  |etwfe(never) - SA| = %.2e\n",
            abs(d2s$est - m$imp$est), abs(etw$est - m$imp$est), abs(etw_never$est - m$sa_att$est)))

if (length(warn_log)) { cat("\nPackage warnings:\n"); print(table(warn_log)) } else cat("\nNo package warnings.\n")

saveRDS(list(did2s = d2s, etwfe = etw, etwfe_never = etw_never,
             compare = list(did2s_minus_bjs = d2s$est - m$imp$est,
                            etwfe_minus_bjs = etw$est - m$imp$est,
                            etwfe_never_minus_sa = etw_never$est - m$sa_att$est),
             versions = c(did2s = as.character(packageVersion("did2s")),
                          etwfe = as.character(packageVersion("etwfe"))),
             warnings = warn_log),
        file.path(code_dir, "res_family.rds"))
cat("Saved res_family.rds\n")
