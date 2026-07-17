# ------------------------------------------------------------------
# 04_modern.R -- Modern staggered-DiD estimators on the Kestrel panel:
#   (1) Sun-Abraham interaction-weighted event study (fixest::sunab)
#   (2) Callaway-Sant'Anna group-time ATTs (did::att_gt) with
#       not-yet-treated controls + simple/dynamic/group aggregation
#   (3) Imputation estimator (didimputation::did_imputation),
#       Borusyak-Jaravel-Spiess.
# ------------------------------------------------------------------
library(data.table)
library(fixest)
library(did)
library(didimputation)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dt    <- readRDS(file.path(code_dir, "kestrel_panel.rds"))
truth <- readRDS(file.path(code_dir, "truth.rds"))

## ---- (1) Sun-Abraham via fixest::sunab ---------------------------
## Never-treated cities need a cohort value outside the sample period.
dt[, g_sunab := fifelse(g == 0, 10000L, as.integer(g))]
m_sa <- feols(lgmv ~ sunab(g_sunab, month) | city_id + month,
              data = dt, cluster = ~city_id)

sa_att <- summary(m_sa, agg = "att")
cat("== Sun-Abraham (sunab): ATT aggregate ==\n")
print(coeftable(sa_att))

ct <- as.data.table(coeftable(m_sa), keep.rownames = "term")
ct[, e := as.integer(gsub(".*::(-?\\d+).*", "\\1", term))]
sa_es <- ct[, .(e, estimate = Estimate, se = `Std. Error`)][order(e)]
cat("\nSun-Abraham event-study coefficients at selected horizons:\n")
print(sa_es[e %in% c(-6, -3, 0, 3, 6, 9), .(e, estimate = round(estimate, 4), se = round(se, 4))])

## ---- (2) Callaway-Sant'Anna via did::att_gt ----------------------
## Fixed seed so the multiplier-bootstrap SEs/CIs reproduce exactly.
set.seed(11)
cs <- att_gt(yname = "lgmv", tname = "month", idname = "city_id",
             gname = "g", data = dt, control_group = "notyettreated",
             clustervars = "city_id", bstrap = TRUE, biters = 1000)

agg_simple  <- aggte(cs, type = "simple")
agg_dynamic <- aggte(cs, type = "dynamic")
agg_group   <- aggte(cs, type = "group")

cat("\n== Callaway-Sant'Anna: simple aggregate ==\n")
print(summary(agg_simple))
cat("\n== Callaway-Sant'Anna: group aggregate ==\n")
print(summary(agg_group))

cs_es <- data.table(e = agg_dynamic$egt, estimate = agg_dynamic$att.egt,
                    se = agg_dynamic$se.egt)
cat("\nCS dynamic (event-study) coefficients at selected horizons:\n")
print(cs_es[e %in% c(-6, -3, 0, 3, 6, 9), .(e, estimate = round(estimate, 4), se = round(se, 4))])

## ---- (3) BJS imputation via didimputation ------------------------
imp <- did_imputation(data = dt, yname = "lgmv", gname = "g",
                      tname = "month", idname = "city_id")
cat("\n== BJS imputation: overall ATT ==\n")
print(imp)

imp_es <- did_imputation(data = dt, yname = "lgmv", gname = "g",
                         tname = "month", idname = "city_id",
                         horizon = TRUE, pretrends = -12:-2)
imp_es_dt <- as.data.table(imp_es)[, .(e = as.integer(term), estimate, se = std.error)][order(e)]
cat("\nImputation event-study coefficients at selected horizons:\n")
print(imp_es_dt[e %in% c(-6, -3, 0, 3, 6, 9), .(e, estimate = round(estimate, 4), se = round(se, 4))])

## ---- Save results ------------------------------------------------
saveRDS(list(
  sa_att   = list(est = coeftable(sa_att)["ATT", "Estimate"],
                  se  = coeftable(sa_att)["ATT", "Std. Error"]),
  sa_es    = sa_es,
  cs_simple = list(est = agg_simple$overall.att, se = agg_simple$overall.se),
  cs_group  = data.table(g = agg_group$egt, att = agg_group$att.egt, se = agg_group$se.egt),
  cs_es     = cs_es,
  imp       = list(est = imp$estimate[1], se = imp$std.error[1]),
  imp_es    = imp_es_dt
), file.path(code_dir, "res_modern.rds"))
cat("\nSaved res_modern.rds\n")
