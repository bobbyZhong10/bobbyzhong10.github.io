# ============================================================
# Chapter 16 - Berry inversion: OLS logit vs IV logit.
# The opening puzzle: OLS understates price sensitivity so badly that the
# implied marginal cost of the median product turns NEGATIVE.
# ============================================================
suppressPackageStartupMessages({library(data.table); library(fixest)})
D <- readRDS("ch20_dgp.rds"); dt <- D$dt; J <- D$J; truth <- D$truth
setkey(dt, t, j)

## Berry (1994) plain-logit inversion: y = ln(s_j) - ln(s0)
dt[, y := log(s) - log(s0)]

## BLP instruments: sum of rival characteristics within market (other firms)
dt[, sumx1_riv := sum(x1) - x1, by = t]
dt[, sumx2_riv := sum(x2) - x2, by = t]
dt[, nriv := .N - 1L, by = t]

## ---- OLS logit ----
ols <- feols(y ~ x1 + x2 + p, data = dt)
a_ols <- -coef(ols)["p"]
cat("--- OLS logit (Berry inversion) ---\n")
cat("alpha_OLS =", round(a_ols,3), " (truth", truth$alpha, ")  beta_x1 =",
    round(coef(ols)["x1"],3), " beta_x2 =", round(coef(ols)["x2"],3), "\n")

## ---- IV logit: instrument p with cost shifter w (+ BLP instruments) ----
iv <- feols(y ~ x1 + x2 | p ~ w + sumx1_riv + sumx2_riv, data = dt)
a_iv <- -coef(iv)["fit_p"]
cat("\n--- IV logit (instrument price with w + BLP instruments) ---\n")
cat("alpha_IV =", round(a_iv,3), " (truth", truth$alpha, ")  beta_x1 =",
    round(coef(iv)["x1"],3), " beta_x2 =", round(coef(iv)["x2"],3), "\n")
# first-stage F
fs <- feols(p ~ x1 + x2 + w + sumx1_riv + sumx2_riv, data = dt)
cat("first-stage coef on w =", round(coef(fs)["w"],3),
    " | partial F(w) approx =", round(fitstat(fs,"f")$f$stat,1), "\n")

## ---- implied median marginal cost under each alpha (single-product logit) ----
## logit own derivative ds_j/dp_j = -alpha s_j (1-s_j); markup = s_j / (alpha s_j (1-s_j))
## = 1 / (alpha (1-s_j)); mc = p - markup
implied_mc <- function(a){ dt[, mk := 1/(a*(1-s))]; dt[, p - mk] }
mc_ols <- implied_mc(a_ols); mc_iv <- implied_mc(a_iv)
cat("\n--- implied marginal cost (single-product logit FOC) ---\n")
cat("OLS: median implied mc =", round(median(mc_ols),2),
    " | share negative =", round(mean(mc_ols<0),3), "\n")
cat("IV : median implied mc =", round(median(mc_iv),2),
    " | share negative =", round(mean(mc_iv<0),3), "\n")
cat("TRUTH: median mc =", round(median(dt$mc),2), "\n")

## ---- own-price elasticity of the median product ----
dt[, e_ols := -a_ols*p*(1-s)]; dt[, e_iv := -a_iv*p*(1-s)]; dt[, e_true := -truth$alpha*p*(1-s)]
cat("\n--- logit own-price elasticity (median) ---\n")
cat("OLS", round(median(dt$e_ols),3), " IV", round(median(dt$e_iv),3),
    " (naive logit-formula truth", round(median(dt$e_true),3), ")\n")

saveRDS(list(a_ols=a_ols, a_iv=a_iv, ols=coef(ols), iv=coef(iv),
             mc_ols=mc_ols, mc_iv=mc_iv,
             Fw=fitstat(fs,"f")$f$stat), "ch20_logit.rds")
cat("saved ch20_logit.rds\n")
