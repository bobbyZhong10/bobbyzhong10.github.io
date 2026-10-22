# ------------------------------------------------------------------
# 05_ht_slopes.R -- blocks P4, P5 and P5b.
#   P4   Hausman-Taylor (plm, random.method = "ht", inst.method = "baltagi"; Wooldridge 2010 §11.3)
#        against RE, FE and pooled OLS for the effect gamma of Nimbus Verified status; the same
#        coefficients from ht_iv() (IV on quasi-demeaned data); the unit-level first stage of W_i on
#        Xbar_1i.
#   P5   random slopes correlated with the within-merchant promotion variability: FE converges to the
#        within-variance-weighted average slope beta_* (truth10 P5$fe_plim), the mean-group estimator to
#        the average slope; tests of correlated random slopes (Wooldridge 2010 §11.7.4).
#   P5b  random slopes correlated with the level Xbar_i only: FE covers the average slope.
# Reads nimbus_P4.rds, nimbus_P5.rds, nimbus_P5b.rds, truth10.rds; writes res_p45.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("panel_fns.R")
suppressPackageStartupMessages({library(plm); library(sandwich)})
tr <- readRDS("truth10.rds")
inside <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])
row_of <- function(est, se, truth) list(est = unname(est), se = unname(se), ci = ci95(unname(est), unname(se)),
                                        covers_truth = inside(ci95(unname(est), unname(se)), truth))

## =================================================================
## P4: Hausman-Taylor
## =================================================================
d4 <- readRDS("nimbus_P4.rds"); T4 <- tr$P4
pd <- pdata.frame(d4, index = c("id", "t"))
ht <- plm(y ~ x1 + x2 + w | x1 | x2, data = pd, model = "random", random.method = "ht", inst.method = "baltagi")
re <- plm(y ~ x1 + x2 + w, data = pd, model = "random")
po <- plm(y ~ x1 + x2 + w, data = pd, model = "pooling")
fe <- plm(y ~ x1 + x2, data = pd, model = "within")
cl <- function(m) vcovHC(m, method = "arellano", type = "HC0", cluster = "group")
Vht <- vcov(ht); Vre <- vcov(re)
hi <- ht_iv(d4$y, d4$x1, d4$x2, d4$w, d4$id, varphi = ercomp(ht)$theta)
u <- d4[d4$t == 1, c("id", "w")]; u$x1bar <- tapply(d4$x1, d4$id, mean)[as.character(u$id)]
u$x2bar <- tapply(d4$x2, d4$id, mean)[as.character(u$id)]
fs <- lm(w ~ x1bar, data = u)
fs_t <- coef(fs)[["x1bar"]] / sqrt(vcovHC(fs, type = "HC0")["x1bar", "x1bar"])
res4 <- list(
  ht = row_of(coef(ht)["w"], sqrt(Vht["w", "w"]), T4$gamma),
  re = row_of(coef(re)["w"], sqrt(Vre["w", "w"]), T4$gamma),
  pool = row_of(coef(po)["w"], sqrt(cl(po)["w", "w"]), T4$gamma),
  re_cl = row_of(coef(re)["w"], sqrt(cl(re)["w", "w"]), T4$gamma),
  ht_beta1 = row_of(coef(ht)["x1"], sqrt(Vht["x1", "x1"]), T4$beta1),
  ht_beta2 = row_of(coef(ht)["x2"], sqrt(Vht["x2", "x2"]), T4$beta2),
  re_beta2 = row_of(coef(re)["x2"], sqrt(Vre["x2", "x2"]), T4$beta2),
  fe_beta1 = row_of(coef(fe)["x1"], sqrt(cl(fe)["x1", "x1"]), T4$beta1),
  fe_beta2 = row_of(coef(fe)["x2"], sqrt(cl(fe)["x2", "x2"]), T4$beta2),
  varphi_ht = unname(ercomp(ht)$theta), varphi_re = unname(ercomp(re)$theta),
  ht_iv_maxdiff = max(abs(hi$coef[names(coef(ht))] - coef(ht))),
  share_certified = mean(u$w), corr_w_x1bar = cor(u$w, u$x1bar), corr_w_alpha = cor(u$w, d4$alpha_true[d4$t == 1]),
  first_stage_coef = unname(coef(fs)[["x1bar"]]), first_stage_F = unname(fs_t^2),
  haus_robust = hausman_robust(mundlak_fit(d4$y, cbind(x1 = d4$x1, x2 = d4$x2), d4$id)),
  N = length(unique(d4$id)), T = max(d4$t))

## =================================================================
## P5 and P5b: random slopes
## =================================================================
slopes_block <- function(d, avg, plim) {
  f <- fe_fit(d$y, d$x, d$id); m <- mg_fit(d$y, d$x, d$id)
  po <- coef(lm.fit(cbind(1, d$x), d$y))[2]
  xw <- within_tf(d$y, d$x, d$id)$X[, 1]
  wi <- tapply(xw^2, d$id, sum); bi <- tapply(d$b_true, d$id, `[`, 1)
  s2hat <- ave(xw^2, d$id) * max(d$t) / (max(d$t) - 1)                 # within variance of X_i, per row
  type <- tapply(d$s_true, d$id, `[`, 1)
  list(fe = row_of(f$b, sqrt(f$V[1, 1]), avg), mg = row_of(m$est, m$se, avg),
       fe_covers_plim = inside(ci95(unname(f$b), sqrt(f$V[1, 1])), plim),
       pooled = unname(po),
       sample_avg_b = mean(bi), sample_weighted_b = sum(wi * bi) / sum(wi),    # sample analogues of beta and beta_*
       test_xbar = slope_test(d$y, d$x, d$id, ave(d$x, d$id))[c("stat", "df", "p")],
       test_s2 = slope_test(d$y, d$x, d$id, s2hat)[c("stat", "df", "p")],
       by_type = data.frame(s = sort(unique(as.vector(type))), n = as.vector(table(type)),
                            mean_b = as.vector(tapply(bi, type, mean)),
                            weight_fe = as.vector(tapply(wi, type, sum)) / sum(wi)),
       unit = data.frame(s = as.vector(type), b = as.vector(bi), b_hat = m$slopes, w = as.vector(wi)))
}
d5 <- readRDS("nimbus_P5.rds"); d5b <- readRDS("nimbus_P5b.rds")
res5 <- slopes_block(d5, tr$P5$avg_slope, tr$P5$fe_plim)
res5b <- slopes_block(d5b, tr$P5b$avg_slope, tr$P5b$fe_plim)
res5$fe_covers_avg <- res5$fe$covers_truth

saveRDS(list(P4 = res4, P5 = res5, P5b = res5b), "res_p45.rds")
cat(sprintf("P4 gamma (truth %.2f): HT %.4f [%.3f, %.3f]  RE %.4f [%.3f, %.3f]  pooled %.4f; ht_iv max diff %.2g; first-stage F %.1f\n",
            T4$gamma, res4$ht$est, res4$ht$ci[1], res4$ht$ci[2], res4$re$est, res4$re$ci[1], res4$re$ci[2],
            res4$pool$est, res4$ht_iv_maxdiff, res4$first_stage_F))
cat(sprintf("P5:  FE %.4f [%.4f, %.4f] (beta_* %.4f, avg %.2f); MG %.4f [%.4f, %.4f]; slope tests p (Xbar) %.3g (s^2) %.3g\n",
            res5$fe$est, res5$fe$ci[1], res5$fe$ci[2], tr$P5$fe_plim, tr$P5$avg_slope, res5$mg$est, res5$mg$ci[1],
            res5$mg$ci[2], res5$test_xbar$p, res5$test_s2$p))
cat(sprintf("P5b: FE %.4f [%.4f, %.4f] (avg %.2f); MG %.4f; slope tests p (Xbar) %.3g (s^2) %.3g\n",
            res5b$fe$est, res5b$fe$ci[1], res5b$fe$ci[2], tr$P5b$avg_slope, res5b$mg$est, res5b$test_xbar$p, res5b$test_s2$p))
