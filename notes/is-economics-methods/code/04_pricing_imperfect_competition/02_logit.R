# ------------------------------------------------------------------
# 02_logit.R -- block C1: pass-through by OLS in the Cirrus-Garnet logit markets under two conducts, the
# homogeneous-product variant (Cournot against Bertrand), and the pass-through curves of the figure.
#   OLS of each firm's price on both firms' marginal costs, market by market, separately for the Bertrand world
#   (kappa = 0) and the joint-maximisation world (kappa = 1); conventional OLS standard errors (Ruling P23). The
#   slope estimates the average pass-through dp_k / dc_j over the markets (truth04.rds: C1$pt_*_mean).
#   Homogeneous variant: OLS of the market price on the two costs under Cournot and under Bertrand.
#   Curves (n = 1 to 10 symmetric firms, common cost shock): Cournot with linear, exponential and constant-
#   elasticity (e = 2) demand; logit (quality 3.5, the mean of the C1 qualities, alpha = 1, cost 2) under Bertrand
#   and under joint ownership of all products.
#   Lerner and conduct: at the joint-maximisation baseline, the marginal cost an analyst would back out by assuming
#   Bertrand, p - 1 / (alpha (1 - s)).
# Reads cirrus_C1.rds, truth04.rds. Writes res_logit.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("pricing_fns.R")
C1 <- readRDS("cirrus_C1.rds")
tr <- readRDS("truth04.rds")
P1 <- tr$params$C1
firms <- c("Cirrus", "Garnet")

## ---- OLS pass-through under the two conducts --------------------------------------------------------
ols_rows <- list(); fits <- list()
for (cond in c("bertrand", "joint")) {
  truth <- tr$C1[[paste0("pt_", cond, "_mean")]]; base <- tr$C1[[paste0("pt_", cond, "_base")]]
  for (k in 1:2) {
    y <- C1[[sprintf("p%d_%s", k, cond)]]
    fit <- lm(y ~ c1 + c2, data = C1)
    fits[[paste(cond, firms[k], sep = "_")]] <- fit
    cf <- summary(fit)$coefficients
    for (j in 1:2) {
      term <- paste0("c", j)
      ols_rows[[length(ols_rows) + 1L]] <- data.frame(
        conduct = cond, price = firms[k], cost = firms[j], estimate = cf[term, "Estimate"], se = cf[term, "Std. Error"],
        truth_mean = truth[k, j], truth_base = base[k, j], r2 = summary(fit)$r.squared, sigma = summary(fit)$sigma)
    }
  }
}
ols <- do.call(rbind, ols_rows)
pick <- function(cond, price, cost, col) ols[[col]][ols$conduct == cond & ols$price == price & ols$cost == cost]

## ---- homogeneous-product variant ------------------------------------------------------------------------
fc <- summary(lm(p_cournot_hom ~ c1 + c2, data = C1))$coefficients
fb <- summary(lm(p_bertrand_hom ~ c1 + c2, data = C1))$coefficients
qc1 <- C1$q1_cournot_hom; qc2 <- C1$q2_cournot_hom
hom <- list(cournot_coef = unname(fc[c("c1", "c2"), "Estimate"]), cournot_se = unname(fc[c("c1", "c2"), "Std. Error"]),
            bertrand_coef = unname(fb[c("c1", "c2"), "Estimate"]), bertrand_se = unname(fb[c("c1", "c2"), "Std. Error"]),
            p_cournot_mean = mean(C1$p_cournot_hom), p_bertrand_mean = mean(C1$p_bertrand_hom),
            markup_cournot_mean = mean(C1$p_cournot_hom - (qc1 * C1$c1 + qc2 * C1$c2) / (qc1 + qc2)),
            markup_bertrand_mean = mean(C1$p_bertrand_hom - pmin(C1$c1, C1$c2)))

## ---- the single-product logit monopoly: PT = 1 - s (Cirrus alone, baseline quality and cost) ----------------
pm1 <- bertrand_logit(P1$delta[["Cirrus"]], P1$alpha, P1$c_base[1])
mono_example <- list(p = pm1, s = logit_share(pm1, P1$delta[["Cirrus"]], P1$alpha),
                     pt = drop(passthrough_logit(P1$delta[["Cirrus"]], P1$alpha, P1$c_base[1], 0)))

## ---- pass-through curves: curvature and conduct ------------------------------------------------------------
ns <- 1:10
curve_par <- list(delta = mean(P1$delta), alpha = P1$alpha, c = P1$c_base[1], ce_e = 2)
common_pt <- function(n, kappa) {
  d <- rep(curve_par$delta, n); cc <- rep(curve_par$c, n)
  sum(passthrough_logit(d, curve_par$alpha, cc, kappa)[1, ])            # every cost moves by the same amount
}
s0_joint <- vapply(ns, function(n) {
  d <- rep(curve_par$delta, n); 1 - sum(logit_share(joint_max_logit(d, curve_par$alpha, rep(curve_par$c, n)), d, curve_par$alpha))
}, numeric(1))
pt_curves <- rbind(
  data.frame(n = ns, model = "cournot_linear", pt = passthrough_cournot(ns)),
  data.frame(n = ns, model = "cournot_exp", pt = rep(1, length(ns))),
  data.frame(n = ns, model = "cournot_ce2", pt = curve_par$ce_e * ns / (curve_par$ce_e * ns - 1)),
  data.frame(n = ns, model = "logit_bertrand", pt = vapply(ns, common_pt, numeric(1), kappa = 0)),
  data.frame(n = ns, model = "logit_joint", pt = vapply(ns, common_pt, numeric(1), kappa = 1)))

## ---- Lerner and conduct at the baseline -------------------------------------------------------------------
pb <- tr$C1$p_bertrand_base; sb <- tr$C1$s_bertrand_base
pj <- tr$C1$p_joint_base; sj <- tr$C1$s_joint_base
lerner <- list(lerner_bertrand = (pb - P1$c_base) / pb, elas_abs_bertrand = P1$alpha * pb * (1 - sb),
               lerner_joint = (pj - P1$c_base) / pj, elas_abs_joint = P1$alpha * pj * (1 - sj),
               mc_if_bertrand_at_joint = pj - 1 / (P1$alpha * (1 - sj)))       # true marginal cost is c_base = 2

res <- list(p_bertrand = cbind(C1$p1_bertrand, C1$p2_bertrand), p_joint = cbind(C1$p1_joint, C1$p2_joint),
            pt_own_bertrand = pick("bertrand", "Cirrus", "Cirrus", "estimate"),
            pt_cross_bertrand = pick("bertrand", "Cirrus", "Garnet", "estimate"),
            pt_own_joint = pick("joint", "Cirrus", "Cirrus", "estimate"),
            pt_cross_joint = pick("joint", "Cirrus", "Garnet", "estimate"),
            se_own_bertrand = pick("bertrand", "Cirrus", "Cirrus", "se"), se_cross_bertrand = pick("bertrand", "Cirrus", "Garnet", "se"),
            se_own_joint = pick("joint", "Cirrus", "Cirrus", "se"), se_cross_joint = pick("joint", "Cirrus", "Garnet", "se"),
            ols = ols, hom = hom, mono_example = mono_example,
            pt_curves = pt_curves, pt_curves_par = curve_par, pt_curves_s0 = s0_joint, lerner = lerner,
            sample = list(p1_bertrand_mean = mean(C1$p1_bertrand), p2_bertrand_mean = mean(C1$p2_bertrand),
                          p1_joint_mean = mean(C1$p1_joint), p2_joint_mean = mean(C1$p2_joint),
                          s1_bertrand_mean = mean(C1$s1_bertrand), s2_bertrand_mean = mean(C1$s2_bertrand),
                          s0_bertrand_mean = mean(1 - C1$s1_bertrand - C1$s2_bertrand),
                          s0_joint_mean = mean(1 - C1$s1_joint - C1$s2_joint),
                          c1_mean = mean(C1$c1), c2_mean = mean(C1$c2), c1_sd = sd(C1$c1), c2_sd = sd(C1$c2)))
saveRDS(res, "res_logit.rds")
print(ols[, c("conduct", "price", "cost", "estimate", "se", "truth_mean")], digits = 4, row.names = FALSE)
cat(sprintf("Homogeneous: Cournot coefficients %.3f, %.3f; Bertrand %.3f, %.3f\n", hom$cournot_coef[1], hom$cournot_coef[2],
            hom$bertrand_coef[1], hom$bertrand_coef[2]))
cat("Saved res_logit.rds\n")
