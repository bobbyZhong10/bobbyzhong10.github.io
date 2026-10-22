# ------------------------------------------------------------------
# 06_unbalanced.R -- block P6: unbalanced panels and attrition (Wooldridge 2010 §19.9).
#   FE on the full balanced panel, on the panel thinned by attrition on alpha_i (s_alpha) and on the
#   panel thinned by attrition on the shock (s_eps); the unbalanced Mundlak identity (Xbar_i over the
#   T_i observed months) and its two-way version (with the unit averages of the month dummies), and the
#   naive two-way Mundlak regression without them; RE on the unbalanced panel (varphi_i varies with T_i)
#   and the robust Hausman test.
# Reads nimbus_P6.rds, truth10.rds; writes res_p6.rds. "se" is cluster-robust by merchant (CR0).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("panel_fns.R")
suppressPackageStartupMessages(library(plm))
tr <- readRDS("truth10.rds"); P6 <- tr$params$P6
d6 <- readRDS("nimbus_P6.rds")
inside <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])

fe_row <- function(d) {
  f <- fe_fit(d$y, d$x, d$id); e <- unname(f$b); se <- sqrt(f$V[1, 1])
  Ti <- tabulate(match(d$id, unique(d$id)))
  list(est = e, se = se, ci = ci95(e, se), covers_truth = inside(ci95(e, se), P6$beta),
       n = nrow(d), G = length(Ti), mean_Ti = mean(Ti), share_Ti_full = mean(Ti == P6$T),
       share_Ti_1 = mean(Ti == 1), keep = nrow(d) / nrow(d6))
}
full <- d6
sa <- d6[d6$s_alpha == 1, ]
se_ <- d6[d6$s_eps == 1, ]
res <- list(full = fe_row(full), attr_alpha = fe_row(sa), attr_eps = fe_row(se_))
res$attr_eps$bias <- res$attr_eps$est - P6$beta
res$attr_eps$z_vs_truth <- (res$attr_eps$est - P6$beta) / res$attr_eps$se
# who leaves: merchant effect of leavers and stayers (month T), GMV of stayers
res$alpha_leavers_alpha <- mean(d6$alpha_true[d6$t == P6$T & d6$s_alpha == 0])
res$alpha_stayers_alpha <- mean(d6$alpha_true[d6$t == P6$T & d6$s_alpha == 1])
res$mean_y_full <- mean(d6$y)
res$share_below_floor <- mean(d6$y[d6$t > 1] < P6$floor)

## ---- unbalanced Mundlak identity, on the panel with attrition on alpha_i ----------------------
pa <- pdata.frame(sa, index = c("id", "t"))
fe_p <- plm(y ~ x, data = pa, model = "within")
tw_p <- plm(y ~ x, data = pa, model = "within", effect = "twoways")
mk <- mundlak_fit(sa$y, sa$x, sa$id)
mk_tw <- mundlak_fit(sa$y, sa$x, sa$id, time_fe = TRUE, t = sa$t)
xbar <- ave(sa$x, sa$id)
naive <- unname(coef(lm(y ~ x + xbar + factor(t), data = data.frame(y = sa$y, x = sa$x, xbar = xbar, t = sa$t)))["x"])
re_p <- plm(y ~ x, data = pa, model = "random")
th <- ercomp(re_p)$theta
res$re_alpha <- list(est = unname(coef(re_p)["x"]),
                     se = sqrt(vcovHC(re_p, method = "arellano", type = "HC0", cluster = "group")["x", "x"]),
                     varphi_min = min(th), varphi_max = max(th))
res$haus_robust_alpha <- hausman_robust(mk)
res$mundlak_se_x <- sqrt(mk$V["x", "x"])

out <- list(P6 = res, fe_unbal = unname(coef(fe_p)), mundlak_unbal = mk, twfe_unbal = unname(coef(tw_p)),
            mundlak_tw_unbal = mk_tw, mundlak_tw_naive_unbal = naive)
saveRDS(out, "res_p6.rds")
cat(sprintf("P6 FE (truth %.2f): full %.4f [%.4f, %.4f]; attrition on alpha %.4f [%.4f, %.4f] (kept %.3f); on the shock %.4f [%.4f, %.4f] (kept %.3f)\n",
            P6$beta, res$full$est, res$full$ci[1], res$full$ci[2], res$attr_alpha$est, res$attr_alpha$ci[1],
            res$attr_alpha$ci[2], res$attr_alpha$keep, res$attr_eps$est, res$attr_eps$ci[1], res$attr_eps$ci[2], res$attr_eps$keep))
cat(sprintf("   unbalanced Mundlak b_x %.10f vs FE %.10f; two-way %.10f vs TWFE %.10f; naive two-way %.10f\n",
            mk$b_x, out$fe_unbal, mk_tw$b_x, out$twfe_unbal, naive))
cat(sprintf("   RE %.4f, varphi_i from %.3f to %.3f; robust Hausman %.1f (p %.2g)\n",
            res$re_alpha$est, res$re_alpha$varphi_min, res$re_alpha$varphi_max, res$haus_robust_alpha$stat, res$haus_robust_alpha$p))
