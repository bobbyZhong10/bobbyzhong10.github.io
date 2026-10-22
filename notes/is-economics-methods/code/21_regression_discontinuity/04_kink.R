# =====================================================================
# 21_regression_discontinuity / 04_kink.R
# Block R3 (aster_R3.rds): regression kink in the commission schedule. The rate b(x) is 12%
# up to USD 1,500 of last quarter's average monthly GMV x and falls 0.6 percentage points
# per USD 100 above it (slope change -0.006 pp per USD, known from the schedule). True kink
# effect -0.05 (log GMV per percentage point of commission); no level jump.
#   rk          headline: rdrobust deriv = 1 with p = 2 (Calonico, Cattaneo and Titiunik 2014
#               kink estimator: local quadratic for the first derivative), slope change of
#               E[y | x] divided by the known slope change of b; robust interval
#   rk_p1       rdrobust deriv = 1, p = 1 (local linear slope change), for comparison
#   rk_hand     rk_ratio() (rd_fns.R, local linear, triangular) at the p = 1 bandwidth: the
#               hand computation of Card, Lee, Pei and Weber (2015); equals rk_p1's
#               conventional estimate, and shows the local-linear slope bias at that h
#   rk_fuzzy    rdrobust deriv = 1, p = 2, fuzzy = b (fuzzy kink: the slope change of b is
#               estimated; with a deterministic schedule it equals the known -0.006)
#   level_jump  rdrobust deriv = 0, p = 1 at 1,500: the placebo jump in the level (truth 0)
# Output: res_kink.rds. Run from this directory:  Rscript 04_kink.R   (a few seconds)
# =====================================================================

suppressWarnings(suppressMessages({
  library(rdrobust)
}))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
source(file.path(here, "rd_fns.R"))
B <- readRDS(file.path(here, "aster_R3.rds"))
k0 <- B$meta$cutoff
db <- B$meta$slope_right - B$meta$slope_left       # known slope change of the schedule (-0.006)
truth <- B$meta$kink_effect

# the schedule in the data is the published one (sharp kink): check before dividing by db
bf <- B$meta$b_fun
stopifnot(isTRUE(all.equal(B$b, bf(B$x), tolerance = 1e-12)),
          isTRUE(all.equal((bf(k0 + 1) - bf(k0)) - (bf(k0) - bf(k0 - 1)), db)))   # the schedule's own kink

# a derivative estimate divided by db (< 0): the interval bounds swap
kink_list <- function(rd, p) {
  ci <- sort(as.numeric(rd$ci["Robust", ]) / db)
  list(est = rd$coef["Conventional", 1] / db, se = rd$se["Conventional", 1] / abs(db),
       est_bc = rd$coef["Robust", 1] / db, se_rob = rd$se["Robust", 1] / abs(db), ci = ci,
       ci_conv = sort(as.numeric(rd$ci["Conventional", ]) / db),
       dslope_y = rd$coef["Conventional", 1], h = rd$bws["h", "left"], b = rd$bws["b", "left"],
       n_eff = as.integer(rd$N_h), p = p)
}

rk2 <- rdrobust(B$y, B$x, c = k0, deriv = 1, p = 2)
rk1 <- rdrobust(B$y, B$x, c = k0, deriv = 1, p = 1)
res <- list(truth = list(kink_effect = truth, dslope_b = db, dslope_y = truth * db, level_jump = 0, cutoff = k0),
            sample = list(N = length(B$x), n_left = sum(B$x < k0), n_right = sum(B$x >= k0),
                          x_min = min(B$x), x_max = max(B$x), floor_from = B$meta$floor_from))
res$rk    <- kink_list(rk2, 2)
res$rk_p1 <- kink_list(rk1, 1)

hand <- rk_ratio(B$y, B$b, B$x, k0, h = res$rk_p1$h)
res$rk_hand <- list(est = hand$est, dslope_y = hand$dslope_y, dslope_b = hand$dslope_b, h = res$rk_p1$h)
stopifnot(abs(hand$dslope_b - db) < 1e-9,                       # the schedule's slope change, recovered exactly
          abs(hand$est - res$rk_p1$est) < 1e-6)                 # same local-linear fit as rdrobust p = 1

fz <- rdrobust(B$y, B$x, c = k0, deriv = 1, p = 2, fuzzy = B$b)
res$rk_fuzzy <- list(est = fz$coef["Conventional", 1], se = fz$se["Conventional", 1],
                     ci = as.numeric(fz$ci["Robust", ]), h = fz$bws["h", "left"])

lv <- rdrobust(B$y, B$x, c = k0)
res$level_jump <- list(est = lv$coef["Conventional", 1], se = lv$se["Conventional", 1],
                       ci = as.numeric(lv$ci["Robust", ]), ci_conv = as.numeric(lv$ci["Conventional", ]),
                       p_rob = lv$pv["Robust", 1], h = lv$bws["h", "left"], n_eff = as.integer(lv$N_h))
res$level_jump$t <- res$level_jump$est / res$level_jump$se

saveRDS(res, file.path(here, "res_kink.rds"))

## ---- report ---------------------------------------------------------------------------------------------
ci_txt <- function(ci, d = 4) sprintf(paste0("[%.", d, "f, %.", d, "f]"), ci[1], ci[2])
cat(sprintf("================ R3: regression kink at USD %d (truth %.3f) ================\n", k0, truth))
cat(sprintf("rdrobust deriv 1, p 2 (headline): %.4f  robust %s  se %.4f  h %.0f  eff. n %d / %d\n", res$rk$est,
            ci_txt(res$rk$ci), res$rk$se, res$rk$h, res$rk$n_eff[1], res$rk$n_eff[2]))
cat(sprintf("rdrobust deriv 1, p 1           : %.4f  robust %s  h %.0f\n", res$rk_p1$est, ci_txt(res$rk_p1$ci), res$rk_p1$h))
cat(sprintf("hand rk_ratio at h %.0f          : %.4f  (dslope_y %.3e / dslope_b %.4f)\n", res$rk_hand$h, res$rk_hand$est,
            res$rk_hand$dslope_y, res$rk_hand$dslope_b))
cat(sprintf("fuzzy kink (fuzzy = b), p 2     : %.4f  robust %s\n", res$rk_fuzzy$est, ci_txt(res$rk_fuzzy$ci)))
cat(sprintf("level jump at the kink          : %.4f  se %.4f  (%.2f se)  robust %s  h %.0f\n", res$level_jump$est,
            res$level_jump$se, res$level_jump$t, ci_txt(res$level_jump$ci), res$level_jump$h))
