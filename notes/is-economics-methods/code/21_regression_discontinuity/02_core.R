# =====================================================================
# 21_regression_discontinuity / 02_core.R
# Block R1 (the inherited Preferred-badge world, aster_R1.rds): every inherited number of
# chapter 21, computed with the legacy calls of 02_rdd_estimate.R (which this script
# replaces), plus two bandwidth checks for the local-polynomial section.
#   (1) naive above-vs-below difference (selection-biased);
#   (2) sharp local-linear RD (rdrobust, triangular kernel, CCT / MSE-optimal bandwidth),
#       conventional and bias-corrected robust inference;
#   (3) bandwidth sensitivity (h / 2, 2 h), CER-optimal bandwidth, local quadratic;
#   (4) global polynomials of order 4, 6, 8 (Gelman and Imbens 2019 warning);
#   (5) rddensity manipulation test;
#   (6) fuzzy RD: first stage, reduced form, Wald ratio, rdrobust fuzzy LATE.
# Output: res_core.rds (fields read by tests/test_reproduce.R and numbers_def.R).
# Run from this directory:  Rscript 02_core.R   (a few seconds)
# =====================================================================

suppressWarnings(suppressMessages({
  library(rdrobust)
  library(rddensity)
}))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
obj  <- readRDS(file.path(here, "aster_R1.rds"))
dat  <- obj$dat
meta <- obj$meta
c_cut <- meta$c_cut

# rdrobust prints nothing, but rddensity and friends can; keep the log clean
quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressWarnings(expr))); out }
ci_of <- function(rd, row) as.numeric(rd$ci[row, ])

res <- list(truth = list(jump = meta$tau_true, first_stage = meta$first_stage_jump_true, cutoff = c_cut),
            sample = list(N = meta$N, n_below = meta$n_below, n_above = meta$n_above,
                          frac_badge_below = meta$frac_badge_below, frac_badge_above = meta$frac_badge_above))

# ---- (1) naive above-vs-below difference (stored by 01_dgp.R, legacy definition) ----
res$naive <- meta$naive_diff

# ---- (2) sharp local-linear RD ------------------------------------------------------
rd <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut,
               kernel = "triangular", p = 1, bwselect = "mserd")
h <- rd$bws["h", "left"]
res$sharp <- list(est = rd$coef["Conventional", 1], se = rd$se["Conventional", 1],
                  est_bc = rd$coef["Robust", 1], se_rob = rd$se["Robust", 1],
                  ci_conv = ci_of(rd, "Conventional"), ci_rob = ci_of(rd, "Robust"),
                  p_conv = rd$pv["Conventional", 1], p_rob = rd$pv["Robust", 1],
                  h = h, b = rd$bws["b", "left"], n_eff = as.integer(rd$N_h),
                  window = c(c_cut - h, c_cut + h))

# ---- (3) bandwidth checks -------------------------------------------------------------
rd_half <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut, kernel = "triangular", p = 1, h = h / 2)
rd_dbl  <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut, kernel = "triangular", p = 1, h = h * 2)
res$bw_sens <- c(rd_half$coef["Conventional", 1], rd_dbl$coef["Conventional", 1])
res$bw_sens_detail <- list(
  half = list(h = h / 2, est = rd_half$coef["Conventional", 1], ci_rob = ci_of(rd_half, "Robust"),
              n_eff = as.integer(rd_half$N_h)),
  dbl  = list(h = h * 2, est = rd_dbl$coef["Conventional", 1], ci_rob = ci_of(rd_dbl, "Robust"),
              n_eff = as.integer(rd_dbl$N_h)))
# CER-optimal bandwidth (Calonico, Cattaneo and Farrell 2020): shorter than the MSE-optimal
# one, chosen for the coverage error of the robust interval rather than for the MSE of the
# point estimate
rd_cer <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut, kernel = "triangular", p = 1, bwselect = "cerrd")
res$cer <- list(h = rd_cer$bws["h", "left"], est = rd_cer$coef["Conventional", 1],
                ci_conv = ci_of(rd_cer, "Conventional"), ci_rob = ci_of(rd_cer, "Robust"),
                n_eff = as.integer(rd_cer$N_h))
# local quadratic at its own MSE-optimal bandwidth
rd_p2 <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut, kernel = "triangular", p = 2, bwselect = "mserd")
res$p2 <- list(h = rd_p2$bws["h", "left"], est = rd_p2$coef["Conventional", 1],
               ci_conv = ci_of(rd_p2, "Conventional"), ci_rob = ci_of(rd_p2, "Robust"),
               n_eff = as.integer(rd_p2$N_h))

# ---- (4) global polynomials (legacy code) ------------------------------------------------
# Rescale the running variable to [-1, 1] and use RAW powers, so that at the cutoff (x = 0)
# every power vanishes and the coefficient on D is exactly the estimated jump.
sc_c <- (dat$score - c_cut) / max(abs(dat$score - c_cut))
poly_jump <- function(k) {
  df <- data.frame(y = dat$y_sharp, D = dat$D_sharp, x = sc_c)
  fit <- lm(y ~ D * poly(x, k, raw = TRUE), data = df)
  as.numeric(coef(fit)["D"])
}
res$poly <- c(poly_jump(4), poly_jump(6), poly_jump(8))
names(res$poly) <- c("p4", "p6", "p8")

# ---- (5) rddensity manipulation test --------------------------------------------------------
dens <- quiet(rddensity(X = dat$score, c = c_cut))
res$density <- list(p = dens$test$p_jk, T = dens$test$t_jk,
                    h_left = dens$h$left, h_right = dens$h$right)

# ---- (6) fuzzy RD ------------------------------------------------------------------------------
rd_fz <- rdrobust(y = dat$y_fuzzy, x = dat$score, c = c_cut, fuzzy = dat$D_fuzzy,
                  kernel = "triangular", p = 1, bwselect = "mserd")
rd_fs <- rdrobust(y = dat$D_fuzzy, x = dat$score, c = c_cut, kernel = "triangular", p = 1, bwselect = "mserd")
rd_rf <- rdrobust(y = dat$y_fuzzy, x = dat$score, c = c_cut, kernel = "triangular", p = 1, bwselect = "mserd")
res$fuzzy <- list(est = rd_fz$coef["Conventional", 1], se = rd_fz$se["Conventional", 1],
                  est_bc = rd_fz$coef["Robust", 1],
                  ci_conv = ci_of(rd_fz, "Conventional"), ci_rob = ci_of(rd_fz, "Robust"),
                  h = rd_fz$bws["h", "left"], n_eff = as.integer(rd_fz$N_h),
                  fs = rd_fs$coef["Conventional", 1], fs_ci_rob = ci_of(rd_fs, "Robust"), fs_h = rd_fs$bws["h", "left"],
                  rf = rd_rf$coef["Conventional", 1], rf_ci_rob = ci_of(rd_rf, "Robust"), rf_h = rd_rf$bws["h", "left"])
res$fuzzy$wald <- res$fuzzy$rf / res$fuzzy$fs      # hand Wald ratio (first stage and reduced form at their own h)

saveRDS(res, file.path(here, "res_core.rds"))

# ---- report (legacy layout and precision) ----------------------------------------------------
f3 <- function(v) sprintf("%.3f", v)
cat("================ R1: inherited RD estimates ================\n")
cat(sprintf("true jump                     : %.3f\n", res$truth$jump))
cat(sprintf("naive above-vs-below diff     : %.3f  (%.1f times the truth)\n", res$naive, res$naive / res$truth$jump))
cat(sprintf("rdrobust conventional est     : %.3f  (se %.3f)\n", res$sharp$est, res$sharp$se))
cat(sprintf("  conventional 95%% CI          : [%.3f, %.3f]\n", res$sharp$ci_conv[1], res$sharp$ci_conv[2]))
cat(sprintf("rdrobust robust 95%% CI        : [%.3f, %.3f]  (p=%.3g)\n", res$sharp$ci_rob[1], res$sharp$ci_rob[2], res$sharp$p_rob))
cat(sprintf("CCT/MSE-optimal bandwidth     : %.3f  (window %.3f to %.3f)\n", h, res$sharp$window[1], res$sharp$window[2]))
cat(sprintf("  eff. n (left/right)          : %d / %d\n", res$sharp$n_eff[1], res$sharp$n_eff[2]))
cat(sprintf("bandwidth sens (h/2, 2h)      : %.3f , %.3f\n", res$bw_sens[1], res$bw_sens[2]))
cat(sprintf("CER-optimal bandwidth         : %.3f  est %.3f robust CI [%.3f, %.3f]\n", res$cer$h, res$cer$est,
            res$cer$ci_rob[1], res$cer$ci_rob[2]))
cat(sprintf("local quadratic (own h %.3f)  : %.3f robust CI [%.3f, %.3f]\n", res$p2$h, res$p2$est, res$p2$ci_rob[1], res$p2$ci_rob[2]))
cat(sprintf("global poly p=4 / 6 / 8       : %s / %s / %s  (unstable)\n", f3(res$poly[1]), f3(res$poly[2]), f3(res$poly[3])))
cat(sprintf("McCrary/rddensity p           : %.3f  (T=%.3f; large p = no manip.)\n", res$density$p, res$density$T))
cat("---------------- fuzzy RD ----------------\n")
cat(sprintf("first-stage jump  est / true  : %.3f / %.3f\n", res$fuzzy$fs, res$truth$first_stage))
cat(sprintf("reduced-form jump est         : %.3f\n", res$fuzzy$rf))
cat(sprintf("hand Wald ratio rf / fs       : %.3f\n", res$fuzzy$wald))
cat(sprintf("fuzzy-RD LATE est             : %.3f  (bw %.3f)\n", res$fuzzy$est, res$fuzzy$h))
cat(sprintf("  fuzzy robust 95%% CI          : [%.3f, %.3f]\n", res$fuzzy$ci_rob[1], res$fuzzy$ci_rob[2]))
