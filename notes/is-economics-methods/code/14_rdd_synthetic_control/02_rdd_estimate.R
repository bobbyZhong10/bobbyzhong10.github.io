# =====================================================================
# 14_rdd_synthetic_control / 02_rdd_estimate.R
# Estimate the Preferred-badge effect four ways and contrast them:
#   (1) naive above-vs-below difference (biased);
#   (2) local-linear RDD via rdrobust (CCT bandwidth, conventional +
#       bias-corrected robust CI);
#   (3) a high-order GLOBAL polynomial (to SHOW overfitting);
#   (4) fuzzy RD (rdrobust with the badge as fuzzy treatment) = the
#       cutoff IV / LATE, tying back to Ch.10.
# Plus a McCrary / rddensity manipulation test.
# =====================================================================

suppressWarnings(suppressMessages({
  library(rdrobust)
  library(rddensity)
}))

in_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(in_dir) == 0 || in_dir == "") in_dir <- "."
obj  <- readRDS(file.path(in_dir, "data_rdd.rds"))
dat  <- obj$dat
meta <- obj$meta
c_cut <- meta$c_cut

res <- list()
res$tau_true   <- meta$tau_true
res$naive_diff <- meta$naive_diff

# ---- (2) local-linear RDD, sharp -----------------------------------
rd <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut,
               kernel = "triangular", p = 1, bwselect = "mserd")
res$bw_cct        <- rd$bws["h", "left"]      # CCT/MSE-optimal bandwidth
res$rd_conv_est   <- rd$coef["Conventional", 1]
res$rd_conv_se    <- rd$se["Conventional", 1]
res$rd_conv_ci    <- rd$ci["Conventional", ]
res$rd_robust_est <- rd$coef["Robust", 1]     # point est same, robust inference
res$rd_robust_ci  <- rd$ci["Robust", ]
res$rd_robust_p   <- rd$pv["Robust", 1]
res$n_eff_left    <- rd$N_h[1]
res$n_eff_right   <- rd$N_h[2]

# ---- bandwidth sensitivity: half and double the CCT bandwidth ------
rd_half <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut,
                    kernel = "triangular", p = 1, h = res$bw_cct / 2)
rd_dbl  <- rdrobust(y = dat$y_sharp, x = dat$score, c = c_cut,
                    kernel = "triangular", p = 1, h = res$bw_cct * 2)
res$rd_bw_half <- rd_half$coef["Conventional", 1]
res$rd_bw_dbl  <- rd_dbl$coef["Conventional", 1]

# ---- (3) high-order GLOBAL polynomial (Gelman-Imbens warning) ------
# Fit y ~ D + poly(score,k) + D:poly(score,k) over the WHOLE support;
# the treatment coefficient is the estimated jump at the cutoff. The
# running variable is rescaled to roughly [-1,1] so the high powers do
# not overflow -- the instability that remains is the genuine Runge /
# boundary overfitting Gelman & Imbens warn about.
# Rescale the running variable to [-1,1] and use RAW powers so that at
# the cutoff (x = 0) every power vanishes -> the coefficient on D is
# exactly the estimated jump at the cutoff.
sc_c <- (dat$score - c_cut) / max(abs(dat$score - c_cut))
poly_jump <- function(k) {
  df <- data.frame(y = dat$y_sharp, D = dat$D_sharp, x = sc_c)
  fit <- lm(y ~ D * poly(x, k, raw = TRUE), data = df)
  coef(fit)["D"]
}
res$poly4_est <- as.numeric(poly_jump(4))
res$poly6_est <- as.numeric(poly_jump(6))
res$poly8_est <- as.numeric(poly_jump(8))

# ---- (4) fuzzy RD: badge is the (imperfect) treatment --------------
rd_fz <- rdrobust(y = dat$y_fuzzy, x = dat$score, c = c_cut,
                  fuzzy = dat$D_fuzzy, kernel = "triangular", p = 1,
                  bwselect = "mserd")
res$fuzzy_est       <- rd_fz$coef["Conventional", 1]
res$fuzzy_robust_ci <- rd_fz$ci["Robust", ]
res$fuzzy_bw        <- rd_fz$bws["h", "left"]

# first-stage jump (compliance) via a sharp RD on the badge itself
rd_fs <- rdrobust(y = dat$D_fuzzy, x = dat$score, c = c_cut,
                  kernel = "triangular", p = 1, bwselect = "mserd")
res$first_stage_est  <- rd_fs$coef["Conventional", 1]
res$first_stage_true <- meta$first_stage_jump_true

# reduced-form jump in the outcome
rd_rf <- rdrobust(y = dat$y_fuzzy, x = dat$score, c = c_cut,
                  kernel = "triangular", p = 1, bwselect = "mserd")
res$reduced_form_est <- rd_rf$coef["Conventional", 1]

# ---- McCrary / rddensity manipulation test -------------------------
dens <- rddensity(X = dat$score, c = c_cut)
res$mccrary_p <- dens$test$p_jk    # jackknife p-value of the density jump
res$mccrary_T <- dens$test$t_jk

saveRDS(res, file.path(in_dir, "res_rdd.rds"))

# ---- report --------------------------------------------------------
cat("================ RDD estimates ================\n")
cat(sprintf("true jump                     : %.3f\n", res$tau_true))
cat(sprintf("naive above-vs-below diff     : %.3f\n", res$naive_diff))
cat(sprintf("rdrobust conventional est     : %.3f  (se %.3f)\n",
            res$rd_conv_est, res$rd_conv_se))
cat(sprintf("  conventional 95%% CI          : [%.3f, %.3f]\n",
            res$rd_conv_ci[1], res$rd_conv_ci[2]))
cat(sprintf("rdrobust robust 95%% CI        : [%.3f, %.3f]  (p=%.3g)\n",
            res$rd_robust_ci[1], res$rd_robust_ci[2], res$rd_robust_p))
cat(sprintf("CCT/MSE-optimal bandwidth     : %.3f\n", res$bw_cct))
cat(sprintf("  eff. n (left/right)          : %d / %d\n",
            res$n_eff_left, res$n_eff_right))
cat(sprintf("bandwidth sens (h/2, 2h)      : %.3f , %.3f\n",
            res$rd_bw_half, res$rd_bw_dbl))
cat(sprintf("global poly p=4 / 6 / 8       : %.3f / %.3f / %.3f  (unstable)\n",
            res$poly4_est, res$poly6_est, res$poly8_est))
cat(sprintf("McCrary/rddensity p           : %.3f  (T=%.3f; large p = no manip.)\n",
            res$mccrary_p, res$mccrary_T))
cat("---------------- fuzzy RD ----------------\n")
cat(sprintf("first-stage jump  est / true  : %.3f / %.3f\n",
            res$first_stage_est, res$first_stage_true))
cat(sprintf("reduced-form jump est         : %.3f\n", res$reduced_form_est))
cat(sprintf("fuzzy-RD LATE est             : %.3f  (bw %.3f)\n",
            res$fuzzy_est, res$fuzzy_bw))
cat(sprintf("  fuzzy robust 95%% CI          : [%.3f, %.3f]\n",
            res$fuzzy_robust_ci[1], res$fuzzy_robust_ci[2]))
