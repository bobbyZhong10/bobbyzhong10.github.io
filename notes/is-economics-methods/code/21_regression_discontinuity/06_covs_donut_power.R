# =====================================================================
# 21_regression_discontinuity / 06_covs_donut_power.R
# Three add-ons of the local-polynomial toolkit, each on its own block.
#   R1c (aster_R1c.rds): covariates in RD (Calonico, Cattaneo, Farrell and Titiunik 2019).
#       nocovs / covs: rdrobust at 80 without and with the predetermined covariates gmv_prev
#       and tenure (each at its own MSE-optimal h); covs_same_h: with covariates at the
#       no-covariate h, so the precision gain is not a bandwidth effect. balance: rdrobust
#       with each covariate as the outcome (a predetermined covariate must not jump).
#   R1h (aster_R1h.rds): heaping at 80.0. heap: rddensity test and the counts at the heap;
#       heap_naive: rdrobust on all sellers; donut: rdrobust after dropping |x - 80| <= 0.5
#       on both sides (Barreca, Guldi, Lindo and Waddell 2011; the headline); donut_heap:
#       dropping only x = 80.0. No true-score benchmark: on x_true the heapers keep d = 1,
#       so the design on the true score is fuzzy (first stage about 0.48), and a sharp
#       rdrobust on x_true does not target the 0.20 jump (Task 14 final review, ruling 25).
#   R1 (aster_R1.rds): power. rdpower (Cattaneo, Titiunik and Vazquez-Bare 2019) of the
#       robust test for a 0.20 jump in the R1 design (the test's number), for a 0.10 jump,
#       and the minimum detectable effect at 80% power (rdmde with init.cond and a timeout:
#       its Newton loop can hang, see the track notes).
# Output: res_extra.rds. Run from this directory:  Rscript 06_covs_donut_power.R  (a few seconds)
# =====================================================================

suppressWarnings(suppressMessages({
  library(rdrobust)
  library(rddensity)
  library(rdpower)
  library(R.utils)
}))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
C  <- readRDS(file.path(here, "aster_R1c.rds"))
H  <- readRDS(file.path(here, "aster_R1h.rds"))
R1o <- readRDS(file.path(here, "aster_R1.rds")); R1 <- R1o$dat
tr <- readRDS(file.path(here, "truth21.rds"))
c_cut <- C$meta$cutoff                                  # 80 in all three blocks
stopifnot(H$meta$cutoff == c_cut, R1o$meta$c_cut == c_cut)
DONUT <- 0.5

quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressWarnings(expr))); out }
rd_list <- function(rd) {
  ci <- as.numeric(rd$ci["Robust", ]); cc <- as.numeric(rd$ci["Conventional", ])
  list(est = rd$coef["Conventional", 1], se = rd$se["Conventional", 1], est_bc = rd$coef["Robust", 1],
       se_rob = rd$se["Robust", 1], ci = ci, ci_len = diff(ci), ci_conv = cc, ci_len_conv = diff(cc),
       p_rob = rd$pv["Robust", 1], h = rd$bws["h", "left"], n_eff = as.integer(rd$N_h))
}
res <- list(truth = list(R1c = tr$R1c$jump, R1h = tr$R1h$jump, R1 = tr$R1$jump))

## ---- R1c: covariates ------------------------------------------------------------------------------
Z <- cbind(gmv_prev = C$gmv_prev, tenure = C$tenure)
r0 <- rdrobust(C$y, C$x, c = c_cut)
r1 <- rdrobust(C$y, C$x, c = c_cut, covs = Z)
r1h <- rdrobust(C$y, C$x, c = c_cut, covs = Z, h = r0$bws["h", "left"], b = r0$bws["b", "left"])
res$nocovs <- rd_list(r0)
res$covs <- rd_list(r1)
# rdrobust returns coef_covs without names, and not in the column order of covs: attach each
# name by matching the value to the same triangular-kernel weighted least squares fit (same h,
# both sides), whose coefficients carry the column names
h1 <- r1$bws["h", "left"]; stopifnot(isTRUE(all.equal(h1, r1$bws["h", "right"])))
s1 <- abs(C$x - c_cut) <= h1
w1 <- pmax(1 - abs(C$x - c_cut) / h1, 0)[s1]
d1 <- as.numeric(C$x >= c_cut)[s1]; x1 <- (C$x - c_cut)[s1]; Z1 <- Z[s1, , drop = FALSE]
wls <- coef(lm(C$y[s1] ~ d1 * x1 + Z1, weights = w1))
ref <- stats::setNames(as.numeric(wls[paste0("Z1", colnames(Z))]), colnames(Z))
stopifnot(abs(wls[["d1"]] - r1$coef["Conventional", 1]) < 1e-8)     # same fit as rdrobust
cc <- as.numeric(r1$coef_covs)
ord <- vapply(cc, function(v) which.min(abs(ref - v)), integer(1))
stopifnot(length(cc) == ncol(Z), !anyDuplicated(ord), max(abs(cc - ref[ord])) < 1e-8)
res$covs$coef_covs <- stats::setNames(cc, colnames(Z)[ord])
res$covs_same_h <- rd_list(r1h)
res$covs_ratio <- list(len_rob = res$covs$ci_len / res$nocovs$ci_len, len_conv = res$covs$ci_len_conv / res$nocovs$ci_len_conv,
                       len_rob_same_h = res$covs_same_h$ci_len / res$nocovs$ci_len,
                       est_change = res$covs$est - res$nocovs$est)
res$balance <- lapply(list(gmv_prev = C$gmv_prev, tenure = C$tenure), function(z) rd_list(rdrobust(z, C$x, c = c_cut)))
# share of the outcome variance near the cutoff that the covariates explain (|x - 80| <= 5,
# cubic in the score on each side): the source of the precision gain
nb <- abs(C$x - c_cut) <= 5
f0 <- lm(C$y[nb] ~ poly(C$x[nb] - c_cut, 3) * I(C$x[nb] >= c_cut))
f1 <- update(f0, . ~ . + Z[nb, ])
res$covs_partial_r2 <- 1 - sum(resid(f1)^2) / sum(resid(f0)^2)

## ---- R1h: heaping and the donut --------------------------------------------------------------------
dens <- quiet(rddensity(X = H$x, c = c_cut))
res$heap <- list(density_p = dens$test$p_jk, density_T = dens$test$t_jk,
                 n_at_80 = sum(H$x == c_cut), n_at_79 = sum(H$x == c_cut - 1),
                 n_at_79p9 = sum(H$x == round(c_cut - 0.1, 1)), n_at_80p1 = sum(H$x == round(c_cut + 0.1, 1)), n_heaped_true = sum(H$heaped_true), share_heaped = mean(H$heaped_true),
                 N = length(H$x), push_max = H$meta$push_max, push_width = H$meta$push_width)
# density test on the true score (not observable): the heap, not the world, makes the test reject
dens_true <- quiet(rddensity(X = H$x_true, c = c_cut))
res$heap$density_p_true <- dens_true$test$p_jk
# the R1h score is reported with one decimal: rdrobust warns about mass points and handles
# them with its default masspoints = "adjust"; quiet() drops the warning
res$heap_naive <- rd_list(quiet(rdrobust(H$y, H$x, c = c_cut)))
keep <- abs(H$x - c_cut) > DONUT
res$donut <- c(rd_list(quiet(rdrobust(H$y[keep], H$x[keep], c = c_cut))), list(radius = DONUT, n_dropped = sum(!keep)))
keep1 <- H$x != c_cut
res$donut_heap <- c(rd_list(quiet(rdrobust(H$y[keep1], H$x[keep1], c = c_cut))), list(n_dropped = sum(!keep1)))
# mean outcome of the heapers against sellers truly at 80 to 80.5: the heapers are worse
# sellers sitting at the cutoff, which pulls the right-hand limit down
res$heap$mean_y_heaped <- mean(H$y[H$heaped_true == 1])
res$heap$mean_y_true_80 <- mean(H$y[H$heaped_true == 0 & H$x >= c_cut & H$x <= c_cut + DONUT])
res$heap$mean_x_true_heaped <- mean(H$x_true[H$heaped_true == 1])

## ---- R1: power --------------------------------------------------------------------------------------
tau_pw <- tr$R1$jump                                    # power to detect the true 0.20 jump, and half of it
pw <- withTimeout(quiet(rdpower(data = cbind(R1$y_sharp, R1$score), cutoff = c_cut, tau = tau_pw)), timeout = 120)
pw1 <- withTimeout(quiet(rdpower(data = cbind(R1$y_sharp, R1$score), cutoff = c_cut, tau = tau_pw / 2)), timeout = 120)
mde <- withTimeout(quiet(rdmde(data = cbind(R1$y_sharp, R1$score), cutoff = c_cut, init.cond = 0.1)),
                   timeout = 60, onTimeout = "error")
res$power <- list(power = pw$power.rbc, power_conv = pw$power.conv, se_rbc = pw$se.rbc, tau = tau_pw,
                  power_tau01 = pw1$power.rbc, power_conv_tau01 = pw1$power.conv,
                  h = pw$samph.l, n_left = pw$Nh.l, n_right = pw$Nh.r,
                  mde = mde$mde, mde_conv = mde$mde.conv, beta = 0.8, alpha = 0.05)

saveRDS(res, file.path(here, "res_extra.rds"))

## ---- report ---------------------------------------------------------------------------------------------
ci_txt <- function(ci) sprintf("[%.3f, %.3f]", ci[1], ci[2])
cat(sprintf("================ R1c: covariates (truth %.3f) ================\n", tr$R1c$jump))
for (nm in c("nocovs", "covs", "covs_same_h")) {
  r <- res[[nm]]
  cat(sprintf("%-12s %.4f robust %s length %.3f (conventional %.3f)  h %.3f  eff. n %d / %d\n", nm, r$est, ci_txt(r$ci),
              r$ci_len, r$ci_len_conv, r$h, r$n_eff[1], r$n_eff[2]))
}
cat(sprintf("length ratio robust %.3f, conventional %.3f, robust at the same h %.3f; estimate change %.4f; partial R2 %.3f\n",
            res$covs_ratio$len_rob, res$covs_ratio$len_conv, res$covs_ratio$len_rob_same_h, res$covs_ratio$est_change,
            res$covs_partial_r2))
for (nm in names(res$balance))
  cat(sprintf("balance %-9s jump %.4f robust %s p %.3f\n", nm, res$balance[[nm]]$est, ci_txt(res$balance[[nm]]$ci),
              res$balance[[nm]]$p_rob))
cat(sprintf("================ R1h: heaping at %.1f (truth %.3f) ================\n", c_cut, tr$R1h$jump))
cat(sprintf("sellers at 80.0: %d (79.9: %d, 80.1: %d, 79.0: %d); heapers %d (%.1f%%)\n", res$heap$n_at_80,
            res$heap$n_at_79p9, res$heap$n_at_80p1, res$heap$n_at_79, res$heap$n_heaped_true, 100 * res$heap$share_heaped))
cat(sprintf("rddensity p %.3g (T %.2f); on the true score p %.3f\n", res$heap$density_p, res$heap$density_T,
            res$heap$density_p_true))
for (nm in c("heap_naive", "donut", "donut_heap"))
  cat(sprintf("%-11s %.4f robust %s  h %.3f\n", nm, res[[nm]]$est, ci_txt(res[[nm]]$ci), res[[nm]]$h))
cat(sprintf("mean y of heapers %.3f vs sellers truly at 80 to 80.5 %.3f (heapers' mean true score %.2f)\n",
            res$heap$mean_y_heaped, res$heap$mean_y_true_80, res$heap$mean_x_true_heaped))
cat("================ R1: power ================\n")
cat(sprintf("rdpower tau %.2f: robust %.4f conventional %.4f (h %.3f, n %d / %d); half of it: robust %.4f\n", tau_pw, res$power$power,
            res$power$power_conv, res$power$h, res$power$n_left, res$power$n_right, res$power$power_tau01))
cat(sprintf("rdmde at 80%% power: robust %.4f conventional %.4f\n", res$power$mde, res$power$mde_conv))
