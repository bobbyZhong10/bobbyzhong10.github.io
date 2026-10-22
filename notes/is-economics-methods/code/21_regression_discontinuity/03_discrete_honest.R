# =====================================================================
# 21_regression_discontinuity / 03_discrete_honest.R
# Block R2 (aster_R2.rds): the running variable is the average of a seller's 10 most recent
# review stars (21 support points 3.0, ..., 5.0); the "4.5 stars & up" filter starts at
# 4.25; true jump 0.15. Five ways to build a 95% interval for the jump, on this draw and in
# a Monte Carlo over fresh draws of the same world (gen_R2):
#   conv      rdrobust local linear, MSE-optimal h, masspoints = "adjust", Conventional row
#   robust    the same fit, bias-corrected robust interval (Calonico, Cattaneo, Titiunik 2014)
#   crv       the conventional fit at the same h with standard errors clustered by support
#             point (Lee and Card 2008; criticised by Kolesar and Rothe 2018)
#   honest    RDHonest (Armstrong and Kolesar 2018, 2020; Kolesar and Rothe 2018): bias-aware
#             interval valid when |second derivative of E[y | x]| <= M on each side, at
#               M = meta$M_bound (7.79, the DGP's true bound: the theorem's premise),
#               M = RDSmoothnessBound lower bound (s = 2; the smallest M the data cannot
#                   reject, a data-driven choice), and M = RDSmoothnessBound point estimate
#   locrand   local randomization in the window rdwinselect picks with the predetermined
#             covariate tenure (Cattaneo, Frandsen and Titiunik 2015): difference in means,
#             randomization p-value and test-inversion interval (rdrandinf), plus the hand
#             lr_diffmeans and the window sensitivity (rdsensitivity)
# Monte Carlo (ruling R11): REPS draws of gen_R2 at the block's n; per-draw seeds are drawn from the
# master seed MASTER_SEED, so the results do not depend on the number of cores. Coverage of
# 0.15 and mean interval length for conv, robust, crv, honest (M_bound; the test's column),
# honest at the per-draw RDSmoothnessBound lower bound and point estimate, and local
# randomization in the window rdwinselect picks on the frozen draw ([4.2, 4.3]), held fixed
# across draws (Neyman / Welch interval; no threshold).
# Output: res_discrete.rds. Run from this directory:  Rscript 03_discrete_honest.R
# (about 70 s on 16 cores; the Monte Carlo is most of it, the rdrandinf interval about 7 s)
# =====================================================================

suppressWarnings(suppressMessages({
  library(rdrobust)
  library(RDHonest)
  library(rdlocrand)
  library(parallel)
}))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
source(file.path(here, "rd_fns.R"))
B <- readRDS(file.path(here, "aster_R2.rds"))
c0 <- B$meta$cutoff
M_bound <- B$meta$M_bound
tau <- B$meta$jump

REPS <- 1000
MASTER_SEED <- 21030
LR_SEED <- 2103           # permutation seed of rdwinselect, rdrandinf, rdsensitivity, lr_diffmeans
# on Windows set CORES <- 1L: mclapply cannot fork there (see 00_run_all.R)
CORES <- max(1L, min(16L, parallel::detectCores() - 2L))

quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressWarnings(expr))); out }
ci_row <- function(rd, row) as.numeric(rd$ci[row, ])

# all interval estimators on one sample (used for the frozen draw and inside the Monte Carlo)
fit_all <- function(y, x, M) {
  rd <- quiet(rdrobust(y, x, c = c0, masspoints = "adjust"))
  h <- rd$bws["h", "left"]
  # same h and pilot b as the conventional fit (b only enters the bias correction, which the
  # Conventional row does not use; with b = h a small h leaves too few support points for it)
  cr <- quiet(rdrobust(y, x, c = c0, masspoints = "adjust", h = h, b = rd$bws["b", "left"], cluster = x, vce = "hc1"))
  df <- data.frame(y = y, x = x)
  ho <- RDHonest(y ~ x, data = df, cutoff = c0, M = M)
  sb <- RDSmoothnessBound(ho, s = 2)
  m_lb <- max(sb$conf.low, 0)
  hl <- RDHonest(y ~ x, data = df, cutoff = c0, M = m_lb)
  he <- RDHonest(y ~ x, data = df, cutoff = c0, M = sb$estimate)
  list(rd = rd, cr = cr, ho = ho$coefficients, hl = hl$coefficients, he = he$coefficients,
       m_lb = m_lb, m_est = sb$estimate, h = h)
}
honest_list <- function(hc) list(est = hc$estimate, se = hc$std.error, max_bias = hc$maximum.bias,
                                 ci = c(hc$conf.low, hc$conf.high), h = hc$bandwidth, eff_obs = hc$eff.obs,
                                 M = hc$M, p = hc$p.value)

## ---- the frozen draw (aster_R2.rds) --------------------------------------------------------------
fa <- fit_all(B$y, B$x, M_bound)
rd <- fa$rd
res <- list(truth = list(jump = tau, M_bound = M_bound, cutoff = c0))
res$sample <- list(N = length(B$x), n_support = length(unique(B$x)), min_mass = min(table(B$x)),
                   n_left = sum(B$x < c0), n_right = sum(B$x >= c0),
                   share_at_5 = mean(B$x == 5), n_support_right = length(unique(B$x[B$x >= c0])))
res$conv   <- list(est = rd$coef["Conventional", 1], se = rd$se["Conventional", 1], ci = ci_row(rd, "Conventional"),
                   h = fa$h, b = rd$bws["b", "left"], n_eff = as.integer(rd$N_h),
                   n_support_in_h = c(left = length(unique(B$x[B$x < c0 & B$x >= c0 - fa$h])),
                                      right = length(unique(B$x[B$x >= c0 & B$x <= c0 + fa$h]))))
res$robust <- list(est = rd$coef["Robust", 1], se = rd$se["Robust", 1], ci = ci_row(rd, "Robust"))
res$crv    <- list(est = fa$cr$coef["Conventional", 1], se = fa$cr$se["Conventional", 1],
                   ci = ci_row(fa$cr, "Conventional"), h = fa$h, n_clusters = length(unique(B$x[abs(B$x - c0) <= fa$h])))
stopifnot(abs(res$crv$est - res$conv$est) < 1e-10)      # same fit, only the standard error differs
res$honest      <- honest_list(fa$ho)                    # M = M_bound (conservative, true bound)
res$honest_lb   <- honest_list(fa$hl)                    # M = RDSmoothnessBound lower bound
res$honest_mest <- honest_list(fa$he)                    # M = RDSmoothnessBound point estimate
res$smooth_bound <- list(est = fa$m_est, lb = fa$m_lb, s = 2)

# local randomization: window from covariate balance (tenure), then inference inside it
ws <- quiet(rdwinselect(B$x, B$tenure, cutoff = c0, wmasspoints = TRUE, reps = 1000, seed = LR_SEED, quietly = TRUE))
wl <- ws$w_left; wr <- ws$w_right
stopifnot(any(B$x == wl), any(B$x == wr), wl < c0, wr >= c0)   # the window ends are support points
grid <- seq(0, 0.3, by = 0.0025)
lr <- quiet(rdrandinf(B$y, B$x, cutoff = c0, wl = wl, wr = wr, ci = c(0.05, grid), reps = 1000, seed = LR_SEED,
                      quietly = TRUE))
lr_ci <- matrix(lr$ci, ncol = 2)
if (nrow(lr_ci) != 1L) stop("03: the rdrandinf confidence set is not one interval")
if (min(lr_ci) <= min(grid) || max(lr_ci) >= max(grid)) stop("03: the rdrandinf interval touches the tau grid edge")
inw <- B$x >= wl & B$x <= wr
hand <- lr_diffmeans(B$y, B$x, c0, w0 = (wr - wl) / 2, reps = 2000, seed = LR_SEED)
stopifnot(abs(hand$est - lr$obs.stat) < 1e-10)
res$locrand <- list(est = lr$obs.stat, ci = as.numeric(lr_ci[1, ]), p = lr$p.value, p_asy = lr$asy.pvalue,
                    window = c(wl, wr), n_left = sum(inw & B$x < c0), n_right = sum(inw & B$x >= c0),
                    p_hand = hand$p, reps = 1000, reps_hand = 2000, grid_step = 0.0025)
res$winselect <- list(window = c(wl, wr), w_right = ws$wlist_right + c0,
                      balance_p = as.numeric(ws$results[, "p-value"]))
# window sensitivity: the first three symmetric windows of rdwinselect ([4.2, 4.3], [4.1, 4.4],
# [4.0, 4.5]); p-values of the sharp nulls tau = 0 and tau = truth, and the difference in
# means in each window
wlist <- round(c0 + ws$wlist_right[1:3], 10)
stopifnot(isTRUE(all.equal(wlist[1], wr)))
sens <- quiet(rdsensitivity(B$y, B$x, cutoff = c0, wlist = wlist, tlist = c(0, tau), reps = 1000, seed = LR_SEED,
                            nodraw = TRUE, quietly = TRUE))
res$sensitivity <- list(w_left = sens$wlist_left, w_right = sens$wlist, tlist = sens$tlist,
                        p = sens$results,                 # rows: tlist, columns: windows
                        est = vapply(wlist, function(w) lr_diffmeans(B$y, B$x, c0, w0 = w - c0, reps = 10,
                                                                    seed = LR_SEED)$est, numeric(1)))

## ---- Monte Carlo -------------------------------------------------------------------------------------
t_mc <- proc.time()[["elapsed"]]
set.seed(MASTER_SEED)
seeds <- sample.int(1e8, REPS)
one_draw <- function(s) {
  D <- gen_R2(n = B$meta$N, seed = s)
  stopifnot(D$meta$jump == tau, D$meta$M_bound == M_bound)
  f <- fit_all(D$y, D$x, M_bound)
  w <- D$x >= wl & D$x <= wr; y1 <- D$y[w & D$x >= c0]; y0 <- D$y[w & D$x < c0]   # the frozen draw's window
  lr_est <- mean(y1) - mean(y0); lr_se <- sqrt(var(y1) / length(y1) + var(y0) / length(y0))
  c(seed = s, est_conv = f$rd$coef["Conventional", 1], h = f$h,
    conv_lo = f$rd$ci["Conventional", 1], conv_hi = f$rd$ci["Conventional", 2],
    rob_lo = f$rd$ci["Robust", 1], rob_hi = f$rd$ci["Robust", 2],
    crv_lo = f$cr$ci["Conventional", 1], crv_hi = f$cr$ci["Conventional", 2],
    est_hon = f$ho$estimate, hon_lo = f$ho$conf.low, hon_hi = f$ho$conf.high,
    m_lb = f$m_lb, hlb_lo = f$hl$conf.low, hlb_hi = f$hl$conf.high,
    m_est = f$m_est, hest_lo = f$he$conf.low, hest_hi = f$he$conf.high,
    est_lr = lr_est, lr_lo = lr_est - qnorm(0.975) * lr_se, lr_hi = lr_est + qnorm(0.975) * lr_se)
}
out <- parallel::mclapply(seeds, one_draw, mc.cores = CORES)
bad <- !vapply(out, is.numeric, logical(1))
if (any(bad)) stop("03: ", sum(bad), " Monte Carlo draws failed, first: ", as.character(out[[which(bad)[1]]]))
draws <- do.call(rbind, out)
if (nrow(draws) != REPS || anyNA(draws)) stop("03: Monte Carlo draws returned NA")
runtime_mc <- proc.time()[["elapsed"]] - t_mc
cov_len <- function(lo, hi) list(cover = mean(draws[, lo] <= tau & tau <= draws[, hi]),
                                 len = mean(draws[, hi] - draws[, lo]))
cl <- list(conv = cov_len("conv_lo", "conv_hi"), robust = cov_len("rob_lo", "rob_hi"), crv = cov_len("crv_lo", "crv_hi"),
           honest = cov_len("hon_lo", "hon_hi"), honest_lb = cov_len("hlb_lo", "hlb_hi"),
           honest_mest = cov_len("hest_lo", "hest_hi"), lr = cov_len("lr_lo", "lr_hi"))
mcse <- function(p) sqrt(p * (1 - p) / REPS)
res$mc <- list(reps = REPS, master_seed = MASTER_SEED, n = B$meta$N, cores = CORES, runtime_sec = runtime_mc,
               cover_conv = cl$conv$cover, cover_robust = cl$robust$cover, cover_crv = cl$crv$cover,
               cover_honest = cl$honest$cover, cover_honest_lb = cl$honest_lb$cover,
               cover_honest_mest = cl$honest_mest$cover, cover_lr = cl$lr$cover,
               len_conv = cl$conv$len, len_robust = cl$robust$len, len_crv = cl$crv$len, len_honest = cl$honest$len,
               len_honest_lb = cl$honest_lb$len, len_honest_mest = cl$honest_mest$len, len_lr = cl$lr$len,
               mcse_conv = mcse(cl$conv$cover), mcse_honest = mcse(cl$honest$cover),
               mean_est_conv = mean(draws[, "est_conv"]), sd_est_conv = sd(draws[, "est_conv"]),
               mean_est_honest = mean(draws[, "est_hon"]), mean_est_lr = mean(draws[, "est_lr"]),
               mean_h_conv = mean(draws[, "h"]), mean_m_lb = mean(draws[, "m_lb"]), median_m_lb = median(draws[, "m_lb"]),
               share_m_lb_zero = mean(draws[, "m_lb"] == 0), mean_m_est = mean(draws[, "m_est"]),
               share_m_est_ge_bound = mean(draws[, "m_est"] >= M_bound),
               lr_window = c(wl, wr), draws = draws)

saveRDS(res, file.path(here, "res_discrete.rds"))

## ---- report ---------------------------------------------------------------------------------------------
ci_txt <- function(ci) sprintf("[%.3f, %.3f]", ci[1], ci[2])
cat(sprintf("================ R2: discrete running variable (truth %.3f) ================\n", tau))
cat(sprintf("support points %d (min mass %d), %d right of %.2f; N %d\n", res$sample$n_support, res$sample$min_mass,
            res$sample$n_support_right, c0, res$sample$N))
cat(sprintf("conventional  %.3f %s  h %.3f (support points in h: %d left, %d right), eff. n %d / %d\n",
            res$conv$est, ci_txt(res$conv$ci), res$conv$h, res$conv$n_support_in_h[1], res$conv$n_support_in_h[2],
            res$conv$n_eff[1], res$conv$n_eff[2]))
cat(sprintf("robust        %.3f %s\n", res$robust$est, ci_txt(res$robust$ci)))
cat(sprintf("CRV by x      %.3f %s  (se %.4f vs conventional %.4f, %d clusters in h)\n", res$crv$est, ci_txt(res$crv$ci),
            res$crv$se, res$conv$se, res$crv$n_clusters))
for (nm in c("honest", "honest_lb", "honest_mest")) {
  hh <- res[[nm]]
  cat(sprintf("%-13s %.3f %s  M %.3f, h %.3f, max bias %.3f, se %.3f\n", nm, hh$est, ci_txt(hh$ci), hh$M, hh$h,
              hh$max_bias, hh$se))
}
cat(sprintf("RDSmoothnessBound (s = 2): estimate %.3f, lower bound %.3f (M_bound %.3f)\n", res$smooth_bound$est,
            res$smooth_bound$lb, M_bound))
cat(sprintf("rdwinselect window [%.1f, %.1f] (balance p %s)\n", wl, wr,
            paste(sprintf("%.3f", head(res$winselect$balance_p, 3)), collapse = ", ")))
cat(sprintf("local rand.   %.4f %s  randomization p %.4f (hand, valid version %.5f), n %d / %d\n", res$locrand$est,
            ci_txt(res$locrand$ci), res$locrand$p, res$locrand$p_hand, res$locrand$n_left, res$locrand$n_right))
cat(sprintf("window sensitivity (diff in means; p of tau = 0, tau = %.2f):\n", tau))
for (k in seq_along(wlist))
  cat(sprintf("  [%.1f, %.1f]  %.4f  p0 %.3f  p_truth %.3f\n", sens$wlist_left[k], wlist[k], res$sensitivity$est[k],
              sens$results[1, k], sens$results[2, k]))
cat(sprintf("Monte Carlo: %d draws, master seed %d, %d cores, %.0f s\n", REPS, MASTER_SEED, CORES, runtime_mc))
for (nm in names(cl)) cat(sprintf("  %-12s coverage %.3f  mean length %.3f\n", nm, cl[[nm]]$cover, cl[[nm]]$len))
cat(sprintf("  mean conventional estimate %.4f (sd %.4f); mean RDSmoothnessBound lower bound %.3f, estimate %.3f\n",
            res$mc$mean_est_conv, res$mc$sd_est_conv, res$mc$mean_m_lb, res$mc$mean_m_est))
