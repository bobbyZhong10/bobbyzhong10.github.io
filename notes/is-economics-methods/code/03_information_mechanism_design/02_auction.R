# ------------------------------------------------------------------
# 02_auction.R -- block V1: revenue equivalence, the optimal reserve and the GPV inversion (Vireo's ad exchange).
#
# 1. Revenues per auction from the bids in vireo_V1.rds under the four formats (first price: the highest bid;
#    second price: the second-highest truthful bid; both again with the reserve r*; unsold: 0), pooled over the
#    10,000 auctions and by bidder count n. The formats are run on the same value draws, so the revenue
#    difference is a paired mean: its standard error is sd(d) / sqrt(N) (the standard error of a mean).
# 2. GPV in quantile form on the first-price bids, separately for each n: v_hat(tau) = b(tau) + tau b'(tau) / (n-1),
#    b(tau) the empirical tau-quantile of the bids, b' a central difference of half-width h = 0.02, at
#    tau = 0.05, 0.06, ..., 0.95; compared with the true value quantile F^{-1}(tau).
# 3. Monte Carlo (seed 3102, Ruling P23): 1,000 replications of the whole 10,000-auction design (first-price
#    bids from fpa_bid_fun, the exact bid interpolated on a grid of 2,001 values); the mean and the 2.5% to
#    97.5% range of each revenue, of the differences, and of the GPV estimate at each (n, tau).
# Output: res_auction.rds. About 30 s.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("mech_fns.R")
tr <- readRDS("truth03.rds"); P1 <- tr$params$V1; seed <- tr$params$seed
V1 <- readRDS("vireo_V1.rds")
D1 <- value_dist(tr$V1$F); rs <- tr$V1$r_star
rng <- function(x) unname(stats::quantile(x, c(0.025, 0.975), type = 7))
se <- function(x) stats::sd(x) / sqrt(length(x))

## ---- 1. revenues by format ----------------------------------------------------------------
auc <- V1[, {
  o <- order(v_true, decreasing = TRUE); v1 <- v_true[o[1]]; v2 <- v_true[o[2]]
  list(n = n[1], v1 = v1, v2 = v2, rev_fp = max(b_fp), rev_sp = v2,
       rev_fp_res = if (v1 < rs) 0 else max(b_fp_res, na.rm = TRUE),
       rev_sp_res = if (v1 < rs) 0 else max(v2, rs))
}, by = auction]
auc[, `:=`(d = rev_fp - rev_sp, d_res = rev_fp_res - rev_sp_res, gain = rev_sp_res - rev_sp)]
shade <- V1[, .(shading = mean((v_true - b_fp) / v_true), bid_over_value = mean(b_fp / v_true)), by = n]
by_n <- auc[, .(N = .N, rev_fp = mean(rev_fp), rev_sp = mean(rev_sp), rev_fp_res = mean(rev_fp_res),
                rev_sp_res = mean(rev_sp_res), se_diff = se(d), se_diff_res = se(d_res),
                gain_reserve = mean(gain), se_gain = se(gain), share_unsold = mean(v1 < rs)), by = n]
by_n <- merge(by_n, shade, by = "n")
by_n[, `:=`(rev_true = tr$V1$rev_true$none[match(n, tr$V1$rev_true$n)],
            rev_true_res = tr$V1$rev_true$reserve[match(n, tr$V1$rev_true$n)],
            p_unsold_true = tr$V1$p_unsold[match(n, P1$n_set)])]
setorder(by_n, n)

## ---- 2. GPV, quantile form, by n ---------------------------------------------------------------
gpv <- rbindlist(lapply(P1$n_set, function(k) {
  g <- gpv_quantile(V1[n == k]$b_fp, k, P1$taus, h = P1$h)
  g[, `:=`(n = k, v_true = D1$quantile(tau))]
  g[, relerr := abs(v_hat - v_true) / v_true]
  g[, .(n, tau, b, v_hat, v_true, relerr)]
}))
gpv_by_n <- gpv[, .(med_relerr = stats::median(relerr), max_relerr = max(relerr), mean_relerr = mean(relerr)), by = n]

## ---- 3. Monte Carlo (seed 3102) ----------------------------------------------------------------
set.seed(seed$V1_mc)
R <- P1$R_mc; ns <- P1$n_set; N <- P1$n_auc_per_n; nt <- length(P1$taus)
bf <- lapply(ns, function(n) fpa_bid_fun(n, D1)); bfr <- lapply(ns, function(n) fpa_bid_fun(n, D1, reserve = rs))
fmt <- c("rev_fp", "rev_sp", "rev_fp_res", "rev_sp_res")
mc_n <- array(NA_real_, c(R, length(ns), length(fmt)), dimnames = list(NULL, ns, fmt))
mc_pool <- matrix(NA_real_, R, 8, dimnames = list(NULL, c(fmt, "se_diff", "se_diff_res", "se_gain", "med_relerr")))
mc_gpv <- array(NA_real_, c(R, length(ns), nt))
for (r in seq_len(R)) {
  out <- vector("list", length(ns)); rel <- vector("list", length(ns))
  for (k in seq_along(ns)) {
    n <- ns[k]
    V <- matrix(D1$quantile(stats::runif(N * n)), N, n, byrow = TRUE)      # one row per auction
    ao <- auction_outcomes(V, rs, bf[[k]], bfr[[k]])
    out[[k]] <- ao
    mc_n[r, k, ] <- colMeans(ao[, ..fmt])
    g <- gpv_quantile(bf[[k]](as.vector(V)), n, P1$taus, h = P1$h)
    mc_gpv[r, k, ] <- g$v_hat
    vt <- D1$quantile(P1$taus); rel[[k]] <- abs(g$v_hat - vt) / vt
  }
  ao <- rbindlist(out)
  mc_pool[r, ] <- c(colMeans(ao[, ..fmt]), se(ao$rev_fp - ao$rev_sp), se(ao$rev_fp_res - ao$rev_sp_res),
                    se(ao$rev_sp_res - ao$rev_sp), stats::median(unlist(rel)))
}
diff_r <- mc_pool[, "rev_fp"] - mc_pool[, "rev_sp"]
diff_res_r <- mc_pool[, "rev_fp_res"] - mc_pool[, "rev_sp_res"]
gain_r <- mc_pool[, "rev_sp_res"] - mc_pool[, "rev_sp"]
mc_by_n <- rbindlist(lapply(seq_along(ns), function(k) {
  x <- mc_n[, k, , drop = TRUE]
  data.table(n = ns[k], format = fmt, mean = colMeans(x),
             lo = apply(x, 2, function(z) rng(z)[1]), hi = apply(x, 2, function(z) rng(z)[2]))
}))
gain_by_n <- rbindlist(lapply(seq_along(ns), function(k) {
  g <- mc_n[, k, "rev_sp_res"] - mc_n[, k, "rev_sp"]; d <- mc_n[, k, "rev_fp"] - mc_n[, k, "rev_sp"]
  data.table(n = ns[k], gain_mean = mean(g), gain_lo = rng(g)[1], gain_hi = rng(g)[2],
             diff_mean = mean(d), diff_lo = rng(d)[1], diff_hi = rng(d)[2])
}))
gpv_band <- rbindlist(lapply(seq_along(ns), function(k) {
  x <- mc_gpv[, k, ]
  data.table(n = ns[k], tau = P1$taus, v_true = D1$quantile(P1$taus), mean = colMeans(x),
             lo = apply(x, 2, function(z) rng(z)[1]), hi = apply(x, 2, function(z) rng(z)[2]))
}))
mc <- list(R = R, diff_mean = mean(diff_r), diff_sd = stats::sd(diff_r), diff_range = rng(diff_r),
           diff_res_mean = mean(diff_res_r), diff_res_sd = stats::sd(diff_res_r), diff_res_range = rng(diff_res_r),
           gain_mean = mean(gain_r), gain_sd = stats::sd(gain_r), gain_range = rng(gain_r),
           rev_fp_mean = mean(mc_pool[, "rev_fp"]), rev_fp_range = rng(mc_pool[, "rev_fp"]),
           rev_sp_mean = mean(mc_pool[, "rev_sp"]), rev_sp_range = rng(mc_pool[, "rev_sp"]),
           rev_fp_res_mean = mean(mc_pool[, "rev_fp_res"]), rev_fp_res_range = rng(mc_pool[, "rev_fp_res"]),
           rev_sp_res_mean = mean(mc_pool[, "rev_sp_res"]), rev_sp_res_range = rng(mc_pool[, "rev_sp_res"]),
           se_diff_mean = mean(mc_pool[, "se_diff"]),
           share_within_3se = mean(abs(diff_r) < 3 * mc_pool[, "se_diff"]),
           share_within_3se_res = mean(abs(diff_res_r) < 3 * mc_pool[, "se_diff_res"]),
           share_gain_positive = mean(gain_r > 0),
           gpv_relerr_mean = mean(mc_pool[, "med_relerr"]), gpv_relerr_range = rng(mc_pool[, "med_relerr"]),
           by_n = mc_by_n, gain_by_n = gain_by_n, gpv_band = gpv_band)

## ---- save ----------------------------------------------------------------------------------------
a1 <- list(N_auc = nrow(auc), r_star = rs,
           rev_fp = mean(auc$rev_fp), rev_sp = mean(auc$rev_sp), rev_fp_reserve = mean(auc$rev_fp_res),
           rev_reserve = mean(auc$rev_sp_res),
           diff = mean(auc$d), se_diff = se(auc$d), diff_res = mean(auc$d_res), se_diff_res = se(auc$d_res),
           gain = mean(auc$gain), se_gain = se(auc$gain), se_rev_sp = se(auc$rev_sp), se_rev_fp = se(auc$rev_fp),
           share_unsold = mean(auc$v1 < rs), by_n = by_n,
           gpv = gpv, gpv_by_n = gpv_by_n, gpv_med_relerr = stats::median(gpv$relerr),
           gpv_med_relerr_by_n = gpv_by_n$med_relerr, gpv_max_relerr = max(gpv$relerr), gpv_mean_relerr = mean(gpv$relerr),
           mc = mc)
saveRDS(a1, "res_auction.rds")

cat(sprintf("revenue per auction (USD per click): first price %.4f, second price %.4f, difference %.4f (se %.4f)\n",
            a1$rev_fp, a1$rev_sp, a1$diff, a1$se_diff))
cat(sprintf("with r* = %.4f: first price %.4f, second price %.4f; reserve gain %.4f (se %.4f); unsold %.4f\n",
            rs, a1$rev_fp_reserve, a1$rev_reserve, a1$gain, a1$se_gain, a1$share_unsold))
for (i in seq_len(nrow(by_n))) with(by_n[i], cat(sprintf(
  "  n = %d: FP %.4f SP %.4f (truth %.4f) | FP+r %.4f SP+r %.4f (truth %.4f) | gain %.4f | shading %.4f\n",
  n, rev_fp, rev_sp, rev_true, rev_fp_res, rev_sp_res, rev_true_res, gain_reserve, shading)))
cat(sprintf("GPV: median relative error %.4f (by n: %s); max %.4f\n", a1$gpv_med_relerr,
            paste(sprintf("%.4f", a1$gpv_med_relerr_by_n), collapse = " "), a1$gpv_max_relerr))
cat(sprintf("Monte Carlo, %d replications: FP - SP mean %.5f, sd %.5f, range [%.4f, %.4f]; se formula %.5f; within 3 se %.3f\n",
            R, mc$diff_mean, mc$diff_sd, mc$diff_range[1], mc$diff_range[2], mc$se_diff_mean, mc$share_within_3se))
cat(sprintf("  reserve gain mean %.4f range [%.4f, %.4f]; GPV median relative error mean %.4f range [%.4f, %.4f]\n",
            mc$gain_mean, mc$gain_range[1], mc$gain_range[2], mc$gpv_relerr_mean, mc$gpv_relerr_range[1], mc$gpv_relerr_range[2]))
cat("Saved res_auction.rds\n")
