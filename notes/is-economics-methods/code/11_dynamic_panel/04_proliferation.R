# ------------------------------------------------------------------
# 04_proliferation.R -- block D2: instrument proliferation (70 key accounts, 16 months, rho = 0.70).
# System GMM (Blundell-Bond: transformation = "ld", month effects: effect = "twoways"), one-step with the
# robust SE, and the Hansen J of the one-step fit (sargan(fit, "twosteps"), as summary.pgmm reports):
#   full        every available lag, lag(y, 2:99): 134 instruments for 70 merchants (r > N)
#   collapsed   the same lags, collapse = TRUE (Roodman 2009): 30 instruments
#   path        lag depth 2:L for L in params$D2$lags, uncollapsed and collapsed (fig_11_proliferation)
# Where proliferation pulls GMM (controller Ruling 41, which corrects Ruling 15; Roodman 2009): with the
# standard one-step weight H = D D', GLS on the differenced equation is exactly within-groups, so DIFFERENCE
# GMM with many instruments is pulled towards FE (0.5494 = FE on this sample, review check D), while SYSTEM GMM
# is pulled mostly towards OLS in levels (a mix with about 0.13 weight on within-groups, check E). D2 is
# system GMM, so its non-instrumented counterpart ols_transformed is OLS on the level equation: pooled OLS
# of Y_it on (1, Y_{i,t-1}, month dummies); FE and OLS in first differences are stored too.
# One-step because with r > N the two-step weight matrix (the inverse of sum_i Z_i' e_i e_i' Z_i, rank <= N)
# is singular: plm then uses a generalised inverse (warning suppressed); for the same reason the Hansen J of
# the full set is mechanically N whatever the data (chapter 7's many-moment problem): J = N exactly when
# r >= N (see below), so the p-value P(chi2_df > N) depends on the instrument count alone.
# Fix round 1: D2's Hansen p > 0.99 is therefore no evidence that the moments are valid. The evidence that
# proliferation destroys the power of the J test is the misspecified MA(1) twin of D2 (gen_D2_ma, seed
# D2_mis, 200 samples): with Y_{i,t-2} an invalid instrument the full set never rejects (J = N), the middle
# set (lags 2:4, r = 68 < N) rejects only sometimes, the collapsed set (r = 30) rejects almost always.
# Monte Carlo: R_mc = 500 samples of the same design (seed D2_mc) for the averages the text quotes, with the
# Hansen J of the full, the middle (lags 2:mid_lag) and the collapsed set; then the 200-sample MA(1) twin.
# Phase 4 fix (audit ch11-1): two J statistics per set in both Monte Carlos. J (keys J_*, rej_*) is the J of the
# one-step fit, the moments and the weight both from the one-step residuals (chapter Theorem 6; exactly N when
# r >= N). J2 (keys J2_*, rej2_*) is the Hansen J of chapter Theorem 4 (iv) (Roodman 2009a, Stata Journal,
# section 2.5): sargan() of the two-step fit, the same weight evaluated at the two-step estimate, which
# minimises that quadratic form, so J2 <= J in every sample and J2 <= N when r >= N. One-step system GMM does
# not use the efficient weight, so its J is not chi2(r - k) and over-rejects (collapsed set: share of p < 0.05
# under the correct model 0.266 with J, 0.056 with J2); the chapter reads size and power from J2. The two-step
# fits draw no random numbers, so every one-step number is unchanged.
# Reads nimbus_D2.rds, truth11.rds; writes res_d2.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("dyn_fns.R")
suppressPackageStartupMessages(library(plm))
tr <- readRDS("truth11.rds"); P <- tr$params$D2
d <- readRDS("nimbus_D2.rds")

sys_fit <- function(pd, L, coll, model = "onestep") {
  f <- stats::as.formula(sprintf("y ~ lag(y, 1) | lag(y, 2:%d)", L))
  suppressWarnings(pgmm(f, pd, effect = "twoways", model = model, transformation = "ld", collapse = coll))
}
## J statistic and p-value of a fit: of a one-step fit, the J of Theorem 6; of a two-step fit, the Hansen J2
J_stat <- function(f) { s <- suppressWarnings(sargan(f)); c(J = unname(s$statistic), p = unname(s$p.value)) }
pd <- pdata.frame(d, index = c("id", "t"))

## ---- the lag-depth path (the last row of each kind is the full or the collapsed fit) ----------------
path <- do.call(rbind, lapply(c(FALSE, TRUE), function(coll) do.call(rbind, lapply(P$lags, function(L) {
  s <- suppressWarnings(gmm_summary(sys_fit(pd, L, coll)))
  data.frame(collapse = coll, max_lag = L, r = s$r, est = s$est, se = s$se, lo = s$ci[1], hi = s$ci[2],
             J = s$J, J_df = s$J_df, J_p = s$J_p)
}))))
full <- suppressWarnings(gmm_summary(sys_fit(pd, 99, FALSE)))
coll <- suppressWarnings(gmm_summary(sys_fit(pd, 99, TRUE)))
ols_lev <- pooled_ar1(d$y, d$id, d$t, time_fe = TRUE)
ols_plain <- pooled_ar1(d$y, d$id, d$t, time_fe = FALSE)
fe <- fe_ar1(d$y, d$id, d$t)
ols_d <- ols_diff_ar1(d$y, d$id, d$t)

## ---- Monte Carlo of the same design --------------------------------------------------------------------
set.seed(tr$params$seed$D2_mc)
mc_m <- t(vapply(seq_len(P$R_mc), function(r) {
  dr <- gen_D2(P); pr <- pdata.frame(dr, index = c("id", "t"))
  gf <- sys_fit(pr, 99, FALSE); gc <- sys_fit(pr, 99, TRUE)
  gm <- sys_fit(pr, P$mid_lag, FALSE)
  j1 <- c(full = J_stat(gf), mid = J_stat(gm), collapsed = J_stat(gc))          # J of the one-step fits
  j2 <- c(full = J_stat(sys_fit(pr, 99, FALSE, "twosteps")), mid = J_stat(sys_fit(pr, P$mid_lag, FALSE, "twosteps")),
          collapsed = J_stat(sys_fit(pr, 99, TRUE, "twosteps")))                # Hansen J2 of the two-step fits
  c(full = coef(gf)[[1]], collapsed = coef(gc)[[1]], ols = pooled_ar1(dr$y, dr$id, dr$t, time_fe = TRUE)$est,
    fe = unname(fe_ar1(dr$y, dr$id, dr$t)$est["rho"]),
    J_full = j1[["full.J"]], J_mid = j1[["mid.J"]], J_collapsed = j1[["collapsed.J"]],
    J_p_full = j1[["full.p"]], J_p_mid = j1[["mid.p"]], J_p_collapsed = j1[["collapsed.p"]],
    J2_full = j2[["full.J"]], J2_mid = j2[["mid.J"]], J2_collapsed = j2[["collapsed.J"]],
    J2_p_full = j2[["full.p"]], J2_p_mid = j2[["mid.p"]], J2_p_collapsed = j2[["collapsed.p"]])
}, numeric(16)))
closer <- abs(mc_m[, "full"] - mc_m[, "ols"]) < abs(mc_m[, "collapsed"] - mc_m[, "ols"])
J1c <- c("J_full", "J_mid", "J_collapsed"); J2c <- c("J2_full", "J2_mid", "J2_collapsed")
mc <- list(R = P$R_mc, mean_full = mean(mc_m[, "full"]), mean_collapsed = mean(mc_m[, "collapsed"]),
           mean_ols = mean(mc_m[, "ols"]), mean_fe = mean(mc_m[, "fe"]),
           sd_full = sd(mc_m[, "full"]), sd_collapsed = sd(mc_m[, "collapsed"]),
           share_J_full_above_099 = mean(mc_m[, "J_p_full"] > 0.99), mean_J_p_full = mean(mc_m[, "J_p_full"]),
           mean_J_p_collapsed = mean(mc_m[, "J_p_collapsed"]), share_J_collapsed_below_005 = mean(mc_m[, "J_p_collapsed"] < 0.05),
           share_J_full_below_005 = mean(mc_m[, "J_p_full"] < 0.05), share_J_mid_below_005 = mean(mc_m[, "J_p_mid"] < 0.05),
           mean_J_p_mid = mean(mc_m[, "J_p_mid"]), J_full_maxdev_from_N = max(abs(mc_m[, "J_full"] - P$N)),
           share_J2_full_above_099 = mean(mc_m[, "J2_p_full"] > 0.99),
           share_J2_full_below_005 = mean(mc_m[, "J2_p_full"] < 0.05), share_J2_mid_below_005 = mean(mc_m[, "J2_p_mid"] < 0.05),
           share_J2_collapsed_below_005 = mean(mc_m[, "J2_p_collapsed"] < 0.05),
           mean_J2_p_full = mean(mc_m[, "J2_p_full"]), mean_J2_p_mid = mean(mc_m[, "J2_p_mid"]),
           mean_J2_p_collapsed = mean(mc_m[, "J2_p_collapsed"]),
           J2_full_max = max(mc_m[, "J2_full"]), J2_minus_J_max = max(mc_m[, J2c] - mc_m[, J1c]),
           share_closer = mean(closer), draws = mc_m)

## ---- misspecified twin: MA(1) errors, Y_{i,t-2} invalid; Hansen J power of the three instrument sets ------
## R_mc = 200 samples of gen_D2_ma (seed D2_mis). Rejection rate = share of samples with Hansen p < 0.05.
## J = N exactly when r >= N (the weight matrix has rank <= N, so J is the squared length of the N-vector of
## ones projected on the span of the unit moments), hence the full set (r = 134 > N = 70) cannot reject.
Pm <- tr$params$D2_ma
mis_fit <- function(pr, L, coll) {
  f <- sys_fit(pr, L, coll); s <- J_stat(f); s2 <- J_stat(sys_fit(pr, L, coll, "twosteps"))
  c(est = coef(f)[[1]], J = s[["J"]], p = s[["p"]], r = n_instruments(f), J2 = s2[["J"]], p2 = s2[["p"]])
}
set.seed(tr$params$seed$D2_mis)
mis_m <- t(vapply(seq_len(Pm$R_mc), function(r) {
  pr <- pdata.frame(gen_D2_ma(Pm), index = c("id", "t"))
  c(full = mis_fit(pr, 99, FALSE), mid = mis_fit(pr, Pm$mid_lag, FALSE), collapsed = mis_fit(pr, 99, TRUE))
}, numeric(18)))
mis <- list(R = Pm$R_mc, ma_coef = Pm$ma_coef, mid_lag = Pm$mid_lag,
            r_full = mis_m[1, "full.r"], r_mid = mis_m[1, "mid.r"], r_collapsed = mis_m[1, "collapsed.r"],
            rej_full = mean(mis_m[, "full.p"] < 0.05), rej_mid = mean(mis_m[, "mid.p"] < 0.05),
            rej_collapsed = mean(mis_m[, "collapsed.p"] < 0.05),
            rej2_full = mean(mis_m[, "full.p2"] < 0.05), rej2_mid = mean(mis_m[, "mid.p2"] < 0.05),
            rej2_collapsed = mean(mis_m[, "collapsed.p2"] < 0.05),
            mean_J2_full = mean(mis_m[, "full.J2"]), mean_J2_mid = mean(mis_m[, "mid.J2"]),
            mean_J2_collapsed = mean(mis_m[, "collapsed.J2"]),
            J2_minus_J_max = max(mis_m[, paste0(c("full", "mid", "collapsed"), ".J2")] -
                                 mis_m[, paste0(c("full", "mid", "collapsed"), ".J")]),
            mean_J_full = mean(mis_m[, "full.J"]), mean_J_mid = mean(mis_m[, "mid.J"]), mean_J_collapsed = mean(mis_m[, "collapsed.J"]),
            mean_p_full = mean(mis_m[, "full.p"]), mean_p_mid = mean(mis_m[, "mid.p"]), mean_p_collapsed = mean(mis_m[, "collapsed.p"]),
            mean_est_full = mean(mis_m[, "full.est"]), mean_est_mid = mean(mis_m[, "mid.est"]),
            mean_est_collapsed = mean(mis_m[, "collapsed.est"]),
            J_full_maxdev_from_N = max(abs(mis_m[, "full.J"] - Pm$N)), draws = mis_m)

res <- list(J_p_full = full$J_p, J_p_collapsed = coll$J_p, J_full = full$J, J_df_full = full$J_df,
            J_collapsed = coll$J, J_df_collapsed = coll$J_df,
            est_full = full$est, se_full = full$se, ci_full = full$ci,
            est_collapsed = coll$est, se_collapsed = coll$se, ci_collapsed = coll$ci,
            m2_p_full = full$m2_p, m2_p_collapsed = coll$m2_p,
            ols_transformed = ols_lev$est, ols_transformed_se = ols_lev$se, ols_levels_no_month = ols_plain$est,
            fe = unname(fe$est["rho"]), fe_se = unname(fe$se["rho"]), ols_diff = ols_d$est,
            n_inst_full = full$r, n_inst_collapsed = coll$r, N = P$N, T = P$T, path = path, mc = mc, mis = mis)
saveRDS(res, "res_d2.rds")

cat(sprintf("=== D2: %d key accounts, %d months, truth rho %.2f ===\n", P$N, P$T, P$rho))
cat(sprintf("system GMM, every lag:  %.4f (SE %.4f), r = %d, Hansen J = %.2f on %d df, p = %.4f\n",
            full$est, full$se, full$r, full$J, full$J_df, full$J_p))
cat(sprintf("system GMM, collapsed:  %.4f (SE %.4f), r = %d, Hansen J = %.2f on %d df, p = %.4f\n",
            coll$est, coll$se, coll$r, coll$J, coll$J_df, coll$J_p))
cat(sprintf("OLS in levels (month dummies) %.4f, without %.4f; FE %.4f; OLS in differences %.4f\n",
            ols_lev$est, ols_plain$est, res$fe, ols_d$est))
print(path, digits = 4, row.names = FALSE)
cat(sprintf("Monte Carlo (R = %d): mean full %.4f, collapsed %.4f, OLS levels %.4f, FE %.4f; Hansen p > 0.99 in %.1f%% (full);\n",
            mc$R, mc$mean_full, mc$mean_collapsed, mc$mean_ols, mc$mean_fe, 100 * mc$share_J_full_above_099))
cat(sprintf("  full closer to OLS than collapsed in %.1f%% of samples; Hansen p < 0.05 (size, correct model): full %.1f%%, lags 2:%d %.1f%%, collapsed %.1f%%\n",
            100 * mc$share_closer, 100 * mc$share_J_full_below_005, P$mid_lag, 100 * mc$share_J_mid_below_005,
            100 * mc$share_J_collapsed_below_005))
cat(sprintf("  max |J - N| of the full set over the %d samples: %.2e\n", mc$R, mc$J_full_maxdev_from_N))
cat(sprintf("  two-step Hansen J2, p < 0.05 (size): full %.1f%%, lags 2:%d %.1f%%, collapsed %.1f%%; full p > 0.99 in %.1f%%; max J2 of the full set %.2f (N = %d); max (J2 - J) %.1e\n",
            100 * mc$share_J2_full_below_005, P$mid_lag, 100 * mc$share_J2_mid_below_005, 100 * mc$share_J2_collapsed_below_005,
            100 * mc$share_J2_full_above_099, mc$J2_full_max, P$N, mc$J2_minus_J_max))
cat(sprintf("=== D2 misspecified twin (MA(1) errors, coefficient %.1f; %d samples): Hansen J rejection rate at 5%% ===\n",
            mis$ma_coef, mis$R))
cat(sprintf("  every lag (r = %d): %.3f (mean J %.2f, N = %d; max |J - N| %.1e); lags 2:%d (r = %d): %.3f (mean J %.2f); collapsed (r = %d): %.3f (mean J %.2f)\n",
            mis$r_full, mis$rej_full, mis$mean_J_full, Pm$N, mis$J_full_maxdev_from_N, mis$mid_lag, mis$r_mid, mis$rej_mid,
            mis$mean_J_mid, mis$r_collapsed, mis$rej_collapsed, mis$mean_J_collapsed))
cat(sprintf("  two-step Hansen J2: every lag %.3f (mean J2 %.2f); lags 2:%d %.3f (mean J2 %.2f); collapsed %.3f (mean J2 %.2f); max (J2 - J) %.1e\n",
            mis$rej2_full, mis$mean_J2_full, mis$mid_lag, mis$rej2_mid, mis$mean_J2_mid, mis$rej2_collapsed, mis$mean_J2_collapsed,
            mis$J2_minus_J_max))
cat("Saved res_d2.rds\n")
