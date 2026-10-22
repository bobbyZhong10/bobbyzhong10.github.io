# ------------------------------------------------------------------
# 02_d1_more.R -- blocks D1 (dyn_main, rho = 0.70) and D1p (dyn_persist, rho = 0.90), the inherited
# samples: the estimators the legacy 03_dynamic_nickell.R does not run.
#   pooled OLS and FE with CR0 SEs; the Hahn-Kuersteiner correction and the half-panel jackknife of FE (Dhaene and
#   Jochmans 2015 eqs. 2.2 and 2.4: the halves weighted by |S|/T, so for T = 7 by 3/7 and 4/7);
#   Anderson-Hsiao with the level and the differenced instrument; one-step Arellano-Bond (robust SE);
#   two-step Arellano-Bond and Blundell-Bond with the conventional and the Windmeijer SE (the legacy
#   fits, refitted: the estimates equal res_dynamic.rds); instrument counts, Hansen J, m1, m2.
# All GMM fits as the legacy code: plm::pgmm, effect = "individual", every available lag lag(y, 2:99).
# Fix round 2 (11.3b I2): a Monte Carlo of the Anderson-Hsiao estimator with the differenced instrument in the D1p
# design (rho = 0.90): the legacy gen_dyn at N = 500, T_obs = 8, rho = 0.90, R_mc = 1,000 samples, seed D1p_ah_mc
# (truth11.rds params); it stores the coverage of the 95% interval, the median and the mean estimate
# (res_d1.rds, D1p$ah_diff_mc). The single D1p sample's interval for this estimator misses 0.90 (0.3116 [-0.2007, 0.8238]).
# Reads nimbus.rds, truth11.rds; writes res_d1.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("dyn_fns.R")
suppressPackageStartupMessages(library(plm))
nb <- readRDS("nimbus.rds")

one_block <- function(d) {
  pd <- pdata.frame(d, index = c("id", "t"))
  gm <- function(model, tf) pgmm(y ~ lag(y, 1) | lag(y, 2:99), pd, effect = "individual", model = model, transformation = tf)
  pool <- pooled_ar1(d$y, d$id, d$t)
  fe <- fe_ar1(d$y, d$id, d$t)
  TT <- fe$T                                                     # regression periods: 8 months - 1 = 7
  sp <- spj_parts(d$y, d$id, d$t)
  fe_l <- list(est = unname(fe$est["rho"]), se = unname(fe$se["rho"]))
  fe_l$ci <- ci95(fe_l$est, fe_l$se)
  hk <- hk_correct(fe_l$est, TT)
  se_hk <- (TT + 1) / TT * fe_l$se                               # hk_correct is linear in rho_hat_FE
  list(pool = list(est = pool$est, se = pool$se, ci = ci95(pool$est, pool$se)),
       fe = fe_l, T = TT, hk = hk, hk_se = se_hk, hk_ci = ci95(hk, se_hk),
       spj = sp$spj, spj_halves = sp$halves, spj_half_mean = sp$half_mean,
       ah = anderson_hsiao(d$y, d$id, d$t, "level"), ah_diff = anderson_hsiao(d$y, d$id, d$t, "diff"),
       ab_onestep = gmm_summary(gm("onestep", "d")), ab = gmm_summary(gm("twosteps", "d")),
       bb = gmm_summary(gm("twosteps", "ld")))
}
D1 <- one_block(nb$dyn_main)
D1p <- one_block(nb$dyn_persist)

## ---- Monte Carlo: Anderson-Hsiao, differenced instrument, D1p design (rho = 0.90) ------------------------------
## Exactly the D1p design (the legacy gen_dyn stored in nimbus.rds, arguments from params$D1p_ah) and exactly one
## estimator, anderson_hsiao(inst = "diff") (the instrument Delta Y_{i,t-2}); set.seed once, then R_mc draws.
## cover = share of samples whose 95% interval est -/+ qnorm(0.975) SE (CR0 by merchant) contains rho = 0.90.
tr <- readRDS("truth11.rds"); Pa <- tr$params$D1p_ah
set.seed(tr$params$seed$D1p_ah_mc)
ah_draws <- t(vapply(seq_len(Pa$R_mc), function(r) {
  dd <- nb$gen_dyn(N = Pa$N, T = Pa$T_obs, rho = Pa$rho, burn = Pa$burn, sigma_a = Pa$sigma_a, sigma_e = Pa$sigma_e)
  a <- anderson_hsiao(dd$y, dd$id, dd$t, "diff")
  c(est = a$est, se = a$se)
}, numeric(2)))
## A just-identified IV estimator has no finite moments (the ratio of two sample moments), so the Monte Carlo mean and
## SD of the estimates are dominated by a few extreme draws and are not estimates of anything: they are stored because
## the text may be asked for them, next to the median, the interquartile range and the central 95% range of the
## estimates, which are what a weak instrument leaves well defined.
ah_ci <- cbind(lo = ah_draws[, "est"] - qnorm(0.975) * ah_draws[, "se"], hi = ah_draws[, "est"] + qnorm(0.975) * ah_draws[, "se"])
ah_cover <- mean(ah_ci[, "lo"] <= Pa$rho & Pa$rho <= ah_ci[, "hi"])
ah_q <- unname(quantile(ah_draws[, "est"], c(0.025, 0.975)))
D1p$ah_diff_mc <- list(R = Pa$R_mc, cover = ah_cover, cover_mcse = sqrt(ah_cover * (1 - ah_cover) / Pa$R_mc),
                       share_ci_below_truth = mean(ah_ci[, "hi"] < Pa$rho), share_ci_above_truth = mean(ah_ci[, "lo"] > Pa$rho),
                       median = median(ah_draws[, "est"]), median_bias = median(ah_draws[, "est"]) - Pa$rho,
                       mean = mean(ah_draws[, "est"]), iqr = IQR(ah_draws[, "est"]), q025 = ah_q[1], q975 = ah_q[2],
                       draws = ah_draws)

res <- list(ah = D1$ah, ah_diff = D1$ah_diff, ab_onestep = D1$ab_onestep, ab_twostep = D1$ab, bb_twostep = D1$bb,
            D1 = D1, D1p = D1p)
saveRDS(res, "res_d1.rds")

pr <- function(lab, s, se = s$se) cat(sprintf("  %-34s %.4f (SE %.4f)\n", lab, s$est, se))
for (b in c("D1", "D1p")) {
  x <- res[[b]]
  cat(sprintf("=== %s (truth rho %.2f; T = %d regression periods) ===\n", b, if (b == "D1") 0.70 else 0.90, x$T))
  pr("pooled OLS", x$pool); pr("FE", x$fe)
  cat(sprintf("  %-34s %.4f (SE %.4f)\n", "FE + Hahn-Kuersteiner", x$hk, x$hk_se))
  cat(sprintf("  %-34s %.4f (halves %s)\n", "half-panel jackknife", x$spj,
              paste(sprintf("%.4f", x$spj_halves), collapse = " ")))
  pr("Anderson-Hsiao, level instrument", x$ah); pr("Anderson-Hsiao, diff instrument", x$ah_diff)
  pr("AB one-step (robust SE)", x$ab_onestep)
  cat(sprintf("  %-34s %.4f (SE %.4f Windmeijer, %.4f conventional); r = %d, Hansen p %.3f, m2 p %.3f\n",
              "AB two-step", x$ab$est, x$ab$se, x$ab$se_conv, x$ab$r, x$ab$J_p, x$ab$m2_p))
  cat(sprintf("  %-34s %.4f (SE %.4f Windmeijer, %.4f conventional); r = %d, Hansen p %.3f\n",
              "BB two-step", x$bb$est, x$bb$se, x$bb$se_conv, x$bb$r, x$bb$J_p))
}
m <- D1p$ah_diff_mc
cat(sprintf("=== D1p Anderson-Hsiao, differenced instrument, Monte Carlo (%d samples of %d merchants, rho = %.2f) ===\n",
            m$R, Pa$N, Pa$rho))
cat(sprintf("  coverage of the 95%% interval %.3f (Monte Carlo SE %.3f; intervals entirely below the truth %.3f, above %.3f)\n",
            m$cover, m$cover_mcse, m$share_ci_below_truth, m$share_ci_above_truth))
cat(sprintf("  median estimate %.4f, mean %.4f (no finite mean: dominated by extreme draws), interquartile range %.3f, central 95%% range [%.3f, %.3f]\n",
            m$median, m$mean, m$iqr, m$q025, m$q975))
cat("Saved res_d1.rds\n")
