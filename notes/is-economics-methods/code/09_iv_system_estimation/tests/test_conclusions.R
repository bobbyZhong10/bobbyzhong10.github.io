# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 9.
# Run from the chapter code directory after 01 to 05 (tests/run_tests.R does this;
# it has sourced chapter_kit.R, so covers() is available).
# The first block is the chapter card, verbatim (n_sellers is the name of the
# coefficient on the number of sellers in the buyer equation). The second block
# checks the controller's ruling 14 and the Monte Carlo behind the single-sample
# statements.
# ------------------------------------------------------------------

## ---- chapter card (verbatim) ------------------------------------------------------
tr <- readRDS("truth09.rds"); s1 <- readRDS("res_system.rds"); s2 <- readRDS("res_sur.rds")
s3 <- readRDS("res_weak.rds"); s4 <- readRDS("res_overid.rds"); b <- tr$S1$buyer[["n_sellers"]]
stopifnot(max(abs(s1$tsls_hand$b - s1$tsls_pkg$b)) < 1e-10, max(abs(s1$three_sls_hand$b - s1$three_sls_pkg$b)) < 1e-6)
stopifnot(!covers(s1$ols$buyer_ci, b), covers(s1$tsls$buyer_ci, b), s1$se$three_sls <= s1$se$tsls)
stopifnot(abs(s1$misspec$three_sls_bias) > abs(s1$misspec$tsls_bias), covers(s1$misspec$tsls_ci, b))   # 一条方程错，3SLS 传染
stopifnot(max(abs(s2$same_X$sur - s2$same_X$ols)) < 1e-10, all(s2$diff_X$se_sur <= s2$diff_X$se_ols))  # Theorem 7.6 与效率增益
stopifnot(s3$cover$wald <= 0.85, s3$cover$ar >= 0.93)                                                  # MC 不少于 1000 次
stopifnot(abs(s4$valid$mc_reject_J - 0.05) < 0.02, s4$invalid$J_p < 0.01)

## ---- ruling 14 and the Monte Carlo behind the card's single-sample statements -------
# S1 is over-identified in the seller equation (two excluded buyer-side shifters), so 3SLS and
# 2SLS differ on the just-identified buyer equation; on the over-identified seller equation they
# coincide exactly (Wooldridge 2010 section 9.3.2). The gain on n_sellers is not a rounding artefact.
stopifnot(abs(s1$diff$buyer_n_sellers) > 1e-3, max(abs(s1$diff$seller)) < 1e-8,
          s1$se$three_sls / s1$se$tsls < 0.95, tr$S1$pop_avar$se_ratio < 0.95)
# without the coupon-fee link in the city plan (independent shifters) 3SLS would not help n_sellers at all
stopifnot(abs(tr$S1$pop_avar$se_ratio_independent - 1) < 1e-8)
# Monte Carlo (R >= 1000 per design): 3SLS has the smaller sampling spread when both equations are right;
# when the seller equation omits coupon x buyers, 3SLS is biased for n_sellers and 2SLS is not
mS1 <- s1$mc$S1; mS1m <- s1$mc$S1m
stopifnot(mS1$R >= 1000, mS1m$R >= 1000, mS1$sd_three_sls < mS1$sd_tsls,
          abs(mS1$mean_tsls - b) < 0.02, abs(mS1$mean_three_sls - b) < 0.02,
          abs(mS1m$mean_tsls - b) < 0.02, abs(mS1m$mean_three_sls - b) > 0.2,
          mS1m$cover_tsls >= 0.93, mS1m$cover_three_sls < 0.5)
# the large-sample limits agree: 3SLS converges to the wrong value in S1m, 2SLS does not
stopifnot(abs(tr$S1m$plim$tsls_buyer - b) < 0.01, abs(tr$S1m$plim$three_sls_buyer - b) > 0.2)
# the contagion depends on the coupon-fee link (02b_contagion_r0.R, one million markets per design, own seeds):
# with corr(coupon, fee) = -0.7 3SLS is off on n_sellers; with independent shifters the n_sellers limits agree
# and the contamination moves to the buyer equation's coupon coefficient; 2SLS is consistent for both coefficients
cr <- readRDS("res_contagion_r0.rds"); bc <- cr$truth[["coupon"]]
stopifnot(cr$N >= 1e6, !any(unlist(cr$seed) %in% unlist(tr$params$seed)), abs(cr$r0$corr_coupon_fee) < 0.01,
          abs(cr$r0$n_sellers_tsls - b) < 0.02, abs(cr$r0$n_sellers_three_sls - b) < 0.02,
          abs(cr$r0$coupon_tsls - bc) < 0.02, cr$r0$coupon_three_sls - bc > 0.2,
          abs(cr$r_minus0p7$n_sellers_tsls - b) < 0.02, b - cr$r_minus0p7$n_sellers_three_sls > 0.2,
          abs(cr$r_minus0p7$coupon_tsls - bc) < 0.02, cr$r_minus0p7$coupon_three_sls - bc > 0.2)
# the 2SLS SE on the n denominator (like-for-like with the 3SLS SE): smaller than the n - k one by sqrt((n - k) / n)
stopifnot(abs(s1$se$tsls_n_denominator / s1$se$tsls - sqrt((s1$n - 5) / s1$n)) < 1e-12,
          s1$se$ratio_n_denominator < 1, s1$se$ratio_n_denominator > s1$se$ratio)
# the misspecified seller equation is detected by its own over-identification test (power about 1;
# single-sample p-values under a true null are reported, not tested)
stopifnot(s1$misspec$J_seller_p < 0.01)
# weak IV: R >= 1000 at the main design and at every grid point; Wald coverage rises to nominal only
# with a strong first stage, AR coverage is nominal at every strength (grid MC standard error 0.007)
stopifnot(s3$cover$R >= 1000, all(s3$grid$R >= 1000),
          max(s3$grid$cover_wald[s3$grid$mu2 <= 1]) < 0.86, s3$grid$cover_wald[which.max(s3$grid$mu2)] > 0.93,
          all(s3$grid$cover_ar >= 0.925))
# overidentification: R >= 1000 valid and 1000 invalid samples; under heteroskedasticity tied to an
# instrument the homoskedastic Sargan statistic rejects a true null more often than the robust Hansen J;
# J detects the invalid referral bonus in (almost) every invalid sample, and 2SLS with it is biased
stopifnot(s4$valid$R >= 1000, s4$invalid$R_power >= 1000, s4$valid$mc_reject_sargan > s4$valid$mc_reject_J,
          s4$invalid$mc_reject_J > 0.95, s4$invalid$mc_bias > 0.1, abs(s4$valid$mc_mean_est - tr$S4$beta[["n_sellers"]]) < 0.02)
