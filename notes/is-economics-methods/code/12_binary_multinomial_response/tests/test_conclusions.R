# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 12 (chapter card, verbatim).
# Run from the chapter code directory after 01 to 08 (tests/run_tests.R does this; covers() comes
# from tools/R/chapter_kit.R, sourced once by the runner).
# ------------------------------------------------------------------
tr <- readRDS("truth12.rds"); x <- readRDS("res_extra.rds")
stopifnot(sign(x$B2$coef_interaction) != sign(x$B2$avg_cross), covers(x$B2$avg_cross_ci, tr$B2$avg_cross))  # Ai-Norton
stopifnot(!covers(x$B3$naive$ape_ci, tr$B3$ape), covers(x$B3$cf$ape_ci, tr$B3$ape))
stopifnot(abs(x$B4$beta) < 0.9 * abs(tr$B4$beta), covers(x$B4$ape_ci, tr$B4$ape))                        # 系数衰减，APE 不变
stopifnot(covers(x$B5$clogit$ci, tr$B5$beta), !covers(x$B5$pooled$ci, tr$B5$beta), abs(x$B5$bc$bias_T4) < abs(x$B5$fe$bias_T4))
# Ruling 29: in the one B6 sample the interval misses its truth (an ordinary sampling miss: beta-hat for
# engagement is 2.6 standard errors above 0.5; the miss is kept, numbers.json records it as
# mult_B6_ape_top_covers_truth = 0). The card's single-sample assertion
#   stopifnot(covers(x$B6$ape_top_ci, tr$B6$ape_top))
# is replaced by its Monte Carlo counterpart (R >= 1000 samples of the B6 design, seeds 1260000 + r): the mean
# estimate is within 0.001 of the target (1% of it) and the nominal 95% interval covers it in >= 93% of samples.
stopifnot(x$B6$mc$R >= 1000, abs(x$B6$mc$mean_est - tr$B6$ape_top) < 0.001, x$B6$mc$cover >= 0.93)
# B7 IIA: the classical Hausman-McFadden test over-rejects at 5% (chapter 12 section 7: 7.0% in 4,000 replications)
stopifnot(x$B7$iia$mc_reject > 0.05 + 2 * sqrt(0.05 * 0.95 / x$B7$iia$mc_R), x$B7$nested$hm_p < 0.01)

# ---- Fix round 1 (review of task 12.1) ----------------------------------------------------------------
# B5, APE of the FE logit: alpaca::getAPEs()$delta already averages over all firm-months, so the all-firm APE of
# the FE MLE and of its bias correction are close to the true APE (the first version scaled delta by the share of
# rows kept and understated it by about a third). Monte Carlo relative biases, both panel lengths.
stopifnot(abs(x$B5$mc$T4$ape_fe_rel_bias) < 0.03, abs(x$B5$mc$T8$ape_fe_rel_bias) < 0.03,
          abs(x$B5$mc$T4$ape_bc_rel_bias) < 0.05, abs(x$B5$mc$T8$ape_bc_rel_bias) < 0.05)
# B5, two designs at T = 4: the one-sample T = 4 panel (months 1 to 4 of the T = 8 panel) is compared with the truth
# of that design (T8first4); the 4-month DGP truths (T4) belong to the Monte Carlo
stopifnot(abs(x$B5$T4$truth$ape - tr$B5$T8first4$ape) < 1e-15, abs(x$B5$T4$truth$pooled_plim - tr$B5$T8first4$pooled_plim) < 1e-15,
          abs(x$B5$T8$truth$ape - tr$B5$ape) < 1e-15, abs(x$B5$mc$T4$ape_true - tr$B5$ape_T4) < 1e-15)
# B7, Hausman and McFadden: the stacked-score (suest-type) statistic is positive semidefinite by construction (no
# negative statistics) and its size under IIA is within the card's tolerance on the same Monte Carlo draws
stopifnot(x$B7$iia$suest$mc_R == x$B7$iia$mc_R, x$B7$iia$suest$mc_neg_share == 0,
          abs(x$B7$iia$suest$mc_reject - 0.05) < 0.03, x$B7$iia$mc_nonpsd_share > 0)

# ---- Fix round 2 (controller ruling 44) ----------------------------------------------------------------------
# B5: with the exact quadrature truths the all-firm APE biases are (reviewer's re-scoring): FE -0.70% (T = 4), -0.22% (T = 8);
# bias-corrected +2.70%, +0.57%
stopifnot(abs(x$B5$mc$T4$ape_fe_rel_bias + 0.0070) < 5e-4, abs(x$B5$mc$T8$ape_fe_rel_bias + 0.0022) < 5e-4,
          abs(x$B5$mc$T4$ape_bc_rel_bias - 0.0270) < 5e-4, abs(x$B5$mc$T8$ape_bc_rel_bias - 0.0057) < 5e-4)
# B5: the error of the all-firm FE APE is an offset of two large pieces: the kept firms' overstated effects (positive) and the
# dropped firms' zero estimated effects (negative); with the sampling difference of the per-panel truth they add up to the bias
for (m in x$B5$mc) {
  stopifnot(m$ape_fe_kept_part > 0.01, m$ape_fe_drop_part < -0.01,
            abs(m$ape_fe_bias - (m$ape_fe_kept_part + m$ape_fe_drop_part + m$ape_fe_sampling_part)) < 1e-12,
            abs(m$ape_fe_kept_part + m$ape_fe_drop_part) < 0.2 * m$ape_fe_kept_part)             # the net error is a small fraction of each piece
}
# the movers-only bias-corrected APE is not a consistent estimator (getAPEs(bc) / share_kept rescales alpaca's per-row bias term):
# neither the code paths nor the keys exist any more; the exact FE movers-only keys stay
nj <- names(jsonlite::read_json("numbers.json"))
stopifnot(!any(grepl("ape_bc_kept", nj)), all(c("panel_T4_ape_fe_kept", "panel_T8_ape_fe_kept", "panel_mc_T4_ape_fe_kept_bias",
          "panel_mc_T8_ape_fe_kept_rel_bias") %in% nj), is.null(x$B5$T4$ape_bc_kept), is.null(x$B5$mc$T4$ape_bc_kept_bias))
