# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 7 (chapter card, verbatim),
# plus one added line (marked) for the algebraic MD = GMM identity.
# Run from the chapter code directory after 01 to 04 (tests/run_tests.R does this; covers() comes from
# tools/R/chapter_kit.R, which the runner sources once, Phase 3 ruling 3).
# ------------------------------------------------------------------
tr <- readRDS("truth07.rds"); md <- readRDS("res_md.rds"); g2 <- readRDS("res_gmm_mc.rds"); g3 <- readRDS("res_weak.rds")
stopifnot(abs(md$G1$theta - md$G1$ivw_mean) < 1e-10, covers(md$G1$ci, tr$G1$theta))  # 对角最优权重的 MD = 逆方差加权
stopifnot(abs(md$G1$mc_reject_J - 0.05) < 0.02, md$G1h$J_p < 0.01)                   # 同质时 size 对，异质时拒绝
stopifnot(abs(md$md_vs_gmm) < 0.1 * md$G1$se)                                       # 最优 MD 与有效 GMM 渐近等价
stopifnot(g2$se_eff <= g2$se_identity, g2$cover$two_step_naive < g2$cover$windmeijer, g2$cover$windmeijer >= 0.92)
stopifnot(g3$cover$wald <= 0.85, g3$cover$s_set >= 0.93)                            # 弱矩：Wald 不足，S 置信集稳健
## added (not on the card): when the GMM weight matrix is built from the unrestricted per-market fits, efficient GMM on
## the pooled micro data and optimal MD on the market slopes are the same estimator, to rounding error
stopifnot(abs(md$md_vs_gmm_unrestricted) < 1e-8)
