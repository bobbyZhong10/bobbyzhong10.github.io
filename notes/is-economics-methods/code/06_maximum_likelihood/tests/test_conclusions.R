# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 6 (chapter card, verbatim).
# Run from the chapter code directory after 01 to 05 (tests/run_tests.R does this).
# The source() line is added so that the file also runs on its own; the six stopifnot() lines
# and the two readRDS() lines are the card's, character for character, except one assertion in
# the third line replaced by controller Ruling 29 (see the comment above it).
# ------------------------------------------------------------------
source("../../tools/R/chapter_kit.R")
tr <- readRDS("truth06.rds"); c1 <- readRDS("res_cens.rds"); c2 <- readRDS("res_qmle.rds")
c3 <- readRDS("res_grouped.rds"); c4 <- readRDS("res_twostep.rds")
stopifnot(!covers(c1$ols_observed$ci, tr$C1$beta), covers(c1$mle$ci, tr$C1$beta))           # 删失下朴素 OLS 错，MLE 对
stopifnot(abs(c1$se_hess / c1$se_opg - 1) < 0.10, abs(c2$se_hess / c2$se_opg - 1) > 0.30)      # 设定对 A = B，设定错分家
# Ruling 29: single-sample coverage is a 5% event; the claim rests on the Monte Carlo (card had covers(c2$qmle$ci, tr$C2$pseudo))
stopifnot(c2$mc$cover_model <= 0.88, c2$mc$cover_sandwich >= 0.93, abs(c2$mc$mean_est - tr$C2$pseudo) < 0.003)  # MC 不少于 1000 次
stopifnot(max(abs(c1$trinity$W - c1$trinity$LR), abs(c1$trinity$LM - c1$trinity$LR)) < 0.10 * c1$trinity$LR)
stopifnot(max(abs(c3$cloglog_pooled - c3$grouped_mle)) < 1e-6)                               # pooled cloglog = grouped MLE
stopifnot(c4$se_corrected >= c4$se_naive, covers(c4$corrected$ci, tr$C4$beta))
# the card's "MC 不少于 1000 次" made checkable (it is a comment on the card)
stopifnot(c2$mc$R >= 1000)
