# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 13 (chapter card, verbatim).
# Run from the chapter code directory after 01 to 06 (tests/run_tests.R does this; covers() comes
# from tools/R/chapter_kit.R, sourced once by the runner).
#
# Cross-chapter check (card): chapter 20's fuzzy_did_fns.R is read here and only here. The chapter-13
# estimation script (06_cic.R) stores the toy data sets and its own CIC counterfactual numerators
# E(Y11) - E(Q_{D10}(Y10)) in res_extra.rds (x$cross$toys, x$cross$num13); this file computes the
# differences from chapter 20's wald_cic (numerator = estimate * first_stage) into x$cross$cic_vs_ch20
# before the card's assertion runs.
# ------------------------------------------------------------------
e20 <- new.env()
sys.source("../20_did_beyond_binary/fuzzy_did_fns.R", envir = e20)

## ---- chapter card (verbatim; the cross-chapter differences are filled in first) ---------------
tr <- readRDS("truth13.rds"); q <- readRDS("res_qte.rds"); x <- readRDS("res_extra.rds")
stopifnot(length(x$cross$toys) >= 3, length(x$cross$toys) == length(x$cross$num13))
x$cross$cic_vs_ch20 <- mapply(function(toy, num13) {
  w <- e20$wald_cic(toy, "y", "d", "g", "t")
  num13 - w$estimate * w$first_stage
}, x$cross$toys, x$cross$num13)
stopifnot(length(x$cross$cic_vs_ch20) == length(x$cross$toys), all(is.finite(x$cross$cic_vs_ch20)))

stopifnot(covers(q$ols$ci, 0), q$qte$hi[q$qte$tau == 0.1] < 0, q$qte$lo[q$qte$tau == 0.9] > 0)   # mean about zero, opposite signs in the tails
stopifnot(covers(x$Q2$rif$ci, tr$Q2$uqe), !covers(x$Q2$cqr$ci, tr$Q2$uqe))
stopifnot(covers(x$Q3$ivqr$ci, tr$Q3$lqte_median), !covers(x$Q3$qr$ci, tr$Q3$lqte_median))   # IV-QR targets beta_D(tau); lqte_median = beta_D(0.5) = 0.10, the one level where it equals the complier LQTE
stopifnot(abs(x$Q4$canay_T20 - tr$Q4$beta_tau) < abs(x$Q4$dummies_T5 - tr$Q4$beta_tau))
stopifnot(all(mapply(covers, x$Q5$cic$ci, tr$Q5$qtt)), !all(mapply(covers, x$Q5$qdid$ci, tr$Q5$qtt)), covers(x$Q5$lqte$ci, tr$Q5$lqte))
stopifnot(max(abs(x$cross$cic_vs_ch20)) < 1e-10)   # on toy data, equal to chapter 20's wald_cic counterfactual numerator (../20_did_beyond_binary/ is read only in tests)
