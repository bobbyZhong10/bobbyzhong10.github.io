# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 5 (chapter card, verbatim).
# Run from the chapter code directory after 01 to 05 (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("../../tools/R/chapter_kit.R"); tr <- readRDS("truth05.rds")
lg <- readRDS("res_logs.rds"); nl <- readRDS("res_nls.rds"); ts <- readRDS("res_twostep.rds"); mc <- readRDS("res_mc.rds")
stopifnot(abs(lg$log1p$usd - lg$log1p$cents) > 0.2 * abs(lg$log1p$usd))            # log(1+y) changes with the unit
stopifnot(abs(lg$ppml$usd$est - lg$ppml$cents$est) < 1e-6, covers(lg$ppml$usd$ci, tr$M1$log_ratio))
stopifnot(max(abs(nl$gn$theta - nl$nls_coef)) < 1e-6, covers(nl$misspec$ci, tr$M2m$pseudo))
stopifnot(ts$se_corrected > ts$se_naive, ts$mc$cover_naive <= 0.88, ts$mc$cover_corrected >= 0.93)
stopifnot(mc$cover_model <= 0.85, mc$cover_sandwich >= 0.93)                         # overdispersion: the sandwich restores coverage
