# ------------------------------------------------------------------
# tests/test_reproduce.R -- the inherited Helix block B1 (legacy code: 01_dgp.R lines 1 to 77 with
# set.seed(512) and set.seed(9001), 02_lpm_naive.R, 03_mle_logit.R with set.seed(778)) reproduces the
# legacy numbers. Chapter card, verbatim. res_helix.rds is written by 04_helix.R; its field names
# organise what the legacy 02_lpm_naive.R prints (lpm) and what 03_mle_logit.R prints and saves (logit).
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
h <- readRDS("res_helix.rds")
stopifnot(abs(h$lpm$slope - 0.2153) < 5e-5, abs(h$lpm$max_phat - 1.2133) < 5e-5,
          abs(h$logit$beta1 - 1.0969) < 5e-5, abs(h$logit$ame - 0.2171) < 5e-5)
