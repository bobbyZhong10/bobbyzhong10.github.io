# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 12 (binary and multinomial response), then run the tests.
#   01_dgp.R            Helix data: inherited block B1 (legacy, set.seed(512) and set.seed(9001)) first,
#                       then B2 to B7 and truth12.rds
#   02_lpm_naive.R      inherited legacy script, block B1: the LPM and its fitted values (prints only)
#   03_mle_logit.R      inherited legacy script, block B1: logit MLE, SEs, AME, consistency simulation
#                       (set.seed(778))                                                    -> res_mle.rds
#   04_helix.R          block B1 under the card's field names, probit comparison, legacy
#                       reconciliation print (the former 04_summary.R)                     -> res_helix.rds
#   05_cross_section.R  blocks B2 (interaction), B3 (control function, 499 bootstrap draws),
#                       B4 (omitted heterogeneity)                                          -> res_cs.rds
#   06_panel.R          block B5: FE logit, bias correction, CMLE, CRE; 1,000-panel Monte Carlo per T;
#                       the one-sample T = 4 panel is months 1 to 4 of the T = 8 panel (own truths),
#                       the Monte Carlo T = 4 a 4-month DGP                                  -> res_panel.rds
#   07_multinomial.R    blocks B6 (ordered probit, 1,000-sample Monte Carlo) and B7 (MNL, Hausman-McFadden
#                       with the stacked-score contrast, 4,000-sample Monte Carlo)          -> res_mult.rds
#   08_collect.R        res_extra.rds (the chapter card's field names)
#   09_figures.R        ../../assets/fig/fig_12_*.svg
#   10_summary.R        numbers.json, sessionInfo.txt
#   tests/run_tests.R   test_fns, test_dgp, test_reproduce, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About two and a half minutes on the author's machine
# (the Monte Carlo runs use up to 8 cores through parallel::mclapply; results do not depend on the number
# of cores because every replication sets its own seed).
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_lpm_naive.R", "03_mle_logit.R", "04_helix.R", "05_cross_section.R", "06_panel.R",
           "07_multinomial.R", "08_collect.R", "09_figures.R", "10_summary.R", "tests/run_tests.R")
t_all <- proc.time()[["elapsed"]]
for (s in steps) {
  cat(sprintf("\n== %s\n", s))
  t0 <- proc.time()[["elapsed"]]
  status <- system2(rscript, shQuote(s))
  cat(sprintf("== %s finished with status %d (%.1f s)\n", s, status, proc.time()[["elapsed"]] - t0))
  if (status != 0) {
    cat(sprintf("STOPPED: %s failed\n", s))
    quit(status = 1)
  }
}
cat(sprintf("\nAll steps done in %.1f s\n", proc.time()[["elapsed"]] - t_all))
