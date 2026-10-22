# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 6 (maximum likelihood), then run the tests.
#   01_dgp.R          Sorrel data (C1 to C4) and truth06.rds (incl. the C2 pseudo-true value on 2e6 draws)
#   02_censoring.R    block C1: naive OLS against the censored-data Weibull MLE; A = B; W, LR, LM -> res_cens.rds
#   03_qmle.R         block C2: Weibull QMLE of a log-logistic mixture; Vuong; Monte Carlo -> res_qmle.rds
#   04_grouped.R      block C3: pooled cloglog = grouped-duration MLE                      -> res_grouped.rds
#   05_twostep.R      block C4: two-step MLE with a generated usage index; Monte Carlo      -> res_twostep.rds
#   06_figures.R      ../../assets/fig/fig_06_*.svg
#   07_summary.R      numbers.json, sessionInfo.txt
#   tests/run_tests.R test_fns, test_dgp, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 30 seconds on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_censoring.R", "03_qmle.R", "04_grouped.R", "05_twostep.R", "06_figures.R", "07_summary.R",
           "tests/run_tests.R")
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
