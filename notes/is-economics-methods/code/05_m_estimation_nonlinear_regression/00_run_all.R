# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 5 (M-estimation and nonlinear regression),
# then run the tests.
#   01_dgp.R          Larkspur data (M1 to M4) and truth05.rds
#   02_logs.R         block M1: log(1 + y) against PPML            -> res_logs.rds
#   03_nls.R          blocks M2, M2m: NLS, pseudo-true value         -> res_nls.rds
#   04_twostep.R      block M3: two-step variance, Monte Carlo       -> res_twostep.rds
#   05_mc_poisson.R   block M4: Poisson QMLE Monte Carlo             -> res_mc.rds
#   06_figures.R      ../../assets/fig/fig_05_*.svg
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
steps <- c("01_dgp.R", "02_logs.R", "03_nls.R", "04_twostep.R", "05_mc_poisson.R", "06_figures.R", "07_summary.R",
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
