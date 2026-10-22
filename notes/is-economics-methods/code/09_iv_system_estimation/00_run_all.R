# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 9 (instrumental variables and system estimation),
# then run the tests.
#   01_dgp.R          Bramble data (S1, S1m, S2, S3, S4, S4v) and truth09.rds
#   02_system.R       blocks S1, S1m: OLS, 2SLS, 3SLS, contagion, Monte Carlo   -> res_system.rds
#   02b_contagion_r0.R block S1m at corr(coupon, fee) = 0 and -0.7, n = 1e6    -> res_contagion_r0.rds
#   03_sur.R          block S2: SUR with the same and with different regressors -> res_sur.rds
#   04_weak.R         block S3: weak instrument, Wald against Anderson-Rubin      -> res_weak.rds
#   05_overid.R       blocks S4, S4v: Hansen J and Sargan, Monte Carlo           -> res_overid.rds
#   06_figures.R      ../../assets/fig/fig_09_*.svg
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
steps <- c("01_dgp.R", "02_system.R", "02b_contagion_r0.R", "03_sur.R", "04_weak.R", "05_overid.R", "06_figures.R", "07_summary.R",
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
