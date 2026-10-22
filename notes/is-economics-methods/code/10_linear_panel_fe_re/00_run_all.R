# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 10 (linear panels: FE, RE, FD and CRE), then run the tests.
#   01_dgp.R              Nimbus data: inherited static block P1 (legacy, set.seed(737)) first, then
#                         P2 to P7 and truth10.rds
#   02_static_fe_re.R     inherited legacy script, block P1: pooled OLS, RE, FE, Hausman  -> res_static.rds
#   03_static_more.R      block P1: between, FD, varphi, Mundlak, Hausman, lead test, MD -> res_p1.rds
#   04_fe_fd_feedback.R   blocks P2 (with a 1,000-sample Monte Carlo) and P3            -> res_p23.rds
#   05_ht_slopes.R        blocks P4, P5, P5b                                             -> res_p45.rds
#   06_unbalanced.R       block P6                                                       -> res_p6.rds
#   07_mc_hausman.R       block P7: 2,000-sample Monte Carlo of the Hausman tests        -> res_p7.rds
#   08_collect.R          res_extra.rds (the chapter card's field names)
#   09_figures.R          ../../assets/fig/fig_10_*.svg
#   10_summary.R          numbers.json, sessionInfo.txt
#   tests/run_tests.R     test_fns, test_dgp, test_reproduce, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 20 seconds on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_static_fe_re.R", "03_static_more.R", "04_fe_fd_feedback.R", "05_ht_slopes.R",
           "06_unbalanced.R", "07_mc_hausman.R", "08_collect.R", "09_figures.R", "10_summary.R", "tests/run_tests.R")
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
