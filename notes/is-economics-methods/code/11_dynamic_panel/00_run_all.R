# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 11 (dynamic panels), then run the tests.
#   01_dgp.R                 Nimbus data: inherited blocks D1, D1p (legacy gen_dyn, seeds 808 and 7372)
#                            first and unchanged, then D2 to D4, the D5 design and truth11.rds
#   02_d1_more.R             blocks D1, D1p: FE corrections, Anderson-Hsiao, one- and two-step GMM, and the
#                            Anderson-Hsiao (differenced instrument) Monte Carlo at rho = 0.90      -> res_d1.rds
#   03_dynamic_nickell.R     inherited legacy script: pooled, FE, AB, BB, Nickell Monte Carlo D1T
#                            (set.seed(303)), rho -> 1                                          -> res_dynamic.rds
#   04_proliferation.R       block D2 (with a 500-sample Monte Carlo)                           -> res_d2.rds
#   05_predetermined.R       block D3                                                           -> res_d3.rds
#   06_initial_conditions.R  block D4                                                           -> res_d4.rds
#   07_bias_correction_mc.R  block D5: 1,000-sample Monte Carlo per T_obs = 8, 16, 32          -> res_d5.rds
#   08_collect.R             res_extra.rds (the chapter card's field names)
#   09_figures.R             ../../assets/fig/fig_11_*.svg
#   10_summary.R             numbers.json, sessionInfo.txt
#   tests/run_tests.R        test_fns, test_dgp, test_reproduce, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 3 minutes on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_d1_more.R", "03_dynamic_nickell.R", "04_proliferation.R", "05_predetermined.R",
           "06_initial_conditions.R", "07_bias_correction_mc.R", "08_collect.R", "09_figures.R", "10_summary.R",
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
