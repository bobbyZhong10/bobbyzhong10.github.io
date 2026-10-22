# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 1 (consumer and producer theory, welfare, uncertainty),
# then run the tests.
#   01_dgp.R            Plover data (W1, W1g, W2 and its Bronars choosers), block W0 by hand, truth01.rds
#   02_welfare.R        blocks W1, W1g: dCS, Willig, Vartia CV, representative consumer,
#                       out-of-support extrapolation, the Willig check (Ruling P13)    -> res_welfare.rds
#   03_revealed_pref.R  block W2: GARP, CCEI, Bronars power                            -> res_rp.rds
#   04_risk.R           block W3: protection choice under CARA, Rabin-type calibration -> res_risk.rds
#   05_figures.R        figures/fig_01_*.svg (Ruling P2)
#   06_summary.R        numbers.json, sessionInfo.txt
#   tests/run_tests.R   test_fns, test_dgp, test_inherit, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 15 seconds on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_welfare.R", "03_revealed_pref.R", "04_risk.R", "05_figures.R", "06_summary.R",
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
