# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 2 (game theory), then run the tests.
#   01_dgp.R            Osprey-Boreal data (G1 to G4), block G0 by hand, truth02.rds
#   02_qre.R            block G1: QRE inversion by OLS on group log-odds (unique, pooled,
#                       by market) and its Monte Carlo coverage                    -> res_qre.rds
#   03_green_porter.R   block G2: on-path price wars, one-shot deviation checks     -> res_gp.rds
#   04_learning.R       block G3: regret matching, fictitious play                  -> res_learn.rds
#   05_global_game.R    block G4: thresholds as noise shrinks, uniqueness bound     -> res_gg.rds
#   06_figures.R        figures/fig_02_*.svg (Ruling P2)
#   07_summary.R        numbers.json, sessionInfo.txt
#   tests/run_tests.R   test_fns, test_dgp, test_inherit, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 25 seconds on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_qre.R", "03_green_porter.R", "04_learning.R", "05_global_game.R", "06_figures.R",
           "07_summary.R", "tests/run_tests.R")
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
