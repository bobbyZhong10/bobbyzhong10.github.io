# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 3 (information economics, contracts, mechanism design),
# then run the tests.
#   01_dgp.R            Vireo data (V1 auctions, V4 market), block V0 by hand, V2 and V3 parameters, truth03.rds
#   02_auction.R        block V1: revenue equivalence, the optimal reserve, quantile GPV, Monte Carlo -> res_auction.rds
#   03_menu.R           block V2: optimal SaaS menus with and without ironing, two types           -> res_menu.rds
#   04_contract.R       block V3: the LEN contract with extra signals, small simulation          -> res_len.rds
#   05_matching.R       block V4: deferred acceptance both ways, manipulation checks             -> res_match.rds
#   06_figures.R        figures/fig_03_*.svg (Ruling P2)
#   07_summary.R        numbers.json, sessionInfo.txt
#   tests/run_tests.R   test_fns, test_dgp, test_inherit, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 30 seconds on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_auction.R", "03_menu.R", "04_contract.R", "05_matching.R", "06_figures.R", "07_summary.R",
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
