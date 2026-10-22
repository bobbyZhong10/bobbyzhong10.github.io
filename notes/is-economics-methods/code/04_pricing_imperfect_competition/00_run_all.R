# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 4 (pricing and imperfect competition), then run the tests.
#   01_dgp.R            Cirrus data (C1 markets, C2 merchants, C3 bundle values), block C0 by hand, truth04.rds
#   02_logit.R          block C1: OLS pass-through under Bertrand and joint maximisation, homogeneous Cournot
#                       against Bertrand, pass-through curves                                    -> res_logit.rds
#   03_personalized.R   block C2: uniform, coarse and fine pricing, BBM, the output test, Monte Carlo -> res_pd.rds
#   04_bundle.R         block C3: pure bundling of J = 1 to 100 goods                          -> res_bundle.rds
#   05_leader.R         block C4: simultaneous against sequential pricing                      -> res_leader.rds
#   06_figures.R        figures/fig_04_*.svg (Ruling P2)
#   07_summary.R        numbers.json, sessionInfo.txt
#   tests/run_tests.R   test_fns, test_dgp, test_inherit, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 5 seconds on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_logit.R", "03_personalized.R", "04_bundle.R", "05_leader.R", "06_figures.R", "07_summary.R",
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
