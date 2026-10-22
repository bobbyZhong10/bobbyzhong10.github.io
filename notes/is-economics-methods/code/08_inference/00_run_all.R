# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 8 (inference: Delta method, bootstrap,
# cluster, testing), then run the tests.
#   01_dgp.R              Halcyon data (WTP, city rollout, daily series)  -> halcyon.rds
#   02_delta_bootstrap.R  WTP ratio: Delta method against the bootstrap   -> res_delta.rds
#   03_cluster.R          city rollout: iid, HC1, CR1, CR2; size; G = 10   -> res_cluster.rds
#   04_hac.R              AR(1) series: Newey-West, coverage               -> res_hac.rds
#   05_testing.R          logit trinity, Wald reparameterization           -> res_test.rds
#   06_figures.R          ../../assets/fig/fig_08_*.svg
#   07_summary.R          numbers.json, sessionInfo.txt
#   tests/run_tests.R     test_reproduce, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 15 seconds on the author's machine.
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_delta_bootstrap.R", "03_cluster.R", "04_hac.R", "05_testing.R", "06_figures.R",
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
