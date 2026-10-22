# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 7 (GMM and minimum distance), then run the tests.
#   01_dgp.R          Tamarind data (G1, G1h, G2, G3) and truth07.rds
#   02_md.R           blocks G1, G1h: minimum distance, J, Monte Carlo; pooled GMM   -> res_md.rds
#   03_gmm_mc.R       block G2: two-step, iterated, CUE, Windmeijer; Monte Carlo     -> res_gmm_mc.rds
#   04_weak.R         block G3: weak moments, Wald against the S set; Monte Carlo     -> res_weak.rds
#   05_figures.R      ../../assets/fig/fig_07_*.svg
#   06_summary.R      numbers.json, sessionInfo.txt
#   tests/run_tests.R test_fns, test_dgp, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 3 minutes on the author's machine
# (03_gmm_mc.R takes most of it).
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_md.R", "03_gmm_mc.R", "04_weak.R", "05_figures.R", "06_summary.R", "tests/run_tests.R")
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
