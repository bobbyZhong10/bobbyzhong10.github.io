# ------------------------------------------------------------------
# 00_run_all.R -- rebuild everything for chapter 13 (quantile regression), then run the tests.
#   01_dgp.R          Marigold data (Q1, Q1r, Q2, Q3, Q4, Q5, Q5f) and truth13.rds
#   02_qte.R          blocks Q1, Q1r: OLS, marginal QTE, conditional QR      -> res_qte.rds
#   03_rif_dr.R       block Q2: RIF regression, conditional QR, DR            -> res_extra.rds (Q2)
#   04_ivqr.R         block Q3: naive QR, IV-QR, complier QTE                 -> res_extra.rds (Q3)
#   05_panel_mc.R     block Q4: panel QR Monte Carlo (R = 1000)               -> res_extra.rds (Q4)
#   06_cic.R          blocks Q5, Q5f: CIC, QDID, switchers' LQTE; toys        -> res_extra.rds (Q5, cross)
#   07_figures.R      ../../assets/fig/fig_13_*.svg
#   08_summary.R      numbers.json, sessionInfo.txt
#   tests/run_tests.R test_fns, test_dgp, test_conclusions, test_numbers_fresh
# Usage: Rscript 00_run_all.R (from any directory). About 1.5 minutes on the author's machine
# (05 uses up to 8 cores through parallel::mclapply).
# Stops at the first script that fails; exits with status 1 when a test fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_qte.R", "03_rif_dr.R", "04_ivqr.R", "05_panel_mc.R", "06_cic.R", "07_figures.R", "08_summary.R",
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
