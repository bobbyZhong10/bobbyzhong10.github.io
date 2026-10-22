# =====================================================================
# 21_regression_discontinuity / 00_run_all.R
# Run the whole chapter 21 pipeline end to end, each script in its own R process, then
# every test. Run from this directory:  Rscript 00_run_all.R   (about 1.5 minutes on 16
# cores; the Monte Carlo of 03_discrete_honest.R takes most of it)
#   01_dgp.R               data blocks R1 (inherited), R1c, R1h, R2, R3, R4, R4g; truth21.rds
#   02_core.R              R1: inherited sharp, fuzzy, bandwidth, polynomial, density numbers
#   03_discrete_honest.R   R2: conventional, robust, clustered, honest, local randomization; MC
#   04_kink.R              R3: regression kink
#   05_multi.R             R4, R4g: multiple cutoffs, extrapolation, geographic RD
#   06_covs_donut_power.R  R1c, R1h, R1: covariates, heaping and donut, power
#   07_figures.R           ../../assets/fig/fig_21_*.svg
#   08_summary.R           numbers.json, numbers_cited.json (if absent), sessionInfo.txt
#   tests/run_tests.R      test_rd_fns, test_dgp, test_reproduce, test_conclusions,
#                          test_numbers_fresh
# Stops at the first script that fails; exits with status 1 when a test fails.
# Paths: the scripts source rd_fns.R from this directory, and 08_summary.R and the tests
# source ../../tools/R/chapter_kit.R, so this directory must stay at code/<chapter>/ inside
# the book workspace. Cores: the Monte Carlo of 03_discrete_honest.R uses parallel::mclapply
# with CORES = detectCores() - 2 (at most 16). mclapply forks, which Windows does not support:
# with more than one core it stops with an error there, and 03 does not fall back on its own,
# so on Windows set CORES <- 1L in 03_discrete_honest.R (a single core takes about 15
# minutes). Each draw has its own seed, so the results do not depend on the number of cores.
# =====================================================================
if (!file.exists("rd_fns.R") || !file.exists("01_dgp.R"))
  stop("00_run_all.R: run it from the chapter 21 code directory (Rscript 00_run_all.R)", call. = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_core.R", "03_discrete_honest.R", "04_kink.R", "05_multi.R",
           "06_covs_donut_power.R", "07_figures.R", "08_summary.R")
t_all <- proc.time()[["elapsed"]]
for (s in steps) {
  cat("\n>>>>>>>>>>>>>>>>>>>>  ", s, "  <<<<<<<<<<<<<<<<<<<<\n")
  t0 <- proc.time()[["elapsed"]]
  st <- system2(rscript, s)
  cat(sprintf("[%s: exit %d, %.1f s]\n", s, st, proc.time()[["elapsed"]] - t0))
  if (st != 0) stop("00_run_all.R: ", s, " failed (exit ", st, ")", call. = FALSE)
}
cat("\n>>>>>>>>>>>>>>>>>>>>   tests/run_tests.R   <<<<<<<<<<<<<<<<<<<<\n")
st <- system2(rscript, file.path("tests", "run_tests.R"))
cat(sprintf("\n00_run_all.R: pipeline and tests finished in %.0f s\n", proc.time()[["elapsed"]] - t_all))
if (st != 0) quit(save = "no", status = 1)
