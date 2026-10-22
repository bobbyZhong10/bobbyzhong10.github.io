# =====================================================================
# 22_synthetic_control_factor_models / 00_run_all.R
# Run the whole chapter 22 pipeline end to end, each script in its own R process, then
# every test. Run from this directory:  Rscript 00_run_all.R
#   01_dgp.R            data blocks S1 (inherited), S2, S3 (gen_S3 in dgp_s3.R), S4; truth22.rds
#   02_core.R           S1: the inherited tidysynth analysis (res_sc_core.rds)
#   03_s1_modern.R      S1: SDID, synthdid's DiD = the 2x2 DiD, ASC, conformal p, scpi
#   04_hull.R           S2: outside the hull; SC, ASC, SDID, penalized SC (LOO-CV penalty)
#   05_staggered.R      S3: TWFE, CS, IFE, MC (default and refit), multisynth, SDID by cohort
#   05b_stag_screen.R   S3: the 24-seed screen (Ruling S6). SLOW (about 2.5 minutes on 8
#                       cores): its result is cached in res_stag_screen.rds and recomputed only
#                       when that file is missing or RUN_SLOW=1 is set:
#                           RUN_SLOW=1 Rscript 00_run_all.R
#                       otherwise the script checks the cached screen against res_stag.rds
#                       and the installed packages (it stops with a RUN_SLOW=1 hint when the
#                       cache is stale) and copies it into res_stag.rds
#   05c_stag_scpi.R     S3: scpi prediction intervals by event month (scdataMulti; Ruling S12),
#                       res_stag_scpi.rds; about 2 minutes
#   06_design.R         S4: synthetic experimental design vs random assignment
#   07_figures.R        ../../assets/fig/fig_22_*.svg
#   08_summary.R        numbers.json, numbers_cited.json (if absent), sessionInfo.txt
#   (helpers: sc_fns.R, dgp_s3.R, stag_fns.R; estimation seeds in est_seeds.R)
#   tests/run_tests.R   test_sc_fns, test_dgp, test_reproduce, test_conclusions,
#                       test_numbers_fresh
# Stops at the first script that fails; exits with status 1 when a test fails.
# Wall time with the cached screen: about 5 minutes, about 2 of them in 05c_stag_scpi.R
# (the exact time depends on the machine and the number of cores).
# 06_design.R prints the wall time of the S4 design search but does not store it, so a
# rerun reproduces res_design.rds byte for byte.
# =====================================================================
if (!file.exists("sc_fns.R") || !file.exists("01_dgp.R"))
  stop("00_run_all.R: run it from the chapter 22 code directory (Rscript 00_run_all.R)", call. = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
steps <- c("01_dgp.R", "02_core.R", "03_s1_modern.R", "04_hull.R", "05_staggered.R", "05b_stag_screen.R",
           "05c_stag_scpi.R", "06_design.R", "07_figures.R", "08_summary.R")
cat(if (!file.exists("res_stag_screen.rds") || identical(Sys.getenv("RUN_SLOW"), "1"))
      "05b_stag_screen.R will RECOMPUTE the 24-seed screen (slow)\n"
    else "05b_stag_screen.R will use the cached screen res_stag_screen.rds (RUN_SLOW=1 recomputes it)\n")
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
