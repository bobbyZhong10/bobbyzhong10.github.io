# tests/run_tests.R
# Run every chapter 22 test, each in its own R process, from the code directory:
#   Rscript tests/run_tests.R
# Order: helpers and generators (Task 9), inherited numbers, Task 10 conclusions, number
# registry freshness. Exits with status 1 when any test fails (all of them are run).
if (!file.exists("sc_fns.R") || !dir.exists("tests"))
  stop("run_tests.R: run it from the chapter 22 code directory (Rscript tests/run_tests.R)", call. = FALSE)

tests <- c("test_sc_fns.R", "test_dgp.R", "test_reproduce.R", "test_conclusions.R", "test_numbers_fresh.R")
rscript <- file.path(R.home("bin"), "Rscript")
status <- integer(0)
for (t in tests) {
  t0 <- proc.time()[["elapsed"]]
  out <- suppressWarnings(system2(rscript, file.path("tests", t), stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); if (is.null(st)) st <- 0L
  status[t] <- st
  cat(sprintf("### %-22s exit %d  (%.1f s)\n", t, st, proc.time()[["elapsed"]] - t0))
  cat(paste0("    ", out), sep = "\n")
}
failed <- names(status)[status != 0L]
if (length(failed)) {
  cat(sprintf("run_tests: %d of %d FAILED: %s\n", length(failed), length(tests), paste(failed, collapse = ", ")))
  quit(save = "no", status = 1)
}
cat(sprintf("run_tests: all %d tests passed\n", length(tests)))
