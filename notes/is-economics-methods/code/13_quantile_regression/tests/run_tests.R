# ------------------------------------------------------------------
# tests/run_tests.R -- run every chapter-13 test file.
#
# Usage (from the chapter code directory):  Rscript tests/run_tests.R
# The script moves to the chapter code directory itself, so it also works when
# called with a full path. Every test file runs with the working directory set to
# the chapter code directory, in its own environment.
#
# Phase 3 ruling 3: tools/R/chapter_kit.R is sourced once here, before the test
# files, so a test file may call covers() without its own source() line.
#
# Each test file is a sequence of stopifnot() calls: it passes when it runs to
# the end. The runner reports PASS or FAIL per file (with the first failing
# assertion) and exits with status 1 when any file fails.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE)))))
})
source("../../tools/R/chapter_kit.R")

known <- c("test_fns.R", "test_dgp.R", "test_conclusions.R", "test_numbers_fresh.R")
found <- sort(list.files("tests", pattern = "^test_.*\\.R$"))
files <- c(intersect(known, found), setdiff(found, known))
missing <- setdiff(known, found)

ok <- setNames(logical(length(files)), files)
for (tf in files) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    source(file.path("tests", tf), local = new.env(parent = globalenv()))
    NULL
  }, error = function(e) conditionMessage(e))
  ok[[tf]] <- is.null(res)
  cat(sprintf("%-24s %s  (%.1f s)%s\n", tf, if (ok[[tf]]) "PASS" else "FAIL",
              proc.time()[["elapsed"]] - t0, if (ok[[tf]]) "" else paste0("\n    ", res)))
}
for (m in missing) cat(sprintf("%-24s MISSING\n", m))
cat(sprintf("%d of %d test files passed\n", sum(ok), length(known) + length(setdiff(found, known))))
if (!all(ok) || length(missing)) quit(status = 1)
