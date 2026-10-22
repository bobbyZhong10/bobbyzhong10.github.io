# tests/test_numbers_fresh.R
# numbers.json must not be older than the res_*.rds files: recompute the number registry
# with collect_numbers21() (numbers_def.R) and compare it with the file 08_summary.R wrote.
# Run from the code directory:  Rscript tests/test_numbers_fresh.R
#
# Part 1 is the Task 7 brief's test, verbatim. Part 2 adds guards.

## ---- Part 1: the brief ------------------------------------------------------------
source("numbers_def.R"); j <- jsonlite::read_json("numbers.json"); n <- collect_numbers21()
stopifnot(setequal(names(j), names(n)), all(abs(unlist(j[names(n)]) - unlist(n)) < 1e-12))

## ---- Part 2: guards ---------------------------------------------------------------
# the registry is not empty, every value is one finite number, and numbers_cited.json exists
# as a JSON object (Task 8 fills it; it starts as {})
stopifnot(length(n) > 100, all(vapply(n, function(v) is.numeric(v) && length(v) == 1L && is.finite(v), logical(1))),
          length(unlist(j[names(n)])) == length(n), file.exists("numbers_cited.json"),
          is.list(jsonlite::read_json("numbers_cited.json")))

cat("test_numbers_fresh: OK\n")
