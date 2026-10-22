# ------------------------------------------------------------------
# tests/test_numbers_fresh.R -- numbers.json must not be older than the res_*.rds files
# (Phase 2 Task 7): its keys and values equal collect_numbers10() to 1e-12.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("numbers_def.R"); j <- jsonlite::read_json("numbers.json"); n <- collect_numbers10()
stopifnot(setequal(names(j), names(n)), all(abs(unlist(j[names(n)]) - unlist(n)) < 1e-12))
# numbers_cited.json is a JSON object (empty for this chapter until Task 10.2 adds cited numbers)
stopifnot(is.list(jsonlite::read_json("numbers_cited.json")))
