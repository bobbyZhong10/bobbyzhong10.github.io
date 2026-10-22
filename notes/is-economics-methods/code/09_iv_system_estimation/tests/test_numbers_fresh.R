# ------------------------------------------------------------------
# tests/test_numbers_fresh.R -- numbers.json must not be older than the res_*.rds files
# (Phase 2 Task 7): its keys and values equal collect_numbers09() to 1e-12.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("numbers_def.R"); j <- jsonlite::read_json("numbers.json"); n <- collect_numbers09()
stopifnot(setequal(names(j), names(n)), all(abs(unlist(j[names(n)]) - unlist(n)) < 1e-12))
# the AR set's shape must not change the key set: bounded flag and number of pieces are always there, the end points
# lo / hi exist only when finite (bounded interval or two rays)
stopifnot(all(c("weak_sample_ar_bounded", "weak_sample_ar_n_pieces", "weak_sample_ar_type") %in% names(n)),
          n$weak_sample_ar_bounded == (n$weak_sample_ar_type == 1), n$weak_sample_ar_n_pieces == c(1, 2, 1)[n$weak_sample_ar_type],
          all(c("weak_sample_ar_lo", "weak_sample_ar_hi") %in% names(n)) == (n$weak_sample_ar_type %in% 1:2))
# numbers_cited.json is a JSON object (empty for this chapter until Task 09.2 adds cited numbers)
stopifnot(is.list(jsonlite::read_json("numbers_cited.json")))
