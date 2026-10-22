# tests/helpers.R
# Helpers shared by the chapter 21 tests (sourced from the code directory, not a test itself).

# Ruling R2 layout of a data block: a list of equal-length, non-empty atomic column vectors
# (at least three of them) plus `meta` (a list). Used by test_rd_fns.R (generator output) and
# test_dgp.R (the stored blocks).
is_block <- function(B) {
  cols <- setdiff(names(B), "meta")
  is.list(B) && !is.data.frame(B) && is.list(B$meta) && length(cols) >= 3 &&
    all(vapply(B[cols], function(v) is.atomic(v) && length(v) > 0, TRUE)) &&
    length(unique(vapply(B[cols], length, 1L))) == 1L
}
