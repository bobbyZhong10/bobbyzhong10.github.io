# tests/test_reproduce.R
# The inherited R1 numbers that chapter 21 prints (sections 1 and 6), checked against
# res_core.rds written by Task 7's 02_core.R from aster_R1.rds. Run from the code
# directory:  Rscript tests/test_reproduce.R
# Task 6 brief, verbatim; fails until 02_core.R exists (Task 7). Comparisons are at
# printed precision (the legacy rds carries about 1e-7 drift), never identical() on doubles.
r <- readRDS("res_core.rds")
stopifnot(abs(r$naive - 1.891) < 5e-4, abs(r$sharp$est - 0.200) < 5e-4, abs(r$sharp$h - 1.876) < 5e-4,
          all(abs(r$sharp$ci_conv - c(0.133, 0.267)) < 5e-4), all(abs(r$sharp$ci_rob - c(0.111, 0.261)) < 5e-4),
          identical(as.integer(r$sharp$n_eff), c(695L, 689L)), all(abs(r$bw_sens - c(0.186, 0.237)) < 5e-4),
          all(abs(r$poly - c(0.556, -0.215, 0.254)) < 5e-4), abs(r$density$p - 0.250) < 5e-4,
          abs(r$fuzzy$fs - 0.749) < 5e-4, abs(r$fuzzy$rf - 0.170) < 5e-4, abs(r$fuzzy$est - 0.241) < 5e-4,
          all(abs(r$fuzzy$ci_rob - c(0.071, 0.236)) < 5e-4), abs(r$fuzzy$h - 3.333) < 5e-4)
cat("test_reproduce: OK\n")
