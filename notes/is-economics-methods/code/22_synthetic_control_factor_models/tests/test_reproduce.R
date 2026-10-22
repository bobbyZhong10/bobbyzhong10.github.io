# tests/test_reproduce.R
# Locks the inherited S1 numbers (legacy chapter 22 synthetic control on the
# Aster fee reform, aster_S1.rds) at their printed precision.
# Run from the code directory:  Rscript tests/test_reproduce.R
# Reads res_sc_core.rds, which 02_core.R writes (fields naive_did, att_post,
# gap_last, pre_rmspe, mspe_ratio, p_perm, intime_gap and w, the SC weights
# named by donor); run it after 02_core.R.
#
# The brief's test, verbatim:
s <- readRDS("res_sc_core.rds")
stopifnot(abs(s$naive_did - 0.244) < 5e-4, abs(s$att_post - 0.086) < 5e-4, abs(s$gap_last - 0.103) < 5e-4,
          abs(s$pre_rmspe - 0.0165) < 5e-5, abs(s$mspe_ratio - 32.777) < 5e-4, abs(s$p_perm - 1/21) < 1e-12,
          abs(s$intime_gap - (-0.0190)) < 5e-5, abs(s$w[["Donor-10"]] - 0.240) < 5e-4)

cat("test_reproduce.R: all checks passed\n")
