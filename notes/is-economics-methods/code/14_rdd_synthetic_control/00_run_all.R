# =====================================================================
# 14_rdd_synthetic_control / 00_run_all.R
# Run the whole Module 12 pipeline end to end.
# Usage:  Rscript 00_run_all.R
# =====================================================================
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
for (s in c("01_dgp_rdd.R", "02_rdd_estimate.R",
            "03_dgp_sc.R", "04_sc_estimate.R",
            "05_figures.R", "06_summary.R")) {
  cat("\n>>>>>>>>>>>>>>>>>>>>  ", s, "  <<<<<<<<<<<<<<<<<<<<\n")
  source(file.path(here, s), local = new.env(), chdir = FALSE)
}
