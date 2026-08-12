# Run all Chapter 22 scripts in order (from this directory).
d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
for (f in c("01_generated_regressor.R", "02_algo_collusion.R", "03_figures.R"))
  system(paste("Rscript", shQuote(file.path(d, f))))
