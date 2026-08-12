# Run all Chapter 21 scripts in order (from this directory).
d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
for (f in c("01_dgp.R", "02_estimate.R", "03_figures.R"))
  system(paste("Rscript", shQuote(file.path(d, f))))
