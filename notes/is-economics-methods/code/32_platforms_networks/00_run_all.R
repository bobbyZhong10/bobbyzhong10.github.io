# Run all Chapter 32 scripts in order (from this directory).
# Rscript writes spaces in the --file= path as "~+~", so decode them, then run
# each script by its relative name from this directory (the workspace path
# contains spaces, which would otherwise break each script's own path lookup).
d <- dirname(normalizePath(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)))
setwd(d)
for (f in c("01_dgp.R", "02_estimate.R", "04_twosided.R", "03_figures.R"))
  system2("Rscript", f)
