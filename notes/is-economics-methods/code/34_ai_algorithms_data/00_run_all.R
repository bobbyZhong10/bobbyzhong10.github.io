# Run all Chapter 34 scripts in order (from this directory).
# Rscript writes spaces in --file= as "~+~"; decode them so paths with spaces work.
f0 <- gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)
d <- dirname(normalizePath(f0))
setwd(d)  # each script is called by its relative name; 04 writes relative to this directory
for (f in c("01_generated_regressor.R", "02_algo_collusion.R", "03_figures.R",
            "04_adoption_complements.R"))
  system(paste("Rscript", shQuote(f)))
