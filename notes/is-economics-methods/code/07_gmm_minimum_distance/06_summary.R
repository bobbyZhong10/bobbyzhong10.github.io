# ------------------------------------------------------------------
# 06_summary.R -- numbers.json and sessionInfo.txt for chapter 7.
# Run LAST, after 01 to 05: numbers.json = collect_numbers07() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers07(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: lm, optim, pchisq; numDeriv: jacobian, grad (gmm_fns.R); sandwich: vcovHC (03, tests);
# gmm: gmm, specTest and plm: pgmm, vcovHC (tests/test_fns.R cross-checks); ggplot2 + svglite: figures (05);
# jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "numDeriv", "sandwich", "gmm", "plm", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
