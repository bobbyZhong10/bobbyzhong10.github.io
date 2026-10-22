# ------------------------------------------------------------------
# 07_summary.R -- numbers.json and sessionInfo.txt for chapter 5.
# Run LAST, after 01 to 06: numbers.json = collect_numbers05() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers05(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: nls, glm, lm; sandwich: sandwich, vcovCL; numDeriv: hessian, jacobian (03, tests);
# fixest: fepois (02); ggplot2 + svglite: figures (06); jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "sandwich", "numDeriv", "fixest", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
