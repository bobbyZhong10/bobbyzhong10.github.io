# ------------------------------------------------------------------
# 07_summary.R -- numbers.json and sessionInfo.txt for chapter 4.
# Run LAST, after 01 to 06: numbers.json = collect_numbers04() (numbers_def.R), written by chapter_kit.R's
# write_numbers_json (full precision, keys sorted); sessionInfo.txt lists every package the chapter's scripts
# and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers04(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: lm, uniroot, optimize, quantile (pricing_fns.R, 01 to 05); nleqslv: logit equilibria and price leadership
# (pricing_fns.R); ggplot2 + svglite: figures (06); jsonlite: numbers.json (chapter_kit, tests);
# numDeriv: grad (tests/test_fns.R)
write_session_info(".", c("stats", "nleqslv", "ggplot2", "svglite", "jsonlite", "numDeriv"))
cat("Saved sessionInfo.txt\n")
