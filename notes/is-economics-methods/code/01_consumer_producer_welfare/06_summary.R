# ------------------------------------------------------------------
# 06_summary.R -- numbers.json and sessionInfo.txt for chapter 1.
# Run LAST, after 01 to 05: numbers.json = collect_numbers01() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers01(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: lm, uniroot, quantile (welfare_fns.R, 01 to 04); ggplot2 + svglite: figures (05);
# jsonlite: numbers.json (chapter_kit, tests); numDeriv: grad (tests/test_fns.R);
# revealedPrefs: checkGarp, a cross-check of garp_check and ccei (tests/test_fns.R, tests/test_conclusions.R)
write_session_info(".", c("stats", "ggplot2", "svglite", "jsonlite", "numDeriv", "revealedPrefs"))
cat("Saved sessionInfo.txt\n")
