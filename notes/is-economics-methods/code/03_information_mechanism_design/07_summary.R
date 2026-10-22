# ------------------------------------------------------------------
# 07_summary.R -- numbers.json and sessionInfo.txt for chapter 3.
# Run LAST, after 01 to 06: numbers.json = collect_numbers03() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers03(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: integrate, uniroot, splinefun, quantile, rnorm (mech_fns.R, 01 to 05); data.table: the V1 bids and
# the menus (mech_fns.R, 01 to 03); ggplot2 + svglite: figures (06); jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "data.table", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
