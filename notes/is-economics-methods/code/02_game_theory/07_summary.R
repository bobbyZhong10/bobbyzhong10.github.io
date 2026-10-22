# ------------------------------------------------------------------
# 07_summary.R -- numbers.json and sessionInfo.txt for chapter 2.
# Run LAST, after 01 to 06: numbers.json = collect_numbers02() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers02(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: uniroot, optimize, qlogis, rbinom (games_fns.R, 01 to 05); nleqslv: logit_qre (games_fns.R);
# data.table: gp_simulate and the G2 panel (games_fns.R, 01, 03); ggplot2 + svglite: figures (06);
# jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "nleqslv", "data.table", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
