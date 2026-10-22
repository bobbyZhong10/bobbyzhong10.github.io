# ------------------------------------------------------------------
# 07_summary.R -- numbers.json and sessionInfo.txt for chapter 9.
# Run LAST, after 01 to 06: numbers.json = collect_numbers09() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers09(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: lm, qr, quantiles (iv_fns, 01 to 05); ivreg: ivreg, summary(diagnostics = TRUE) (02, 05, tests);
# systemfit: 2SLS, 3SLS, SUR (02, 03, tests); sandwich: vcovHC, vcovCL (tests); ivmodel: AR.test, CLR, LIML,
# Fuller (04, tests); gmm: two-step GMM J (tests); ggplot2 + svglite: figures (06); jsonlite: numbers.json
write_session_info(".", c("stats", "ivreg", "systemfit", "sandwich", "ivmodel", "gmm", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
