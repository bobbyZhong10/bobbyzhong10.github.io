# ------------------------------------------------------------------
# 08_summary.R -- numbers.json and sessionInfo.txt for chapter 13.
# Run LAST, after 01 to 07: numbers.json = collect_numbers13() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers13(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: lm, glm.fit, density, quantile; quantreg: rq, rq.fit, summary.rq, rq.fit.sfn (02 to 05, tests);
# SparseM: matrix.csr for fe_dummy_qr (05); sandwich: vcovHC (02, 03); parallel: mclapply (05);
# ggplot2 + svglite: figures (07); jsonlite: numbers.json (chapter_kit, tests);
# qte, rifreg: cross-checks in tests/test_fns.R only
write_session_info(".", c("stats", "parallel", "quantreg", "SparseM", "sandwich", "ggplot2", "svglite", "jsonlite",
                          "qte", "rifreg"))
cat("Saved sessionInfo.txt\n")
