# ------------------------------------------------------------------
# 07_summary.R -- numbers.json and sessionInfo.txt for chapter 6.
# Run LAST, after 01 to 06: numbers.json = collect_numbers06() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers06(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: glm, glm.fit, optim, integrate (tests); survival: survreg, survfit (02, 03, 06, tests);
# maxLik: maxLik, maxLik::hessian (ml_fns.R fit_weibull, tests); numDeriv: hessian, jacobian, grad
# (info_hat, trinity, pseudo_true_weibull, tests); sandwich: vcovCL (04), sandwich, vcovHC (tests);
# lmtest: lrtest (02); ggplot2 + patchwork + svglite: figures (06); jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "survival", "maxLik", "numDeriv", "sandwich", "lmtest", "ggplot2", "patchwork",
                          "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
