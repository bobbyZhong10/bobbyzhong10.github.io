# ------------------------------------------------------------------
# 10_summary.R -- numbers.json and sessionInfo.txt for chapter 12.
# Run LAST, after 01 to 09: numbers.json = collect_numbers12() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists
# every package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

nums <- collect_numbers12(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: glm, lm, glm.fit, optim, integrate; sandwich: vcovHC, vcovCL (03, 04, 05, 06, tests); lmtest (legacy 03);
# numDeriv: grad, hessian, jacobian (bin_fns.R); parallel: mclapply (05, 06, 07);
# marginaleffects: avg_slopes checks (05, 07, tests); survival: clogit (06, tests); alpaca: feglm, biasCorr,
# getAPEs (06); lme4: glmer (06); ordinal: clm (07, tests); MASS: polr (07, tests); nnet: multinom (07, bin_fns.R);
# mlogit and dfidx: hmftest check (07, tests); ggplot2 + svglite: figures (09); jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "parallel", "sandwich", "lmtest", "numDeriv", "marginaleffects", "survival", "alpaca",
                          "lme4", "ordinal", "MASS", "nnet", "mlogit", "dfidx", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
