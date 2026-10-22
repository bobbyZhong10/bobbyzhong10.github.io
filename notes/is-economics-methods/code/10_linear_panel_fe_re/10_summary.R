# ------------------------------------------------------------------
# 10_summary.R -- numbers.json and sessionInfo.txt for chapter 10 (replaces the legacy 05_summary.R,
# whose reconciliation lines for the static panel are printed below).
# Run LAST, after 01 to 09: numbers.json = collect_numbers10() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists every
# package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

## ---- the legacy reconciliation of the static panel (block P1) ------------------------------------
st <- readRDS("res_static.rds")
L <- function(...) cat(sprintf(...), "\n")
cat("--- static panel (beta, truth 0.50) ---\n")
L("pooled OLS = %.4f (SE %.4f)", st$pool, st$se_pool)
L("random effects = %.4f", st$re)
L("fixed effects = %.4f (SE %.4f; cluster-by-merchant %.4f)", st$fe, st$se_fe, st$se_fe_cl)
L("Hausman chi2 = %.1f  p = %.3g", st$haus_stat, st$haus_p)

nums <- collect_numbers10(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: lm, lm.fit, qr; plm: plm, pdata.frame, phtest, ercomp, vcovHC, pwartest, pwfdtest, pvcm (02 to 07,
# tests); sandwich: vcovHC, vcovCL (05, tests); fixest: feols (tests); clubSandwich: vcovCR (tests);
# ggplot2 + svglite: figures (09); jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "plm", "sandwich", "fixest", "clubSandwich", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
