# ------------------------------------------------------------------
# 10_summary.R -- numbers.json and sessionInfo.txt for chapter 11 (replaces the legacy 05_summary.R,
# whose reconciliation lines for the dynamic panel are printed below).
# Run LAST, after 01 to 09: numbers.json = collect_numbers11() (numbers_def.R), written by
# chapter_kit.R's write_numbers_json (full precision, keys sorted); sessionInfo.txt lists every
# package the chapter's scripts and tests use.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("../../tools/R/chapter_kit.R")
source("numbers_def.R")

## ---- the legacy reconciliation of the dynamic panel (blocks D1, D1p, D1T) --------------------------
dy <- readRDS("res_dynamic.rds")
L <- function(...) cat(sprintf(...), "\n")
cat("--- dynamic panel (rho, truth 0.70) ---\n")
L("pooled OLS = %.4f (up)", dy$pool)
L("within/FE  = %.4f (Nickell, down)", dy$fe)
L("Arellano-Bond (diff GMM)  = %.4f", dy$ab)
L("Blundell-Bond (sys GMM)   = %.4f", dy$bb)
# dy$sargan_p, the legacy "Sargan p 0.290", is the two-step Hansen J p-value (plm::sargan(fit, "twosteps") of the
# two-step AB fit, the "Sargan test" line of summary.pgmm); it equals d1_ab2_J_p in numbers.json. It is not a
# one-step Sargan statistic (plm's one-step version is unscaled and is not stored: dyn_fns.R, gmm_summary).
L("AB tests: Sargan p (two-step Hansen J) = %.3f | AR(1) p = %.4f | AR(2) p = %.3f", dy$sargan_p, dy$ar1_p, dy$ar2_p)
cat("--- Nickell bias vs T_obs (FE / AB / legacy approx rho - (1 + rho)/(T_obs - 1); truth 0.70) ---\n")
print(round(dy$nick, 4))
cat("--- rho near 1 (truth 0.90) ---\n")
L("difference GMM = %.4f (SE %.4f)", dy$ab9, dy$se_ab9)
L("system GMM     = %.4f (SE %.4f)", dy$bb9, dy$se_bb9)

nums <- collect_numbers11(".")
write_numbers_json(nums, ".")
cat(sprintf("Saved numbers.json (%d numbers)\n", length(nums)))

# stats: lm, model.matrix, rnorm; plm: plm, pgmm, pdata.frame, sargan, mtest, vcovHC (01 to 07, tests);
# sandwich: vcovCL (tests); AER: ivreg (tests); tools: md5sum (tests); ggplot2 + svglite: figures (09);
# jsonlite: numbers.json (chapter_kit, tests)
write_session_info(".", c("stats", "tools", "plm", "sandwich", "AER", "ggplot2", "svglite", "jsonlite"))
cat("Saved sessionInfo.txt\n")
