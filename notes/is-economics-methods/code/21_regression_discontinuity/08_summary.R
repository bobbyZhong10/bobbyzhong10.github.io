# =====================================================================
# 21_regression_discontinuity / 08_summary.R
# Write the number registry and the session record of chapter 21, and print a summary:
#   numbers.json        write_numbers_json(collect_numbers21(), ".") (numbers_def.R), every
#                       number the chapter may print, at full precision
#   numbers_cited.json  numbers quoted from papers (filled by Task 8); created as the empty
#                       object {} only when it does not exist, never overwritten
#   sessionInfo.txt     write_session_info() after attaching every package the scripts use
# Run LAST (after 01 to 07), from this directory:  Rscript 08_summary.R
# =====================================================================

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
source(file.path(here, "..", "..", "tools", "R", "chapter_kit.R"))
source(file.path(here, "numbers_def.R"))

nums <- collect_numbers21(here)
write_numbers_json(nums, here)
cited <- file.path(here, "numbers_cited.json")
if (!file.exists(cited)) writeLines("{}", cited)
pkgs <- c("rdrobust", "rddensity", "lpdensity", "rdlocrand", "rdmulti", "RDHonest", "rdpower", "R.utils",
          "ggplot2", "dplyr", "svglite", "patchwork", "jsonlite", "parallel")
write_session_info(here, pkgs)

f3 <- function(v) sprintf("%.3f", v)
n <- nums
cat("\n########################  CHAPTER 21 SUMMARY  ########################\n")
cat(sprintf("numbers.json: %d keys; numbers_cited.json %s; sessionInfo.txt written\n", length(n),
            if (identical(readLines(cited, warn = FALSE), "{}")) "is {}" else "kept as is"))
cat("\n--- R1 (inherited): Preferred badge, true jump 0.200 ---\n")
cat("  naive", f3(n$r1_naive), "| sharp", f3(n$r1_sharp_est), "robust CI [", f3(n$r1_sharp_ci_rob_lo), ",",
    f3(n$r1_sharp_ci_rob_hi), "] h", f3(n$r1_sharp_h), "| fuzzy", f3(n$r1_fuzzy_est), "| density p", f3(n$r1_density_p), "\n")
cat("  power for 0.20:", sprintf("%.4f", n$r1_power_power), "| MDE", f3(n$r1_power_mde), "\n")
cat("--- R1c covariates: CI length ratio", f3(n$r1c_ratio_len_rob), "estimate change", sprintf("%.4f", n$r1c_ratio_est_change), "\n")
cat("--- R1h heaping: density p", sprintf("%.2g", n$r1h_heap_density_p), "naive", f3(n$r1h_heap_naive_est), "donut",
    f3(n$r1h_donut_est), "\n")
cat("--- R2 discrete (truth 0.150): conv", f3(n$r2_conv_est), "honest(M_bound) [", f3(n$r2_honest_mbound_ci_lo), ",",
    f3(n$r2_honest_mbound_ci_hi), "] locrand [", f3(n$r2_locrand_ci_lo), ",", f3(n$r2_locrand_ci_hi), "]\n")
cat(sprintf("    MC (%d draws): coverage conv %.3f robust %.3f crv %.3f honest %.3f honest_lb %.3f lr %.3f\n", n$r2_mc_reps,
            n$r2_mc_cover_conv, n$r2_mc_cover_robust, n$r2_mc_cover_crv, n$r2_mc_cover_honest, n$r2_mc_cover_honest_lb,
            n$r2_mc_cover_lr))
cat("--- R3 kink (truth -0.050):", sprintf("%.4f", n$r3_rk_est), "robust CI [", sprintf("%.4f", n$r3_rk_ci_lo), ",",
    sprintf("%.4f", n$r3_rk_ci_hi), "] | hand local linear", sprintf("%.4f", n$r3_rk_hand_est), "| level jump",
    sprintf("%.4f (%.2f se)", n$r3_level_est, n$r3_level_t), "\n")
cat("--- R4 multi-cutoff: by cutoff", f3(n$r4_c70_est), f3(n$r4_c80_est), f3(n$r4_c90_est), "| weights",
    f3(n$r4_c70_weight), f3(n$r4_c80_weight), f3(n$r4_c90_weight), "| pooled", f3(n$r4_pooled_est), "normalized",
    f3(n$r4_pooled_norm_est), "| extrap", f3(n$r4_extrap_est), "(truth 0.325)\n")
cat("--- R4g geographic: share of boundary CIs covering", f3(n$r4g_share_covered), "\n")
cat("\n####################################################################\n")
