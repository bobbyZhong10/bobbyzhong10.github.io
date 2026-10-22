# =====================================================================
# 22_synthetic_control_factor_models / 08_summary.R
# Write the number registry and the session record of chapter 22, and print a summary:
#   numbers.json        write_numbers_json(collect_numbers22(), ".") (numbers_def.R), every
#                       number the chapter may print, at full precision
#   numbers_cited.json  numbers quoted from papers (filled by Task 11); created as the empty
#                       object {} only when it does not exist, never overwritten
#   sessionInfo.txt     write_session_info() after attaching every package the scripts and
#                       tests use
# Run LAST (after 01 to 07), from this directory:  Rscript 08_summary.R
# =====================================================================

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
source(file.path(here, "..", "..", "tools", "R", "chapter_kit.R"))
source(file.path(here, "numbers_def.R"))

nums <- collect_numbers22(here)
write_numbers_json(nums, here)
cited <- file.path(here, "numbers_cited.json")
if (!file.exists(cited)) writeLines("{}", cited)
# every package a script or test of this directory loads (02-07, sc_fns.R, stag_fns.R, tests),
# plus the solvers tidysynth calls through :: (kernlab::ipop and kernlab::primal, optimx::optimx,
# rgenoud::genoud), which library(tidysynth) does not attach
pkgs <- c("tidysynth", "dplyr", "synthdid", "augsynth", "scpi", "fixest", "did", "fect",
          "quadprog", "lpSolve", "pensynth", "clarabel", "kernlab", "optimx", "rgenoud",
          "ggplot2", "svglite", "patchwork", "jsonlite", "parallel")
write_session_info(here, pkgs)

f3 <- function(v) sprintf("%.3f", v); f4 <- function(v) sprintf("%.4f", v)
n <- nums
cat("\n########################  CHAPTER 22 SUMMARY  ########################\n")
cat(sprintf("numbers.json: %d keys; numbers_cited.json %s; sessionInfo.txt written\n", length(n),
            if (identical(readLines(cited, warn = FALSE), "{}")) "is {}" else "kept as is"))
cat("\n--- S1 (inherited, tidysynth): fee reform in Aster-Metro, truth", f3(n$truth_s1_att_post), "---\n")
cat("  naive DiD", f3(n$s1_core_naive_did), "| SC", f3(n$s1_core_att_post), "pre-RMSPE", f4(n$s1_core_pre_rmspe),
    "| MSPE ratio", f3(n$s1_core_mspe_ratio), "p", f3(n$s1_core_p_perm), "| in-time placebo", f4(n$s1_core_intime_gap), "\n")
cat("--- S1 modern: did_estimate", f4(n$s1_did_uniform), "(= naive DiD) | SDID", f4(n$s1_sdid_est), "[",
    f4(n$s1_sdid_ci_lo), ",", f4(n$s1_sdid_ci_hi), "] | ASC", f4(n$s1_asc_est), "| SC (V = I)", f4(n$s1_sc_hand_est),
    "pre-RMSPE", f4(n$s1_sc_pre_rmspe), "\n")
cat("    conformal p: iid", f4(n$s1_conformal_p), "moving block", f4(n$s1_conformal_p_moving_block),
    "| scpi covers", n$s1_scpi_n_covered, "of 12 months\n")
cat("--- S2 (outside the hull), truth", f4(n$truth_s2_att_post), ": SC", f4(n$s2_sc_est), "pre-RMSPE",
    f4(n$s2_sc_pre_rmspe), "| ASC", f4(n$s2_asc_est), "(", n$s2_asc_n_neg, "negative weights ) | SDID", f4(n$s2_sdid_est),
    "| penalized SC pen", format(n$s2_pen_pen, digits = 3), f4(n$s2_pen_est), "with", n$s2_pen_n_nonzero, "donors\n")
cat("--- S3 (staggered), truth", f4(n$truth_s3_att), ": TWFE", f4(n$s3_twfe_est), "| CS", f4(n$s3_cs_est), "| IFE",
    f4(n$s3_ife_est), "| MC", f4(n$s3_mc_est), "(default", f4(n$s3_mc_default_est), ") | multisynth", f4(n$s3_multisynth_est),
    "| SDID cohorts", f4(n$s3_sdid_cohort_est), "\n")
cat("    screen (24 draws), intervals covering 0.06: TWFE", n$s3_screen_covered_twfe, "| CS", n$s3_screen_covered_cs,
    "| IFE", n$s3_screen_covered_ife, "| MC default", n$s3_screen_covered_mc_default, "| MC refit", n$s3_screen_covered_mc,
    "| multisynth", n$s3_screen_covered_multisynth, "| MC between truth and TWFE", n$s3_screen_mc_between_mc, "\n")
cat("    scpi (scdataMulti, event months 0 to", n$s3_scpi_n_months - 1, "): mean effect", f4(n$s3_scpi_est),
    "| intervals covering 0.06:", n$s3_scpi_n_covered, "of", n$s3_scpi_n_months, "| bounds", f4(n$s3_scpi_ci_lo_min),
    "to", f4(n$s3_scpi_ci_hi_max), "\n")
cat("--- S4 (pilot design), tau", f3(n$truth_s4_tau), ": RMSE synthetic", f4(n$s4_rmse_synthetic_design), "| random + DiM",
    f4(n$s4_rmse_random_design), "| random + DiD", f4(n$s4_rmse_random_did), "\n")
cat("\n####################################################################\n")
