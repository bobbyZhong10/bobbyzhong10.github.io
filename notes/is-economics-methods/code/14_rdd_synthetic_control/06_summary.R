# =====================================================================
# 14_rdd_synthetic_control / 06_summary.R
# Collect every quotable number for Module 12 and print to console.
# Run 01..05 first (or run 00_run_all.R).
# =====================================================================

in_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(in_dir) == 0 || in_dir == "") in_dir <- "."

rdd  <- readRDS(file.path(in_dir, "res_rdd.rds"))
sc   <- readRDS(file.path(in_dir, "res_sc.rds"))

f <- function(x) formatC(x, format = "f", digits = 3)

cat("\n########################  MODULE 12 SUMMARY  ########################\n")
cat("\n--- RDD: Preferred badge on log GMV (sharp) ---\n")
cat("  true cutoff jump              :", f(rdd$tau_true), "\n")
cat("  naive above-vs-below diff     :", f(rdd$naive_diff), "\n")
cat("  rdrobust conventional est     :", f(rdd$rd_conv_est),
    "  95% CI [", f(rdd$rd_conv_ci[1]), ",", f(rdd$rd_conv_ci[2]), "]\n")
cat("  rdrobust robust  95% CI       : [", f(rdd$rd_robust_ci[1]), ",",
    f(rdd$rd_robust_ci[2]), "]  (p =", formatC(rdd$rd_robust_p, format="e", digits=2), ")\n")
cat("  CCT/MSE-optimal bandwidth     :", f(rdd$bw_cct),
    " (eff. n L/R =", rdd$n_eff_left, "/", rdd$n_eff_right, ")\n")
cat("  bandwidth sensitivity h/2, 2h :", f(rdd$rd_bw_half), ",", f(rdd$rd_bw_dbl), "\n")
cat("  global poly p=4 / 6 / 8       :", f(rdd$poly4_est), "/",
    f(rdd$poly6_est), "/", f(rdd$poly8_est), "  (unstable)\n")
cat("  McCrary/rddensity p-value     :", f(rdd$mccrary_p),
    " (T =", f(rdd$mccrary_T), "; no manipulation)\n")

cat("\n--- RDD: fuzzy variant (cutoff IV / LATE, ties to Ch.10) ---\n")
cat("  first-stage jump est / true   :", f(rdd$first_stage_est), "/",
    f(rdd$first_stage_true), "\n")
cat("  reduced-form jump est         :", f(rdd$reduced_form_est), "\n")
cat("  fuzzy-RD LATE est             :", f(rdd$fuzzy_est),
    "  robust 95% CI [", f(rdd$fuzzy_robust_ci[1]), ",",
    f(rdd$fuzzy_robust_ci[2]), "]\n")

cat("\n--- Synthetic control: Aster fee reform on metro log GMV ---\n")
cat("  true post-period avg effect   :", f(sc$true_post_avg), "\n")
cat("  true effect at last month     :", f(sc$true_last), "\n")
cat("  naive 2x2 DiD (FAILS)         :", f(sc$did_est), "\n")
cat("  SC post-period avg gap (est)  :", f(sc$sc_post_avg), "\n")
cat("  SC gap at last month          :", f(sc$sc_post_last), "\n")
cat("  pre-period RMSPE              :", formatC(sc$pre_rmspe, format="f", digits=4), "\n")
cat("  in-time placebo avg gap       :", formatC(sc$intime_avg, format="f", digits=4),
    " (~0, as it should be)\n")
cat("  permutation p-value          :", f(sc$perm_p),
    " (rank 1 of", sc$n_perm, "; post/pre RMSPE ratio", f(sc$rmspe_ratio), ")\n")
cat("\n  SC donor weights (top):\n")
print(sc$top_weights)
cat("  (true generating anchors were donors", paste(sc$anchors, collapse=","),
    "=", paste(sc$w_true, collapse=","), ";\n",
    "   SC found a DIFFERENT convex combo -> weights are not unique when\n",
    "   donors outnumber factors; pre-fit + counterfactual are what matter.)\n")
cat("\n####################################################################\n")
