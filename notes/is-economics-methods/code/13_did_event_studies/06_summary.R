# ------------------------------------------------------------------
# 06_summary.R -- Collect all estimates against the true simple ATT.
# ------------------------------------------------------------------
library(data.table)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
truth  <- readRDS(file.path(code_dir, "truth.rds"))
r_twfe <- readRDS(file.path(code_dir, "res_twfe.rds"))
r_bac  <- readRDS(file.path(code_dir, "res_bacon.rds"))
r_mod  <- readRDS(file.path(code_dir, "res_modern.rds"))

tab <- data.table(
  estimator = c("True simple ATT (DGP)",
                "Static TWFE (feols)",
                "CS simple ATT (att_gt, not-yet-treated)",
                "SA average ATT (sunab)",
                "BJS imputation (did_imputation)"),
  estimate = c(truth$true_att_simple,
               r_twfe$static$est,
               r_mod$cs_simple$est,
               r_mod$sa_att$est,
               r_mod$imp$est),
  se = c(NA_real_,
         r_twfe$static$se,
         r_mod$cs_simple$se,
         r_mod$sa_att$se,
         r_mod$imp$se)
)
tab[, `:=`(estimate = round(estimate, 4), se = round(se, 4))]
tab[, share_of_truth := round(estimate / truth$true_att_simple, 3)]

cat("== Module 11 summary: estimators vs truth ==\n")
print(tab)

cat("\n== Goodman-Bacon weights by comparison type ==\n")
print(r_bac$summary[, .(type, weight = round(weight, 4),
                        avg_estimate = round(avg_estimate, 4))])

cat("\n== CS group ATTs vs truth by cohort ==\n")
cmp <- merge(r_mod$cs_group, truth$true_att_by_g, by = "g")
print(cmp[, .(g, cs_att = round(att.x, 4), cs_se = round(se, 4),
              true_att = round(att.y, 4))])

fwrite(tab, file.path(code_dir, "summary_table.csv"))
cat("\nSaved summary_table.csv\n")
