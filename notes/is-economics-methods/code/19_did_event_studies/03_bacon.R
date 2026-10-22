# ------------------------------------------------------------------
# 03_bacon.R -- Goodman-Bacon (2021) decomposition of the static
# TWFE estimate into its 2x2 DiD components.
# ------------------------------------------------------------------
library(data.table)
library(bacondecomp)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dt <- readRDS(file.path(code_dir, "kestrel_panel.rds"))

df <- as.data.frame(dt[, .(city_id, month, lgmv, D)])
bd <- bacon(lgmv ~ D, data = df, id_var = "city_id", time_var = "month")

cat("== Goodman-Bacon decomposition: all 2x2 components ==\n")
print(cbind(bd[, c("treated", "untreated", "type")],
            round(bd[, c("weight", "estimate")], 4)))

summ <- as.data.table(bd)[, .(weight = sum(weight),
                              avg_estimate = weighted.mean(estimate, weight)),
                          by = type]
summ[, contribution := weight * avg_estimate]
cat("\n== Summary by comparison type ==\n")
print(summ[, .(type, weight = round(weight, 4),
               avg_estimate = round(avg_estimate, 4),
               contribution = round(contribution, 4))])

cat(sprintf("\nWeighted sum over all 2x2s (= TWFE beta): %.4f\n",
            sum(bd$weight * bd$estimate)))

saveRDS(list(components = as.data.table(bd), summary = summ),
        file.path(code_dir, "res_bacon.rds"))
cat("Saved res_bacon.rds\n")
