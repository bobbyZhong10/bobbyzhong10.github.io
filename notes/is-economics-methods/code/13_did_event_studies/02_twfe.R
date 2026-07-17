# ------------------------------------------------------------------
# 02_twfe.R -- Static TWFE and the naive dynamic TWFE event study
# on the Kestrel Marketplace panel. Clustered by city.
# ------------------------------------------------------------------
library(data.table)
library(fixest)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dt    <- readRDS(file.path(code_dir, "kestrel_panel.rds"))
truth <- readRDS(file.path(code_dir, "truth.rds"))

## ---- Static TWFE: Y_it = alpha_i + lambda_t + beta D_it + eps ----
m_static <- feols(lgmv ~ D | city_id + month, data = dt, cluster = ~city_id)
cat("== Static TWFE ==\n")
print(coeftable(m_static))
cat(sprintf("\nStatic TWFE beta = %.4f (SE %.4f) vs true simple ATT = %.4f\n",
            coef(m_static)[["D"]], se(m_static)[["D"]], truth$true_att_simple))

## ---- Naive dynamic TWFE event study ------------------------------
## Event-time dummies relative to first treatment, e = -1 omitted.
## Never-treated cities get event time -1 so that all their dummies
## are zero (they sit in the reference category). Endpoints binned at
## e <= -12 and e >= 18.
dt[, ev_naive := fifelse(g > 0, pmin(pmax(month - g, -12L), 18L), -1L)]
m_naive <- feols(lgmv ~ i(ev_naive, ref = -1) | city_id + month,
                 data = dt, cluster = ~city_id)
cat("\n== Naive dynamic TWFE event study (ref e = -1) ==\n")
print(coeftable(m_naive))

## Tidy coefficients for later figures/tables
ct <- as.data.table(coeftable(m_naive), keep.rownames = "term")
ct[, e := as.integer(gsub(".*::(-?\\d+)$", "\\1", term))]
naive_es <- ct[, .(e, estimate = Estimate, se = `Std. Error`)][order(e)]

saveRDS(list(static = list(est = coef(m_static)[["D"]], se = se(m_static)[["D"]]),
             naive_es = naive_es),
        file.path(code_dir, "res_twfe.rds"))

cat("\nNaive event-study coefficients at selected horizons:\n")
print(naive_es[e %in% c(-6, -3, 0, 3, 6, 9), .(e, estimate = round(estimate, 4), se = round(se, 4))])
cat("Saved res_twfe.rds\n")
