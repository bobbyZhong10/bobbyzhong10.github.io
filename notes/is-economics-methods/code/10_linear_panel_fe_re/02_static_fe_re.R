# ============================================================================
# 02_static_fe_re.R  --  static panel: pooled OLS / RE are biased when the
# merchant effect correlates with the regressor; FE (within) is consistent;
# the Hausman test detects the difference.
# ============================================================================

suppressMessages(library(plm))
N <- readRDS("nimbus.rds"); truth <- N$truth
d  <- pdata.frame(N$static, index = c("id", "t"))

pool <- plm(y ~ x, data = d, model = "pooling")
re   <- plm(y ~ x, data = d, model = "random")
fe   <- plm(y ~ x, data = d, model = "within")

haus <- phtest(fe, re)                                  # Hausman FE vs RE

cat("=== static panel: beta_hat under three estimators (truth 0.50) ===\n")
cat(sprintf("pooled OLS beta = %.4f (SE %.4f)\n",
            coef(pool)["x"], sqrt(vcov(pool)["x","x"])))
cat(sprintf("random effects  = %.4f (SE %.4f)\n",
            coef(re)["x"], sqrt(vcov(re)["x","x"])))
cat(sprintf("fixed effects   = %.4f (SE %.4f)\n",
            coef(fe)["x"], sqrt(vcov(fe)["x","x"])))
cat(sprintf("\nHausman chi2 = %.2f, df = %d, p = %.3g\n",
            haus$statistic, haus$parameter, haus$p.value))
cat("(RE strongly rejected: alpha_i is correlated with x, so use FE)\n")

## clustered-by-merchant SE for the FE estimator (panel-robust)
se_cl <- sqrt(vcovHC(fe, method = "arellano", cluster = "group")["x","x"])
cat(sprintf("\nFE beta SE: default %.4f vs cluster-by-merchant %.4f\n",
            sqrt(vcov(fe)["x","x"]), se_cl))

saveRDS(list(pool = coef(pool)["x"], re = coef(re)["x"], fe = coef(fe)["x"],
             se_pool = sqrt(vcov(pool)["x","x"]), se_fe = sqrt(vcov(fe)["x","x"]),
             se_fe_cl = se_cl, haus_stat = unname(haus$statistic),
             haus_p = haus$p.value), "res_static.rds")
