# ------------------------------------------------------------------
# 05_predetermined.R -- block D3: lagged GMV and a predetermined promotion budget with feedback
# (4,000 merchants, 8 months; Y_it = 0.6 Y_{i,t-1} + 0.5 X_it + alpha_i + eps_it, X_it responds to eps_{i,t-1}).
#   FE (within) of Y_it on (Y_{i,t-1}, X_it), CR0 SEs: inconsistent for rho (Nickell) and for beta.
#   Difference GMM (Arellano-Bond, effect = "individual") with the instrument sets that match three
#   assumptions on X (Wooldridge 2010 §11.6.1; s is the lag order of E[Z_{i,t-s} Delta eps_it] = 0):
#     predetermined  lag(y, 2:99) + lag(x, 1:99)   X_{i,t-s}, s >= 1: valid here (the main estimate)
#     endogenous     lag(y, 2:99) + lag(x, 2:99)   s >= 2: valid, fewer instruments
#     exogenous      lag(y, 2:99) + lag(x, 0:99)   adds X_it, which the feedback makes invalid
#   One-step with the robust SE (main); the two-step with the Windmeijer SE alongside.
# Reads nimbus_D3.rds, truth11.rds; writes res_d3.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("dyn_fns.R")
suppressPackageStartupMessages(library(plm))
tr <- readRDS("truth11.rds"); P <- tr$params$D3
d <- readRDS("nimbus_D3.rds"); pd <- pdata.frame(d, index = c("id", "t"))

fe <- fe_ar1(d$y, d$id, d$t, X = cbind(x = d$x))
fe_beta <- list(est = unname(fe$est["x"]), se = unname(fe$se["x"])); fe_beta$ci <- ci95(fe_beta$est, fe_beta$se)
fe_rho <- list(est = unname(fe$est["rho"]), se = unname(fe$se["rho"])); fe_rho$ci <- ci95(fe_rho$est, fe_rho$se)

ab_fit <- function(xlags, model = "onestep") {
  f <- stats::as.formula(sprintf("y ~ lag(y, 1) + x | lag(y, 2:99) + lag(x, %s)", xlags))
  pgmm(f, pd, effect = "individual", model = model, transformation = "d")
}
g_pre <- ab_fit("1:99"); g_pre2 <- ab_fit("1:99", "twosteps"); g_end <- ab_fit("2:99"); g_exo <- ab_fit("0:99")
res <- list(fe = fe_beta, fe_rho = fe_rho,
            ab = gmm_summary(g_pre, "x"), ab_rho = gmm_summary(g_pre, "lag(y, 1)"),
            ab_twostep = gmm_summary(g_pre2, "x"), ab_twostep_rho = gmm_summary(g_pre2, "lag(y, 1)"),
            ab_endog = gmm_summary(g_end, "x"), ab_endog_rho = gmm_summary(g_end, "lag(y, 1)"),
            ab_exog = gmm_summary(g_exo, "x"), ab_exog_rho = gmm_summary(g_exo, "lag(y, 1)"),
            N = P$N, T = P$T)
saveRDS(res, "res_d3.rds")

cat(sprintf("=== D3: %d merchants, %d months; truth rho %.2f, beta %.2f ===\n", P$N, P$T, P$rho, P$beta))
cat(sprintf("FE:                     rho %.4f (SE %.4f)  beta %.4f (SE %.4f) [%.4f, %.4f]\n",
            fe_rho$est, fe_rho$se, fe_beta$est, fe_beta$se, fe_beta$ci[1], fe_beta$ci[2]))
for (k in c("ab", "ab_twostep", "ab_endog", "ab_exog")) {
  b <- res[[k]]; r <- res[[paste0(k, "_rho")]]
  cat(sprintf("%-22s  rho %.4f (SE %.4f)  beta %.4f (SE %.4f) [%.4f, %.4f]; r = %d, Hansen p %.4f, m2 p %.3f\n",
              k, r$est, r$se, b$est, b$se, b$ci[1], b$ci[2], b$r, b$J_p, b$m2_p))
}
cat("Saved res_d3.rds\n")
