# ============================================================================
# 04_hac.R  --  serial correlation in a daily platform time series: naive OLS
# SEs understate, Newey-West HAC corrects. Same lesson as clustering (positive
# correlation makes observations count for less than n independent draws), now
# along time instead of within groups.
# ============================================================================

suppressMessages({ library(sandwich); library(lmtest) })
H <- readRDS("halcyon.rds"); truth <- H$truth
ts <- H$ts_dat

m <- lm(yt ~ xt, data = ts); b <- unname(coef(m)["xt"])
se_iid <- sqrt(vcov(m)["xt","xt"])
# Newey-West with Bartlett kernel; rule-of-thumb lag ~ 0.75 T^(1/3)
L <- floor(0.75 * nrow(ts)^(1/3))
se_nw <- sqrt(NeweyWest(m, lag = L, prewhite = FALSE, adjust = TRUE)["xt","xt"])

cat("=== daily platform series, AR(1) errors (T=400) ===\n")
cat(sprintf("beta_hat = %.4f   (truth %.2f)\n", b, truth$beta_ts))
cat(sprintf("naive OLS  SE = %.4f   t = %5.2f\n", se_iid, b/se_iid))
cat(sprintf("Newey-West SE = %.4f   t = %5.2f   (lag L=%d)\n", se_nw, b/se_nw, L))
cat(sprintf("HAC / naive SE ratio = %.2f\n", se_nw/se_iid))

## ---- Monte-Carlo: CI coverage of the true beta ----------------------------
set.seed(414)
covr <- rowMeans(replicate(2000, {
  xt <- as.numeric(arima.sim(list(ar = 0.5), n = 400))
  u  <- as.numeric(arima.sim(list(ar = 0.6), n = 400, sd = 1))
  yt <- 1.0 + truth$beta_ts * xt + u
  mm <- lm(yt ~ xt)
  bb <- unname(coef(mm)["xt"])
  s_i <- sqrt(vcov(mm)["xt","xt"])
  s_n <- sqrt(NeweyWest(mm, lag = L, prewhite = FALSE, adjust = TRUE)["xt","xt"])
  c(iid = abs(bb - truth$beta_ts) <= 1.96*s_i,
    nw  = abs(bb - truth$beta_ts) <= 1.96*s_n)
}))
cat(sprintf("\n95%% CI coverage over 2000 series: naive OLS = %.3f | Newey-West = %.3f\n",
            covr["iid"], covr["nw"]))

saveRDS(list(b=b, se_iid=se_iid, se_nw=se_nw, L=L, cov_iid=covr["iid"],
             cov_nw=covr["nw"]), "res_hac.rds")
