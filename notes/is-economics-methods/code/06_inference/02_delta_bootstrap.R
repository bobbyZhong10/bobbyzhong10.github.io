# ============================================================================
# 02_delta_bootstrap.R  --  WTP = b_time/b_fee: Delta method vs the bootstrap.
# Main sample (fee varies a lot -> strong denominator): Delta and bootstrap
# agree. Weak-denominator variant (fee barely varies): the ratio's sampling
# distribution is skewed/heavy-tailed, the symmetric Delta CI misleads, and the
# percentile-t bootstrap is more honest.
# ============================================================================

H <- readRDS("halcyon.rds")
truth <- H$truth

## WTP point estimate and Delta-method SE from one fitted model
wtp_delta <- function(dat) {
  m  <- lm(sat ~ time + fee, data = dat)
  b  <- coef(m); V <- vcov(m)
  bt <- b["time"]; bf <- b["fee"]
  wtp <- bt / bf
  grad <- c(0, 1/bf, -bt/bf^2)                 # d(wtp)/d(b0,b_time,b_fee)
  se  <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  c(wtp = as.numeric(wtp), se = se, bt = as.numeric(bt), bf = as.numeric(bf),
    se_bf = sqrt(V["fee","fee"]))
}

## nonparametric pairs bootstrap: SE, percentile CI, percentile-t CI
boot_wtp <- function(dat, B = 4000, seed = 1) {
  set.seed(seed); n <- nrow(dat)
  base <- wtp_delta(dat); wtp_hat <- base["wtp"]; se_hat <- base["se"]
  out <- t(sapply(1:B, function(b){
    d <- dat[sample(n, n, replace = TRUE), ]
    wd <- wtp_delta(d); c(wd["wtp"], wd["se"])
  }))
  wstar <- out[,1]; sstar <- out[,2]
  se_boot <- sd(wstar)
  ci_perc <- quantile(wstar, c(.025, .975))
  tstar   <- (wstar - wtp_hat) / sstar          # studentized
  ci_pt   <- wtp_hat - quantile(tstar, c(.975, .025)) * se_hat
  list(wtp_hat = as.numeric(wtp_hat), se_delta = as.numeric(se_hat),
       se_boot = se_boot, ci_delta = wtp_hat + c(-1,1)*1.96*se_hat,
       ci_perc = ci_perc, ci_pt = ci_pt, wstar = wstar,
       skew = mean((wstar-mean(wstar))^3)/sd(wstar)^3)
}

cat("=== WTP = b_time / b_fee  (truth =", truth$WTP, ") ===\n\n")

cat("--- MAIN sample (fee varies a lot: strong denominator) ---\n")
bm <- boot_wtp(H$wtp_main, seed = 101)
cat(sprintf("point WTP        = %.4f\n", bm$wtp_hat))
cat(sprintf("Delta SE         = %.4f\n", bm$se_delta))
cat(sprintf("bootstrap SE     = %.4f\n", bm$se_boot))
cat(sprintf("Delta 95%% CI     = [%.4f, %.4f]\n", bm$ci_delta[1], bm$ci_delta[2]))
cat(sprintf("percentile 95%% CI= [%.4f, %.4f]\n", bm$ci_perc[1], bm$ci_perc[2]))
cat(sprintf("percentile-t 95%% CI=[%.4f, %.4f]\n", bm$ci_pt[1], bm$ci_pt[2]))
cat(sprintf("bootstrap skewness = %.3f\n", bm$skew))

cat("\n--- WEAK-denominator variant (fee barely varies: b_fee imprecise) ---\n")
bw <- boot_wtp(H$wtp_weak, seed = 202)
base_w <- wtp_delta(H$wtp_weak)
cat(sprintf("b_fee estimate   = %.4f (SE %.4f)  <- imprecise denominator\n",
            base_w["bf"], base_w["se_bf"]))
cat(sprintf("point WTP        = %.4f\n", bw$wtp_hat))
cat(sprintf("Delta SE         = %.4f\n", bw$se_delta))
cat(sprintf("bootstrap SE     = %.4f\n", bw$se_boot))
cat(sprintf("Delta 95%% CI     = [%.4f, %.4f]  (symmetric, misleading)\n",
            bw$ci_delta[1], bw$ci_delta[2]))
cat(sprintf("percentile 95%% CI= [%.4f, %.4f]\n", bw$ci_perc[1], bw$ci_perc[2]))
cat(sprintf("percentile-t 95%% CI=[%.4f, %.4f]  (asymmetric, honest)\n",
            bw$ci_pt[1], bw$ci_pt[2]))
cat(sprintf("bootstrap skewness = %.3f\n", bw$skew))

saveRDS(list(main = bm, weak = bw, base_w = base_w), "res_delta.rds")
