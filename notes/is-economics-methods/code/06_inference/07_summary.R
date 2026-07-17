# ============================================================================
# 07_summary.R  --  master reconciliation of every number chapter 6 cites.
# ============================================================================
H <- readRDS("halcyon.rds"); truth <- H$truth
dl <- readRDS("res_delta.rds"); cl <- readRDS("res_cluster.rds")
hc <- readRDS("res_hac.rds"); te <- readRDS("res_test.rds")
L <- function(...) cat(sprintf(...), "\n")
cat("================ CHAPTER 6 NUMBER RECONCILIATION ================\n\n")

cat("--- Delta vs bootstrap (WTP = b_time/b_fee, truth 0.20) ---\n")
L("MAIN: WTP=%.4f  Delta SE=%.4f  boot SE=%.4f  skew=%.2f",
  dl$main$wtp_hat, dl$main$se_delta, dl$main$se_boot, dl$main$skew)
L("  Delta CI [%.3f,%.3f]  perc CI [%.3f,%.3f]  perc-t CI [%.3f,%.3f]",
  dl$main$ci_delta[1], dl$main$ci_delta[2], dl$main$ci_perc[1], dl$main$ci_perc[2],
  dl$main$ci_pt[1], dl$main$ci_pt[2])
L("WEAK: b_fee=%.3f (SE %.3f)  WTP=%.3f  Delta SE=%.3f  boot SE=%.3f  skew=%.1f",
  dl$base_w["bf"], dl$base_w["se_bf"], dl$weak$wtp_hat, dl$weak$se_delta,
  dl$weak$se_boot, dl$weak$skew)
L("  Delta CI [%.2f,%.2f]  perc CI [%.2f,%.2f]  perc-t CI [%.2f,%.2f]",
  dl$weak$ci_delta[1], dl$weak$ci_delta[2], dl$weak$ci_perc[1], dl$weak$ci_perc[2],
  dl$weak$ci_pt[1], dl$weak$ci_pt[2])

cat("\n--- clustering (G=40 rollout, truth tau=0.25) ---\n")
L("tau_hat=%.4f  naive SE=%.4f (t=%.2f)  HC1 SE=%.4f  CR1 SE=%.4f (t=%.2f)  CR2 SE=%.4f",
  cl$tau_hat, cl$se_iid, cl$tau_hat/cl$se_iid, cl$se_hc1, cl$se_cr1,
  cl$tau_hat/cl$se_cr1, cl$se_cr2)
L("SE inflation CR1/naive=%.2f  Moulton predicts %.2f", cl$se_cr1/cl$se_iid, cl$moulton)
L("false-positive rate @5%% G=40: naive=%.3f cluster-z=%.3f cluster-t=%.3f",
  cl$s40["iid"], cl$s40["cr1_z"], cl$s40["cr1_t"])
L("false-positive rate @5%% G=10: naive=%.3f cluster-z=%.3f cluster-t=%.3f",
  cl$s10["iid"], cl$s10["cr1_z"], cl$s10["cr1_t"])
L("few-cluster p (t=%.2f): normal-z=%.4f  t(9)=%.4f  wild-cluster boot=%.4f",
  cl$t10, cl$p10_z, cl$p10_t, cl$p10_wcb)

cat("\n--- HAC (T=400, AR(1) errors) ---\n")
L("beta=%.4f  naive SE=%.4f (t=%.2f)  Newey-West SE=%.4f (t=%.2f, lag %d)  ratio=%.2f",
  hc$b, hc$se_iid, hc$b/hc$se_iid, hc$se_nw, hc$b/hc$se_nw, hc$L, hc$se_nw/hc$se_iid)
L("CI coverage: naive=%.3f  Newey-West=%.3f", hc$cov_iid, hc$cov_nw)

cat("\n--- testing trinity (H0: g1=0, chi^2_1) ---\n")
L("g1_hat=%.4f  W=%.2f (p=%.5f)  LR=%.2f (p=%.5f)  LM=%.2f (p=%.5f)",
  te$g1_hat, te$W, te$pW, te$LR, te$pLR, te$LM, te$pLM)
L("Wald non-invariance (H0: g1=0.5): on g1 W=%.2f (p=%.4f); on log(g1) W=%.2f (p=%.4f)",
  te$W_direct, te$pWd, te$W_repar, te$pWr)
cat("\n================================================================\n")
