# ============================================================================
# 08_summary.R  --  master reconciliation table of every number the chapter
# cites, each against its oracle / Monte-Carlo target.
# ============================================================================
oracle <- readRDS("oracle.rds")
mle <- readRDS("res_mle.rds"); qm <- readRDS("res_qmle.rds")
gm  <- readRDS("res_gmm.rds"); wk <- readRDS("res_weak.rds")

line <- function(...) cat(sprintf(...), "\n")
cat("================ CHAPTER 5 NUMBER RECONCILIATION ================\n\n")

cat("--- adoption / logit MLE ---\n")
line("logit MLE  b0 = %.4f  b1 = %.4f   (truth -0.40 / 1.10)", mle$coef[1], mle$coef[2])
line("AME (logit) = %.4f   (oracle %.4f)", mle$ame, oracle$ame_logit)
line("b1 Fisher SE = %.4f   sandwich SE = %.4f  (agree: correct spec)",
     mle$se_fisher[2], mle$se_sand[2])
cat("consistency of b1_hat (mean / sd) across n = 100..102400:\n")
print(round(mle$consistency[c("mean","sd"),], 4))

cat("\n--- usage / Poisson QMLE (overdispersed) ---\n")
line("Poisson QMLE b0 = %.4f  b1 = %.4f   (truth 1.00 / 0.50; consistent)",
     qm$coef[1], qm$coef[2])
line("b1 model/Fisher SE = %.4f   sandwich SE = %.4f   (ratio %.2f)",
     qm$se_model[2], qm$se_sand[2], qm$se_sand[2]/qm$se_model[2])
line("Pearson dispersion phi = %.3f", qm$phi)
line("MC truth SD(b1) = %.4f ; avg Fisher SE = %.4f ; avg sandwich SE = %.4f",
     qm$sd_true, qm$mean_se_mod, qm$mean_se_sand)
line("95%% coverage: naive Fisher = %.3f   sandwich = %.3f", qm$cov_model, qm$cov_sand)

cat("\n--- GMM ---\n")
line("just-ID  h=(1,x)      b1 = %.4f   (= Poisson QMLE)", gm$b_just[2])
line("over-ID  identity W   b1 = %.4f  (SE %.4f)", gm$b_id[2], gm$se_id[2])
line("over-ID  efficient W  b1 = %.4f  (SE %.4f)  eff/id SE = %.3f",
     gm$b_eff[2], gm$se_eff[2], gm$se_eff[2]/gm$se_id[2])
line("Hansen J (correct spec)      = %.4f  p = %.4f", gm$Jstat, gm$Jp)
line("Hansen J (misspecified mean) = %.4f  p = %.6f", gm$Jmis, gm$Jmisp)

cat("\n--- weak identification ---\n")
line("SD(b1) strong = %.4f   weak = %.4f   (%.1fx)",
     wk$strong_sd, wk$weak_sd, wk$weak_sd/wk$strong_sd)
line("criterion curvature strong = %.4f   weak = %.4f   (%.0fx)",
     wk$strong_curv, wk$weak_curv, wk$strong_curv/wk$weak_curv)
line("weak: median b1 = %.4f  avg SE = %.4f  95%% cover = %.3f",
     wk$weak_median, wk$weak_meanse, wk$weak_cov)
cat("\n================================================================\n")
