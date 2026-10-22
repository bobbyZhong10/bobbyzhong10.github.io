# ============================================================================
# 02_experiment.R  —  Study A (randomized pilot)
#   difference-in-means ATE, analytic (Neyman) SE, and a permutation p-value.
# ============================================================================
d   <- readRDS("northwind.rds")
rct <- d$rct

y1 <- rct$Y[rct$D == 1]; y0 <- rct$Y[rct$D == 0]
n1 <- length(y1);        n0 <- length(y0)

dim_est <- mean(y1) - mean(y0)
se_neyman <- sqrt(var(y1)/n1 + var(y0)/n0)          # Neyman conservative SE
t_stat <- dim_est / se_neyman
ci <- dim_est + c(-1, 1) * qnorm(0.975) * se_neyman

# --- randomization (permutation) inference: sharp null of no effect ----------
set.seed(101)
B <- 5000L
perm <- replicate(B, {
  Dp <- sample(rct$D)
  mean(rct$Y[Dp == 1]) - mean(rct$Y[Dp == 0])
})
p_perm <- mean(abs(perm) >= abs(dim_est))

cat("==== Study A : randomized pilot ====\n")
cat(sprintf("n treated / control : %d / %d\n", n1, n0))
cat(sprintf("diff-in-means (ATE) : %.4f\n", dim_est))
cat(sprintf("Neyman SE           : %.4f\n", se_neyman))
cat(sprintf("t / normal p        : %.3f / %.3g\n", t_stat, 2*pnorm(-abs(t_stat))))
cat(sprintf("95%% CI              : [%.4f, %.4f]\n", ci[1], ci[2]))
cat(sprintf("permutation p (B=%d): %.4f\n", B, p_perm))
cat(sprintf("true ATE (RCT)      : %.4f\n", d$truths$true_ATE_rct))

saveRDS(list(dim_est=dim_est, se=se_neyman, ci=ci,
             p_perm=p_perm, n1=n1, n0=n0),
        "res_experiment.rds")
