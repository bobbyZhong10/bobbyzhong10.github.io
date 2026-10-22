dat <- readRDS("lyra_labels.rds")
X <- model.matrix(~ treatment + x, dat)
coef_ols <- function(y, X) drop(solve(crossprod(X), crossprod(X, y)))

b_oracle <- coef_ols(dat$true_satisfaction, X)
b_naive <- coef_ols(dat$llm_satisfaction, X)
g <- dat$gold == 1
b_gold <- coef_ols(dat$true_satisfaction[g], X[g, , drop = FALSE])

# PPI rectification: prediction estimate on the large sample plus the
# gold-sample mean discrepancy in the same estimating equation.
b_pred_gold <- coef_ols(dat$llm_satisfaction[g], X[g, , drop = FALSE])
b_ppi <- b_naive + (b_gold - b_pred_gold)

set.seed(1112)
B <- 300
boot <- matrix(NA_real_, B, ncol(X))
for (b in seq_len(B)) {
  iu <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
  ig <- sample(which(g), sum(g), replace = TRUE)
  bu <- coef_ols(dat$llm_satisfaction[iu], X[iu, , drop = FALSE])
  bg <- coef_ols(dat$true_satisfaction[ig], X[ig, , drop = FALSE])
  bpg <- coef_ols(dat$llm_satisfaction[ig], X[ig, , drop = FALSE])
  boot[b, ] <- bu + bg - bpg
}
se_ppi <- apply(boot, 2, sd)

out <- data.frame(
  estimator = c("oracle", "LLM label naive", "gold only", "PPI rectified"),
  treatment = c(b_oracle[2], b_naive[2], b_gold[2], b_ppi[2]),
  se = c(summary(lm(true_satisfaction ~ treatment + x, dat))$coef[2, 2],
         summary(lm(llm_satisfaction ~ treatment + x, dat))$coef[2, 2],
         summary(lm(true_satisfaction ~ treatment + x, dat, subset = g))$coef[2, 2],
         se_ppi[2])
)
saveRDS(out, "lyra_results.rds")
print(transform(out, treatment = round(treatment, 3), se = round(se, 3)), row.names = FALSE)
