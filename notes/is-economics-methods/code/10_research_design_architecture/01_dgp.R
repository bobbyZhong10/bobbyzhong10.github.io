set.seed(1010)

N <- 12000
seller_id <- seq_len(N)
baseline_quality <- rnorm(N)
region <- sample(1:120, N, replace = TRUE)
high_need <- as.integer(baseline_quality < 0)
u_service <- rnorm(N)

# The platform randomizes early access within region strata.
early_access <- ave(runif(N), region, FUN = function(x) as.integer(x <= median(x)))

# Engagement is measured after assignment. It is both affected by treatment
# and related to latent service capability, so it is not a baseline control.
engagement <- 0.50 * early_access + 0.30 * baseline_quality +
  0.60 * u_service + rnorm(N, sd = 0.70)

direct_effect <- 0.08
engagement_effect <- 0.12
tau_total <- direct_effect + 0.50 * engagement_effect

log_sales <- 2 + 0.40 * baseline_quality + direct_effect * early_access +
  engagement_effect * engagement + 0.50 * u_service + rnorm(N, sd = 0.60)

aurora <- data.frame(
  seller_id, region, baseline_quality, high_need, early_access,
  engagement, log_sales
)

saveRDS(
  list(data = aurora, truth = c(total = tau_total, direct = direct_effect)),
  "aurora_design.rds"
)

cat(sprintf("N = %d | treated share = %.3f\n", N, mean(early_access)))
cat(sprintf("true total effect = %.3f | structural direct effect = %.3f\n",
            tau_total, direct_effect))
