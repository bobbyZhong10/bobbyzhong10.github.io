D <- readRDS("aurora_design.rds")
dat <- D$data

fit_diff <- lm(log_sales ~ early_access, data = dat)
fit_pre <- lm(log_sales ~ early_access + baseline_quality + factor(region), data = dat)
fit_bad <- lm(log_sales ~ early_access + baseline_quality + engagement + factor(region),
              data = dat)

estimate <- c(
  difference_in_means = coef(fit_diff)["early_access"],
  baseline_adjusted = coef(fit_pre)["early_access"],
  post_treatment_control = coef(fit_bad)["early_access"]
)
se <- c(
  difference_in_means = coef(summary(fit_diff))["early_access", "Std. Error"],
  baseline_adjusted = coef(summary(fit_pre))["early_access", "Std. Error"],
  post_treatment_control = coef(summary(fit_bad))["early_access", "Std. Error"]
)

out <- data.frame(specification = names(estimate), estimate = unname(estimate), se = unname(se))
saveRDS(list(table = out, truth = D$truth), "aurora_results.rds")
print(transform(out, estimate = round(estimate, 3), se = round(se, 3)), row.names = FALSE)
