D <- readRDS("nova_sensitivity.rds"); dat <- D$data
fit <- lm(log_sales ~ adopt + x, dat)
fit_oracle <- lm(log_sales ~ adopt + x + u_management, dat)
fit_nc <- lm(pre_policy_complaints ~ adopt + x, dat)
b <- coef(summary(fit))["adopt", "Estimate"]
se <- coef(summary(fit))["adopt", "Std. Error"]
tval <- abs(b/se); df <- df.residual(fit)

# Cinelli-Hazlett robustness value for reducing the point estimate to zero.
fq <- tval/sqrt(df)
rv <- 0.5*(sqrt(fq^4 + 4*fq^2) - fq^2)

# Transparent bias contour under bounds on two partial R-squared values.
r2d <- seq(0, 0.30, length.out=81)
r2y <- seq(0, 0.30, length.out=81)
grid <- expand.grid(r2_d_u_x=r2d, r2_y_u_dx=r2y)
grid$bias <- se*sqrt(df)*sqrt(grid$r2_y_u_dx*grid$r2_d_u_x/(1-grid$r2_d_u_x))
grid$adjusted <- b-grid$bias

# A simple identified region indexed by an interpretable maximum bias B.
bias_cap <- seq(0, 0.45, by=0.01)
region <- data.frame(bias_cap=bias_cap, lower=b-bias_cap, upper=b+bias_cap)

out <- list(
  naive=c(est=b, se=se), oracle=coef(summary(fit_oracle))["adopt", c("Estimate","Std. Error")],
  negative_control=coef(summary(fit_nc))["adopt", c("Estimate","Std. Error")],
  rv=rv, grid=grid, region=region, truth=D$truth
)
saveRDS(out, "nova_results.rds")
cat(sprintf("naive %.3f (%.3f) | oracle %.3f (%.3f) | negative control %.3f (%.3f) | RV %.3f\n",
            out$naive[1],out$naive[2],out$oracle[1],out$oracle[2],
            out$negative_control[1],out$negative_control[2],rv))
