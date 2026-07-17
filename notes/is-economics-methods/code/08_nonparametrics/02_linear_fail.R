# ============================================================================
# 02_linear_fail.R  --  the "one number that lies": a constant-elasticity
# (log-demand on log-price) fit reports a single elasticity and hides that the
# true price sensitivity swings from near-zero to strongly elastic across the
# price range. Its residuals show obvious leftover structure.
# ============================================================================

S <- readRDS("solstice.rds"); truth <- S$truth
d <- S$demand

## constant-elasticity model: log-demand y on log-price
cel <- lm(y ~ log(x), data = d)
elas_hat <- unname(coef(cel)["log(x)"])            # the single reported elasticity
## linear-in-price model, for the implied local slope
lin <- lm(y ~ x, data = d)
slope_lin <- unname(coef(lin)["x"])

cat("=== the constant-elasticity fit and what it hides ===\n")
cat(sprintf("reported constant elasticity (log-log) = %.3f\n", elas_hat))
cat(sprintf("reported constant slope (linear)       = %.3f\n", slope_lin))
cat("\ntrue LOCAL elasticity x*g'(x) at the three prices:\n")
cat(sprintf("  price 2.0 : %.3f   (nearly inelastic -- room to raise price)\n", truth$elas[1]))
cat(sprintf("  price 5.5 : %.3f   (very elastic -- the sensitive zone)\n", truth$elas[2]))
cat(sprintf("  price 8.0 : %.3f\n", truth$elas[3]))

## the linear fit misfits: R^2 looks fine, but residuals are systematically
## structured (curvature the line cannot follow) -- a RESET-style check
rs <- resid(lin)
reset <- lm(rs ~ poly(fitted(lin), 3))             # residuals on fitted^2,^3
Fp <- anova(lm(rs ~ 1), reset)$`Pr(>F)`[2]
cat(sprintf("\nlinear fit R^2 = %.3f  (looks respectable)\n", summary(lin)$r.squared))
cat(sprintf("but residuals carry curvature: RESET-type F p = %.2e (reject 'linear is enough')\n", Fp))

saveRDS(list(elas_hat = elas_hat, slope_lin = slope_lin,
             r2 = summary(lin)$r.squared, reset_p = Fp), "res_linear.rds")
