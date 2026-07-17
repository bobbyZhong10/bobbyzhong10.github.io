D <- readRDS("meridian_marketplace.rds")

est_individual <- function(dat){
  fit <- lm(y ~ d + factor(market) + factor(period),dat)
  c(est=coef(fit)["d"],se=coef(summary(fit))["d","Std. Error"])
}

# Cluster and switchback assignments occur at market-period level. Aggregate
# first so inference uses 2,400 randomized cells, not 192,000 order rows.
est_market_period <- function(dat){
  agg <- aggregate(cbind(y,d) ~ market + period, dat, mean)
  fit <- lm(y ~ d + factor(market) + factor(period), agg)
  c(est=coef(fit)["d"],se=coef(summary(fit))["d","Std. Error"])
}
r_ind <- est_individual(D$individual)
r_clu <- est_market_period(D$cluster)
r_sw <- est_market_period(D$switchback)

# Two-stage saturation estimands from the known response surface.
two_stage <- c(
  direct_at_half=D$truth["direct"],
  spillover_on_control=D$truth["spillover"],
  total_all_vs_none=D$truth["overall"]
)
out <- data.frame(design=c("individual","market-period cluster","switchback"),
                  estimate=c(r_ind[1],r_clu[1],r_sw[1]),se=c(r_ind[2],r_clu[2],r_sw[2]))
saveRDS(list(table=out,two_stage=two_stage,truth=D$truth),"meridian_results.rds")
print(transform(out,estimate=round(estimate,3),se=round(se,3)),row.names=FALSE)
print(round(two_stage,3))
