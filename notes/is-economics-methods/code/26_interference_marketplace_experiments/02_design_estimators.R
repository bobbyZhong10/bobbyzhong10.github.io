D <- readRDS("meridian_marketplace.rds")

est_individual <- function(dat){
  fit <- lm(y ~ d + factor(market) + factor(period),dat)
  c(est=coef(fit)["d"],se=coef(summary(fit))["d","Std. Error"])
}

# Cluster and switchback outcomes are first aggregated to the 2,400 market-period
# cells. The SE then follows the unit of randomization: in the cluster design each
# cell is randomized, so the classical SE on cells is used; in the switchback only
# the start phase of each market is random (120 draws), so the SE is clustered by
# market (CR1).
est_market_period <- function(dat, cluster=NULL){
  agg <- aggregate(cbind(y,d) ~ market + period, dat, mean)
  fit <- lm(y ~ d + factor(market) + factor(period), agg)
  se <- coef(summary(fit))["d","Std. Error"]
  if(!is.null(cluster)){
    X <- model.matrix(fit); cl <- agg[[cluster]]
    bread <- solve(crossprod(X)); meat <- crossprod(rowsum(X*resid(fit), cl))
    n_cl <- length(unique(cl)); n_obs <- nrow(X); k <- ncol(X)
    se <- sqrt((bread %*% meat %*% bread)["d","d"]*n_cl/(n_cl-1)*(n_obs-1)/(n_obs-k))
  }
  c(est=coef(fit)["d"],se=se)
}
r_ind <- est_individual(D$individual)
r_clu <- est_market_period(D$cluster)
r_sw <- est_market_period(D$switchback, cluster="market")

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
