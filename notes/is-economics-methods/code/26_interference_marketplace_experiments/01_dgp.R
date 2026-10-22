set.seed(1818)
M <- 120; Tm <- 20; n <- 80
direct <- 0.10; congestion <- -0.06
market_fe <- rnorm(M, sd=0.20); time_fe <- rnorm(Tm, sd=0.08)
market_time_shock <- matrix(rnorm(M*Tm, sd=0.12), M, Tm)
base <- expand.grid(market=1:M, period=1:Tm, unit=1:n)

make_y <- function(d, dat) {
  sat <- ave(d, dat$market, dat$period, FUN=mean)
  2 + market_fe[dat$market] + time_fe[dat$period] + direct*d +
    congestion*sat + market_time_shock[cbind(dat$market, dat$period)] +
    rnorm(nrow(dat), sd=0.45)
}

# Individual randomization within every market-period.
individual <- base
individual$d <- rbinom(nrow(individual),1,0.5)
individual$y <- make_y(individual$d, individual)

# Cluster randomization at the market-period level.
cluster <- base
assign_mt <- matrix(rbinom(M*Tm,1,0.5),M,Tm)
cluster$d <- assign_mt[cbind(cluster$market,cluster$period)]
cluster$y <- make_y(cluster$d,cluster)

# Switchback: each market alternates treatment in randomized two-period blocks.
switchback <- base
sw <- matrix(0L,M,Tm)
for(m in 1:M){
  first <- rbinom(1,1,0.5)
  sw[m,] <- rep(rep(c(first,1-first),each=2),length.out=Tm)
}
switchback$d <- sw[cbind(switchback$market,switchback$period)]
switchback$y <- make_y(switchback$d,switchback)

saveRDS(list(individual=individual,cluster=cluster,switchback=switchback,
             truth=c(direct=direct,overall=direct+congestion,spillover=0.5*congestion)),
        "meridian_marketplace.rds")
cat(sprintf("observations per design = %d | direct %.3f | global rollout %.3f\n",
            nrow(base),direct,direct+congestion))
