# ============================================================
# Chapter 18 - Entry games & partial identification.
# Running case: "Volt" ride-hailing platform brands deciding whether to enter
# each of many small local markets. Complete-information simultaneous entry
# game with strategic interaction (each entrant lowers rivals' profit) and a
# market-level unobservable xi_m common to all firms in a market.
# Part 1: DGP + a naive probit (the opening puzzle) + Bresnahan-Reiss ordered model.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(18)

## ---- structure ----
Mkt <- 600L          # local markets
Fpot<- 3L            # potential entrant brands (A, B, C)

## ---- true parameters (firm-market profit) ----
# pi_fm = b0 + bS*size_m + a*Z_fm - delta*N_{-f,m} + xi_m + eps_fm ; enter iff pi >= 0
b0    <- -0.3
bS    <-  1.6
a_Z   <-  0.9
delta <-  1.1        # competitive effect: each additional rival cuts profit
sd_xi <-  0.8        # market-level unobservable (common to firms in a market)
sd_eps<-  0.8

## ---- covariates ----
dt_m <- data.table(m = 1:Mkt, size = rnorm(Mkt, 0, 1))
Z   <- matrix(rnorm(Mkt*Fpot, 0, 1), Mkt, Fpot)
xi  <- rnorm(Mkt, 0, sd_xi)
eps <- matrix(rnorm(Mkt*Fpot, 0, sd_eps), Mkt, Fpot)

## ---- equilibrium of the entry game per market (sequential-profitability selection) ----
solve_entry <- function(size, Zrow, epsrow, xim){
  base <- b0 + bS*size + a_Z*Zrow + xim + epsrow    # profit if alone
  ord  <- order(base, decreasing = TRUE)
  entered <- rep(FALSE, Fpot); nin <- 0L
  for (f in ord){
    if (base[f] - delta*nin >= 0){ entered[f] <- TRUE; nin <- nin + 1L }
  }
  entered
}
config <- t(sapply(1:Mkt, function(m) solve_entry(dt_m$size[m], Z[m,], eps[m,], xi[m])))
Nentry <- rowSums(config)
dt_m[, N := Nentry]

cat("=== DGP: entry configuration ===\n")
print(table(Nentry))
cat("share markets with >=1 entrant:", round(mean(Nentry>=1),3),
    " | entry rate A/B/C:", round(colMeans(config),3), "\n")
cat("corr(size, N):", round(cor(dt_m$size, Nentry),3), "\n")

## build firm-market long data for the naive probit
long <- data.table(m=rep(1:Mkt,each=Fpot), f=rep(1:Fpot,Mkt),
                   enter=as.integer(t(config)), Z=as.vector(t(Z)),
                   size=rep(dt_m$size,each=Fpot))
long[, Nriv := rep(Nentry,each=Fpot) - enter]     # number of OTHER entrants

## ==========================================================
## OPENING PUZZLE: naive probit of entry on size, own Z, and number of rivals
## treating the number of rivals as exogenous. N_{-f} is endogenous (rivals enter
## more where xi is high => positively correlated with own profit shock), so the
## estimated competitive effect is biased toward zero / wrong sign.
## ==========================================================
naive <- glm(enter ~ size + Z + Nriv, data=long, family=binomial(link="probit"))
cat("\n--- Naive probit (treats rivals' entry as exogenous) ---\n")
cat("coef on Nriv (should be NEGATIVE = competition hurts):",
    round(coef(naive)["Nriv"],3), "  (true delta = -", delta, ")\n", sep="")
cat("  size coef:", round(coef(naive)["size"],3),
    "  Z coef:", round(coef(naive)["Z"],3), "\n")

## ==========================================================
## Bresnahan-Reiss ordered probit for the NUMBER of entrants (unique despite
## multiplicity). Recovers market-size effect and the entry thresholds; widening
## per-firm thresholds => competition toughens with each entrant.
## ==========================================================
ord_probit_ll <- function(par, y, X){
  K <- max(y); beta <- par[1:ncol(X)]
  cuts <- c(-Inf, cumsum(c(par[ncol(X)+1], exp(par[(ncol(X)+2):(ncol(X)+K)]))), Inf)
  xb <- as.vector(X %*% beta); ll <- 0
  for (n in 0:K){ idx<-which(y==n)
    p <- pnorm(cuts[n+2]-xb[idx]) - pnorm(cuts[n+1]-xb[idx]); ll <- ll+sum(log(pmax(p,1e-12))) }
  -ll
}
X <- cbind(size = dt_m$size); K <- max(Nentry)
fit_br <- optim(c(0.8,-0.3,rep(log(1),K-1)), ord_probit_ll, y=Nentry, X=X,
                method="BFGS", control=list(maxit=800,reltol=1e-10), hessian=TRUE)
p <- fit_br$par; cuts <- cumsum(c(p[2], exp(p[3:(1+K)])))
gS <- p[1]
# per-firm entry thresholds in size units: s_N = (cut_N / gamma_S) / N
sthr <- (cuts/gS)/(1:K)
cat("\n--- Bresnahan-Reiss ordered probit (N on market size) ---\n")
cat("size coefficient gamma_S:", round(gS,3), "\n")
cat("cut points (size index):", round(cuts,3), "\n")
cat("per-firm entry thresholds s_N:", round(sthr,3), "\n")
cat("threshold ratios s_2/s_1:", round(sthr[2]/sthr[1],3),
    " s_3/s_2:", round(sthr[3]/sthr[2],3), " (>1 => competition toughens)\n")

saveRDS(list(dt_m=dt_m, Z=Z, eps=eps, xi=xi, config=config, Nentry=Nentry, long=long,
             truth=list(b0=b0,bS=bS,a_Z=a_Z,delta=delta,sd_xi=sd_xi,sd_eps=sd_eps),
             naive=coef(naive), br=list(gS=gS,cuts=cuts,sthr=sthr)),
        "ch22_dgp.rds")
cat("saved ch22_dgp.rds\n")
