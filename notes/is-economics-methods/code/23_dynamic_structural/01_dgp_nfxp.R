# ============================================================
# Chapter 19 - Dynamic structural: single-agent dynamic discrete choice (Rust).
# Running case: a merchant on the "Nimbus" platform decides each period whether
# to UPGRADE their storefront toolkit. Staleness x accumulates (raising a
# friction/lost-sales cost); upgrading pays a fixed cost and resets x to 0.
# This is Rust's optimal-stopping / renewal problem.
# Part 1: solve the model, simulate a panel, estimate by NFXP; a myopic foil.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(19)

## ---- state space and primitives ----
K   <- 20L               # staleness bins x = 0..19
xs  <- 0:(K-1)
beta<- 0.95              # discount factor (fixed; not separately identified, see Magnac-Thesmar)
## flow utility (negative costs): keep = -theta1 * x ; upgrade = -RC
theta1_true <- 0.15      # per-bin friction/lost-sales cost of staleness
RC_true     <- 4.0       # fixed upgrade cost
## staleness increments under KEEP: j in {0,1,2} with probs; UPGRADE resets to 0 then increments
pj <- c(0.30, 0.50, 0.20)

# transition matrices (K x K). Under keep from x: x' = min(x+j, K-1). Under upgrade: from 0.
Tmat <- function(){
  Pk <- matrix(0,K,K); P0 <- matrix(0,K,K)
  for(x in 1:K){ for(j in 0:2){ xn<-min((x-1)+j, K-1)+1; Pk[x,xn]<-Pk[x,xn]+pj[j+1] } }
  for(j in 0:2){ xn<-min(0+j,K-1)+1; P0[,xn]<-P0[,xn]+pj[j+1] }   # upgrade: reset to 0 then increment
  list(Pk=Pk, P0=P0)
}
Tr <- Tmat()

## ---- solve the value function (EV contraction, Type-I EV) ----
solve_V <- function(theta1, RC, beta, tol=1e-12, maxit=5000){
  u_keep <- -theta1*xs           # K
  u_up   <- rep(-RC, K)
  V <- rep(0, K)
  for(it in 1:maxit){
    v0 <- u_keep + beta*(Tr$Pk %*% V)     # keep
    v1 <- u_up   + beta*(Tr$P0 %*% V)     # upgrade (continuation from reset)
    m  <- pmax(v0, v1)
    Vn <- as.vector(m + log(exp(v0-m)+exp(v1-m)))   # log-sum (Euler const dropped)
    if(max(abs(Vn-V))<tol){ V<-Vn; break }
    V <- Vn
  }
  v0 <- u_keep + beta*(Tr$Pk %*% V); v1 <- u_up + beta*(Tr$P0 %*% V)
  Pup <- as.vector(1/(1+exp(v0-v1)))               # P(upgrade | x)
  list(V=V, v0=as.vector(v0), v1=as.vector(v1), Pup=Pup)
}
sol <- solve_V(theta1_true, RC_true, beta)
cat("=== true model solution ===\n")
cat("P(upgrade|x) at x=0,5,10,15,19:", round(sol$Pup[c(1,6,11,16,20)],3), "\n")

## ---- simulate a panel of merchants ----
Nm <- 1000L; Tp <- 30L                    # merchants, periods
sim_one <- function(){
  x<-1L; out<-integer(Tp); xr<-integer(Tp)   # x index in 1..K
  for(t in 1:Tp){
    xr[t]<-x-1L
    up <- runif(1) < sol$Pup[x]
    out[t]<-as.integer(up)
    if(up){ x<- min(0+sample(0:2,1,prob=pj), K-1)+1L }
    else  { x<- min((x-1)+sample(0:2,1,prob=pj), K-1)+1L }
  }
  data.table(x=xr, i=out)
}
dat <- rbindlist(lapply(1:Nm, function(m){d<-sim_one(); d[,mid:=m]; d}))
cat("panel:",nrow(dat),"merchant-periods | upgrade rate:",round(mean(dat$i),3),
    "| mean staleness:",round(mean(dat$x),2),"\n")

## ==========================================================
## NFXP: outer loop MLE over (theta1, RC); inner loop solves the fixed point.
## ==========================================================
nll <- function(par, beta_use){
  th1<-par[1]; RC<-par[2]
  s<-solve_V(th1, RC, beta_use)
  Pup<-pmin(pmax(s$Pup,1e-8),1-1e-8)
  ll<-sum(ifelse(dat$i==1, log(Pup[dat$x+1]), log(1-Pup[dat$x+1])))
  -ll
}
fit <- optim(c(0.1,3), nll, beta_use=beta, method="BFGS",
             control=list(maxit=300,reltol=1e-10), hessian=TRUE)
se <- sqrt(diag(solve(fit$hessian)))
cat("\n--- NFXP (dynamic, beta=0.95) ---\n")
cat("theta1 =",round(fit$par[1],3),"(SE",round(se[1],3),", truth",theta1_true,")\n")
cat("RC     =",round(fit$par[2],3),"(SE",round(se[2],3),", truth",RC_true,")\n")

## ==========================================================
## OPENING PUZZLE: a myopic model (beta=0) ignores that merchants upgrade in
## anticipation of future friction. It fits choices but recovers the wrong
## structural costs (and, sec 2, the wrong counterfactual).
## ==========================================================
fit_my <- optim(c(0.1,3), nll, beta_use=0, method="BFGS",
                control=list(maxit=300,reltol=1e-10))
cat("\n--- Myopic (beta=0) ---\n")
cat("theta1 =",round(fit_my$par[1],3),"  RC =",round(fit_my$par[2],3),
    "  (truth 0.15, 4.0)\n")

saveRDS(list(dat=dat, sol=sol, Tr=Tr, xs=xs, K=K, beta=beta, pj=pj,
             truth=list(theta1=theta1_true, RC=RC_true),
             nfxp=list(par=fit$par, se=se), myopic=fit_my$par,
             solve_V_args=list()), "ch23_dgp.rds")
cat("saved ch23_dgp.rds\n")
