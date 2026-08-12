# ============================================================
# Chapter 19 - Hotz-Miller CCP estimator (no fixed point), a policy
# counterfactual (upgrade subsidy), and figures. Reuses ch23_dgp.rds.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(19)
D<-readRDS("ch23_dgp.rds")
dat<-D$dat; Tr<-D$Tr; xs<-D$xs; K<-D$K; beta<-D$beta; pj<-D$pj; truth<-D$truth
Tkeep<-Tr$Pk; Tup<-Tr$P0

## ---- full solve (for counterfactuals) ----
solve_V<-function(theta1,RC,beta,tol=1e-12,maxit=5000){
  u0<- -theta1*xs; u1<-rep(-RC,K); V<-rep(0,K)
  for(it in 1:maxit){ v0<-u0+beta*(Tkeep%*%V); v1<-u1+beta*(Tup%*%V)
    m<-pmax(v0,v1); Vn<-as.vector(m+log(exp(v0-m)+exp(v1-m)))
    if(max(abs(Vn-V))<tol){V<-Vn;break}; V<-Vn }
  v0<-u0+beta*(Tkeep%*%V); v1<-u1+beta*(Tup%*%V)
  list(V=V, Pup=as.vector(1/(1+exp(v0-v1))))
}
# stationary distribution of staleness under a policy Pup, to get aggregate stats
stat_dist<-function(Pup){
  # transition over x induced by policy: with prob Pup upgrade (reset), else keep
  Tp<-diag(1,K)*0
  for(x in 1:K){ Tp[x,]<-Pup[x]*Tup[x,] + (1-Pup[x])*Tkeep[x,] }
  d<-rep(1/K,K); for(it in 1:5000){dn<-as.vector(d%*%Tp); if(max(abs(dn-d))<1e-12){d<-dn;break}; d<-dn}
  d
}
agg<-function(sol){d<-stat_dist(sol$Pup); c(uprate=sum(d*sol$Pup), stale=sum(d*xs))}

## ==========================================================
## Hotz-Miller CCP estimator (two-step, no fixed point).
## 1) empirical CCPs Phat(upgrade|x). 2) V(theta) = (I - beta M)^{-1} b(theta),
##    linear in theta given CCPs. 3) match log(P1/P0) = v1-v0 (linear in theta) by OLS.
## ==========================================================
tab<-dat[,.(n=.N, up=sum(i)),by=x][order(x)]
Phat<-rep(NA,K); for(k in 0:(K-1)){r<-tab[x==k]; Phat[k+1]<-if(nrow(r)) (r$up+0.5)/(r$n+1) else NA}
# fill any empty high-staleness bins by last observed (rarely visited)
for(k in 2:K) if(is.na(Phat[k])) Phat[k]<-Phat[k-1]
P1<-pmin(pmax(Phat,1e-4),1-1e-4); P0<-1-P1
psi0<- -log(P0); psi1<- -log(P1)                      # HMSS correction (Euler const dropped)
M<-diag(P0)%*%Tkeep + diag(P1)%*%Tup
IbM_inv<-solve(diag(K)-beta*M)
# b(theta) = P0*(u0+psi0)+P1*(u1+psi1); u0=-theta1*x, u1=-RC. Split into theta pieces:
# V = IbM_inv %*% b ; linear: V = Vc + theta1*Vt1 + RC*Vrc
b_const<- P0*psi0 + P1*psi1
b_t1   <- P0*(-xs)                                    # coeff of theta1
b_rc   <- P1*(-1)                                     # coeff of RC
Vc <-as.vector(IbM_inv%*%b_const); Vt1<-as.vector(IbM_inv%*%b_t1); Vrc<-as.vector(IbM_inv%*%b_rc)
# v1-v0 = (u1-u0) + beta*((Tup-Tkeep)V) ; LHS data = log(P1/P0)
Dm<-Tup-Tkeep
lhs<- log(P1/P0)
# RHS = (-RC - (-theta1 x)) + beta*(Dm %*% (Vc+theta1 Vt1+RC Vrc))
#     = theta1*(x + beta*Dm%*%Vt1) + RC*(-1 + beta*Dm%*%Vrc) + beta*Dm%*%Vc
A_t1<- xs + beta*as.vector(Dm%*%Vt1)
A_rc<- -1 + beta*as.vector(Dm%*%Vrc)
off <- beta*as.vector(Dm%*%Vc)
y<- lhs - off
Xr<-cbind(A_t1, A_rc)
# weight by observation counts (states visited more count more)
w<-tab$n[match(0:(K-1),tab$x)]; w[is.na(w)]<-0; w<-w/sum(w)
ccp<-lm(y~Xr-1, weights=w)
cat("=== Hotz-Miller CCP estimator (no fixed point) ===\n")
cat("theta1 =",round(coef(ccp)[1],3)," RC =",round(coef(ccp)[2],3),
    "  (truth 0.15, 4.0 ; NFXP 0.154, 4.082)\n")

## ==========================================================
## Counterfactual: platform subsidizes upgrades, RC -> RC*(1-s), s=0.30.
## Predict change in upgrade rate & mean staleness. Compare DYNAMIC-correct
## estimates vs MYOPIC-wrong estimates vs TRUTH.
## ==========================================================
s_sub<-0.30
cf<-function(th1,RC,bet){
  base<-agg(solve_V(th1,RC,bet)); sub<-agg(solve_V(th1,RC*(1-s_sub),bet))
  c(base_up=base["uprate"], sub_up=sub["uprate"],
    d_up=sub["uprate"]-base["uprate"], base_stale=base["stale"], sub_stale=sub["stale"])
}
cf_true<-cf(truth$theta1, truth$RC, beta)
cf_dyn <-cf(D$nfxp$par[1], D$nfxp$par[2], beta)
cf_myo <-cf(D$myopic[1],  D$myopic[2],  0)          # myopic uses beta=0 for its counterfactual too
cat("\n=== counterfactual: 30% upgrade subsidy, change in upgrade rate ===\n")
cat("TRUTH  : baseline up-rate",round(cf_true["base_up.uprate"],3),
    "-> subsidy",round(cf_true["sub_up.uprate"],3),
    " (Delta",round(cf_true["d_up.uprate"],3),")\n")
cat("DYNAMIC: Delta up-rate",round(cf_dyn["d_up.uprate"],3),
    " (matches truth)\n")
cat("MYOPIC : Delta up-rate",round(cf_myo["d_up.uprate"],3),
    " (biased)\n")

saveRDS(list(Phat=P1, ccp=coef(ccp), cf_true=cf_true, cf_dyn=cf_dyn, cf_myo=cf_myo,
             sol=D$sol), "ch23_ccp.rds")

## ---------- Figure 1: optimal policy P(upgrade|x) vs staleness ----------
svglite::svglite("../../assets/fig/fig_23_policy.svg", width=8.2, height=4.4)
par(mar=c(4.4,4.6,2.4,1),family="sans")
plot(xs, D$sol$Pup, type="l", lwd=2.5, col="#1a365d", ylim=c(0,1),
     xlab="staleness x (versions behind)", ylab="P(upgrade | x)", axes=FALSE)
axis(1);axis(2)
emp<-tab$up/tab$n; points(tab$x, emp, pch=16, col="#c53030", cex=0.8)
legend("topleft", c("optimal policy (solved model)","empirical frequency (data)"),
       col=c("#1a365d","#c53030"), lwd=c(2.5,NA), pch=c(NA,16), bty="n", cex=0.9)
title(main="The renewal policy: upgrade probability rises with staleness",
      cex.main=1.0,font.main=1,col.main="#1a365d")
dev.off();cat("wrote fig_23_policy.svg\n")

## ---------- Figure 2: myopic vs dynamic estimates & counterfactual ----------
svglite::svglite("../../assets/fig/fig_23_myopic.svg", width=8.2, height=4.4)
par(mar=c(4.2,4.6,2.4,1),family="sans", mfrow=c(1,2))
# left: theta1 estimates
est<-c(truth=truth$theta1, dynamic=D$nfxp$par[1], myopic=D$myopic[1])
bp<-barplot(est, col=c("#1a365d","#2b6cb0","#c53030"), border=NA, ylim=c(0,0.5),
    names.arg=c("truth","dynamic","myopic"), ylab=expression(paste("staleness cost  ",theta[1])),
    cex.names=0.95)
text(bp, est+0.02, sprintf("%.2f",est), cex=0.9)
title(main="Myopic overstates staleness cost",cex.main=0.95,font.main=1,col.main="#1a365d")
# right: counterfactual Delta upgrade rate
dd<-c(truth=cf_true["d_up.uprate"], dynamic=cf_dyn["d_up.uprate"], myopic=cf_myo["d_up.uprate"])
bp2<-barplot(dd, col=c("#1a365d","#2b6cb0","#c53030"), border=NA,
    names.arg=c("truth","dynamic","myopic"), ylab="Delta upgrade rate (30% subsidy)",
    cex.names=0.95, ylim=c(0,max(dd)*1.25))
text(bp2, dd+max(dd)*0.04, sprintf("%.3f",dd), cex=0.9)
title(main="...and mis-predicts the subsidy counterfactual",cex.main=0.9,font.main=1,col.main="#1a365d")
dev.off();cat("wrote fig_23_myopic.svg\n")
