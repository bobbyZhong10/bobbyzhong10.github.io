# ============================================================
# Chapter 17 - Supply side: recover marginal costs from the Chapter 16 demand,
# simulate the merger of firm 1 and firm 2, compute GUPPI and welfare.
# Done THREE ways (true / BLP-estimated / IV-logit demand) to show that
# getting demand right changes the merger verdict.
# Reuses ../20_demand_blp/ch20_dgp.rds and ch20_blp.rds.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(17)
D <- readRDS("../20_demand_blp/ch20_dgp.rds")
B <- readRDS("../20_demand_blp/ch20_blp.rds")
L <- readRDS("../20_demand_blp/ch20_logit.rds")
dt<-D$dt; J<-D$J; Tm<-D$Tm; truth<-D$truth
R<-D$R; nu_x1<-D$nu_x1; nu_p<-D$nu_p
setkey(dt,t,j)
X1<-matrix(dt$x1,nrow=J); X2<-matrix(dt$x2,nrow=J); Pobs<-matrix(dt$p,nrow=J)
Sobs<-matrix(dt$s,nrow=J); MC<-matrix(dt$mc,nrow=J)

## ---- ownership matrices ----
H0<-diag(J)                                   # pre-merger: 6 single-product firms
H1<-diag(J); H1[1,2]<-1; H1[2,1]<-1           # post-merger: firms 1 & 2 combined

## ---- a demand "model" = list(a = non-price mean utility JxTm, alpha, sig1, sigp) ----
## individual shares in market t given prices p (length J)
mk_shares<-function(a_t, p, alpha, sig1, sigp){
  mu<-outer(sig1*nu_x1, X1_t) + 0    # placeholder; X1_t set by caller via closure
  NULL
}
# to keep X1 per-market available, write market-indexed functions
shares_t<-function(a_t, x1_t, p, alpha, sig1, sigp){
  mu<-outer(sig1*nu_x1, x1_t) + outer(sigp*nu_p, p)     # R x J
  V <- sweep(mu,2,a_t - alpha*p,"+")                     # a_j - alpha p_j + mu
  eV<-exp(V); sij<-eV/(1+rowSums(eV))
  list(s=colMeans(sij), sij=sij)
}
# Delta_jk = -ds_j/dp_k (Convention A, positive diagonal); a_i = alpha - sigp*nu_p
Delta_t<-function(sij, alpha, sigp){
  ai<-alpha - sigp*nu_p                                  # R (price sensitivity)
  Del<-matrix(0,J,J)
  for(j in 1:J)for(k in 1:J)
    Del[j,k]<- if(j==k) mean(ai*sij[,j]*(1-sij[,j])) else -mean(ai*sij[,j]*sij[,k])
  Del
}
# recover mc in market t under ownership H, demand model
mc_t<-function(a_t,x1_t,alpha,sig1,sigp,H,p){
  sh<-shares_t(a_t,x1_t,p,alpha,sig1,sigp)
  Del<-Delta_t(sh$sij,alpha,sigp)
  p - solve(H*Del, sh$s)
}
# Morrow-Skerlos zeta-map: solve post-merger prices given mc, ownership H
solve_merger_t<-function(a_t,x1_t,alpha,sig1,sigp,H,mc,p0){
  p<-p0
  for(it in 1:300){
    sh<-shares_t(a_t,x1_t,p,alpha,sig1,sigp); s<-sh$s; sij<-sh$sij
    ai<-alpha - sigp*nu_p
    Lam<-sapply(1:J,function(j) mean(ai*sij[,j]))        # Lambda_jj (aggregate price margin)
    Gam<-matrix(0,J,J)
    for(j in 1:J)for(k in 1:J) Gam[j,k]<-mean(ai*sij[,j]*sij[,k])
    eta<-p-mc
    eta_new<-(s + (H*Gam)%*%eta)/Lam                     # zeta map: eta = Lam^{-1}[s + (H.Gam)eta]
    p_new<-mc+as.numeric(eta_new)
    if(max(abs(p_new-p))<1e-11){p<-p_new;break}
    p<-0.5*p+0.5*p_new
  }
  p
}
# consumer surplus per consumer in market t (log-sum / alpha_i), Small-Rosen
CS_t<-function(a_t,x1_t,p,alpha,sig1,sigp){
  mu<-outer(sig1*nu_x1,x1_t)+outer(sigp*nu_p,p)
  V<-sweep(mu,2,a_t-alpha*p,"+")
  ai<-alpha - sigp*nu_p
  ls<-log(1+rowSums(exp(V)))
  mean(ls/ai)                                            # avg over consumers, dollars
}

## ---- build the three demand models' non-price utilities a_jt ----
# TRUE: a = beta0 + b1 x1 + b2 x2 + xi
a_true<-matrix(truth$beta0+truth$beta_x1*dt$x1+truth$beta_x2*dt$x2+dt$xi,nrow=J)
# BLP: a = delta_hat + alpha_blp * p_obs
a_blp <- B$delta_hat + B$alpha*Pobs
# IV logit: a = (ln s - ln s0) + alpha_iv * p_obs ; sig=0 (plain logit)
a_log <- (log(Sobs)-log(matrix(dt$s0,nrow=J))) + L$a_iv*Pobs

models<-list(
  true = list(a=a_true, alpha=truth$alpha, sig1=truth$sigma_x1, sigp=truth$sigma_p),
  blp  = list(a=a_blp,  alpha=B$alpha,     sig1=B$sigma_x1,     sigp=B$sigma_p),
  logit= list(a=a_log,  alpha=L$a_iv,      sig1=0,              sigp=0)
)

run_merger<-function(m){
  dP1<-dP2<-dPr<-numeric(Tm); mcm<-numeric(Tm)
  CSpre<-CSpost<-PSpre<-PSpost<-numeric(Tm)
  guppi1<-guppi2<-numeric(Tm)
  for(t in 1:Tm){
    a_t<-m$a[,t]; x1_t<-X1[,t]; p_obs<-Pobs[,t]
    mc<-mc_t(a_t,x1_t,m$alpha,m$sig1,m$sigp,H0,p_obs)    # recover mc (pre-merger conduct)
    mcm[t]<-median(mc)
    p1<-solve_merger_t(a_t,x1_t,m$alpha,m$sig1,m$sigp,H1,mc,p_obs)  # post-merger prices
    dP1[t]<-p1[1]/p_obs[1]-1; dP2[t]<-p1[2]/p_obs[2]-1
    dPr[t]<-mean(p1[3:J]/p_obs[3:J]-1)
    # GUPPI for brand 1 and 2 (pre-merger diversion & margins)
    sh<-shares_t(a_t,x1_t,p_obs,m$alpha,m$sig1,m$sigp); Del<-Delta_t(sh$sij,m$alpha,m$sigp)
    D12<- -Del[2,1]/Del[1,1]; D21<- -Del[1,2]/Del[2,2]  # wait: Del=-ds/dp; D_{1->2}=ds2/dp1 / -ds1/dp1
    D12<- (-Del[2,1])/Del[1,1]; D21<- (-Del[1,2])/Del[2,2]
    guppi1[t]<-D12*(p_obs[2]-mc[2])/p_obs[1]
    guppi2[t]<-D21*(p_obs[1]-mc[1])/p_obs[2]
    # welfare: CS per consumer pre/post; PS = sum (p-mc) s
    CSpre[t]<-CS_t(a_t,x1_t,p_obs,m$alpha,m$sig1,m$sigp)
    CSpost[t]<-CS_t(a_t,x1_t,p1,m$alpha,m$sig1,m$sigp)
    s_post<-shares_t(a_t,x1_t,p1,m$alpha,m$sig1,m$sigp)$s
    PSpre[t]<-sum((p_obs-mc)*sh$s); PSpost[t]<-sum((p1-mc)*s_post)
  }
  list(dP1=median(dP1),dP2=median(dP2),dPr=median(dPr),
       guppi1=median(guppi1),guppi2=median(guppi2),
       dCS=mean(CSpost-CSpre), dPS=mean(PSpost-PSpre),
       CSpre=mean(CSpre))
}

cat("=== Merger of firm 1 & firm 2 (brands 1,2 = closest substitutes) ===\n")
res<-lapply(models,run_merger)
cat(sprintf("%-7s  dP1     dP2     dP_rivals  GUPPI1  GUPPI2   dCS/cons   dPS\n",""))
for(nm in names(res)){r<-res[[nm]]
  cat(sprintf("%-7s  %+.1f%%  %+.1f%%   %+.1f%%     %.3f   %.3f    %+.4f   %+.4f\n",
      nm, 100*r$dP1,100*r$dP2,100*r$dPr,r$guppi1,r$guppi2,r$dCS,r$dPS))}

saveRDS(list(res=res, H0=H0, H1=H1), "ch21_merger.rds")
cat("saved ch21_merger.rds\n")
