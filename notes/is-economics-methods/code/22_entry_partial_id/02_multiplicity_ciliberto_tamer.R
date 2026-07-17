# ============================================================
# Chapter 18 - multiplicity of equilibria and Ciliberto-Tamer partial ID.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(18)
D<-readRDS("ch22_dgp.rds")
dt_m<-D$dt_m; Zm<-D$Z; truth<-D$truth; config<-D$config; Nentry<-D$Nentry
Fpot<-3L; Mkt<-nrow(dt_m)
b0<-truth$b0;bS<-truth$bS;a_Z<-truth$a_Z;delta_true<-truth$delta;sd_xi<-truth$sd_xi;sd_eps<-truth$sd_eps

## all 2^F configurations (rows = subsets of {1,2,3})
subsets<-as.matrix(expand.grid(rep(list(0:1),Fpot)))   # 8 x 3
## is config e (0/1 vector) a pure-strategy Nash eq given standalone profits 'base'?
is_nash<-function(e, base, delta){
  n<-sum(e)
  ok<-TRUE
  for(f in 1:Fpot){
    if(e[f]==1){ if(base[f]-delta*(n-1) < 0){ok<-FALSE;break} }
    else       { if(base[f]-delta*n     >= 0){ok<-FALSE;break} }
  }
  ok
}
nash_set<-function(base, delta){
  which(apply(subsets,1,function(e) is_nash(e,base,delta)))
}

## ---- multiplicity in the realized DGP ----
mult<-integer(Mkt)
for(m in 1:Mkt){
  base<-b0+bS*dt_m$size[m]+a_Z*Zm[m,]+D$xi[m]+D$eps[m,]
  mult[m]<-length(nash_set(base,delta_true))
}
cat("=== multiplicity of equilibria (realized DGP) ===\n")
cat("distribution of # equilibria per market:\n"); print(table(mult))
cat("share of markets with MULTIPLE equilibria:", round(mean(mult>1),3), "\n")
# among 2-entrant outcomes, how often is the IDENTITY indeterminate?
cat("(the NUMBER of entrants is unique in", round(mean(sapply(1:Mkt,function(m){
  base<-b0+bS*dt_m$size[m]+a_Z*Zm[m,]+D$xi[m]+D$eps[m,]
  ns<-nash_set(base,delta_true); length(unique(rowSums(subsets[ns,,drop=FALSE])))==1})),3),
  "of markets even when identity is not)\n")

## ==========================================================
## Ciliberto-Tamer moment inequalities: identified set for delta.
## For each candidate delta (others fixed at truth), integrate over (xi,eps):
##   P(c is unique eq | X_m)  <=  P(observe c)  <=  P(c is a possible eq | X_m)
## Unconditional moments: freq_c in [ mean_m P_unique_c , mean_m P_possible_c ].
## Q(delta) = sum_c [ max(0, meanUnique-freq)^2 + max(0, freq-meanPossible)^2 ].
## ==========================================================
S<-80L                                          # simulation draws per market
XI <-matrix(rnorm(Mkt*S,0,sd_xi),Mkt,S)
EPS<-array(rnorm(Mkt*Fpot*S,0,sd_eps),c(Mkt,Fpot,S))
# code each config as an integer 1..8 (row index in subsets)
obs_code<-apply(config,1,function(e) which(apply(subsets,1,function(s) all(s==e))))
freq<-tabulate(obs_code, nbins=nrow(subsets))/Mkt

ct_Q<-function(delta){
  Puniq<-matrix(0,Mkt,nrow(subsets)); Ppos<-matrix(0,Mkt,nrow(subsets))
  for(m in 1:Mkt){
    cu<-numeric(nrow(subsets)); cp<-numeric(nrow(subsets))
    for(s in 1:S){
      base<-b0+bS*dt_m$size[m]+a_Z*Zm[m,]+XI[m,s]+EPS[m,,s]
      ns<-nash_set(base,delta)
      if(length(ns)==1) cu[ns]<-cu[ns]+1
      cp[ns]<-cp[ns]+1
    }
    Puniq[m,]<-cu/S; Ppos[m,]<-cp/S
  }
  mu<-colMeans(Puniq); mp<-colMeans(Ppos)
  sum(pmax(0,mu-freq)^2 + pmax(0,freq-mp)^2)
}
grid<-seq(0.3,2.2,by=0.1)
Q<-sapply(grid, ct_Q)
cat("\n=== Ciliberto-Tamer identified set for delta ===\n")
cat("delta grid  :",grid,"\n")
cat("Q(delta)*1e3:",round(Q*1e3,2),"\n")
# identified set = region where the moment inequalities are (approximately) satisfied,
# i.e. Q at its floor. Use a small absolute threshold on Q*1e3.
inset<-grid[Q*1e3 <= 0.15]
cat("identified set (Q*1e3 <= 0.15): [", min(inset), ",", max(inset), "]  (truth =",delta_true,")\n")

## ==========================================================
## Opening-puzzle contrast: a point estimate that ASSUMES a selection rule.
## Simulated MLE of delta (others fixed at truth) under two selection rules.
## Rule A = sequential by profitability (the DGP's rule, correct).
## Rule B = fixed brand priority A>B>C (a different, equally arbitrary rule).
## Both "fit"; they give different delta -> a point estimate is rule-dependent.
## ==========================================================
sel_profit<-function(base,delta){                # most profitable enters first
  ord<-order(base,decreasing=TRUE); e<-rep(0L,Fpot); n<-0L
  for(f in ord) if(base[f]-delta*n>=0){e[f]<-1L;n<-n+1L}; e}
sel_brand<-function(base,delta){                 # fixed brand priority 1>2>3
  e<-rep(0L,Fpot); n<-0L
  for(f in 1:Fpot) if(base[f]-delta*n>=0){e[f]<-1L;n<-n+1L}; e}
sim_nll<-function(delta, rule){
  P<-matrix(0,Mkt,nrow(subsets))
  for(m in 1:Mkt){cnt<-numeric(nrow(subsets))
    for(s in 1:S){base<-b0+bS*dt_m$size[m]+a_Z*Zm[m,]+XI[m,s]+EPS[m,,s]
      e<-rule(base,delta); code<-which(apply(subsets,1,function(z) all(z==e))); cnt[code]<-cnt[code]+1}
    P[m,]<-cnt/S}
  -sum(log(pmax(P[cbind(1:Mkt,obs_code)],1e-6)))
}
gridB<-seq(0.6,1.8,by=0.1)
nllA<-sapply(gridB,function(d) sim_nll(d,sel_profit))
nllB<-sapply(gridB,function(d) sim_nll(d,sel_brand))
cat("\n=== point estimate depends on the assumed selection rule ===\n")
cat("rule A (profitability, correct): delta_hat =",gridB[which.min(nllA)],"\n")
cat("rule B (brand priority, wrong) : delta_hat =",gridB[which.min(nllB)],"\n")

saveRDS(list(mult=mult, freq=freq, grid=grid, Q=Q, inset=inset,
             gridB=gridB, nllA=nllA, nllB=nllB,
             dA=gridB[which.min(nllA)], dB=gridB[which.min(nllB)],
             mult_share=mean(mult>1)), "ch22_ct.rds")
cat("saved ch22_ct.rds\n")
