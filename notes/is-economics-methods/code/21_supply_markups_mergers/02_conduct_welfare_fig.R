# ============================================================
# Chapter 17 - testing conduct (Bertrand vs collusion), the efficiency offset,
# and figures. Reuses Chapter 16 demand (BLP) and DGP.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(17)
D <- readRDS("../20_demand_blp/ch20_dgp.rds")
B <- readRDS("../20_demand_blp/ch20_blp.rds")
M <- readRDS("ch21_merger.rds")
dt<-D$dt; J<-D$J; Tm<-D$Tm; truth<-D$truth; R<-D$R; nu_x1<-D$nu_x1; nu_p<-D$nu_p
setkey(dt,t,j)
X1<-matrix(dt$x1,nrow=J); Pobs<-matrix(dt$p,nrow=J); Sobs<-matrix(dt$s,nrow=J)

## demand helpers (use estimated BLP demand: the analyst's object)
alpha<-B$alpha; sig1<-B$sigma_x1; sigp<-B$sigma_p
a_blp<-B$delta_hat + alpha*Pobs
shares_t<-function(a_t,x1_t,p){mu<-outer(sig1*nu_x1,x1_t)+outer(sigp*nu_p,p)
  V<-sweep(mu,2,a_t-alpha*p,"+");eV<-exp(V);sij<-eV/(1+rowSums(eV));list(s=colMeans(sij),sij=sij)}
Delta_t<-function(sij){ai<-alpha-sigp*nu_p;Del<-matrix(0,J,J)
  for(j in 1:J)for(k in 1:J)Del[j,k]<-if(j==k)mean(ai*sij[,j]*(1-sij[,j]))else -mean(ai*sij[,j]*sij[,k]);Del}

## ---- recover mc under two conduct assumptions ----
H_bert<-diag(J); H_coll<-matrix(1,J,J)     # Bertrand vs perfect collusion (joint monopoly)
mc_bert<-mc_coll<-matrix(0,J,Tm)
for(t in 1:Tm){a_t<-a_blp[,t];x1_t<-X1[,t];p<-Pobs[,t]
  sh<-shares_t(a_t,x1_t,p);Del<-Delta_t(sh$sij)
  mc_bert[,t]<-p-solve(H_bert*Del,sh$s)
  mc_coll[,t]<-p-solve(H_coll*Del,sh$s)}
cat("=== recovered marginal cost under two conduct assumptions ===\n")
cat("Bertrand : median",round(median(mc_bert),2)," share<=0",round(mean(mc_bert<=0),3),"\n")
cat("Collusion: median",round(median(mc_coll),2)," share<=0",round(mean(mc_coll<=0),3),"\n")
cat("TRUTH mc : median",round(median(dt$mc),2),"\n")

## conduct test: which recovered mc is better explained by a cost model?
## supply moment: mc = g0 + gw*w + gx1*x1 + omega ; instruments z = [w, x1, rival w]
dtm<-copy(dt); dtm[,mc_b:=as.vector(mc_bert)]; dtm[,mc_c:=as.vector(mc_coll)]
dtm[,sumw_riv:=sum(w)-w,by=t]
Z<-cbind(1,dtm$w,dtm$x1,dtm$sumw_riv)
Xs<-cbind(1,dtm$w,dtm$x1); PZ<-Z%*%solve(crossprod(Z))%*%t(Z)
gmmQ<-function(y){b<-solve(t(Xs)%*%PZ%*%Xs,t(Xs)%*%PZ%*%y);om<-y-Xs%*%b
  g<-t(Z)%*%om; as.numeric(t(g)%*%solve(crossprod(Z))%*%g)/nrow(Z)}
Qb<-gmmQ(dtm$mc_b); Qc<-gmmQ(dtm$mc_c)
cat("\n=== conduct test (supply-moment GMM objective; lower = better fit) ===\n")
cat("Bertrand  Q =",format(Qb,digits=4),"\n")
cat("Collusion Q =",format(Qc,digits=4),"\n")
cat("preferred conduct:", ifelse(Qb<Qc,"Bertrand (correct)","Collusion"),
    " ; ratio Qc/Qb =",round(Qc/Qb,1),"\n")

## ---- efficiency offset: mc reduction on merged brands that neutralizes brand-2 price rise ----
H1<-diag(J);H1[1,2]<-1;H1[2,1]<-1
solve_merger_t<-function(a_t,x1_t,H,mc,p0){p<-p0
  for(it in 1:300){sh<-shares_t(a_t,x1_t,p);s<-sh$s;sij<-sh$sij;ai<-alpha-sigp*nu_p
    Lam<-sapply(1:J,function(j)mean(ai*sij[,j]));Gam<-matrix(0,J,J)
    for(j in 1:J)for(k in 1:J)Gam[j,k]<-mean(ai*sij[,j]*sij[,k])
    eta<-p-mc;etan<-(s+(H*Gam)%*%eta)/Lam;pn<-mc+as.numeric(etan)
    if(max(abs(pn-p))<1e-11){p<-pn;break};p<-0.5*p+0.5*pn};p}
# find efficiency e in {0..0.3} making median brand-2 post price ~ pre
effs<-seq(0,0.30,by=0.02); dP2e<-numeric(length(effs))
for(ei in seq_along(effs)){e<-effs[ei];d2<-numeric(Tm)
  for(t in 1:Tm){a_t<-a_blp[,t];x1_t<-X1[,t];p<-Pobs[,t]
    sh<-shares_t(a_t,x1_t,p);mc<-p-solve(diag(J)*Delta_t(sh$sij),sh$s)
    mc[1:2]<-mc[1:2]*(1-e)
    p1<-solve_merger_t(a_t,x1_t,H1,mc,p);d2[t]<-p1[2]/p[2]-1}
  dP2e[ei]<-median(d2)}
e_star<-approx(dP2e,effs,xout=0)$y
cat("\n=== efficiency offset ===\n")
cat("mc reduction on merged brands needed to neutralize brand-2 price rise: e* =",
    round(100*e_star,1),"%\n")

saveRDS(list(mc_bert=mc_bert,mc_coll=mc_coll,Qb=Qb,Qc=Qc,e_star=e_star,
             effs=effs,dP2e=dP2e), "ch21_conduct.rds")

## ---------- Figure 1: merger price effects, true vs BLP vs logit ----------
res<-M$res
mat<-rbind(`brand 1`=c(res$true$dP1,res$blp$dP1,res$logit$dP1),
           `brand 2`=c(res$true$dP2,res$blp$dP2,res$logit$dP2),
           rivals   =c(res$true$dPr,res$blp$dPr,res$logit$dPr))*100
svglite::svglite("../../assets/fig/fig_21_merger.svg", width=8.2, height=4.4)
par(mar=c(4,4.4,2.6,1),family="sans")
cols<-c("#1a365d","#2b6cb0","#c53030")
bp<-barplot(t(mat),beside=TRUE,col=cols,border=NA,ylim=c(0,max(mat)*1.2),
   ylab="post-merger price change (%)",names.arg=rownames(mat),cex.names=1.0)
legend("topright",c("true demand","BLP estimate","IV-logit"),fill=cols,border=NA,bty="n",cex=0.95)
title(main="Merger price effects: IV-logit under-predicts; BLP tracks the truth",
      cex.main=0.98,font.main=1,col.main="#1a365d")
dev.off();cat("wrote fig_21_merger.svg\n")

## ---------- Figure 2: conduct test - recovered mc vs cost shifter w ----------
svglite::svglite("../../assets/fig/fig_21_conduct.svg", width=8.2, height=4.4)
par(mar=c(4.2,4.4,2.6,1),family="sans")
ss<-sample(1:nrow(dtm),600)
plot(dtm$w[ss],dtm$mc_b[ss],pch=16,col="#2b6cb066",cex=0.6,
     xlab="cost shifter w",ylab="recovered marginal cost",ylim=range(c(dtm$mc_b,dtm$mc_c)))
points(dtm$w[ss],dtm$mc_c[ss],pch=16,col="#c5303066",cex=0.6)
abline(lm(mc_b~w,dtm),col="#2b6cb0",lwd=2.5); abline(lm(mc_c~w,dtm),col="#c53030",lwd=2.5)
abline(h=0,lty=3,col="#718096")
legend("topleft",c("Bertrand (correct conduct)","perfect collusion (wrong)"),
       col=c("#2b6cb0","#c53030"),pch=16,bty="n",cex=0.9)
title(main="Testing conduct: Bertrand mc tracks the cost shifter; collusion mc is distorted (many below zero)",
      cex.main=0.9,font.main=1,col.main="#1a365d")
dev.off();cat("wrote fig_21_conduct.svg\n")
