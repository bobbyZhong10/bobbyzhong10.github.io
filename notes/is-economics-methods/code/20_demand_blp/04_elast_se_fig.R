# ============================================================
# Chapter 16 - elasticity matrices (logit IIA vs BLP vs truth), diversion,
# GMM standard errors, and figures.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(16)
D<-readRDS("ch20_dgp.rds"); L<-readRDS("ch20_logit.rds"); B<-readRDS("ch20_blp.rds")
dt<-D$dt; J<-D$J; Tm<-D$Tm; truth<-D$truth; R<-D$R; nu_x1<-D$nu_x1; nu_p<-D$nu_p
setkey(dt,t,j)
X1<-matrix(dt$x1,nrow=J); X2<-matrix(dt$x2,nrow=J); P<-matrix(dt$p,nrow=J); Sobs<-matrix(dt$s,nrow=J)

## ---- average elasticity matrix across markets ----
# BLP / truth: random coefficients; logit: homogeneous with alpha
elast_rc<-function(delta,sig1,sigp,alpha){
  E<-matrix(0,J,J)
  for(t in 1:Tm){
    mu<-outer(sig1*nu_x1,X1[,t])+outer(sigp*nu_p,P[,t])
    V<-sweep(mu,2,delta[,t],"+");eV<-exp(V);sij<-eV/(1+rowSums(eV));s<-colMeans(sij)
    ai<-alpha-sigp*nu_p
    for(j in 1:J)for(k in 1:J){
      d<-if(j==k) -mean(ai*sij[,j]*(1-sij[,j])) else mean(ai*sij[,j]*sij[,k])
      E[j,k]<-E[j,k]+ d*(P[k,t]/s[j])
    }
  }
  E/Tm
}
elast_logit<-function(alpha){    # homogeneous logit, IIA
  E<-matrix(0,J,J)
  for(t in 1:Tm){ s<-Sobs[,t]
    for(j in 1:J)for(k in 1:J){
      d<-if(j==k) -alpha*P[j,t]*(1-s[j]) else alpha*P[k,t]*s[k]
      E[j,k]<-E[j,k]+d } }
  E/Tm
}
delta_true<-matrix(truth$beta0+truth$beta_x1*dt$x1+truth$beta_x2*dt$x2-truth$alpha*dt$p+dt$xi,nrow=J)
E_true<-elast_rc(delta_true,truth$sigma_x1,truth$sigma_p,truth$alpha)
E_blp <-elast_rc(B$delta_hat,B$sigma_x1,B$sigma_p,B$alpha)
E_log <-elast_logit(L$a_iv)

## diversion from product 1: D_{1->k} = -(ds_k/dp_1)/(ds_1/dp_1) = -E_{k1}*s_k/(E_{11}*s_1)...
## simpler: use aggregate ds/dp. Recompute ds matrices averaged.
dmat<-function(delta,sig1,sigp,alpha){
  Dsum<-matrix(0,J,J);ssum<-numeric(J)
  for(t in 1:Tm){mu<-outer(sig1*nu_x1,X1[,t])+outer(sigp*nu_p,P[,t])
    V<-sweep(mu,2,delta[,t],"+");eV<-exp(V);sij<-eV/(1+rowSums(eV));s<-colMeans(sij);ai<-alpha-sigp*nu_p
    Dd<-matrix(0,J,J)
    for(j in 1:J)for(k in 1:J) Dd[j,k]<-if(j==k) mean(ai*sij[,j]*(1-sij[,j])) else -mean(ai*sij[,j]*sij[,k])
    Dsum<-Dsum+Dd;ssum<-ssum+s}
  list(D=Dsum/Tm,s=ssum/Tm)   # D[j,k]=-ds_j/dp_k arranged; ds_1/dp_1=-D[1,1]; ds_k/dp_1=+D[k,1]
}
div_from<-function(dd,src=1){ D<-dd$D; d<- -D[-src,src]/D[src,src]; setNames(d, paste0("to_",(1:J)[-src])) }
DD_true<-dmat(delta_true,truth$sigma_x1,truth$sigma_p,truth$alpha)
DD_blp <-dmat(B$delta_hat,B$sigma_x1,B$sigma_p,B$alpha)
# logit IIA diversion from product1: s_k/(1-s_1)
slog<-rowMeans(Sobs); div_log<-slog[-1]/(1-slog[1])
cat("=== diversion from product 1 (high quality, nest-mate = product 2) ===\n")
cat("        to2    to3    to4    to5    to6\n")
cat("truth ", round(div_from(DD_true),3),"\n")
cat("BLP   ", round(div_from(DD_blp),3),"\n")
cat("logit ", round(div_log,3),"\n")

cat("\n=== own-price elasticities (diag), median product ===\n")
cat("truth diag:",round(diag(E_true),2),"\n")
cat("BLP   diag:",round(diag(E_blp),2),"\n")
cat("logit diag:",round(diag(E_log),2),"\n")

## ---- GMM sandwich SEs (numerical Jacobian for sigma) ----
dt[, sumx1_riv:=sum(x1)-x1,by=t];dt[,sumx2_riv:=sum(x2)-x2,by=t];dt[,sumw_riv:=sum(w)-w,by=t]
dt[,dif_x1:={v<-x1;sapply(seq_along(v),function(i)sum((v[i]-v[-i])^2))},by=t]
dt[,dif_x2:={v<-x2;sapply(seq_along(v),function(i)sum((v[i]-v[-i])^2))},by=t]
dt[,dif_w:={v<-w;sapply(seq_along(v),function(i)sum((v[i]-v[-i])^2))},by=t]
Xlin<-cbind(1,dt$x1,dt$x2,dt$p)
Zmat<-cbind(1,dt$x1,dt$x2,dt$w,dt$sumx1_riv,dt$sumx2_riv,dt$sumw_riv,dt$dif_x1,dt$dif_x2,dt$dif_w)
N<-nrow(Xlin)
xi_hat<-as.vector(B$delta_hat)-Xlin%*%c(B$beta["beta0"],B$beta["beta_x1"],B$beta["beta_x2"],B$beta["p"])
# moment jacobian: columns for theta1 = -Z'X/N ; for sigma numeric via delta(sigma)
pred<-function(d,mu){V<-sweep(mu,2,d,"+");eV<-exp(V);colMeans(eV/(1+rowSums(eV)))}
contract<-function(mu,so,d0){d<-d0;for(it in 1:3000){sp<-pred(d,mu);dn<-d+log(so)-log(sp);if(max(abs(dn-d))<1e-12){d<-dn;break};d<-dn};d}
delta_of<-function(s1,sp){mu<-lapply(1:Tm,function(t)outer(s1*nu_x1,X1[,t])+outer(sp*nu_p,P[,t]))
  dl<-matrix(0,J,Tm);d0<-B$delta_hat;for(t in 1:Tm)dl[,t]<-contract(mu[[t]],Sobs[,t],d0[,t]);as.vector(dl)}
h<-1e-4
dd_s1<-(delta_of(B$sigma_x1+h,B$sigma_p)-delta_of(B$sigma_x1-h,B$sigma_p))/(2*h)
dd_sp<-(delta_of(B$sigma_x1,B$sigma_p+h)-delta_of(B$sigma_x1,B$sigma_p-h))/(2*h)
# xi = delta - X theta1 ; dxi/dtheta1 = -X ; dxi/dsigma = ddelta/dsigma
Gcols<-cbind(-Xlin, dd_s1, dd_sp)          # N x 6 (beta0,bx1,bx2,p, s1, sp)
G<-t(Zmat)%*%Gcols/N
Wm<-solve(crossprod(Zmat)/N)
S<-crossprod(Zmat*as.vector(xi_hat))/N
bread<-solve(t(G)%*%Wm%*%G)
V<-bread%*%(t(G)%*%Wm%*%S%*%Wm%*%G)%*%bread/N
se<-sqrt(diag(V)); names(se)<-c("beta0","beta_x1","beta_x2","alpha_negp","sigma_x1","sigma_p")
se["alpha_negp"]<-se["alpha_negp"]  # note: coef on p is -alpha; se same
cat("\n=== BLP GMM standard errors ===\n")
cat("alpha  est",round(B$alpha,3)," se",round(se["alpha_negp"],3),"\n")
cat("beta_x1est",round(B$beta["beta_x1"],3)," se",round(se["beta_x1"],3),"\n")
cat("beta_x2est",round(B$beta["beta_x2"],3)," se",round(se["beta_x2"],3),"\n")
cat("sig_x1 est",round(B$sigma_x1,3)," se",round(se["sigma_x1"],3),"\n")
cat("sig_p  est",round(B$sigma_p,3)," se",round(se["sigma_p"],3),"\n")

saveRDS(list(E_true=E_true,E_blp=E_blp,E_log=E_log,
             div_true=div_from(DD_true),div_blp=div_from(DD_blp),div_log=div_log,
             se=se, DD_true=DD_true, DD_blp=DD_blp), "ch20_elast.rds")

## ---------- Figure 1: recovery of the price coefficient / demand slope ----------
svglite::svglite("../../assets/fig/fig_20_recovery.svg", width=8.2, height=4.4)
par(mar=c(4.4,4.4,2.6,1), family="sans")
est<-c(OLS=L$a_ols, `IV logit`=L$a_iv, BLP=B$alpha, truth=truth$alpha)
cols<-c("#c53030","#dd8452","#2b6cb0","#1a365d")
bp<-barplot(est, col=cols, border=NA, ylim=c(min(est)-0.15,0.85),
    ylab=expression(paste("price coefficient  ", alpha, "  (per $)")),
    names.arg=c("OLS logit","IV logit","full BLP","truth"), cex.names=1.0)
abline(h=0,col="#718096"); abline(h=truth$alpha,lty=3,col="#1a365d")
text(bp, est+ifelse(est<0,-0.06,0.04), sprintf("%.2f",est), cex=0.95)
title(main="Recovering the price coefficient: OLS gets the sign wrong; BLP hits the truth",
      cex.main=0.98,font.main=1,col.main="#1a365d")
dev.off(); cat("wrote fig_20_recovery.svg\n")

## ---------- Figure 2: diversion from product 1 (substitution) ----------
M<-rbind(truth=div_from(DD_true), BLP=div_from(DD_blp), `logit (IIA)`=div_log)
svglite::svglite("../../assets/fig/fig_20_diversion.svg", width=8.2, height=4.4)
par(mar=c(4.2,4.4,2.6,1), family="sans")
cols2<-c("#1a365d","#2b6cb0","#c53030")
barplot(M, beside=TRUE, col=cols2, border=NA, ylim=c(0,max(M)*1.25),
   names.arg=c("to 2 (nest-mate)","to 3","to 4","to 5","to 6"),
   ylab="diversion ratio from product 1", cex.names=0.95)
legend("topright", rownames(M), fill=cols2, border=NA, bty="n", cex=0.95)
title(main="Diversion from product 1: BLP recovers the strong nest-mate substitution logit misses",
      cex.main=0.95,font.main=1,col.main="#1a365d")
dev.off(); cat("wrote fig_20_diversion.svg\n")
