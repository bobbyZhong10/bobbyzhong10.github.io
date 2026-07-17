# ============================================================
# Chapter 15 - Nested logit and PANEL mixed logit (SML), with control function.
# Then own/cross elasticities and diversion ratios vs the DGP truth.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(15)
D <- readRDS("ch19_data.rds")
long <- D$long; J <- D$J; nest <- D$nest; truth <- D$truth
setkey(long, oid, alt)
Nocc <- uniqueN(long$oid)

mkmat <- function(col) matrix(long[[col]], ncol=J+1L, byrow=TRUE)   # Nocc x (J+1)
RATING<-mkmat("rating"); PRICE<-mkmat("price"); CF<-mkmat("cf")
NEST1 <- (mkmat("nestj")==1L)*1.0
DA<-mkmat("d1");DB<-mkmat("d2");DC<-mkmat("d3");DD<-mkmat("d4")
CHOSEN<-mkmat("chosen")
nestrow <- mkmat("nestj")[1,]              # nest id of each column (same all rows)
in1 <- which(nestrow==1L); in2 <- which(nestrow==2L)
uid_occ <- long[alt==0L, uid]              # user id per occasion (ordered by oid)

baseV <- function(dA,dB,dC,dD,br,rho) dA*DA+dB*DB+dC*DC+dD*DD+br*RATING+rho*CF

## ---- Nested logit MLE (common lambda inside; outside lambda=1) ----
negll_nl <- function(par){
  V <- baseV(par[1],par[2],par[3],par[4],par[5],par[7]) + par[6]*PRICE
  lam <- plogis(par[8])
  Vs <- V; Vs[,c(in1,in2)] <- V[,c(in1,in2)]/lam
  eV <- exp(Vs)
  Dout<-eV[,1]; D1<-rowSums(eV[,in1,drop=FALSE]); D2<-rowSums(eV[,in2,drop=FALSE])
  den <- Dout + D1^lam + D2^lam
  P <- matrix(0,Nocc,J+1L); P[,1]<-Dout/den
  P[,in1]<-eV[,in1,drop=FALSE]*(D1^(lam-1))/den
  P[,in2]<-eV[,in2,drop=FALSE]*(D2^(lam-1))/den
  -sum(log(pmax(rowSums(P*CHOSEN),1e-12)))
}
fit_nl <- optim(c(0.5,0.2,0.3,-0.05,0.69,-0.245,0.30,qlogis(0.7)), negll_nl,
                method="BFGS", control=list(maxit=800,reltol=1e-10), hessian=TRUE)
se_nl<-sqrt(diag(solve(fit_nl$hessian))); pnl<-fit_nl$par; lam_hat<-plogis(pnl[8])
se_lam<-dlogis(pnl[8])*se_nl[8]
cat("--- Nested logit MLE (with control function) ---\n")
cat("alpha =",round(-pnl[6],3)," br =",round(pnl[5],3)," rho(cf) =",round(pnl[7],3),"\n")
cat("lambda =",round(lam_hat,3)," (SE ",round(se_lam,3),
    ")  within-nest corr 1-lambda^2 =",round(1-lam_hat^2,3),"\n",sep="")
cat("WTP per star = $",round(pnl[5]/(-pnl[6]),2),"\n",sep="")

## ---- PANEL mixed logit via SML ----
## random coef on price (alpha_i) and nest1 (lam_i), fixed within user across occasions.
## person likelihood: P_i = (1/R) sum_r prod_t p_{it}^r
R <- 100L
Nu <- D$Nusr
NUd  <- matrix(rnorm(Nu*R), Nu, R)     # per-user price-coef draws
ETAd <- matrix(rnorm(Nu*R), Nu, R)     # per-user nest-taste draws
# expand user draws to occasions
NU_occ  <- NUd[uid_occ,,drop=FALSE]    # Nocc x R
ETA_occ <- ETAd[uid_occ,,drop=FALSE]
usr_index <- split(seq_len(Nocc), uid_occ)   # occasion rows per user

negll_mxl <- function(par){
  dA<-par[1];dB<-par[2];dC<-par[3];dD<-par[4];br<-par[5]
  ab<-par[6]; sa<-exp(par[7]); sl<-exp(par[8]); rho<-par[9]
  B <- baseV(dA,dB,dC,dD,br,rho)
  # accumulate per-user L_i^r then average over r
  logLi_r <- matrix(0, Nu, R)
  for (r in 1:R){
    alpha_r <- ab + sa*NU_occ[,r]        # Nocc
    lam_r   <- sl*ETA_occ[,r]
    V <- B - alpha_r*PRICE + lam_r*NEST1
    V <- V - apply(V,1,max)
    eV <- exp(V); P <- eV/rowSums(eV)
    pc <- pmax(rowSums(P*CHOSEN),1e-12)  # chosen prob per occasion
    lpc <- log(pc)
    # sum log chosen-prob within user -> logLi_r
    logLi_r[,r] <- as.numeric(tapply(lpc, uid_occ, sum))
  }
  mx <- apply(logLi_r,1,max)
  Pi <- mx + log(rowMeans(exp(logLi_r - mx)))   # log( mean_r exp(logLi_r) )
  -sum(Pi)
}
st_mxl <- c(0.57,0.21,0.26,-0.06,0.69,0.245,log(0.15),log(0.7),0.30)
fit_mxl <- optim(st_mxl, negll_mxl, method="BFGS",
                 control=list(maxit=600,reltol=1e-7), hessian=TRUE)
se_mxl <- sqrt(diag(solve(fit_mxl$hessian)))
pm <- fit_mxl$par
# transform sigma back and delta-method SE
sa_hat<-exp(pm[7]); sl_hat<-exp(pm[8])
se_sa<-sa_hat*se_mxl[7]; se_sl<-sl_hat*se_mxl[8]
est_m <- c(pm[1:6], sigma_a=sa_hat, sigma_l=sl_hat, cf=pm[9])
se_m  <- c(se_mxl[1:6], se_sa, se_sl, se_mxl[9])
names(est_m)<-c("dA","dB","dC","dD","rating","alpha_bar","sigma_a","sigma_l","cf")
cat("\n--- PANEL mixed logit MLE (SML, R=100) ---\n")
print(round(rbind(est=est_m, se=se_m),3))
cat("WTP per star (at mean) = $",round(est_m["rating"]/est_m["alpha_bar"],2),"\n",sep="")
cat("TRUTH: alpha_bar=",truth$alpha_bar," sigma_a=",truth$sigma_a,
    " sigma_l=",truth$sigma_l," br=",truth$beta_r,"\n",sep="")

## ---- Substitution: raise product A price by $1; diversion vs TRUTH ----
grid<-D$grid
rbar<-sapply(1:4,function(k) mean(grid[j==k]$rating))
pbar<-sapply(1:4,function(k) mean(grid[j==k]$price))
shares_plain <- function(pv,par){V<-c(0,par[1:4]+par[5]*rbar+par[6]*pv);e<-exp(V-max(V));e/sum(e)}
shares_nested<- function(pv,par,lam){
  V<-c(0,par[1:4]+par[5]*rbar+par[6]*pv); Vs<-V; Vs[2:5]<-V[2:5]/lam
  eV<-exp(Vs-max(Vs)); i1<-c(F,T,T,F,F);i2<-c(F,F,F,T,T)
  D1<-sum(eV[i1]);D2<-sum(eV[i2]);Do<-eV[1]; den<-Do+D1^lam+D2^lam
  P<-numeric(5);P[1]<-Do/den;P[i1]<-eV[i1]*D1^(lam-1)/den;P[i2]<-eV[i2]*D2^(lam-1)/den;P
}
Rb<-6000L
shares_mixed<-function(pv,par){
  set.seed(202);nu<-rnorm(Rb);eta<-rnorm(Rb)
  base<-c(0,par["dA"],par["dB"],par["dC"],par["dD"])+c(0,par["rating"]*rbar);nst<-c(0,1,1,0,0)
  S<-numeric(5)
  for(r in 1:Rb){a<-par["alpha_bar"]+par["sigma_a"]*nu[r];l<-par["sigma_l"]*eta[r]
    V<-base-a*c(0,pv)+l*nst;e<-exp(V-max(V));S<-S+e/sum(e)}
  S/Rb
}
shares_truth<-function(pv){
  set.seed(202);nu<-rnorm(Rb);eta<-rnorm(Rb)
  base<-c(0,truth$delta)+c(0,truth$beta_r*rbar);nst<-c(0,1,1,0,0)
  S<-numeric(5)
  for(r in 1:Rb){a<-truth$alpha_bar+truth$sigma_a*nu[r];l<-truth$sigma_l*eta[r]
    V<-base-a*c(0,pv)+l*nst;e<-exp(V-max(V));S<-S+e/sum(e)}
  S/Rb
}
diversion<-function(sf,...){
  s0<-sf(pbar,...);p1<-pbar;p1[1]<-p1[1]+1;s1<-sf(p1,...);ds<-s1-s0
  own_elas<-(ds[2]/s0[2])/(1/pbar[1]); div<- -ds[-2]/ds[2]
  names(div)<-c("outside","B","C","D"); list(s0=s0,own_elas=own_elas,div=div)
}
tr<-diversion(shares_truth)
pl<-diversion(shares_plain, par=c(D$est$cf[1:5], D$est$cf[6]))  # CF-corrected plain logit
nl<-diversion(shares_nested,par=pnl[1:6],lam=lam_hat)
mx<-diversion(shares_mixed, par=est_m)
cat("\n--- Own-price elasticity of A and diversion from A (price +$1) ---\n")
cat("             own_elas  D->outside  D->B    D->C    D->D\n")
prow<-function(nm,o) cat(sprintf("%-12s %8.3f  %8.3f %8.3f %8.3f %8.3f\n",
   nm,o$own_elas,o$div["outside"],o$div["B"],o$div["C"],o$div["D"]))
prow("TRUTH",tr);prow("plain logit",pl);prow("nested",nl);prow("mixed",mx)
cat("\nrepresentative shares (0,A,B,C,D) TRUTH:",round(tr$s0,3),"\n")

saveRDS(list(pnl=pnl,se_nl=se_nl,lam_hat=lam_hat,se_lam=se_lam,
             est_m=est_m,se_m=se_m,
             div=list(truth=tr,plain=pl,nested=nl,mixed=mx),rbar=rbar,pbar=pbar),
        "ch19_estimates.rds")
cat("saved ch19_estimates.rds\n")
