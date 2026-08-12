# ============================================================
# Chapter 16 - Full BLP: contraction mapping + GMM over the random-coefficient
# parameters, concentrating out the linear parameters by 2SLS.
# Recovers alpha, beta, sigma; then true vs estimated own-elasticities and markups.
# Saves demand derivatives for Chapter 17.
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(16)
D <- readRDS("ch20_dgp.rds"); dt <- D$dt; J <- D$J; Tm <- D$Tm; truth <- D$truth
R <- D$R; nu_x1 <- D$nu_x1; nu_p <- D$nu_p
setkey(dt, t, j)

## observed shares and design matrices, organized by market
Sobs <- matrix(dt$s,  nrow=J)             # J x Tm
X1   <- matrix(dt$x1, nrow=J); X2 <- matrix(dt$x2, nrow=J); P <- matrix(dt$p, nrow=J)
# linear design X = [1, x1, x2, p]
Xlin <- cbind(1, dt$x1, dt$x2, dt$p)
# instruments: own (x1,x2,w); BLP rival sums (x1,x2,w); Gandhi-Houde differentiation IVs
dt[, sumx1_riv := sum(x1)-x1, by=t]; dt[, sumx2_riv := sum(x2)-x2, by=t]
dt[, sumw_riv  := sum(w)-w,  by=t]                       # rival cost shifters (identify price RC)
dt[, dif_x1 := { v<-x1; sapply(seq_along(v), function(i) sum((v[i]-v[-i])^2)) }, by=t]
dt[, dif_x2 := { v<-x2; sapply(seq_along(v), function(i) sum((v[i]-v[-i])^2)) }, by=t]
dt[, dif_w  := { v<-w;  sapply(seq_along(v), function(i) sum((v[i]-v[-i])^2)) }, by=t]
Zmat <- cbind(1, dt$x1, dt$x2, dt$w, dt$sumx1_riv, dt$sumx2_riv, dt$sumw_riv,
              dt$dif_x1, dt$dif_x2, dt$dif_w)
PZ   <- Zmat %*% solve(crossprod(Zmat)) %*% t(Zmat)   # projection for 2SLS

## market-level mu given theta2 = (sig_x1, sig_p): list of R x J matrices
make_mu <- function(sig1, sigp){
  lapply(1:Tm, function(t) outer(sig1*nu_x1, X1[,t]) + outer(sigp*nu_p, P[,t]))
}
## predicted inside shares in market t given delta_t (length J) and mu_t (RxJ)
pred_share <- function(delta_t, mu_t){
  V <- sweep(mu_t, 2, delta_t, "+"); eV <- exp(V)
  colMeans(eV / (1 + rowSums(eV)))
}
## contraction for market t
contract_t <- function(mu_t, sobs_t, d0, tol=1e-12, maxit=2000){
  d <- d0
  for (it in 1:maxit){
    sp <- pred_share(d, mu_t)
    dn <- d + log(sobs_t) - log(sp)
    if (max(abs(dn-d)) < tol){ d <- dn; break }
    d <- dn
  }
  d
}
## full delta(theta2) across markets
delta_theta <- function(sig1, sigp, dstart){
  mu <- make_mu(sig1, sigp)
  out <- matrix(0, J, Tm)
  for (t in 1:Tm) out[,t] <- contract_t(mu[[t]], Sobs[,t], dstart[,t])
  out
}

## GMM objective over theta2 (concentrate out linear params by 2SLS)
dstart <- log(Sobs) - log(matrix(dt$s0, nrow=J))   # logit delta as warm start
gmm_obj <- function(par){
  sig1 <- abs(par[1]); sigp <- abs(par[2])
  delta <- delta_theta(sig1, sigp, dstart)
  dvec <- as.vector(delta)                          # stacked j-fastest, matches dt order? dt keyed t,j
  # dt is keyed (t,j): row order is t=1 j=1..J, t=2..; as.vector(J x Tm) is column-major => j fastest within t. matches.
  th1 <- solve(t(Xlin)%*%PZ%*%Xlin, t(Xlin)%*%PZ%*%dvec)   # 2SLS
  xi  <- dvec - Xlin%*%th1
  g   <- t(Zmat)%*%xi
  as.numeric(t(g) %*% solve(crossprod(Zmat)) %*% g)
}

cat("optimizing BLP GMM over (sigma_x1, sigma_p)...\n")
fit <- optim(c(0.8,0.4), gmm_obj, method="Nelder-Mead",
             control=list(reltol=1e-9, maxit=400))
sig1_hat <- abs(fit$par[1]); sigp_hat <- abs(fit$par[2])
# recover linear params at optimum
delta_hat <- delta_theta(sig1_hat, sigp_hat, dstart)
dvec <- as.vector(delta_hat)
th1 <- as.numeric(solve(t(Xlin)%*%PZ%*%Xlin, t(Xlin)%*%PZ%*%dvec))
names(th1) <- c("beta0","beta_x1","beta_x2","p")
alpha_hat <- -th1["p"]
cat("\n--- BLP estimates ---\n")
cat("alpha    =", round(alpha_hat,3), "  (truth", truth$alpha, ")\n")
cat("beta_x1  =", round(th1["beta_x1"],3), "  (truth", truth$beta_x1, ")\n")
cat("beta_x2  =", round(th1["beta_x2"],3), "  (truth", truth$beta_x2, ")\n")
cat("beta0    =", round(th1["beta0"],3), "  (truth", truth$beta0, ")\n")
cat("sigma_x1 =", round(sig1_hat,3), "  (truth", truth$sigma_x1, ")\n")
cat("sigma_p  =", round(sigp_hat,3), "  (truth", truth$sigma_p, ")\n")

## ---- elasticities & markups: TRUE vs BLP-estimated, per market ----
## demand derivative matrix in market t: Delta_jk = -ds_j/dp_k (as arranged for markup)
deriv_mkt <- function(delta_t, X1t, Pt, sig1, sigp, alpha){
  mu <- outer(sig1*nu_x1, X1t) + outer(sigp*nu_p, Pt)
  V  <- sweep(mu, 2, delta_t, "+"); eV <- exp(V); sij <- eV/(1+rowSums(eV))
  s  <- colMeans(sij); alpha_i <- alpha + sigp*nu_p
  Del <- matrix(0,J,J)
  for (j in 1:J) for (k in 1:J){
    if (j==k) Del[j,k] <-  mean(alpha_i*sij[,j]*(1-sij[,j]))
    else      Del[j,k] <- -mean(alpha_i*sij[,j]*sij[,k])
  }
  list(s=s, Del=Del)   # Del is +ds_j/dp arranged so ds_j/dp_j>0 sign handled: own = -Del_jj
}
# own elasticity e_jj = (ds_j/dp_j)(p_j/s_j) ; ds_j/dp_j = -Del_jj (Del_jj is |own deriv|)
own_elas <- function(delta,X1m,Pm,sig1,sigp,alpha){
  e <- numeric(Tm*J); idx<-1
  for (t in 1:Tm){ o<-deriv_mkt(delta[,t],X1m[,t],Pm[,t],sig1,sigp,alpha)
    for (j in 1:J){ e[idx] <- (-o$Del[j,j])*(Pm[j,t]/o$s[j]); idx<-idx+1 } }
  e
}
# true delta uses true xi; reconstruct true delta from DGP
delta_true <- matrix(truth$beta0 + truth$beta_x1*dt$x1 + truth$beta_x2*dt$x2 -
                     truth$alpha*dt$p + dt$xi, nrow=J)
e_true <- own_elas(delta_true, X1, P, truth$sigma_x1, truth$sigma_p, truth$alpha)
e_blp  <- own_elas(delta_hat,  X1, P, sig1_hat, sigp_hat, alpha_hat)
cat("\n--- own-price elasticity (median across products) ---\n")
cat("TRUE BLP", round(median(e_true),3), " | BLP est", round(median(e_blp),3), "\n")

## implied marginal cost under BLP (single-product FOC, full derivative matrix)
markup_blp <- numeric(Tm*J); idx<-1
for (t in 1:Tm){ o<-deriv_mkt(delta_hat[,t],X1[,t],P[,t],sig1_hat,sigp_hat,alpha_hat)
  H<-diag(J); mk<-solve(H*o$Del, o$s)        # single-product ownership
  for(j in 1:J){ markup_blp[idx]<-mk[j]; idx<-idx+1 } }
mc_blp <- dt$p - markup_blp
cat("implied mc (BLP): median", round(median(mc_blp),2),
    " share<0", round(mean(mc_blp<0),3), " (truth median", round(median(dt$mc),2),")\n")

saveRDS(list(alpha=alpha_hat, beta=th1, sigma_x1=sig1_hat, sigma_p=sigp_hat,
             delta_hat=delta_hat, e_true=e_true, e_blp=e_blp,
             mc_blp=mc_blp, markup_blp=markup_blp,
             gmm_val=fit$value), "ch20_blp.rds")
cat("saved ch20_blp.rds\n")
