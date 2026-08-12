# ============================================================
# Chapter 15 - Discrete choice: DGP (PANEL), plain logit MLE, control function
# Running case: "Selva" food-delivery platform, order-level choice.
# Each user places T orders over the window (panel), which is what identifies
# the mixing distribution in the mixed logit (Chapter 15, sec 4.4).
# Truth = random-coefficients logit: heterogeneous price sensitivity alpha_i
#         and a within-nest (cuisine-family) taste lam_i, plus price
#         endogeneity (price correlated with unobserved quality xi).
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(15)

## ---- design ------------------------------------------------
M   <- 300L          # markets (city-weeks): distinct product-characteristic bundles
J   <- 4L            # inside restaurants A,B,C,D  (outside good = 0)
nest <- c(1L,1L,2L,2L)   # A,B in nest 1 (comfort);  C,D in nest 2 (asian)
Nusr <- 1200L        # users
Tocc <- 8L           # orders per user (panel length)

# true utility parameters (utils)
delta_true <- c(0.80, 0.50, 0.60, 0.30)   # brand intercepts (outside = 0)
beta_r     <- 0.90                          # taste for rating (per star)
alpha_bar  <- 0.30                          # mean disutility of price (per $)
sigma_a    <- 0.12                          # sd of random price coefficient
sigma_l    <- 0.90                          # sd of nest-1 taste (within-nest corr)
sigma_xi   <- 0.45                          # sd of unobserved quality shock

## ---- product characteristics by market --------------------
base_rating <- c(4.4, 4.0, 4.2, 3.8)
grid <- CJ(m = 1:M, j = 1:J)
grid[, nestj := nest[j]]
grid[, rating := pmin(pmax(base_rating[j] + rnorm(.N, 0, 0.15), 3.0), 5.0)]
grid[, w  := rnorm(.N, 0, 1)]                       # cost instrument
grid[, xi := rnorm(.N, 0, sigma_xi)]               # unobserved quality
price_base <- c(16, 13, 15, 12)
grid[, price := price_base[j] + 1.6*w + 2.2*xi + 0.8*(rating - 4.0) + rnorm(.N, 0, 0.6)]
cat("corr(price, xi) =", round(cor(grid$price, grid$xi), 3),
    " | price range =", round(range(grid$price),1), "\n")

## ---- users and their random coefficients ------------------
usr <- data.table(uid = 1:Nusr)
usr[, alpha_i := alpha_bar + sigma_a*rnorm(Nusr)]
usr[, lam_i   := sigma_l*rnorm(Nusr)]

## ---- occasions: each user draws T markets, chooses once each --
occ <- CJ(uid = 1:Nusr, t = 1:Tocc)
occ[, oid := .I]
occ[, m := sample(1:M, .N, replace = TRUE)]
occ <- merge(occ, usr, by = "uid")
setkey(grid, m, j)

# simulate the chosen alternative from the TRUE random-coefficients model
Nocc <- nrow(occ)
choose_one <- function(mk, ai, li){
  g <- grid[.(mk), on = "m"]
  V <- delta_true[g$j] + beta_r*g$rating - ai*g$price + li*(g$nestj==1L) + g$xi
  V <- c(0, V)                                   # outside good utility 0
  u <- V - log(-log(runif(J+1)))                 # Type-I EV draws
  which.max(u) - 1L
}
occ[, choice := mapply(choose_one, m, alpha_i, lam_i)]
cat("\naggregate choice shares (0=outside):\n"); print(round(prop.table(table(occ$choice)),3))

## ---- long data: one row per occasion x alternative --------
long <- occ[, .(uid, oid, m, choice)][, .(alt = 0:J), by = .(uid, oid, m, choice)]
long <- merge(long, grid[, .(m, j, rating, price, w, xi)],
              by.x = c("m","alt"), by.y = c("m","j"), all.x = TRUE)
long[alt == 0L, `:=`(rating = 0, price = 0, w = 0, xi = 0)]
long[, nestj := ifelse(alt==0L, 0L, nest[pmax(alt,1L)])]
long[, chosen := as.integer(choice == alt)]
for (jj in 1:J) long[, paste0("d", jj) := as.integer(alt == jj)]
setkey(long, oid, alt)

## ==========================================================
## Plain conditional logit MLE (ignores xi -> price endogenous)
## ==========================================================
Xnames <- c("d1","d2","d3","d4","rating","price")
Xmat   <- as.matrix(long[, ..Xnames])
chosenv <- long$chosen
negll_logit <- function(par){
  V  <- Xmat %*% par
  Vm <- matrix(V, ncol=J+1L, byrow=TRUE); Vm <- Vm - apply(Vm,1,max)
  eV <- exp(Vm); P <- eV/rowSums(eV)
  ch <- matrix(chosenv, ncol=J+1L, byrow=TRUE)
  -sum(log(pmax(rowSums(P*ch),1e-12)))
}
fit_logit <- optim(c(0,0,0,0,0.5,-0.2), negll_logit, method="BFGS",
                   control=list(maxit=500,reltol=1e-10), hessian=TRUE)
b <- fit_logit$par; se_logit <- sqrt(diag(solve(fit_logit$hessian)))
names(b) <- c("dA","dB","dC","dD","rating","price_coef")
cat("\n--- Plain logit MLE (biased: ignores price endogeneity) ---\n")
print(round(rbind(est=b, se=se_logit),3))
alpha_naive <- -b["price_coef"]; wtp_naive <- b["rating"]/alpha_naive
cat("naive alpha =", round(alpha_naive,3),
    " | WTP per star = $", round(wtp_naive,2), "\n", sep="")

## ==========================================================
## Control function: first-stage price on instrument w, include residual
## ==========================================================
gin <- grid[, .(m,j,rating,price,w)]
fs  <- lm(price ~ w + rating + factor(j), data = gin)
Fw  <- summary(fs)$coef["w","t value"]^2
cat("\n--- First stage: price ~ w + rating + brand --- coef_w =",
    round(coef(fs)["w"],3), " partial-F(w) =", round(Fw,1), "\n")
gin[, cf := residuals(fs)]
long <- merge(long, gin[, .(m,j,cf)], by.x=c("m","alt"), by.y=c("m","j"), all.x=TRUE)
long[alt==0L, cf := 0]; setkey(long, oid, alt)
Xnames2 <- c("d1","d2","d3","d4","rating","price","cf")
Xmat2   <- as.matrix(long[, ..Xnames2])
negll_cf <- function(par){
  V  <- Xmat2 %*% par
  Vm <- matrix(V, ncol=J+1L, byrow=TRUE); Vm <- Vm - apply(Vm,1,max)
  eV <- exp(Vm); P <- eV/rowSums(eV)
  ch <- matrix(chosenv, ncol=J+1L, byrow=TRUE)
  -sum(log(pmax(rowSums(P*ch),1e-12)))
}
fit_cf <- optim(c(0,0,0,0,0.5,-0.2,0), negll_cf, method="BFGS",
                control=list(maxit=500,reltol=1e-10), hessian=TRUE)
bcf <- fit_cf$par; se_cf <- sqrt(diag(solve(fit_cf$hessian)))
names(bcf) <- c("dA","dB","dC","dD","rating","price_coef","cf")
cat("\n--- Logit with control function ---\n"); print(round(rbind(est=bcf,se=se_cf),3))
alpha_cf <- -bcf["price_coef"]; wtp_cf <- bcf["rating"]/alpha_cf
cat("cf alpha =", round(alpha_cf,3), " | WTP per star = $", round(wtp_cf,2), "\n", sep="")
cat("TRUTH: alpha_bar =", alpha_bar, " WTP per star = $", round(beta_r/alpha_bar,2), "\n", sep="")

saveRDS(list(occ=occ, usr=usr, grid=grid, long=long, nest=nest, J=J, Tocc=Tocc, Nusr=Nusr,
             est=list(logit=b, se_logit=se_logit, cf=bcf, se_cf=se_cf, Fw=Fw,
                      alpha_naive=alpha_naive, wtp_naive=wtp_naive,
                      alpha_cf=alpha_cf, wtp_cf=wtp_cf),
             truth=list(delta=delta_true,beta_r=beta_r,alpha_bar=alpha_bar,
                        sigma_a=sigma_a,sigma_l=sigma_l,sigma_xi=sigma_xi)),
        "ch19_data.rds")
cat("\nsaved ch19_data.rds  (Nocc =", Nocc, ")\n")
