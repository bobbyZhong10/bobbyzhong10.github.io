# ============================================================
# Chapter 16 - Aggregate demand (BLP): DGP with a supply side.
# Running case: "Selva" platform, market-level data on J meal-kit brands
# across T local markets. True demand = random-coefficients logit (BLP);
# true supply = multiproduct Bertrand-Nash, so equilibrium prices are
# correlated with the unobserved quality xi (=> price endogeneity).
# This DGP is shared with Chapter 17 (marginal-cost recovery + merger).
# ============================================================
suppressPackageStartupMessages({library(data.table)})
set.seed(16)

## ---- structure --------------------------------------------
Tm <- 300L          # markets (city-weeks)
J  <- 6L            # brands, present in every market; + outside good
# baseline ownership: 6 single-product firms (each brand its own firm).
# Chapter 17 merges firm 1 and firm 2.
firm <- 1:J

## ---- true parameters --------------------------------------
beta0   <-  1.5     # intercept (with outside good normalized 0)
beta_x1 <-  1.6     # taste for quality x1
beta_x2 <-  0.8     # taste for attribute x2 (breadth)
alpha   <-  0.7     # mean price disutility (per $ ; prices ~ $6-13)
sigma_x1<-  1.1     # random coef sd on x1 (drives substitution toward similar quality)
sigma_p <-  0.15    # small random coef sd on price (kept modest: a normal price RC with
                    # large sd flips sign for a mass of consumers and destabilizes pricing)
# cost: mc = gamma0 + gamma_w w + gamma_x1 x1 + omega
gamma0  <-  2.2
gamma_w <-  0.8
gamma_x1<-  0.7
sd_xi   <-  1.3
sd_omega<-  0.30

## ---- product characteristics ------------------------------
brand_q  <- c(7.0, 6.6, 5.2, 5.6, 4.4, 4.0)/2   # brand base quality x1 (nest-ish: 1&2 high, 5&6 low)
brand_x2 <- c(1.0, 0.4, 0.8, 0.2, 0.9, 0.3)     # attribute (breadth)
dt <- CJ(t = 1:Tm, j = 1:J)
dt[, x1 := brand_q[j]  + rnorm(.N, 0, 0.50)]           # more cross-market variation => identifies sigma
dt[, x2 := brand_x2[j] + rnorm(.N, 0, 0.30)]
dt[, w  := 0.8 + 0.5*brand_q[j] + rnorm(.N, 0, 0.6)]   # cost shifter (input cost)
dt[, xi := rnorm(.N, 0, sd_xi)]
dt[, omega := rnorm(.N, 0, sd_omega)]
dt[, mc := gamma0 + gamma_w*w + gamma_x1*x1 + omega]
dt[, firm := firm[j]]
setkey(dt, t, j)

## ---- simulation draws for the TRUE shares (shared across markets) ----
R <- 400L
nu_x1 <- rnorm(R); nu_p <- rnorm(R)     # random-coefficient draws

# predicted shares in one market given delta (vector length J) and prices
shares_mkt <- function(delta, x1, p){
  # mu_ij = sigma_x1*nu_x1*x1 + sigma_p*nu_p*p   (R x J)
  mu <- outer(sigma_x1*nu_x1, x1) + outer(sigma_p*nu_p, p)    # R x J
  V  <- sweep(mu, 2, delta, "+")                              # R x J, add delta_j
  eV <- exp(V); denom <- 1 + rowSums(eV)                      # outside = 1
  sij <- eV / denom                                           # R x J
  list(s = colMeans(sij), sij = sij, denom = denom)
}

# multiproduct-Bertrand equilibrium prices in one market, given cost & demand
# FOC: p = mc + markup, markup = (H .* Delta)^{-1} s ; Delta_jk = -ds_j/dp_k
solve_prices <- function(mkt){
  x1 <- mkt$x1; x2 <- mkt$x2; mc <- mkt$mc; xi <- mkt$xi; own <- outer(mkt$firm, mkt$firm, "==")*1
  delta_of_p <- function(p) beta0 + beta_x1*x1 + beta_x2*x2 - alpha*p + xi
  p <- mc + 1.0                       # start
  for (it in 1:200){
    delta <- delta_of_p(p)
    sm <- shares_mkt(delta, x1, p); s <- sm$s; sij <- sm$sij
    # demand derivatives with random coef on price:
    # ds_j/dp_k = -(1/R) sum_r alpha_ir sij_r,j (1{j=k} - sij_r,k), alpha_ir = alpha + sigma_p nu_p
    alpha_i <- alpha - sigma_p*nu_p                  # R
    Delta <- matrix(0, J, J)
    for (j in 1:J) for (k in 1:J){
      if (j==k) Delta[j,k] <- mean(alpha_i * sij[,j]*(1 - sij[,j]))
      else      Delta[j,k] <- -mean(alpha_i * sij[,j]*sij[,k])
    }
    # note Delta here = -ds_j/dp_k arranged so markup = (own .* Delta)^{-1} s
    HDelta <- own * Delta
    markup <- solve(HDelta, s)
    p_new <- mc + markup
    if (max(abs(p_new - p)) < 1e-10){ p <- p_new; break }
    p <- 0.5*p + 0.5*p_new
  }
  list(p = p, s = s, markup = markup)
}

## ---- solve equilibrium market by market -------------------
dt[, p := NA_real_]; dt[, s := NA_real_]; dt[, markup := NA_real_]
for (tt in 1:Tm){
  mkt <- dt[t==tt]
  sol <- solve_prices(mkt)
  dt[t==tt, p := sol$p]; dt[t==tt, s := sol$s]; dt[t==tt, markup := sol$markup]
}
dt[, s0 := 1 - sum(s), by = t]        # outside share
dt[, mc_over_p := mc/p]

cat("=== DGP summary ===\n")
cat("price  range:", round(range(dt$p),2), " median", round(median(dt$p),2), "\n")
cat("mc     range:", round(range(dt$mc),2), "\n")
cat("markup (p-mc) median:", round(median(dt$markup),3),
    " | Lerner (p-mc)/p median:", round(median(dt$markup/dt$p),3), "\n")
cat("inside share median:", round(median(dt$s),4),
    " | outside share median:", round(median(dt$s0),3), "\n")
cat("corr(price, xi):", round(cor(dt$p, dt$xi),3),
    " | corr(price, w):", round(cor(dt$p, dt$w),3), "\n")
cat("any negative markup? ", any(dt$markup<=0), "\n")

saveRDS(list(dt=dt, J=J, Tm=Tm, firm=firm, R=R, nu_x1=nu_x1, nu_p=nu_p,
             truth=list(beta0=beta0,beta_x1=beta_x1,beta_x2=beta_x2,alpha=alpha,
                        sigma_x1=sigma_x1,sigma_p=sigma_p,
                        gamma0=gamma0,gamma_w=gamma_w,gamma_x1=gamma_x1)),
        "ch20_dgp.rds")
cat("saved ch20_dgp.rds\n")
