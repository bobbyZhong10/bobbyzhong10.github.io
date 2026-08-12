# ============================================================================
# 03_cluster.R  --  the "significant" result that isn't: naive i.i.d. SEs on
# a city-clustered rollout vastly overstate precision. Cluster-robust SEs fix
# it; with FEW clusters even they over-reject, and the wild cluster bootstrap /
# randomization inference restore honest size.
# ============================================================================

suppressMessages({ library(sandwich); library(lmtest) })
H <- readRDS("halcyon.rds"); truth <- H$truth; mkclus <- H$mkclus

d <- H$clus_main                                    # G=40 cities, n_c=200
m <- lm(y ~ D, data = d)
tau_hat <- coef(m)["D"]

se_iid  <- sqrt(vcov(m)["D","D"])
se_hc1  <- sqrt(vcovHC(m, type = "HC1")["D","D"])
se_cr1  <- sqrt(vcovCL(m, cluster = ~city, type = "HC1")["D","D"])
se_cr2  <- sqrt(vcovCL(m, cluster = ~city, type = "HC2")["D","D"])

cat("=== one clustered rollout (G=40 cities, 200 sessions each) ===\n")
cat(sprintf("tau_hat = %.4f   (truth %.2f)\n", tau_hat, truth$tau))
cat(sprintf("naive i.i.d. SE = %.4f   t = %5.2f   p = %.4g\n",
            se_iid, tau_hat/se_iid, 2*pnorm(-abs(tau_hat/se_iid))))
cat(sprintf("HC1 robust   SE = %.4f   t = %5.2f\n", se_hc1, tau_hat/se_hc1))
cat(sprintf("cluster CR1  SE = %.4f   t = %5.2f   p = %.4g\n",
            se_cr1, tau_hat/se_cr1, 2*pt(-abs(tau_hat/se_cr1), df = truth$G-1)))
cat(sprintf("cluster CR2  SE = %.4f   t = %5.2f\n", se_cr2, tau_hat/se_cr2))
cat(sprintf("SE inflation cluster/naive = %.2f   (Moulton predicts %.2f)\n",
            se_cr1/se_iid, truth$moulton))

## ---- wild cluster bootstrap (Rademacher), imposing the null tau=0 ----------
wild_cluster_p <- function(dat, B = 999, seed = 7) {
  set.seed(seed)
  m1 <- lm(y ~ D, data = dat); t_obs <- coef(m1)["D"]/sqrt(vcovCL(m1, ~city, type="HC1")["D","D"])
  m0 <- lm(y ~ 1, data = dat); u0 <- residuals(m0); fit0 <- fitted(m0)  # null-imposed
  cl <- unique(dat$city)
  tb <- replicate(B, {
    w <- sample(c(-1,1), length(cl), replace = TRUE)[match(dat$city, cl)]
    yb <- fit0 + u0 * w
    mb <- lm(yb ~ dat$D)
    coef(mb)[2] / sqrt(vcovCL(mb, ~dat$city, type="HC1")[2,2])
  })
  mean(abs(tb) >= abs(t_obs))
}
cat(sprintf("\nwild cluster bootstrap p (G=40) = %.4f\n", wild_cluster_p(d)))

## ---- Monte-Carlo size under the TRUE null tau=0 ----------------------------
size_study <- function(G, R = 1000, seed = 11) {
  set.seed(seed)
  rej <- t(replicate(R, {
    dd <- mkclus(G, truth$n_c, tau = 0)
    mm <- lm(y ~ D, data = dd); b <- unname(coef(mm)["D"])
    se_i <- sqrt(vcov(mm)["D","D"]); se_c <- sqrt(vcovCL(mm, ~city, type="HC1")["D","D"])
    c(iid   = as.numeric(abs(b/se_i) > 1.96),
      cr1_z = as.numeric(abs(b/se_c) > 1.96),
      cr1_t = as.numeric(abs(b/se_c) > qt(.975, G-1)))
  }))
  colMeans(rej)
}
cat("\n=== rejection rate of a TRUE null tau=0 at nominal 5% (1000 reps) ===\n")
s40 <- size_study(40); s10 <- size_study(10, seed = 12)
cat(sprintf("G=40: naive iid = %.3f | cluster(z) = %.3f | cluster t(G-1) = %.3f\n",
            s40["iid"], s40["cr1_z"], s40["cr1_t"]))
cat(sprintf("G=10: naive iid = %.3f | cluster(z) = %.3f | cluster t(G-1) = %.3f\n",
            s10["iid"], s10["cr1_z"], s10["cr1_t"]))

## ---- few clusters (G=10): one dataset with a real effect, compare p-values -
set.seed(313)
d10 <- mkclus(10, truth$n_c, tau = 0.25)
m10 <- lm(y ~ D, data = d10); b10 <- unname(coef(m10)["D"])
se10_cr <- sqrt(vcovCL(m10, ~city, type = "HC1")["D","D"])
t10 <- b10 / se10_cr
p10_z <- 2*pnorm(-abs(t10)); p10_t <- 2*pt(-abs(t10), df = 9)
p10_wcb <- wild_cluster_p(d10, B = 1999, seed = 99)
cat("\n=== few clusters (G=10), a real effect: three p-values for the same t ===\n")
cat(sprintf("tau_hat = %.3f  cluster t = %.2f\n", b10, t10))
cat(sprintf("normal p (z)          = %.4f  (over-optimistic with few clusters)\n", p10_z))
cat(sprintf("t(G-1=9) p            = %.4f\n", p10_t))
cat(sprintf("wild cluster bootstrap = %.4f  (the trustworthy one)\n", p10_wcb))

saveRDS(list(tau_hat=unname(tau_hat), se_iid=se_iid, se_hc1=se_hc1, se_cr1=se_cr1,
             se_cr2=se_cr2, moulton=truth$moulton, s40=s40, s10=s10,
             t10=t10, p10_z=p10_z, p10_t=p10_t, p10_wcb=p10_wcb), "res_cluster.rds")
