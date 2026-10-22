# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of panel_fns.R on small synthetic panels.
# Run from the chapter code directory (tests/run_tests.R does this).
# Every hand-written estimator is checked against a package (plm, fixest, sandwich,
# clubSandwich) or against an algebraic identity; the card's five functions (within_tf,
# mundlak_fit, hausman_classic, hausman_robust, lead_test) come first.
# ------------------------------------------------------------------
source("panel_fns.R")
suppressPackageStartupMessages({library(plm); library(fixest); library(sandwich); library(clubSandwich)})
rel <- function(a, b) max(abs(unname(a) - unname(b))) / max(1, max(abs(unname(b))))
wald <- function(b, V) drop(crossprod(b, solve(V, b)))

## ---- synthetic panels: balanced (N = 150, T = 5) and an unbalanced subsample ----------
set.seed(10901)
N <- 150; TT <- 5
id <- rep(1:N, each = TT); t <- rep(1:TT, times = N)
a <- rnorm(N)
x1 <- 0.6 * a[id] + rnorm(N * TT)
x2 <- rnorm(N * TT) + 0.3 * t
y <- 1 + 0.5 * x1 - 0.3 * x2 + a[id] + 0.2 * t + rnorm(N * TT) * (1 + 0.5 * abs(x1))
X <- cbind(x1 = x1, x2 = x2)
bal <- data.frame(id, t, y, x1, x2)
keep <- t == 1 | runif(N * TT) > 0.2 + 0.2 * (a[id] > 0)        # gaps, and selection on alpha
unb <- bal[keep, ]
Xu <- as.matrix(unb[, c("x1", "x2")])
pb <- pdata.frame(bal, index = c("id", "t")); pu <- pdata.frame(unb, index = c("id", "t"))

## ---- within_tf: Y_it - Ybar_i and X_it - Xbar_i over the T_i observed periods ----------
w <- within_tf(unb$y, Xu, unb$id)
stopifnot(identical(names(w), c("y", "X")), is.matrix(w$X), identical(colnames(w$X), c("x1", "x2")),
          max(abs(rowsum(cbind(w$y, w$X), unb$id))) < 1e-10,
          max(abs(w$y - (unb$y - ave(unb$y, unb$id)))) < 1e-12)
stopifnot(identical(colnames(within_tf(bal$y, bal$x1, bal$id)$X), "x"))         # a vector X becomes one column "x"

## ---- fe_fit = plm within (coefficients, classical and cluster vcov), balanced and unbalanced ----
for (d in list(list(dat = bal, pd = pb, X = X), list(dat = unb, pd = pu, X = Xu))) {
  f <- fe_fit(d$dat$y, d$X, d$dat$id)
  p <- plm(y ~ x1 + x2, data = d$pd, model = "within")
  stopifnot(rel(f$b, coef(p)) < 1e-10, rel(f$V_classic, vcov(p)) < 1e-10,
            rel(f$V, vcovHC(p, method = "arellano", type = "HC0", cluster = "group")) < 1e-10)
  # least-squares dummy-variable regression gives the same slopes (Wooldridge 2010 §10.5.3)
  lsdv <- lm(y ~ x1 + x2 + factor(id), data = d$dat)
  stopifnot(rel(coef(lsdv)[c("x1", "x2")], f$b) < 1e-10)
}

## ---- re_fit = plm random (Swamy-Arora), be_fit = plm between, fd_fit = plm fd (balanced) ----
re <- re_fit(bal$y, X, bal$id)
pre <- plm(y ~ x1 + x2, data = pb, model = "random")
stopifnot(rel(re$b, coef(pre)) < 1e-10, rel(re$V_classic, vcov(pre)) < 1e-10,
          abs(re$varphi - ercomp(pre)$theta) < 1e-10,
          abs(re$sigma2_eps - ercomp(pre)$sigma2[["idios"]]) < 1e-10,
          abs(re$sigma2_alpha - ercomp(pre)$sigma2[["id"]]) < 1e-10)
stopifnot(inherits(tryCatch(re_fit(unb$y, Xu, unb$id), error = function(e) e), "error"))   # balanced panels only
be <- be_fit(bal$y, X, bal$id)
pbe <- plm(y ~ x1 + x2, data = pb, model = "between")
stopifnot(rel(be$b, coef(pbe)) < 1e-10, rel(be$V_classic, vcov(pbe)) < 1e-10)
fd <- fd_fit(bal$y, X, bal$id, bal$t)
pfd <- plm(y ~ x1 + x2, data = pb, model = "fd")
stopifnot(rel(fd$b, coef(pfd)) < 1e-10, rel(fd$V, vcovHC(pfd, method = "arellano", type = "HC0", cluster = "group")) < 1e-10)
# RE is a matrix-weighted average of between and within (balanced, Wooldridge 2010 §10.7.2)
fe <- fe_fit(bal$y, X, bal$id)
xw <- within_tf(bal$y, X, bal$id)$X
xb <- rowsum(X, bal$id) / TT; xb <- sweep(xb, 2, colMeans(xb))
SW <- crossprod(xw); SB <- TT * crossprod(xb); Wb <- solve(SW + (1 - re$varphi)^2 * SB) %*% ((1 - re$varphi)^2 * SB)
stopifnot(rel(re$b[c("x1", "x2")], drop(Wb %*% be$b[c("x1", "x2")] + (diag(2) - Wb) %*% fe$b)) < 1e-10)

## ---- mundlak_fit: b_x = FE exactly (balanced and unbalanced), b_xbar = BE - FE (balanced) ----
mb <- mundlak_fit(bal$y, X, bal$id)
mu <- mundlak_fit(unb$y, Xu, unb$id)
stopifnot(identical(names(mb)[1:3], c("b_x", "b_xbar", "V")),
          max(abs(mb$b_x - fe$b)) < 1e-8, max(abs(mu$b_x - fe_fit(unb$y, Xu, unb$id)$b)) < 1e-8,
          max(abs(mb$b_xbar - (be$b[c("x1", "x2")] - fe$b))) < 1e-8)
# V = cluster-robust (by unit, no small-sample factor) = sandwich::vcovCL(HC0) = clubSandwich CR0
xbar_u <- apply(Xu, 2, function(v) ave(v, unb$id))
du <- data.frame(y = unb$y, x1 = Xu[, 1], x2 = Xu[, 2], xbar_x1 = xbar_u[, 1], xbar_x2 = xbar_u[, 2])
lmu <- lm(y ~ x1 + x2 + xbar_x1 + xbar_x2, data = du)
nm <- names(coef(lmu))
stopifnot(identical(sort(rownames(mu$V)), sort(nm)),
          rel(mu$V[nm, nm], vcovCL(lmu, cluster = unb$id, type = "HC0", cadjust = FALSE)) < 1e-10,
          rel(mu$V[nm, nm], as.matrix(vcovCR(lmu, cluster = unb$id, type = "CR0"))) < 1e-10)
# hausman_robust = cluster-robust Wald of xi = 0
hr <- hausman_robust(mu)
Vx <- mu$V[c("xbar_x1", "xbar_x2"), c("xbar_x1", "xbar_x2")]
stopifnot(identical(names(hr), c("stat", "df", "p")), hr$df == 2,
          abs(hr$stat - wald(mu$b_xbar, Vx)) < 1e-8, abs(hr$p - pchisq(hr$stat, 2, lower.tail = FALSE)) < 1e-12)

## ---- mundlak_fit(time_fe = TRUE) = two-way FE, balanced and unbalanced ---------------------
for (d in list(list(dat = bal, pd = pb, X = X), list(dat = unb, pd = pu, X = Xu))) {
  m2 <- mundlak_fit(d$dat$y, d$X, d$dat$id, time_fe = TRUE, t = d$dat$t)
  tw <- coef(feols(y ~ x1 + x2 | id + t, data = d$dat, notes = FALSE))
  stopifnot(max(abs(m2$b_x - tw)) < 1e-8,
            max(abs(m2$b_x - coef(plm(y ~ x1 + x2, data = d$pd, model = "within", effect = "twoways")))) < 1e-8)
}
# without the unit averages of the time dummies the identity fails in the unbalanced panel
naive <- coef(lm(y ~ x1 + x2 + xbar_x1 + xbar_x2 + factor(unb$t), data = du))[c("x1", "x2")]
stopifnot(max(abs(naive - coef(feols(y ~ x1 + x2 | id + t, data = unb, notes = FALSE)))) > 1e-6)
stopifnot(inherits(tryCatch(mundlak_fit(bal$y, X, bal$id, time_fe = TRUE), error = function(e) e), "error"))

## ---- hausman_classic = plm::phtest(fe, re) ------------------------------------------------
pfe <- plm(y ~ x1 + x2, data = pb, model = "within")
hc <- hausman_classic(coef(pfe), vcov(pfe), coef(pre)[c("x1", "x2")], vcov(pre)[c("x1", "x2"), c("x1", "x2")])
ph <- phtest(pfe, pre)
stopifnot(identical(names(hc), c("stat", "df", "p")), abs(hc$stat - unname(ph$statistic)) < 1e-8,
          hc$df == unname(ph$parameter), abs(hc$p - ph$p.value) < 1e-12)
# the same statistic from the hand-written fits
hc2 <- hausman_classic(fe$b, fe$V_classic, re$b[c("x1", "x2")], re$V_classic[c("x1", "x2"), c("x1", "x2")])
stopifnot(abs(hc2$stat - hc$stat) < 1e-8)

## ---- lead_test: within regression with X_{i,t+1}, cluster-robust Wald of delta = 0 ----------
lp <- lead_test(unb$y, Xu, unb$id, unb$t)
o <- order(unb$id, unb$t); du2 <- unb[o, ]
nxt <- c(du2$id[-1] == du2$id[-nrow(du2)] & du2$t[-1] == du2$t[-nrow(du2)] + 1, FALSE)
du2$x1_lead <- c(du2$x1[-1], NA); du2$x2_lead <- c(du2$x2[-1], NA)
dl <- du2[nxt, ]
pl <- plm(y ~ x1 + x2 + x1_lead + x2_lead, data = pdata.frame(dl, index = c("id", "t")), model = "within")
Vl <- vcovHC(pl, method = "arellano", type = "HC0", cluster = "group")[c("x1_lead", "x2_lead"), c("x1_lead", "x2_lead")]
st <- wald(coef(pl)[c("x1_lead", "x2_lead")], Vl)
stopifnot(is.numeric(lp), length(lp) == 1, abs(attr(lp, "stat") - st) < 1e-8, attr(lp, "df") == 2,
          abs(lp - pchisq(st, 2, lower.tail = FALSE)) < 1e-12)
# feedback from last month's shock to this month's promotions: the lead test rejects
set.seed(10902)
P3t <- list(N = 1500, T = 5, alpha_sd = 1, x_mean = 5, x_alpha = 0.4, feedback = -0.6, x_sd = 1, beta = 0.5,
            mu0 = 10, sigma_e = 1)
d3 <- gen_P3(P3t)
stopifnot(lead_test(d3$y, d3$x, d3$id, d3$t) < 1e-6)

## ---- fd_ar1: AR(1) coefficient of the FD residuals (Wooldridge 2010 §10.6.3) ----------------
set.seed(10903)
P2t <- list(N = 3000, T = 6, alpha_sd = 1, x_mean = 5, x_alpha = 0.4, rho_x = 0.5, beta = 0.5, mu0 = 10, sigma_e = 1)
d2 <- gen_P2(P2t)
ri <- fd_ar1(fd_fit(d2$y_iid, d2$x, d2$id, d2$t)); rr <- fd_ar1(fd_fit(d2$y_rw, d2$x, d2$id, d2$t))
stopifnot(abs(ri$rho + 0.5) < 0.03, abs(rr$rho) < 0.03, ri$p_zero < 1e-6, rr$p_half < 1e-6)

## ---- chamberlain_md: minimum distance on the reduced form (Wooldridge 2010 §11.1.2) ---------
set.seed(10904)
Nc <- 20000; Tc <- 3
idc <- rep(1:Nc, each = Tc); tc <- rep(1:Tc, Nc)
xc <- rnorm(Nc * Tc) + rep(rnorm(Nc), each = Tc)
xi <- c(0.3, 0.2, -0.1)
alc <- 0.4 + drop(matrix(xc, Nc, Tc, byrow = TRUE) %*% xi) + rnorm(Nc, 0, 0.5)
yc <- 0.7 * xc + alc[idc] + rnorm(Nc * Tc, 0, 0.3)
cm <- chamberlain_md(yc, xc, idc, tc)
stopifnot(abs(cm$beta - 0.7) < 4 * cm$se, cm$se < 0.01, max(abs(cm$xi - xi)) < 0.03,
          cm$df == (Tc - 1) * (1 + Tc * 1) - 1, cm$p > 1e-4)
stopifnot(chamberlain_md(d3$y, d3$x, d3$id, d3$t)$p < 1e-6)      # strict exogeneity fails under feedback

## ---- ht_iv: Hausman-Taylor = IV on quasi-demeaned data (plm, random.method = "ht") ----------
set.seed(10905)
P4t <- list(N = 400, T = 4, alpha_sd = 2, x1_mean = 4, x1_mu = 1.2, x1_sd = 0.8, x2_mean = 5, x2_alpha = 0.3,
            x2_sd = 1, w_mu = 0.9, w_alpha = 0.35, w_cut = 0.3, beta1 = 0.3, beta2 = 0.5, gamma = 1.5,
            mu0 = 8, sigma_e = 1.5)
d4 <- gen_P4(P4t)
ht <- plm(y ~ x1 + x2 + w | x1 | x2, data = pdata.frame(d4, index = c("id", "t")),
          model = "random", random.method = "ht", inst.method = "baltagi")
hi <- ht_iv(d4$y, d4$x1, d4$x2, d4$w, d4$id, varphi = ercomp(ht)$theta)
stopifnot(rel(hi$coef[names(coef(ht))], coef(ht)) < 1e-8)

## ---- mg_fit = mean of the unit-by-unit OLS slopes (plm::pvcm, model = "within") -------------
pv <- pvcm(y ~ x1, data = pb, model = "within")
mg <- mg_fit(bal$y, bal$x1, bal$id)
stopifnot(abs(mg$est - mean(coef(pv)[, "x1"])) < 1e-10, abs(mg$se - sd(coef(pv)[, "x1"]) / sqrt(N)) < 1e-10)

## ---- slope_test: FE with X_it x (Z_i - mean Z), cluster-robust Wald (Wooldridge 2010 §11.7.4) ----
zc <- ave(bal$x1, bal$id); zc <- zc - mean(zc)
dd <- data.frame(id = bal$id, t = bal$t, y = bal$y, x1 = bal$x1, xz = bal$x1 * zc)
ps <- plm(y ~ x1 + xz, data = pdata.frame(dd, index = c("id", "t")), model = "within")
Vs <- vcovHC(ps, method = "arellano", type = "HC0", cluster = "group")
stz <- coef(ps)[["xz"]]^2 / Vs["xz", "xz"]
s1 <- slope_test(bal$y, bal$x1, bal$id, ave(bal$x1, bal$id))
stopifnot(abs(s1$stat - stz) < 1e-8, s1$df == 1)
set.seed(10906)
P5t <- list(N = 1000, T = 6, alpha_sd = 1, mu_mean = 5, mu_alpha = 0.4, mu_sd = 0.5, s_types = c(0.5, 1, 1.5, 2),
            s_probs = rep(0.25, 4), avg_slope = 0.5, kappa = 0.5, eta_sd = 0.05, design = "level", mu0 = 8, sigma_e = 1)
d5 <- gen_P5(P5t)
stopifnot(slope_test(d5$y, d5$x, d5$id, ave(d5$x, d5$id))$p < 1e-6)

## ---- ci95 ----------------------------------------------------------------------------------
stopifnot(isTRUE(all.equal(ci95(1, 0.5), c(1 - qnorm(0.975) * 0.5, 1 + qnorm(0.975) * 0.5))))
