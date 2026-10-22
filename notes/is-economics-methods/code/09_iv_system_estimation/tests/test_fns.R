# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of iv_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# A toy two-equation simultaneous system (not the Bramble data) checks every function against a
# package or an algebraic identity:
#   tsls               ivreg::ivreg, sandwich::vcovHC(type = "HC0"); OLS when Z = X; (Z'X)^-1 Z'y when just identified
#   sur_fgls           systemfit(method = "SUR", methodResidCov = "noDfCor"); Theorems 7.5 and 7.6
#   three_sls          systemfit(method = "3SLS", methodResidCov = "noDfCor"); robust V = sandwich::vcovCL(residFit = FALSE);
#                      section 9.3.2 identities (just identified, diagonal Omega, one over-identified equation)
#   order/rank         truth tables
#   ar_test, ar_ci     ivmodel::AR.test (statistic, p-value, interval or whole line); a grid check for two rays
#   hansen_j, sargan   gmm::gmm (two-step, uncentered MDS) specTest; ivreg diagnostics Sargan
#   eff_F, first_stage_F   robust t^2 with one instrument; homoskedastic version = the usual F
#   cf_test            the control-function coefficient on the endogenous regressor equals 2SLS
#   reduced_form, pop_avar_sys   equilibrium identity; 3SLS avar <= 2SLS avar, equal under diagonal Omega
# ------------------------------------------------------------------
source("iv_fns.R")
suppressPackageStartupMessages({
  library(ivreg); library(systemfit); library(sandwich); library(ivmodel); library(gmm)
})
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

## ---- toy system --------------------------------------------------------------------------
## y1 = 1 + 0.8 y2 + 0.5 z1 + 0.3 w + e1           (excluded: z2, z3; over-identified by one)
## y2 = -1 + 0.4 y1 + 0.7 z2 - 0.5 z3 + 0.2 w + e2 (excluded: z1; just identified)
set.seed(901)
n <- 700
w <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n); z3 <- rnorm(n) + 0.4 * z2
E <- matrix(rnorm(2 * n), n) %*% chol(matrix(c(1, 0.6, 0.6, 1.5), 2))
E[, 1] <- E[, 1] * exp(0.3 * z1)                                       # heteroskedastic e1 (for the robust checks)
c1 <- 1 + 0.5 * z1 + 0.3 * w + E[, 1]; c2 <- -1 + 0.7 * z2 - 0.5 * z3 + 0.2 * w + E[, 2]
y1 <- (c1 + 0.8 * c2) / (1 - 0.32); y2 <- c2 + 0.4 * y1
d <- data.frame(y1, y2, z1, z2, z3, w)
Z <- cbind(const = 1, z1 = z1, z2 = z2, z3 = z3, w = w)
X1 <- cbind(const = 1, y2 = y2, z1 = z1, w = w)
X2 <- cbind(const = 1, y1 = y1, z2 = z2, z3 = z3, w = w)
eqs <- list(e1 = y1 ~ y2 + z1 + w, e2 = y2 ~ y1 + z2 + z3 + w)

## ---- chapter card: return values -------------------------------------------------------------
t1 <- tsls(y1, X1, Z)
stopifnot(is.list(t1), all(c("b", "V") %in% names(t1)), length(t1$b) == 4, all(dim(t1$V) == c(4, 4)),
          identical(names(t1$b), colnames(X1)))
s <- sur_fgls(list(y1, y2), list(X1, X2))
stopifnot(all(c("b", "V", "Sigma") %in% names(s)), length(s$b) == 9, all(dim(s$Sigma) == c(2, 2)))
th <- three_sls(list(y1, y2), list(X1, X2), Z)
stopifnot(all(c("b", "V") %in% names(th)), length(th$b) == 9, all(dim(th$V) == c(9, 9)))
stopifnot(isTRUE(order_condition(2, 1)), isTRUE(order_condition(1, 1)), identical(order_condition(0, 1), FALSE))
rc <- rank_condition(matrix(c(0.5, -0.2), 2, 1))
stopifnot(identical(names(rc), c("rank", "ok")), rc$rank == 1, isTRUE(rc$ok))

## ---- tsls: ivreg, HC0 sandwich, homoskedastic vcov, special cases --------------------------------
iv1 <- ivreg(y1 ~ y2 + z1 + w | z1 + z2 + z3 + w, data = d)
stopifnot(max(abs(t1$b - coef(iv1))) < 1e-10,
          rel(t1$V, unname(sandwich::vcovHC(iv1, type = "HC0"))) < 1e-8,
          rel(tsls(y1, X1, Z, vcov = "homo")$V, unname(vcov(iv1))) < 1e-8)
ols1 <- lm(y1 ~ y2 + z1 + w)
stopifnot(max(abs(tsls(y1, X1, X1)$b - coef(ols1))) < 1e-10,                       # Z = X: 2SLS is OLS
          rel(tsls(y1, X1, X1)$V, unname(sandwich::vcovHC(ols1, type = "HC0"))) < 1e-8,
          max(abs(ols_fit(y1, X1)$b - coef(ols1))) < 1e-10,
          rel(ols_fit(y1, X1, vcov = "homo")$V, unname(vcov(ols1))) < 1e-8)
Zj <- cbind(1, z1, z2, w)                                                           # just identified: IV = (Z'X)^-1 Z'y
stopifnot(max(abs(tsls(y1, X1, Zj)$b - drop(solve(crossprod(Zj, X1), crossprod(Zj, y1))))) < 1e-10)
stopifnot(inherits(try(tsls(y1, X1, Z[, 1:3]), silent = TRUE), "try-error"))      # fewer instruments than regressors

## ---- sur_fgls: systemfit SUR, Theorems 7.5 and 7.6 ----------------------------------------------
sf <- systemfit(list(e1 = y1 ~ z1 + w, e2 = y2 ~ z2 + z3 + w), method = "SUR", data = d, methodResidCov = "noDfCor")
s0 <- sur_fgls(list(y1, y2), list(cbind(1, z1, w), cbind(1, z2, z3, w)))
stopifnot(max(abs(s0$b - coef(sf))) < 1e-10, rel(s0$V, unname(vcov(sf))) < 1e-8,
          rel(s0$Sigma, unname(sf$residCovEst)) < 1e-10)
Xs <- cbind(1, z1, z2, w)                                                           # Theorem 7.6: identical regressors
sx <- sur_fgls(list(y1, y2), list(Xs, Xs))
stopifnot(max(abs(sx$b - c(qr.coef(qr(Xs), y1), qr.coef(qr(Xs), y2)))) < 1e-10)
sd0 <- sur_fgls(list(y1, y2), list(cbind(1, z1, w), cbind(1, z2, z3, w)), Sigma = diag(c(2, 3)))   # Theorem 7.5
stopifnot(max(abs(sd0$b - c(qr.coef(qr(cbind(1, z1, w)), y1), qr.coef(qr(cbind(1, z2, z3, w)), y2)))) < 1e-10)
# with the same Omega-hat, FGLS standard errors never exceed equation-by-equation OLS ones (Gauss-Markov)
se_ols0 <- sqrt(c(s0$Sigma[1, 1] * diag(solve(crossprod(cbind(1, z1, w)))),
                  s0$Sigma[2, 2] * diag(solve(crossprod(cbind(1, z2, z3, w))))))
stopifnot(all(sqrt(diag(s0$V)) <= se_ols0))

## ---- three_sls: systemfit 3SLS, robust sandwich, section 9.3.2 identities ------------------------------
f3 <- systemfit(eqs, method = "3SLS", inst = ~ z1 + z2 + z3 + w, data = d, methodResidCov = "noDfCor")
stopifnot(max(abs(th$b - coef(f3))) < 1e-8, rel(th$V, unname(vcov(f3))) < 1e-8)
thr <- three_sls(list(y1, y2), list(X1, X2), Z, vcov = "robust")
Vcl <- sandwich::vcovCL(f3, cluster = rep(seq_len(n), 2), type = "HC0", cadjust = FALSE, residFit = FALSE)
stopifnot(max(abs(thr$b - th$b)) < 1e-12, rel(thr$V, unname(Vcl)) < 1e-6)
t2 <- tsls(y2, X2, Z)
f2 <- systemfit(eqs, method = "2SLS", inst = ~ z1 + z2 + z3 + w, data = d)
stopifnot(max(abs(c(t1$b, t2$b) - coef(f2))) < 1e-10)
# eq 2 just identified, eq 1 over-identified: 3SLS = 2SLS for eq 1 exactly, not for eq 2
stopifnot(max(abs(th$b[1:4] - t1$b)) < 1e-10, max(abs(th$b[5:9] - t2$b)) > 1e-4)
# every equation just identified: 3SLS = 2SLS
Zjj <- cbind(1, z1, z2, w); X2j <- cbind(1, y1, z2, w)
thj <- three_sls(list(y1, y2), list(X1, X2j), Zjj)
stopifnot(max(abs(thj$b - c(tsls(y1, X1, Zjj)$b, tsls(y2, X2j, Zjj)$b))) < 1e-10)
# diagonal Omega-hat: 3SLS = 2SLS
thd <- three_sls(list(y1, y2), list(X1, X2), Z, Sigma = diag(c(1, 2)))
stopifnot(max(abs(thd$b - c(t1$b, t2$b))) < 1e-10)

## ---- rank_condition: rank deficiency, zero and wide matrices -------------------------------------
stopifnot(rank_condition(matrix(0, 1, 1))$rank == 0, !rank_condition(matrix(0, 1, 1))$ok,
          rank_condition(cbind(c(1, 2, 3), c(2, 4, 6)))$rank == 1, !rank_condition(cbind(c(1, 2, 3), c(2, 4, 6)))$ok,
          rank_condition(cbind(c(1, 0, 3), c(2, 4, 6)))$ok,
          !rank_condition(matrix(1, 1, 2))$ok,                                           # fewer rows than columns
          rank_condition(matrix(1e-10, 1, 1))$rank == 0, rank_condition(matrix(1e-6, 1, 1))$rank == 1)
stopifnot(inherits(try(order_condition(-1, 1), silent = TRUE), "try-error"))

## ---- reduced form: y' = z' Pi + v', Pi = -Delta Gamma^-1, v' = -eps' Gamma^-1 ---------------------
## toy structure in the form y'Gamma + z'Delta + eps' = 0 (column m = equation m)
Gam <- rbind(c(-1, 0.4), c(0.8, -1))
Del <- rbind(c(1, -1), c(0.5, 0), c(0, 0.7), c(0, -0.5), c(0.3, 0.2))
Pi <- reduced_form(Gam, Del)
Ypred <- Z %*% Pi - E %*% solve(Gam)
stopifnot(max(abs(Ypred - cbind(y1, y2))) < 1e-10)

## ---- Anderson-Rubin: statistic, p-value and confidence set against ivmodel ----------------------
W <- cbind(w)
ivm <- ivmodel(Y = y1, D = y2, Z = cbind(z2, z3), X = W)
for (b0 in c(0, 0.8, 2)) {
  a_pkg <- AR.test(ivm, beta0 = b0)
  a <- ar_test(y1, y2, cbind(z2, z3), cbind(1, W), b0)
  stopifnot(abs(a$stat - a_pkg$Fstat) < 1e-8 * max(1, a_pkg$Fstat), abs(a$p - a_pkg$p.value) < 1e-10,
            all(c(a$df1, a$df2) == a_pkg$df))
}
ci <- ar_ci(y1, y2, cbind(z2, z3), cbind(1, W))
ci_pkg <- AR.test(ivm)$ci
stopifnot(ci$type == "bounded", abs(ci$lo - ci_pkg[1, "lower"]) < 1e-6, abs(ci$hi - ci_pkg[1, "upper"]) < 1e-6,
          ar_covers(ci, 0.8) == (ci$lo <= 0.8 && 0.8 <= ci$hi))
set.seed(902)                                                                       # an irrelevant instrument: the whole line
zz <- rnorm(n); ivw <- ivmodel(Y = y1, D = y2, Z = cbind(zz), X = W)
ciw <- ar_ci(y1, y2, cbind(zz), cbind(1, W))
stopifnot(grepl("Whole Real Line", AR.test(ivw)$ci.info) == (ciw$type == "real_line"), isTRUE(ar_covers(ciw, 123)))
# a weak instrument whose AR set is two rays: points inside satisfy AR <= critical value, points outside do not
set.seed(903)
zw <- rnorm(300); vw <- rnorm(300); ew <- 0.9 * vw + sqrt(1 - 0.81) * rnorm(300)
xw <- 0.04 * zw + vw; yw <- xw + ew
cr <- ar_ci(yw, xw, cbind(zw), cbind(rep(1, 300)))
crit <- qf(0.95, 1, 298)
stats_at <- function(b0) ar_test(yw, xw, cbind(zw), cbind(rep(1, 300)), b0)$stat
stopifnot(cr$type == "two_rays",
          stats_at(cr$lo - 1) <= crit, stats_at(cr$hi + 1) <= crit, stats_at((cr$lo + cr$hi) / 2) > crit,
          !ar_covers(cr, (cr$lo + cr$hi) / 2), ar_covers(cr, cr$hi + 1))
for (b0 in seq(-5, 5, by = 0.25)) stopifnot(ar_covers(cr, b0) == (stats_at(b0) <= crit))   # the set is exactly {AR <= c}

## ---- Hansen J (robust, two-step GMM) and Sargan (homoskedastic) -------------------------------
J <- hansen_j(y1, X1, Z)
g2 <- gmm::gmm(y1 ~ y2 + z1 + w, ~ z1 + z2 + z3 + w, data = d, type = "twoStep", vcov = "MDS",
               centeredVcov = FALSE, crit = 1e-12)
jt <- gmm::specTest(g2)$test
stopifnot(J$df == 1, abs(J$J - jt[1, 1]) < 1e-6, abs(J$p - jt[1, 2]) < 1e-6, max(abs(J$b - coef(g2))) < 1e-6)
sg <- sargan(y1, X1, Z)
dg <- summary(iv1, diagnostics = TRUE)$diagnostics
stopifnot(sg$df == 1, abs(sg$stat - dg["Sargan", "statistic"]) < 1e-8, abs(sg$p - dg["Sargan", "p-value"]) < 1e-8)
stopifnot(is.na(sargan(y2, X2, Z)$stat) || sargan(y2, X2, Z)$df == 0)                # just identified: no test

## ---- first-stage F and Montiel Olea-Pflueger effective F ---------------------------------------
fs <- first_stage_F(y1, cbind(z1), cbind(1, z2, z3, w))
lmf <- lm(y1 ~ z2 + z3 + w + z1)
stopifnot(abs(fs$F - summary(lmf)$coefficients["z1", 3]^2) < 1e-8, fs$df1 == 1, fs$df2 == n - 5)
stopifnot(abs(eff_F(y1, cbind(z1), cbind(1, z2, z3, w)) -
                (coef(lmf)[["z1"]] / sqrt(sandwich::vcovHC(lmf, type = "HC0")["z1", "z1"]))^2) < 1e-8)
fs2 <- first_stage_F(y2, cbind(z2, z3), cbind(1, z1, w))                         # = ivreg's "Weak instruments" F
stopifnot(abs(fs2$F - dg["Weak instruments", "statistic"]) < 1e-8, fs2$df1 == 2,
          abs(eff_F(y2, cbind(z2, z3), cbind(1, z1, w), vcov = "homo") - fs2$F) < 1e-8)   # homoskedastic effective F = usual F

## ---- control-function (regression) form of the endogeneity test (Wooldridge 2010 eq. 6.15) -------
cf <- cf_test(y1, X1, Z, endog = "y2")
stopifnot(max(abs(cf$b[colnames(X1)] - t1$b)) < 1e-10, is.finite(cf$t), cf$p >= 0, cf$p <= 1)

## ---- population asymptotic variances of 2SLS and 3SLS -----------------------------------------
Q <- crossprod(Z) / n                                                       # sample moments as a stand-in population
Pi1 <- solve(crossprod(Z), crossprod(Z, X1)); Pi2 <- solve(crossprod(Z), crossprod(Z, X2))
Om <- matrix(c(1, 0.6, 0.6, 1.5), 2)
pa <- pop_avar_sys(list(Pi1, Pi2), Q, Om)
stopifnot(all(diag(pa$V_3sls)[1:4] <= diag(pa$V_2sls[[1]]) + 1e-12),
          all(eigen(pa$V_2sls_joint - pa$V_3sls, symmetric = TRUE)$values > -1e-10),
          max(abs(pa$V_2sls_joint[1:4, 1:4] - pa$V_2sls[[1]])) < 1e-12)
pd <- pop_avar_sys(list(Pi1, Pi2), Q, diag(c(1, 1.5)))
stopifnot(max(abs(pd$V_3sls - pd$V_2sls_joint)) < 1e-10, max(abs(pd$V_2sls_joint[1:4, 5:9])) == 0)   # diagonal Omega: no gain
# the formula is the probability limit of n times the sample 3SLS variance (same Pi, Q, Omega)
th_om <- three_sls(list(y1, y2), list(X1, X2), Z, Sigma = Om)
stopifnot(rel(th_om$V * n, pa$V_3sls) < 1e-8)

## ---- ci95 ----
stopifnot(isTRUE(all.equal(ci95(1, 0.5), c(1 - qnorm(0.975) * 0.5, 1 + qnorm(0.975) * 0.5))))
