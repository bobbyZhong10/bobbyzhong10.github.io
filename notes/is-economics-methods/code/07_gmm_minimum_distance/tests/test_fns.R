# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of gmm_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# The chapter card has no function-level assertions of its own; this file checks the four card
# functions (md_fit, gmm_fit, windmeijer_vcov, s_set) against closed forms, against the R
# packages gmm (two-step, iterative, CUE; vcov = "MDS" with the matching centeredVcov, never the
# default "HAC") and plm (pgmm's Windmeijer-corrected two-step variance), and against numerical
# derivatives (numDeriv).
# ------------------------------------------------------------------
source("gmm_fns.R")
suppressPackageStartupMessages({ library(gmm); library(plm) })
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

## =================================================================
## md_fit
## =================================================================
## ---- homogeneous h(theta) = theta * iota_r with diagonal V_pi: MD = inverse-variance weighted mean ----
set.seed(701)
r <- 8; pi_hat <- rnorm(r, -1.5, 0.3); se <- runif(r, 0.1, 0.5); w <- 1 / se^2
md <- md_fit(pi_hat, diag(se^2), function(th) rep(th, r), theta0 = 0)
stopifnot(identical(names(md)[1:4], c("theta", "vcov", "J", "df")), md$df == r - 1, md$converged,
          abs(md$theta - sum(w * pi_hat) / sum(w)) < 1e-10, abs(md$vcov - 1 / sum(w)) < 1e-12,
          abs(md$J - sum(w * (pi_hat - md$theta)^2)) < 1e-10,                  # J = Cochran's Q
          abs(md$p - pchisq(md$J, r - 1, lower.tail = FALSE)) < 1e-14)
## ---- linear h(theta) = H theta with a full V_pi: optimal MD is GLS -----------------------------------------
H <- cbind(1, seq_len(r) / r); A <- matrix(rnorm(r * r), r); V <- crossprod(A) / r + 0.1 * diag(r)
Vi <- solve(V)
md2 <- md_fit(pi_hat, V, function(th) drop(H %*% th), c(0, 0))
b_gls <- solve(t(H) %*% Vi %*% H, t(H) %*% Vi %*% pi_hat)
stopifnot(max(abs(md2$theta - b_gls)) < 1e-9, rel(md2$vcov, solve(t(H) %*% Vi %*% H)) < 1e-8,
          abs(md2$J - drop(t(pi_hat - H %*% b_gls) %*% Vi %*% (pi_hat - H %*% b_gls))) < 1e-8,
          max(abs(md2$H - H)) < 1e-8)
## ---- a non-optimal W gives the sandwich (H'WH)^-1 H'W V W H (H'WH)^-1 ---------------------------------------
md3 <- md_fit(pi_hat, V, function(th) drop(H %*% th), c(0, 0), W = diag(r))
B <- solve(crossprod(H))
stopifnot(max(abs(md3$theta - B %*% crossprod(H, pi_hat))) < 1e-9, rel(md3$vcov, B %*% t(H) %*% V %*% H %*% B) < 1e-8)
## ---- nonlinear h: same minimiser and criterion as optim; vcov (H'V^-1 H)^-1 at the estimate ------------------
h_nl <- function(th) th[1] * exp(th[2] * seq_len(r) / r)
V_nl <- V / 100
pi_nl <- h_nl(c(-1.2, 0.4)) + drop(t(chol(V_nl)) %*% rnorm(r))
md4 <- md_fit(pi_nl, V_nl, h_nl, c(-1, 0))
obj <- function(th) { e <- pi_nl - h_nl(th); drop(t(e) %*% solve(V_nl, e)) }
op <- optim(md4$theta + c(0.05, -0.05), obj, method = "BFGS", control = list(reltol = 1e-15, maxit = 2000))
Hn <- numDeriv::jacobian(h_nl, md4$theta)
stopifnot(md4$converged, md4$df == r - 2, max(abs(md4$theta - op$par)) < 1e-5, abs(md4$J - obj(md4$theta)) < 1e-10,
          md4$J <= op$value + 1e-10, rel(md4$vcov, solve(t(Hn) %*% solve(V_nl) %*% Hn)) < 1e-6)

## =================================================================
## gmm_fit against the gmm package (linear IV, heteroskedastic, overidentified)
## =================================================================
set.seed(702); n <- 400
Zs <- matrix(rnorm(n * 3), n); e <- rnorm(n); v <- 0.6 * e + rnorm(n)
p <- drop(Zs %*% c(0.5, 0.4, 0.3)) + v; y <- 1 + 0.8 * p + e * (1 + 0.6 * abs(Zs[, 1]))
dl <- list(y = y, X = cbind(1, p), Z = cbind(1, Zs))
dd <- data.frame(y = y, p = p, z1 = Zs[, 1], z2 = Zs[, 2], z3 = Zs[, 3])
W1 <- solve(crossprod(dl$Z) / n)                                                # first step = 2SLS, as in gmm's formula interface
tsls <- solve(t(dl$X) %*% dl$Z %*% W1 %*% t(dl$Z) %*% dl$X, t(dl$X) %*% dl$Z %*% W1 %*% t(dl$Z) %*% y)
D_lin <- -crossprod(dl$Z, dl$X) / n
for (cen in c(FALSE, TRUE)) {
  f2 <- gmm_fit(g_iv, c(0, 0), dl, "two_step", W1 = W1, centered = cen)
  ref2 <- gmm(y ~ p, ~ z1 + z2 + z3, data = dd, type = "twoStep", vcov = "MDS", centeredVcov = cen)
  S_fin <- moment_var(g_iv(f2$theta, dl), cen)
  stopifnot(identical(names(f2)[1:4], c("theta", "vcov", "J", "df")), f2$df == 2, f2$converged,
            max(abs(f2$theta1 - tsls)) < 1e-9, max(abs(f2$theta - coef(ref2))) < 1e-8,
            max(abs(f2$D - D_lin)) < 1e-8,                                     # D_0-hat = dgbar/dtheta' (here -Z'X/n)
            abs(f2$J - as.numeric(specTest(ref2)$test[1])) < 1e-6,
            # gmm reports (D'S(theta2)^-1 D)^-1 / n; gmm_fit reports the weight it used, (D'S(theta1)^-1 D)^-1 / n
            rel(solve(t(f2$D) %*% solve(S_fin) %*% f2$D) / n, unname(vcov(ref2))) < 1e-6,
            rel(f2$vcov, solve(t(f2$D) %*% f2$W %*% f2$D) / n) < 1e-12,
            rel(f2$W, solve(moment_var(g_iv(f2$theta1, dl), cen))) < 1e-12)
  fi <- gmm_fit(g_iv, c(0, 0), dl, "iterated", W1 = W1, centered = cen)
  refi <- gmm(y ~ p, ~ z1 + z2 + z3, data = dd, type = "iterative", vcov = "MDS", centeredVcov = cen,
              crit = 1e-13, itermax = 1000)
  stopifnot(fi$converged, max(abs(fi$theta - coef(refi))) < 1e-8, abs(fi$J - as.numeric(specTest(refi)$test[1])) < 1e-6,
            rel(fi$vcov, unname(vcov(refi))) < 1e-6)
  fc <- gmm_fit(g_iv, c(0, 0), dl, "cue", W1 = W1, centered = cen)
  refc <- gmm(y ~ p, ~ z1 + z2 + z3, data = dd, type = "cue", vcov = "MDS", centeredVcov = cen,
              t0 = f2$theta, method = "BFGS", control = list(reltol = 1e-15, maxit = 5000))
  stopifnot(fc$converged, max(abs(fc$theta - coef(refc))) < 1e-5,
            abs(fc$J - n * drop(crossprod(colMeans(g_iv(fc$theta, dl)), solve(moment_var(g_iv(fc$theta, dl), cen), colMeans(g_iv(fc$theta, dl)))))) < 1e-10)
  # no other value of theta gives a smaller CUE criterion than the optimum (local check)
  crit <- function(th) { G <- g_iv(th, dl); n * drop(crossprod(colMeans(G), solve(moment_var(G, cen), colMeans(G)))) }
  stopifnot(all(vapply(list(c(1e-4, 0), c(-1e-4, 0), c(0, 1e-4), c(0, -1e-4)), function(s) crit(fc$theta + s), 1) >= fc$J - 1e-12))
  assign(paste0("fc_", cen), fc)
}
## the CUE minimiser does not depend on centring (Sherman-Morrison: J_c = J_u / (1 - J_u / n))
stopifnot(max(abs(fc_TRUE$theta - fc_FALSE$theta)) < 1e-6, abs(fc_TRUE$J - fc_FALSE$J / (1 - fc_FALSE$J / n)) < 1e-6)
## identity weight (W1 = NULL): one step with W = I_r, robust sandwich via gmm_sandwich
f_id <- gmm_min(g_iv, c(0, 0), dl, diag(4))
ref_id <- gmm(y ~ p, ~ z1 + z2 + z3, data = dd, wmatrix = "ident", vcov = "MDS", centeredVcov = FALSE)
stopifnot(max(abs(f_id$theta - coef(ref_id))) < 1e-8,
          rel(gmm_sandwich(f_id$D, diag(4), moment_var(g_iv(f_id$theta, dl)), n), unname(vcov(ref_id))) < 1e-6)
## exactly identified: every type returns the IV estimate, J = 0, df = 0
dl_ji <- list(y = y, X = cbind(1, p), Z = cbind(1, Zs[, 1]))
f_ji <- gmm_fit(g_iv, c(0, 0), dl_ji, "two_step")
stopifnot(f_ji$df == 0, abs(f_ji$J) < 1e-12,
          max(abs(f_ji$theta - solve(crossprod(dl_ji$Z, dl_ji$X), crossprod(dl_ji$Z, y)))) < 1e-9)
## a nonlinear moment (exponential mean with instruments) against gmm's function interface, identity first step
set.seed(703)
xn <- rnorm(n); zn <- cbind(1, xn, xn^2); yn <- exp(0.3 + 0.5 * xn) * rexp(n)
dn <- list(y = yn, x = xn, Z = zn)
g_exp <- function(th, d) d$Z * (d$y - exp(th[1] + th[2] * d$x))
fn2 <- gmm_fit(g_exp, c(0, 0), dn, "two_step")
refn <- gmm(function(th, x) g_exp(th, dn), x = dn$Z, t0 = c(0, 0), type = "twoStep", vcov = "MDS", centeredVcov = FALSE,
            method = "BFGS", control = list(reltol = 1e-15, maxit = 5000))
stopifnot(fn2$converged, max(abs(fn2$theta - coef(refn))) < 1e-5)

## =================================================================
## windmeijer_vcov
## =================================================================
## ---- its derivative matrix is d theta2 / d theta1': the two-step estimate as a function of the point at which
## ---- the weight S(theta)^-1 is evaluated, differentiated numerically ----------------------------------------
for (cen in c(FALSE, TRUE)) {
  f2 <- gmm_fit(g_iv, c(0, 0), dl, "two_step", W1 = W1, centered = cen)
  wv <- windmeijer_vcov(g_iv, f2$theta1, f2$theta, dl, W1, centered = cen)
  theta2_of <- function(tb) {                                   # closed-form linear GMM with the weight S(tb)^-1
    W <- solve(moment_var(g_iv(tb, dl), cen)); XZ <- crossprod(dl$X, dl$Z)
    drop(solve(XZ %*% W %*% t(XZ), XZ %*% W %*% crossprod(dl$Z, dl$y)))
  }
  D_num <- numDeriv::jacobian(theta2_of, f2$theta1)
  stopifnot(max(abs(theta2_of(f2$theta1) - f2$theta)) < 1e-9, all(dim(wv) == c(2, 2)), rel(attr(wv, "dtheta2_dtheta1"), D_num) < 1e-6,
            rel(attr(wv, "naive"), f2$vcov) < 1e-9,           # numerical Jacobians: ~1e-11 apart
            rel(attr(wv, "V1"), gmm_sandwich(f2$D, W1, moment_var(g_iv(f2$theta1, dl), cen), n)) < 1e-9,
            isSymmetric(unname(wv), tol = 1e-10), all(eigen(wv, only.values = TRUE)$values > 0))
}
## ---- plm::pgmm (Arellano-Bond, two steps) reports the Windmeijer (2005) variance; rebuild it from pgmm's
## ---- own pieces: unit-level moments W_i'(dy_i - dX_i theta), first-step weight A1 ---------------------------
set.seed(704); N <- 150; Tt <- 7
id <- rep(seq_len(N), each = Tt); tt <- rep(seq_len(Tt), N)
a_i <- rnorm(N); xp <- rnorm(N * Tt); yp <- numeric(N * Tt)
for (i in seq_len(N)) {
  ix <- (i - 1) * Tt + seq_len(Tt)
  yp[ix[1]] <- a_i[i] + rnorm(1)
  for (s in 2:Tt) yp[ix[s]] <- 0.5 * yp[ix[s - 1]] + 0.4 * xp[ix[s]] + a_i[i] + rnorm(1) * (1 + 0.5 * abs(xp[ix[s]]))
}
pd <- pdata.frame(data.frame(id, tt, y = yp, x = xp), index = c("id", "tt"))
pf <- pgmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = pd, effect = "individual", model = "twosteps", transformation = "d")
dp <- list(W = pf$W, M = pf$model)
g_ab <- function(th, d) t(mapply(function(Wi, Mi) drop(crossprod(Wi, Mi[, 1] - Mi[, -1, drop = FALSE] %*% th)),
                                 d$W, d$M))
t1 <- unname(pf$coefficients[[1]]); t2 <- unname(pf$coefficients[[2]])
f_ab <- gmm_fit(g_ab, c(0, 0), dp, "two_step", W1 = pf$A1)
wv_ab <- windmeijer_vcov(g_ab, t1, t2, dp, pf$A1)
stopifnot(max(abs(f_ab$theta1 - t1)) < 1e-8, max(abs(f_ab$theta - t2)) < 1e-8,
          rel(f_ab$vcov, unname(pf$vcov)) < 1e-6,                              # naive two-step variance
          rel(unname(wv_ab), unname(vcovHC(pf))) < 1e-6,                        # Windmeijer-corrected variance
          all(diag(wv_ab) > diag(f_ab$vcov)))

## =================================================================
## s_set
## =================================================================
## concentrated scalar moments z_i (y~_i - theta p~_i) with mean-zero shocks z (the form used in blocks G2 and G3)
set.seed(705); n3 <- 600
z3 <- sapply(1:3, function(j) sample(rep(c(-1, 1), n3 / 2)))
u3 <- rnorm(n3); p3 <- drop(z3 %*% c(0.10, 0.06, 0.03)) + 0.5 * u3 + rnorm(n3); y3 <- 2 - 1.5 * p3 + u3
dc <- list(Zs = z3, yc = y3 - mean(y3), pc = p3 - mean(p3))
grid <- seq(-3, 0, by = 0.01)
ss <- s_set(g_conc, grid, dc)
S_direct <- vapply(grid, function(th) { G <- g_conc(th, dc); nrow(G) * drop(crossprod(colMeans(G), solve(crossprod(G) / nrow(G), colMeans(G)))) }, 1)
stopifnot(is.data.frame(ss), identical(names(ss), c("theta", "S", "p", "in_set")), nrow(ss) == length(grid),
          max(abs(ss$S - S_direct)) < 1e-9, max(abs(ss$p - pchisq(S_direct, 3, lower.tail = FALSE))) < 1e-12,
          identical(ss$in_set, S_direct <= qchisq(0.95, 3)),
          identical(s_set(g_conc, grid, dc, level = 0.90)$in_set, S_direct <= qchisq(0.90, 3)),
          any(ss$in_set), !all(ss$in_set), abs(s_stat(g_conc, -1.5, dc) - s_set(g_conc, -1.5, dc)$S) < 1e-12)
## CUE with the uncentred S_n(theta) minimises the S statistic: its J is the minimum of S over theta
fcc <- gmm_fit(g_conc, 0, dc, "cue")
stopifnot(abs(fcc$J - s_stat(g_conc, fcc$theta, dc)) < 1e-10, all(ss$S >= fcc$J - 1e-10), fcc$df == 2)

## =================================================================
## data helpers
## =================================================================
## market_slopes: per-market OLS slope of y on log price with the HC0 standard error (= sandwich::vcovHC "HC0")
set.seed(706)
dm <- data.frame(market = factor(rep(c("a", "b"), c(300, 450))), log_price = log(10) + 0.08 * rbinom(750, 1, 0.5))
dm$y <- 3 - 1.6 * dm$log_price + rnorm(750, 0, 0.8)
ms <- market_slopes(dm)
for (m in c("a", "b")) {
  lmf <- lm(y ~ log_price, data = dm[dm$market == m, ])
  stopifnot(abs(ms$pi_hat[ms$market == m] - coef(lmf)[2]) < 1e-10,
            abs(ms$se[ms$market == m] - sqrt(sandwich::vcovHC(lmf, type = "HC0")[2, 2])) < 1e-10)
}
## g_pool: the pooled moments (market dummies and dummies x log price) of block G1, linear in (a_1..a_r, theta)
dlp <- pool_data(dm)
b0 <- c(1, 2, -1); Gp <- g_pool(b0, dlp)
stopifnot(all(dim(Gp) == c(750, 4)), max(abs(colSums(Gp[, 1:2]) - tapply(dm$y - b0[as.integer(dm$market)] - b0[3] * dm$log_price, dm$market, sum))) < 1e-9)
