# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of ml_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# The card states no assertions for ml_fns.R; these check every card function and helper
# against packages (survival::survreg, maxLik, sandwich, stats::glm), against numerical
# derivatives (numDeriv) and against closed forms. Toy data use seeds 1 to 9 (not the
# chapter's block seeds, which are in 01_dgp.R).
# ------------------------------------------------------------------
source("ml_fns.R")
suppressPackageStartupMessages({ library(survival); library(maxLik) })
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

## ---- toy right-censored Weibull data ---------------------------------------------------
set.seed(1); n <- 3000; x <- rnorm(n); X <- cbind(1, x)
ts <- (rexp(n) / exp(-1 + 0.5 * x))^(1 / 1.4); cc <- runif(n, 0.5, 4)
t <- pmin(ts, cc); d <- as.numeric(ts <= cc)

## ---- nll_weibull_cens: the log likelihood survreg maximises (AFT -> PH conversion) -------
sv <- survreg(Surv(t, d) ~ x, dist = "weibull", control = survreg.control(rel.tolerance = 1e-12, maxiter = 100))
ph <- survreg_to_ph(sv)
stopifnot(length(ph$par) == 3, all(dim(ph$V) == c(3, 3)))
stopifnot(abs(-nll_weibull_cens(ph$par, t, d, X) - sv$loglik[2]) < 1e-6)          # same log likelihood (time scale)
qi <- nll_weibull_cens(ph$par, t, d, X, by_obs = TRUE)
stopifnot(length(qi) == n, abs(sum(qi) - nll_weibull_cens(ph$par, t, d, X)) < 1e-8)
# by hand for one censored and one uncensored observation: q = -[delta (x'b + log k + (k - 1) log t) - exp(x'b) t^k]
p0 <- c(-1, 0.5, log(1.4)); i1 <- which(d == 1)[1]; i0 <- which(d == 0)[1]
q_hand <- function(i) -(d[i] * (sum(X[i, ] * p0[1:2]) + p0[3] + (1.4 - 1) * log(t[i])) - exp(sum(X[i, ] * p0[1:2])) * t[i]^1.4)
stopifnot(abs(nll_weibull_cens(p0, t, d, X, by_obs = TRUE)[i1] - q_hand(i1)) < 1e-12,
          abs(nll_weibull_cens(p0, t, d, X, by_obs = TRUE)[i0] - q_hand(i0)) < 1e-12)

## ---- weibull_parts: analytic scores and Hessian = numDeriv ---------------------------------
wp <- weibull_parts(p0, t, d, X)
stopifnot(all(dim(wp$S) == c(n, 3)), all(dim(wp$H) == c(3, 3)), max(abs(wp$q - nll_weibull_cens(p0, t, d, X, by_obs = TRUE))) < 1e-12)
stopifnot(rel(wp$S, numDeriv::jacobian(function(p) nll_weibull_cens(p, t, d, X, by_obs = TRUE), p0)) < 1e-6,
          rel(wp$H, numDeriv::hessian(function(p) nll_weibull_cens(p, t, d, X), p0)) < 1e-6)

## ---- fit_weibull: maxLik and the Newton engine reach the survreg optimum ---------------------
fm <- fit_weibull(t, d, X)                                   # maxLik, method "NR"
fn <- fit_weibull(t, d, X, engine = "newton")
stopifnot(fm$converged, fn$converged, max(abs(fm$par - ph$par)) < 1e-5, max(abs(fn$par - fm$par)) < 1e-8,
          fm$grad_max < 1e-6, fn$grad_max < 1e-6)
stopifnot(rel(fm$V_hess, ph$V) < 1e-3)                        # survreg's vcov, carried to the PH scale by the delta method
# the maxLik object: maxLik::hessian (written with the package prefix, numDeriv masks it) and sandwich
stopifnot(inherits(fm$fit, "maxLik"), rel(-maxLik::hessian(fm$fit), weibull_parts(fm$par, t, d, X)$H) < 1e-6,
          rel(sandwich::sandwich(fm$fit), fm$V_sand) < 1e-6, rel(vcov(fm$fit), fm$V_hess) < 1e-6)
# shape fixed at kappa = 1: the exponential model, as survreg(dist = "exponential") fits it
fe <- fit_weibull(t, d, X, log_kappa = 0)
se_exp <- survreg(Surv(t, d) ~ x, dist = "exponential", control = survreg.control(rel.tolerance = 1e-12))
stopifnot(fe$converged, fe$par[3] == 0, max(abs(fe$par[1:2] - (-coef(se_exp)))) < 1e-5)

## ---- info_hat (numDeriv) = analytic A and B; sandwich pieces ---------------------------------
ih <- info_hat(fm$par, nll_weibull_cens, t = t, delta = d, X = X)
stopifnot(identical(names(ih), c("A", "B")), all(dim(ih$A) == c(3, 3)), all(dim(ih$B) == c(3, 3)))
stopifnot(rel(ih$A, fm$A) < 1e-6, rel(ih$B, fm$B) < 1e-6,
          rel(fm$A, weibull_parts(fm$par, t, d, X)$H / n) < 1e-12,
          rel(fm$B, crossprod(weibull_parts(fm$par, t, d, X)$S) / n) < 1e-12,
          rel(fm$V_hess, solve(fm$A) / n) < 1e-10, rel(fm$V_opg, solve(fm$B) / n) < 1e-10,
          rel(fm$V_sand, solve(fm$A) %*% fm$B %*% solve(fm$A) / n) < 1e-10)
# correctly specified model: information matrix equality A = B up to sampling error
stopifnot(max(abs(sqrt(diag(fm$V_hess) / diag(fm$V_opg)) - 1)) < 0.1)

## ---- trinity -----------------------------------------------------------------------------
# (a) Gaussian linear model with sigma = 1 known: the log likelihood is exactly quadratic, so
#     W = LR = LM = RSS_r - RSS_u (Wooldridge 2010 section 12.6: the three coincide for a quadratic objective)
set.seed(2); nl <- 400; Xl <- cbind(1, rnorm(nl), rnorm(nl)); yl <- drop(Xl %*% c(1, 0.3, 0.15)) + rnorm(nl)
nll_lin <- function(par, y, X, by_obs = FALSE) { q <- (y - drop(X %*% par))^2 / 2; if (by_obs) q else sum(q) }
bu <- qr.coef(qr(Xl), yl); br <- c(qr.coef(qr(Xl[, 1:2]), yl), 0)
rss <- function(b) sum((yl - drop(Xl %*% b))^2)
tl <- trinity(bu, br, nll_lin, function(p) p[3], y = yl, X = Xl)
# (W and LM go through numDeriv's Hessian and gradient: equal to LR up to numerical-derivative error)
stopifnot(identical(names(tl), c("W", "LR", "LM", "df")), tl$df == 1,
          abs(tl$LR - (rss(br) - rss(bu))) < 1e-8, abs(tl$W / tl$LR - 1) < 1e-5, abs(tl$LM / tl$LR - 1) < 1e-5)
# two restrictions: df = 2
br2 <- c(mean(yl), 0, 0)
tl2 <- trinity(bu, br2, nll_lin, function(p) p[2:3], y = yl, X = Xl)
stopifnot(tl2$df == 2, abs(tl2$W / tl2$LR - 1) < 1e-5, abs(tl2$LM / tl2$LR - 1) < 1e-5, abs(tl2$LR - (rss(br2) - rss(bu))) < 1e-8)
# (b) Weibull, H0: kappa = 1 (log kappa = 0): the three statistics by hand from the analytic pieces
tw <- trinity(fm$par, fe$par, nll_weibull_cens, function(p) p[3], t = t, delta = d, X = X)
g_r <- colSums(weibull_parts(fe$par, t, d, X)$S); H_r <- weibull_parts(fe$par, t, d, X)$H
stopifnot(abs(tw$LR - 2 * (fe$nll - fm$nll)) < 1e-8,
          abs(tw$W / (fm$par[3]^2 / solve(weibull_parts(fm$par, t, d, X)$H)[3, 3]) - 1) < 1e-5,
          abs(tw$LM / drop(t(g_r) %*% solve(H_r, g_r)) - 1) < 1e-5, tw$df == 1)

## ---- person_period and the grouped-duration likelihood ----------------------------------------
# toy: exit in (0, 1]; exit in (2, 3]; censored at 3; beyond the last cut (censored at 4); censored at 2.5 (off a cut)
pp <- person_period(c(0.5, 2.3, 3, 5.5, 2.5), c(1, 1, 0, 1, 0), cbind(fee_usd = c(10, 20, 30, 40, 50)), 0:4)
stopifnot(is.data.frame(pp), identical(names(pp), c("id", "m", "y", "fee_usd")),
          identical(pp$id, c(1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 5L, 5L)),
          identical(pp$m, c(1L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 1L, 2L)),
          identical(pp$y, c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)), all(pp$fee_usd == c(10, 20, 20, 20, 30, 30, 30, 40, 40, 40, 40, 50, 50)))
gd <- grouped_data(c(0.5, 2.3, 3, 5.5, 2.5), c(1, 1, 0, 1, 0), 0:4)
stopifnot(identical(gd$mi, c(1L, 3L, 3L, 4L, 2L)), identical(gd$yi, c(1, 1, 0, 0, 0)))
# grouped_nll (Wooldridge 2010 eq. 22.45 with eq. 22.47) = minus the Bernoulli cloglog log likelihood of the
# person-period rows, at any parameter value (the algebraic identity behind "pooled cloglog = grouped MLE")
set.seed(3); ng <- 1500; xg <- cbind(fee_usd = sample(c(40, 60, 80), ng, TRUE), referred = rbinom(ng, 1, 0.3))
tg <- (rexp(ng) / exp(-2.5 + 0.02 * xg[, 1] - 0.4 * xg[, 2]))^(1 / 0.9); cg <- sample(3:8, ng, TRUE)
dg <- as.numeric(tg <= cg); tg <- pmin(tg, cg)
gdat <- grouped_data(tg, dg, 0:6); ppg <- person_period(tg, dg, xg, 0:6)
for (k in 1:3) {
  par <- c(rnorm(6, -1.5, 0.3), rnorm(2, 0, 0.02))
  eta <- par[ppg$m] + drop(as.matrix(ppg[, c("fee_usd", "referred")]) %*% par[7:8])
  ll_rows <- sum(dbinom(ppg$y, 1, -expm1(-exp(eta)), log = TRUE))
  stopifnot(abs(grouped_nll(par, gdat$mi, gdat$yi, xg) + ll_rows) < 1e-8)
}
# analytic gradient and Hessian of grouped_nll = numDeriv
gm <- grouped_mle(gdat$mi, gdat$yi, xg, M = 6)
stopifnot(gm$converged, length(gm$par) == 8, max(abs(gm$grad)) < 1e-8,
          max(abs(grouped_grad(gm$par, gdat$mi, gdat$yi, xg) - gm$grad)) < 1e-12)
pt <- gm$par + 0.05
stopifnot(rel(grouped_grad(pt, gdat$mi, gdat$yi, xg), numDeriv::grad(function(p) grouped_nll(p, gdat$mi, gdat$yi, xg), pt)) < 1e-6,
          rel(grouped_hess(pt, gdat$mi, gdat$yi, xg), numDeriv::hessian(function(p) grouped_nll(p, gdat$mi, gdat$yi, xg), pt)) < 1e-6)
# the grouped MLE = pooled cloglog with interval dummies (glm)
gl <- glm(y ~ 0 + factor(m) + fee_usd + referred, family = binomial(link = "cloglog"), data = ppg,
          control = glm.control(epsilon = 1e-14, maxit = 500))
stopifnot(gl$converged, max(abs(unname(coef(gl)) - gm$par)) < 1e-6)

## ---- ols_weibull: the log-linear regression implied by the Weibull model ---------------------------
# complete data: log t = -x'beta / kappa + eps / kappa, eps standard minimum-Gumbel (mean -Euler, var pi^2 / 6)
set.seed(4); no <- 200000; xo <- cbind(1, sample(c(45, 60, 75, 90), no, TRUE), rbinom(no, 1, 0.3))
to <- (rexp(no) / exp(drop(xo %*% c(-3.4, 0.02, -0.4))))^(1 / 0.9)
ow <- ols_weibull(log(to), xo)
stopifnot(abs(ow$kappa - 0.9) < 0.01, max(abs(ow$beta - c(-3.4, 0.02, -0.4)) / c(1, 0.02, 0.4)) < 0.05,
          all(dim(ow$V_beta) == c(3, 3)), length(ow$se_beta) == 3)
# delta-method variance = numerical Jacobian of the map (b, s2) -> beta applied to the stacked sandwich
map <- function(th) { k <- pi / sqrt(6 * th[4]); c(-k * th[1] - 0.5772156649015329, -k * th[2:3]) }
Jm <- numDeriv::jacobian(map, c(ow$b, ow$s2))
stopifnot(rel(ow$V_beta, Jm %*% ow$V_bs2 %*% t(Jm)) < 1e-6)
# stacked sandwich for (b, s2): its b block is the HC0 OLS variance
stopifnot(rel(ow$V_bs2[1:3, 1:3], unname(sandwich::vcovHC(lm(log(to) ~ xo[, 2] + xo[, 3]), type = "HC0"))) < 1e-8)

## ---- two-step: two_step_vcov and twostep_weibull ---------------------------------------------------
set.seed(5); S <- matrix(rnorm(600), 200); A <- crossprod(matrix(rnorm(9), 3)) + diag(3)
stopifnot(rel(two_step_vcov(S, A, matrix(0, 3, 2), matrix(rnorm(400), 200)), solve(A) %*% crossprod(S) %*% solve(A) / 200^2) < 1e-12)
P4t <- list(months = 12, growth = 1, fees = c(45, 60, 75, 90), p_fee = c(0.3, 0.35, 0.2, 0.15),
            hh_sizes = 1:5, p_hh = c(0.15, 0.35, 0.2, 0.2, 0.1), gamma = c(-1, 0.25, 0.35),
            beta = c(-3.7, 0.02, -0.8), kappa = 0.9)
set.seed(6); d4 <- gen_C4(P4t, n = 1500)
f4 <- twostep_weibull(d4)
stopifnot(f4$converged, length(f4$gamma) == 3, length(f4$par) == 4, all(dim(f4$F2) == c(4, 3)), all(dim(f4$R1) == c(1500, 3)))
# first step = glm Poisson; second step = fit_weibull on the generated index
g4 <- glm(usage_wk1 ~ hh_size + quiz_std, family = poisson(), data = d4, control = glm.control(epsilon = 1e-12))
stopifnot(max(abs(f4$gamma - unname(coef(g4)))) < 1e-7)
# stacked (joint) M-estimator: psi_i = (first-step Poisson scores, second-step Weibull scores); theta block of the
# joint sandwich G^-1 Omega G^-1' / n = two_step_vcov with A2 = G22, F2 = G21, R1 = -s1 G11^-1'
Z4 <- cbind(1, d4$hh_size, d4$quiz_std); n4 <- nrow(Z4)
psi <- function(phi) {
  g <- phi[1:3]; th <- phi[4:7]; mu <- exp(drop(Z4 %*% g))
  cbind(-Z4 * (d4$usage_wk1 - mu), weibull_parts(th, d4$t_months, d4$delta, cbind(1, d4$fee_usd, drop(Z4 %*% g)))$S)
}
phi_hat <- c(f4$gamma, f4$par); P <- psi(phi_hat)
stopifnot(max(abs(colMeans(P))) < 1e-8)
G <- numDeriv::jacobian(function(p) colMeans(psi(p)), phi_hat)
stopifnot(max(abs(G[1:3, 4:7])) < 1e-8)
Om <- crossprod(P) / n4; Gi <- solve(G); V_joint <- Gi %*% Om %*% t(Gi) / n4
stopifnot(rel(f4$F2, G[4:7, 1:3]) < 1e-6, rel(f4$R1, -P[, 1:3] %*% t(solve(G[1:3, 1:3]))) < 1e-6,
          rel(f4$V_corrected, V_joint[4:7, 4:7]) < 1e-5,
          rel(f4$V_naive, solve(G[4:7, 4:7]) / n4) < 1e-5)
stopifnot(f4$V_corrected[3, 3] > f4$V_naive[3, 3])                # beta_u != 0: the first step adds variance

## ---- block C2 helpers: log-logistic mixture and the quadrature behind the pseudo-true value -----------
P2t <- list(months = 36, growth = 1.04, fees = c(45, 60, 75, 90), p_fee = c(0.3, 0.35, 0.2, 0.15), p_ref = 0.3,
            g0 = -0.5, g_ref = -1, med_deal = 1.5, p_deal = 0.5, med_reg = 12, p_reg = 5, theta_fee = 0.02, fee_ref = 60)
# survival, density and hazard of the mixture agree with each other
tt <- c(0.1, 1, 5, 12, 30)
stopifnot(max(abs(c2_dens(tt, 60, 0, P2t) - sapply(tt, function(s0) numDeriv::grad(function(s) -c2_surv(s, 60, 0, P2t), s0)))) < 1e-8,
          max(abs(c2_hazard(tt, 75, 1, P2t) - c2_dens(tt, 75, 1, P2t) / c2_surv(tt, 75, 1, P2t))) < 1e-12)
cobj <- c2_cells(P2t)
cl <- cobj$cells
stopifnot(nrow(cl) == 4 * 2 * 36, abs(sum(cl$w) - 1) < 1e-12)
i7 <- which(cl$fee == 60 & cl$ref == 0 & cl$c == 7)
g0u <- function(u) exp(u) * c2_dens(exp(u), 60, 0, P2t)          # density of u = log t* (fee 60, not referred)
stopifnot(abs(cl$F[i7] - (1 - c2_surv(7, 60, 0, P2t))) < 1e-10,
          abs(cl$L[i7] - integrate(function(u) u * g0u(u), -300, log(7), rel.tol = 1e-12)$value) < 1e-8)
pm <- c2_power_moments(cobj, 0.8)
stopifnot(abs(pm$P[i7] - integrate(function(u) exp(0.8 * u) * g0u(u), -300, log(7), rel.tol = 1e-12)$value) < 1e-8,
          abs(pm$Plog[i7] - integrate(function(u) exp(0.8 * u) * u * g0u(u), -300, log(7), rel.tol = 1e-12)$value) < 1e-8)
# analytic gradient of the population objective = numDeriv
pz <- c(-3, 0.018, -0.17, log(0.76))
stopifnot(rel(pop_nll_weibull_grad(pz, cobj), numDeriv::grad(function(p) pop_nll_weibull(p, cobj), pz)) < 1e-6)

## ---- vuong and the log-logistic log likelihood by observation ------------------------------------------
sl <- survreg(Surv(t, d) ~ x, dist = "loglogistic", control = survreg.control(rel.tolerance = 1e-12))
lli <- loglogistic_ll_obs(sl, t, d, X)
stopifnot(length(lli) == n, abs(sum(lli) - sl$loglik[2]) < 1e-6)
dl <- -nll_weibull_cens(fm$par, t, d, X, by_obs = TRUE) - lli
vu <- vuong(-nll_weibull_cens(fm$par, t, d, X, by_obs = TRUE), lli)
stopifnot(abs(vu$stat - sqrt(n) * mean(dl) / sqrt(mean((dl - mean(dl))^2))) < 1e-10,
          abs(vu$p_two_sided - 2 * pnorm(-abs(vu$stat))) < 1e-14)

## ---- ci95 ----------------------------------------------------------------------------------------
stopifnot(isTRUE(all.equal(ci95(1, 0.5), c(1 - qnorm(0.975) * 0.5, 1 + qnorm(0.975) * 0.5))))
