# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of bin_fns.R against packages (glm, marginaleffects, sandwich,
# survival::clogit, mlogit::hmftest, ordinal::clm, MASS::polr) and against numerical derivatives.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("bin_fns.R")
suppressPackageStartupMessages({ library(marginaleffects); library(survival) })
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

## ---- links -----------------------------------------------------------------------------------------
u <- seq(-4, 4, by = 0.5)
for (lk in c("probit", "logit")) {
  L <- link_fns(lk)
  stopifnot(max(abs(L$Fp(u) - numDeriv::grad(L$F, u))) < 1e-8, max(abs(L$Fpp(u) - numDeriv::grad(L$Fp, u))) < 1e-8)
}

## ---- ape_index: = marginaleffects::avg_slopes and = the average numerical derivative -------------------
set.seed(1201); n <- 1500
x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); X <- cbind("(Intercept)" = 1, x1 = x1, x2 = x2)
for (lk in c("probit", "logit")) {
  y <- rbinom(n, 1, link_fns(lk)$F(-0.3 + 0.8 * x1 + 0.5 * x2))
  g <- glm(y ~ x1 + x2, family = binomial(lk), data = data.frame(y, x1, x2), control = glm.control(epsilon = 1e-12, maxit = 100))
  a <- ape_index(coef(g), X, "x1", lk)
  ms <- avg_slopes(g, variables = "x1")
  stopifnot(abs(a - ms$estimate) < 1e-6, abs(a - ape_index(coef(g), X, 2, lk)) < 1e-15)
  num <- mean(numDeriv::grad(function(z) link_fns(lk)$F(drop(cbind(1, z, x2) %*% coef(g))), x1))
  stopifnot(abs(a - num) < 1e-8)
  # binary_mle = glm (coefficients, model-based variance, HC0 sandwich)
  bm <- binary_mle(y, X, lk)
  stopifnot(max(abs(bm$coef - coef(g))) < 1e-8, rel(bm$vcov_model, unname(vcov(g))) < 1e-6,
            rel(bm$vcov_sand, unname(sandwich::vcovHC(g, type = "HC0"))) < 1e-6,
            abs(bm$loglik - as.numeric(logLik(g))) < 1e-6)
  # avg_effect_inference: its delta-method part is marginaleffects' standard error
  ai <- avg_effect_inference(function(t) link_fns(lk)$Fp(drop(X %*% t)) * t[2], bm$coef, bm$vcov_model)
  stopifnot(abs(ai$est - a) < 1e-12, abs(ai$se_delta / ms$std.error - 1) < 1e-3,
            abs(ai$se_sample - sd(link_fns(lk)$Fp(drop(X %*% bm$coef)) * bm$coef[2]) / sqrt(n)) < 1e-12,
            abs(ai$se^2 - ai$se_delta^2 - ai$se_sample^2) < 1e-15, isTRUE(all.equal(ai$ci, ci95(ai$est, ai$se))))
}

## ---- cross_partial: = numerical cross derivative of the response probability ---------------------------
z1 <- rnorm(n); z2 <- 0.3 * z1 + rnorm(n)
Xi <- cbind("(Intercept)" = 1, z1 = z1, z2 = z2, "z1:z2" = z1 * z2)
b <- c(0.5, 0.9, 0.6, 0.25)
stopifnot(interaction_col(Xi, 2, 3) == 4, is.na(interaction_col(Xi[, 1:3], 2, 3)))
for (lk in c("probit", "logit")) {
  Fl <- link_fns(lk)$F
  P <- function(a1, a2) Fl(b[1] + b[2] * a1 + b[3] * a2 + b[4] * a1 * a2)
  h <- 1e-4
  num <- (P(z1 + h, z2 + h) - P(z1 + h, z2 - h) - P(z1 - h, z2 + h) + P(z1 - h, z2 - h)) / (4 * h^2)
  ci <- cross_partial_i(b, Xi, 2, 3, lk)
  stopifnot(max(abs(ci - num)) < 1e-6, abs(cross_partial(b, Xi, "z1", "z2", lk) - mean(num)) < 1e-7,
            abs(cross_partial(b, Xi, 2, 3, lk) - cross_partial(b, Xi, 3, 2, lk)) < 1e-15)
  # no interaction column: the cross-partial is beta_j beta_j' F''(u), not zero
  b3 <- b[1:3]
  stopifnot(abs(cross_partial(b3, Xi[, 1:3], 2, 3, lk) - mean(b3[2] * b3[3] * link_fns(lk)$Fpp(drop(Xi[, 1:3] %*% b3)))) < 1e-15,
            abs(cross_partial(b3, Xi[, 1:3], 2, 3, lk)) > 1e-3)
}

## ---- rivers_vuong: two steps by hand, the ASF double average, reproducible bootstrap --------------------
set.seed(1211); n3 <- 600
X3 <- rnorm(n3); Z3 <- rbinom(n3, 1, 0.5); v3 <- rnorm(n3)
D3 <- 0.2 * X3 + 0.9 * Z3 + v3; Y3 <- as.integer(-0.2 + 0.3 * X3 + 0.5 * D3 + 0.5 * v3 + sqrt(0.75) * rnorm(n3) > 0)
Xe <- cbind(1, X3)
f1 <- lm(D3 ~ X3 + Z3); vh <- resid(f1)
g2 <- glm(Y3 ~ X3 + D3 + vh, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
rv <- rivers_vuong(Y3, Xe, D3, Z3, B = 20, seed = 77)
cb <- coef(g2)
ape_loop <- 0
for (i in seq_len(n3)) ape_loop <- ape_loop + mean(cb[3] * dnorm(cb[1] + cb[2] * X3[i] + cb[3] * D3[i] + cb[4] * vh))
ape_loop <- ape_loop / n3
stopifnot(max(abs(unname(rv$beta) - unname(cb))) < 1e-8, max(abs(rv$pi - unname(coef(f1)))) < 1e-10,
          abs(rv$ape - ape_loop) < 1e-10, abs(rv$rho_t - summary(g2)$coefficients["vh", "z value"]) < 1e-4,
          rv$B == 20, nrow(rv$boot) == 20, abs(rv$ape_se - sd(rv$boot[, "ape"])) < 1e-15,
          isTRUE(all.equal(rv$ape_ci, ci95(rv$ape, rv$ape_se))))
rv2 <- rivers_vuong(Y3, Xe, D3, Z3, B = 20, seed = 77, cores = 2L)
stopifnot(identical(rv$boot, rv2$boot))                                   # same seed, same draws (also in parallel)
# the bootstrap resamples whole rows for both steps: draw 1 = the two steps on the resampled rows
set.seed(77); i1 <- sample.int(n3, n3 * 20, replace = TRUE)[1:n3]
stopifnot(max(abs(rv$boot[1, ] - c(rv_fit(Y3, Xe, D3, cbind(Z3), i1)$coef, rv_fit(Y3, Xe, D3, cbind(Z3), i1)$ape))) < 1e-12)

## ---- ivprobit_mle: first-order conditions and the structural scale ----------------------------------
iv <- ivprobit_mle(Y3, Xe, D3, Z3)
stopifnot(iv$convergence == 0, all(is.finite(iv$se_beta)), abs(iv$sigma_v - sd(vh)) < 0.05,
          abs(iv$r) < 1, iv$ape > 0)
# the scaled coefficients beta / sigma_e are close to the two-step ones (both consistent for them)
se_e <- sqrt(1 - iv$r^2)
stopifnot(max(abs(iv$beta / se_e - unname(cb[1:3]))) < 0.1)

## ---- cre_probit = pooled probit on (1, X, Xbar) with cluster-robust variance ------------------------------
set.seed(1212); N <- 300; T <- 5
id <- rep(seq_len(N), each = T); ci_ <- rnorm(N); Xp <- rep(ci_, each = T) + rnorm(N * T)
Yp <- as.integer(0.5 * rep(ci_, each = T) + 0.7 * Xp + rnorm(N * T) > 0)
cre <- cre_probit(Yp, Xp, id)
xb <- ave(Xp, id)
gp <- glm(Yp ~ Xp + xb, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
stopifnot(max(abs(unname(cre$beta) - unname(coef(gp)))) < 1e-8,
          abs(cre$ape - mean(dnorm(predict(gp, type = "link"))) * coef(gp)[2]) < 1e-8,
          rel(cre$vcov, unname(sandwich::vcovCL(gp, cluster = id, type = "HC0", cadjust = FALSE))) < 1e-6,
          names(cre$ape) == "x1", length(cre$ape_ci$x1) == 2)

## ---- fe_logit_cmle = survival::clogit (exact conditional likelihood); T = 2 is a logit on Delta X ---------
set.seed(1213)
al <- rnorm(N); Xf <- cbind(rnorm(N * T) + rep(al, each = T), rnorm(N * T))
Yf <- rbinom(N * T, 1, plogis(rep(al, each = T) + Xf %*% c(0.8, -0.4)))
cm <- fe_logit_cmle(Yf, Xf, id)
cl <- clogit(Yf ~ Xf[, 1] + Xf[, 2] + strata(id))
stopifnot(max(abs(unname(cm$coef) - unname(coef(cl)))) < 1e-6, rel(cm$vcov_model, unname(vcov(cl))) < 1e-5,
          all(eigen(cm$vcov_robust)$values > 0))
s2 <- id %in% which(tapply(Yf, id, sum) %in% 1:(T - 1))
stopifnot(cm$n_units_used == length(unique(id[s2])))
id2 <- rep(seq_len(N), each = 2); y2 <- matrix(Yf, T)[1:2, ]; x2m <- matrix(Xf[, 1], T)[1:2, ]
cm2 <- fe_logit_cmle(as.vector(y2), cbind(as.vector(x2m)), id2)
sw <- colSums(y2) == 1
g2t <- glm(y2[2, sw] ~ 0 + I(x2m[2, sw] - x2m[1, sw]), family = binomial(), control = glm.control(epsilon = 1e-12))
stopifnot(abs(cm2$coef - coef(g2t)) < 1e-6)                                  # Chamberlain (1980): T = 2

## ---- MNL and nested logit probabilities, APE, draws -------------------------------------------------
W <- cbind(1, rnorm(400), rnorm(400)); Bm <- rbind(c(0.3, 1, 0.2), c(-0.5, 0.4, 0.9))
Pm <- mnl_probs(Bm, W)
stopifnot(max(abs(rowSums(Pm) - 1)) < 1e-12, max(abs(nl_probs(Bm, W, 1) - Pm)) < 1e-12,
          max(abs(rowSums(nl_probs(Bm, W, 0.4)) - 1)) < 1e-12)
num <- colMeans(numDeriv::jacobian(function(z) as.vector(mnl_probs(Bm, cbind(1, z, W[, 3]))), W[, 2]) |>
                  (\(J) sapply(1:3, function(m) rowSums(J[(m - 1) * 400 + 1:400, , drop = FALSE])))())
stopifnot(max(abs(mnl_ape(Bm, W, 2) - num)) < 1e-7, abs(sum(mnl_ape(Bm, W, 2))) < 1e-14)
# nested logit: within the paid nest the odds are exp((V_Pro - V_Ent) / lambda)
Pn <- nl_probs(Bm, W, 0.4)
stopifnot(max(abs(log(Pn[, 2] / Pn[, 3]) - drop(W %*% (Bm[1, ] - Bm[2, ])) / 0.4)) < 1e-10)
set.seed(1214); Pd <- matrix(c(0.2, 0.5, 0.3), 2e5, 3, byrow = TRUE); dc <- draw_choice(Pd)
stopifnot(max(abs(tabulate(dc, 3) / 2e5 - c(0.2, 0.5, 0.3))) < 0.005)

## ---- hm_test = mlogit::hmftest (both the non-base and the base category dropped) -------------------------
set.seed(1215); n7 <- 1500
d7 <- data.frame(x = rnorm(n7), s = rnorm(n7))
P7 <- mnl_probs(rbind(c(0.5, 1, 0.3), c(0, 0.4, 1)), cbind(1, d7$x, d7$s))
d7$choice <- factor(c("Free", "Pro", "Enterprise")[draw_choice(P7)], levels = c("Free", "Pro", "Enterprise"))
long <- do.call(rbind, lapply(levels(d7$choice), function(a)
  data.frame(id = seq_len(n7), alt = a, chosen = d7$choice == a, x = d7$x, s = d7$s)))
long <- long[order(long$id, match(long$alt, levels(d7$choice))), ]
dl <- dfidx::dfidx(long, idx = c("id", "alt"), choice = "chosen")
suppressPackageStartupMessages(library(mlogit))
mf <- mlogit(chosen ~ 0 | x + s, data = dl, reflevel = "Free")
hE <- hm_test(choice ~ x + s, d7, "Enterprise")
# the restricted fits are passed as models: hmftest(mf, <alternatives>) re-evaluates the data by name
hm_ml <- hmftest(mf, mlogit(chosen ~ 0 | x + s, data = dl, reflevel = "Free", alt.subset = c("Free", "Pro")))
stopifnot(abs(hE$stat - unname(hm_ml$statistic)) < 1e-3 * max(1, abs(hE$stat)), hE$df == 3, hE$base == "Free")
mfP <- mlogit(chosen ~ 0 | x + s, data = dl, reflevel = "Pro")
hF <- hm_test(choice ~ x + s, d7, "Free")
hm_mlF <- hmftest(mfP, mlogit(chosen ~ 0 | x + s, data = dl, reflevel = "Pro", alt.subset = c("Pro", "Enterprise")))
stopifnot(abs(hF$stat - unname(hm_mlF$statistic)) < 1e-3 * max(1, abs(hF$stat)), hF$base == "Pro",
          abs(hF$p - pchisq(hF$stat, 3, lower.tail = FALSE)) < 1e-15)
# the multinom fit is the mlogit fit
stopifnot(max(abs(coef(hE$full)["Pro", ] - coef(mf)[c("(Intercept):Pro", "x:Pro", "s:Pro")])) < 1e-4)

## ---- ordered probit: probabilities, APE on the top category = marginaleffects on ordinal::clm -----------
suppressPackageStartupMessages(library(ordinal))
set.seed(1216); n6 <- 2000
p6 <- list(n = n6, beta = c(0.5, -0.3), delta_G = 0.2, sigma = c(1, 1.5), cuts = c(-1.9, -1.1, -0.2, 0.8),
           bill_median_kusd = 1.2, sd_logbill = 0.5, p_managed = 0.3)
d6 <- gen_B6(p6)
f6 <- clm(rating ~ x1 + bill_kusd + G, scale = ~ G, data = d6, link = "probit")
cf6 <- coef(f6)
idx6 <- drop(cbind(d6$x1, d6$bill_kusd, d6$G) %*% cf6[5:7]); sg6 <- exp(cf6[8] * d6$G)
pr6 <- oprobit_probs(idx6, cf6[1:4], sg6)
stopifnot(max(abs(rowSums(pr6) - 1)) < 1e-12, max(abs(pr6 - unname(predict(f6, newdata = d6[, c("x1", "bill_kusd", "G")], type = "prob")$fit))) < 1e-8)
ms6 <- avg_slopes(f6, variables = "x1")
stopifnot(abs(mean(oprobit_top_i(cf6["x1"], idx6, cf6["4|5"], sg6)) - ms6$estimate[ms6$group == "5"]) < 1e-5)
f6h <- clm(rating ~ x1 + bill_kusd + G, data = d6, link = "probit")
p6h <- MASS::polr(rating ~ x1 + bill_kusd + G, data = d6, method = "probit", Hess = TRUE)
stopifnot(max(abs(coef(f6h)[5:7] - coef(p6h))) < 1e-4, max(abs(coef(f6h)[1:4] - p6h$zeta)) < 1e-4)

## =====================================================================================================
## Fix round 1 (review of task 12.1): FE logit APE scaling, population-APE variance, stacked-score HM test
## =====================================================================================================

## ---- alpaca::getAPEs()$delta averages over ALL firm-months (nobs.full), never/always adopters included ----
suppressPackageStartupMessages(library(alpaca))
set.seed(1219); Nf <- 250; Tf <- 4
pf <- list(N = Nf, T = Tf, sigma_c = 0.7, sigma_a = 1, xi0 = -0.3, xi = 1, beta = 0.8)
dpn <- gen_B5(pf, N = Nf, T = Tf); dpn$fid <- factor(dpn$id)
fe5 <- alpaca::feglm(Y ~ X | fid, data = dpn, family = binomial("logit"))
bc5 <- alpaca::biasCorr(fe5)
a_fe <- alpaca::getAPEs(fe5)
nfull <- fe5$nobs[["nobs.full"]]
stopifnot(nfull == nrow(dpn), fe5$nobs[["nobs"]] < nfull)                        # some firms never or always adopt and are dropped
# the brief's identity: getAPEs()$delta = sum(beta-hat F'(eta-hat)) / nobs.full, to 1e-10
stopifnot(abs(unname(a_fe$delta) - sum(unname(coef(fe5)) * dlogis(fe5$eta)) / nfull) < 1e-10)
# ... hence it is NOT the mean over the kept rows: the movers-only APE is delta / share_kept
kept <- !(ave(dpn$Y, dpn$id) %in% c(0, 1))
stopifnot(sum(kept) == fe5$nobs[["nobs"]],
          abs(mean(unname(coef(fe5)) * dlogis(fe5$eta)) - unname(a_fe$delta) * nfull / fe5$nobs[["nobs"]]) < 1e-10)

## ---- fe_logit_ape: all-firm APE = delta, movers-only = delta / share_kept, population-APE standard error --
fa <- fe_logit_ape(fe5, bc5, id = dpn$id, y = dpn$Y)
a_bc <- alpaca::getAPEs(bc5)
stopifnot(abs(fa$ape_fe - unname(a_fe$delta)) < 1e-12, abs(fa$ape_bc - unname(a_bc$delta)) < 1e-12,
          abs(fa$share_kept - fe5$nobs[["nobs"]] / nfull) < 1e-15,
          abs(fa$ape_fe_movers - fa$ape_fe / fa$share_kept) < 1e-12,
          is.null(fa$ape_bc_movers),                                         # no movers-only bias-corrected APE (fix round 2)
          abs(fa$ape_fe_movers - mean(unname(coef(fe5)) * dlogis(fe5$eta))) < 1e-10,
          # the individual effects: zero for the dropped firms, beta-hat F'(alpha-hat_i + beta-hat X_it) for the others
          all(fa$g[!kept] == 0), abs(sum(fa$g) / nfull - fa$ape_fe) < 1e-10,
          max(abs(unname(coef(fe5)) * dlogis(alpaca::getFEs(fe5)$fid[as.character(dpn$id[kept])] + unname(coef(fe5)) * dpn$X[kept]) - fa$g[kept])) < 1e-10)
# standard errors: getAPEs' variance is conditional on the covariates; the population APE adds the firm-clustered
# sampling variance of the individual effects (the pattern of avg_effect_inference)
tot <- tapply(fa$g - mean(fa$g), dpn$id, sum)
stopifnot(abs(fa$se_sample - sqrt(sum(tot^2)) / nfull) < 1e-12,
          abs(fa$se_cond_fe - sqrt(a_fe$vcov[1, 1])) < 1e-12, abs(fa$se_cond_bc - sqrt(a_bc$vcov[1, 1])) < 1e-12,
          abs(fa$se_fe^2 - fa$se_cond_fe^2 - fa$se_sample^2) < 1e-15, abs(fa$se_bc^2 - fa$se_cond_bc^2 - fa$se_sample^2) < 1e-15,
          fa$se_fe > fa$se_cond_fe, fa$se_sample > 0, isTRUE(all.equal(fa$ci_bc, ci95(fa$ape_bc, fa$se_bc))))
# sampling_se_clustered is the sampling part of avg_effect_inference (with and without clusters)
gi <- rexp(240); cl_id <- rep(1:60, each = 4)
ai_c <- avg_effect_inference(function(t) gi * t, 1, matrix(0.01), cluster = cl_id)
ai_u <- avg_effect_inference(function(t) gi * t, 1, matrix(0.01))
stopifnot(abs(sampling_se_clustered(gi, cl_id) - ai_c$se_sample) < 1e-12, abs(sampling_se_clustered(gi) - ai_u$se_sample) < 1e-12)

## ---- hm_test_stacked: stacked-score (suest-type) Hausman statistic; PSD by construction -------------------
suppressPackageStartupMessages(library(sandwich))
hs <- hm_test_stacked(hE)                                                         # hE: drop "Enterprise", base Free
stopifnot(identical(hs$df, hE$df), hs$stat >= 0, hs$min_eigen > 0, isTRUE(all.equal(hs$p, pchisq(hs$stat, hs$df, lower.tail = FALSE))),
          max(abs(colSums(hs$scores$full))) < 1e-3, max(abs(colSums(hs$scores$restricted))) < 1e-3)
# scores = sandwich::estfun of nnet::multinom (full model) and of the restricted fit embedded in n rows
mfx <- nnet::multinom(choice ~ x + s, data = d7, trace = FALSE, reltol = 1e-12, maxit = 1000, Hess = TRUE)
rows_r <- which(d7$choice != "Enterprise")
dsx <- droplevels(d7[rows_r, ])
mrx <- nnet::multinom(choice ~ x + s, data = dsx, trace = FALSE, reltol = 1e-12, maxit = 1000, Hess = TRUE)
Ef <- sandwich::estfun(mfx); Er <- matrix(0, n7, 3, dimnames = list(NULL, colnames(sandwich::estfun(mrx)))); Er[rows_r, ] <- sandwich::estfun(mrx)
stopifnot(max(abs(hs$scores$full[, colnames(Ef)] - Ef)) < 1e-5, max(abs(hs$scores$restricted - Er)) < 1e-5)
# the statistic by an independent route (estfun, bread = n vcov): (b_r - b_f)' [Var(IF_r - IF_f)]^{-1} (b_r - b_f)
Vrx <- vcov(mrx); dimnames(Vrx) <- list(colnames(Er), colnames(Er)); Vfx <- vcov(mfx)
IFd <- Er %*% Vrx - (Ef %*% Vfx)[, colnames(Er)]
difx <- unname(coef(mrx)) - unname(coef(mfx)["Pro", ])
statx <- drop(t(difx) %*% solve(crossprod(IFd)) %*% difx)
stopifnot(abs(hs$stat - statx) < 1e-4 * max(1, statx))
# sandwich of the full model close to its model-based variance (information-matrix equality holds under IIA)
stopifnot(max(abs(diag(crossprod(Ef %*% Vfx)) / diag(Vfx) - 1)) < 0.25)
# per-observation analytic scores = numerical derivatives of the log likelihood (rows 1 to 5)
Wm <- cbind(1, d7$x, d7$s); Bm7 <- coef(mfx); lv <- levels(d7$choice)
for (i in 1:5) {
  ll_i <- function(b) log(mnl_probs(matrix(b, 2, 3, byrow = TRUE), Wm[i, , drop = FALSE])[match(as.character(d7$choice[i]), lv)])
  stopifnot(max(abs(numDeriv::grad(ll_i, as.vector(t(Bm7))) - hs$scores$full[i, colnames(Ef)])) < 1e-6)
}
# dropping the base category (re-levelled before both fits): same algebra
hsF <- hm_test_stacked(hF)
stopifnot(hsF$stat >= 0, hsF$min_eigen > 0, identical(hsF$df, hF$df))

## =====================================================================================================
## Fix round 2: exact B5 truths by two-dimensional Gauss-Hermite quadrature (controller ruling 44)
## =====================================================================================================
## (X_it, alpha_i) of gen_B5 is bivariate normal with a law that does not depend on t:
##   Var X = sigma_c^2 + 1, Cov(X, alpha) = xi (sigma_c^2 + 1/T), Var alpha = xi^2 (sigma_c^2 + 1/T) + sigma_a^2.
P5q <- list(sigma_c = 0.7, sigma_a = 1, xi0 = -0.3, xi = 1, beta = 0.8)
gh <- gauss_hermite_normal(40)
stopifnot(abs(sum(gh$w) - 1) < 1e-12, abs(sum(gh$w * gh$z)) < 1e-12, abs(sum(gh$w * gh$z^2) - 1) < 1e-12,
          abs(sum(gh$w * gh$z^4) - 3) < 1e-10)                                   # exact for polynomials of degree <= 79
q8 <- b5_quadrature(P5q, 8, 120); q8b <- b5_quadrature(P5q, 8, 60); q4 <- b5_quadrature(P5q, 4, 120)
stopifnot(max(abs(unlist(q8) - unlist(q8b))) < 1e-9,                              # stable in the number of nodes
          q8$pooled_plim > P5q$beta, q4$pooled_plim > q8$pooled_plim, q4$ape < q8$ape)
# independent route: trapezoid grid (normal_grid) in two independent standard normals and a Cholesky factor of the covariance matrix
ng2 <- normal_grid(301L, 8); g2 <- expand.grid(i = seq_along(ng2$z), j = seq_along(ng2$z)); w2 <- ng2$w[g2$i] * ng2$w[g2$j]
for (Tm in c(8, 4)) {
  Sg <- matrix(c(P5q$sigma_c^2 + 1, P5q$xi * (P5q$sigma_c^2 + 1 / Tm), P5q$xi * (P5q$sigma_c^2 + 1 / Tm),
                 P5q$xi^2 * (P5q$sigma_c^2 + 1 / Tm) + P5q$sigma_a^2), 2)
  Zs <- cbind(ng2$z[g2$i], ng2$z[g2$j]) %*% chol(Sg)                              # rows (X, alpha - xi0)
  X2 <- Zs[, 1]; al2 <- P5q$xi0 + Zs[, 2]
  pp <- plogis(al2 + P5q$beta * X2)
  fit <- optim(c(0, 1), function(b) -sum(w2 * (pp * plogis(b[1] + b[2] * X2, log.p = TRUE) + (1 - pp) * plogis(-(b[1] + b[2] * X2), log.p = TRUE))),
               method = "BFGS", control = list(reltol = 1e-15, maxit = 1000))
  qq <- if (Tm == 8) q8 else q4
  stopifnot(abs(qq$ape - sum(w2 * P5q$beta * dlogis(al2 + P5q$beta * X2))) < 1e-8, abs(qq$rate - sum(w2 * pp)) < 1e-8,
            abs(qq$pooled_plim - fit$par[2]) < 1e-5, abs(qq$pooled_int - fit$par[1]) < 1e-5,
            abs(qq$pooled_ape_plim - sum(w2 * dlogis(fit$par[1] + fit$par[2] * X2)) * fit$par[2]) < 1e-5)
}
