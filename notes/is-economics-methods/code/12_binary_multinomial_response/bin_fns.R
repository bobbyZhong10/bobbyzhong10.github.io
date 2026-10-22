# ------------------------------------------------------------------
# bin_fns.R -- functions for chapter 12 (binary and multinomial response models).
#
# Notation (notation_v3.md 1.10; the code keeps short legacy names, the comments give the symbols):
#   y         binary response Y_i (Helix adoption; the legacy data call it D)
#   X         n x k design matrix with rows X_i' (constant included where the model has one)
#   beta      index coefficients beta; X %*% beta is the index X_i'beta
#   link      F(.): "probit" (Phi) or "logit" (Lambda(.)); its density is F'(.), its second derivative F''(.)
#   j, k      regressor positions j and j' (the card's argument k is j' in the text)
#   x_endog   the endogenous regressor D_i of section 15.7.2 (coefficient beta_D)
#   Z         excluded instruments Z_i; first stage D_i = X_i'pi_1 + Z_i'pi_2 + v_i
#   id        unit index i of a balanced panel with T periods
#   B (MNL)   M x k matrix whose row m is beta_(m)' (category 0 is the base and has no row)
#   lambda    nested-logit dissimilarity parameter of the paid nest (lambda = 1: MNL, IIA holds)
#   V         an estimated variance matrix of the parameters (Avar-hat / n)
# Wooldridge (2010) chapters 15 and 16 have no numbered theorems: comments cite equation numbers.
# ------------------------------------------------------------------

## ---- links: F, F' and F'' ---------------------------------------------------------------------
link_fns <- function(link = c("probit", "logit")) {
  link <- match.arg(link)
  if (link == "probit") {
    list(F = pnorm, Fp = dnorm, Fpp = function(u) -u * dnorm(u))
  } else {
    list(F = plogis, Fp = dlogis, Fpp = function(u) { p <- plogis(u); p * (1 - p) * (1 - 2 * p) })
  }
}

col_index <- function(X, j) {
  if (is.character(j)) {
    w <- match(j, colnames(X))
    if (anyNA(w)) stop("column ", j[is.na(w)][1], " not in X")
    return(w)
  }
  j
}

#' APE of regressor j in the index model P(Y = 1 | X) = F(X'beta) (Wooldridge 2010 eq. 15.32):
#' the sample average of F'(X_i'beta) beta_j. Valid when x_j enters the index only through
#' beta_j x_j (no interaction or power of x_j in X; use cross_partial's ingredients otherwise).
#' beta: beta; X: rows X_i'; j: position or column name of x_j.
ape_index <- function(beta, X, j, link = c("probit", "logit")) {
  L <- link_fns(link)
  X <- as.matrix(X); j <- col_index(X, j)
  u <- drop(X %*% beta)
  mean(L$Fp(u)) * unname(beta[j])
}

#' Column of X that equals x_j * x_j' (the interaction regressor of beta_jj'), or NA.
interaction_col <- function(X, j, k) {
  prod <- X[, j] * X[, k]
  hit <- which(vapply(seq_len(ncol(X)), function(c)
    c != j && c != k && max(abs(X[, c] - prod)) <= 1e-12 * max(1, max(abs(prod))), logical(1)))
  if (length(hit) == 0L) NA_integer_ else hit[1]
}

#' Per-observation cross-partial effect d^2 P(x) / (d x_j d x_j') at x = X_i (Ai and Norton 2003;
#' the formula as printed in Norton, Wang and Ai 2004, p.155):
#'   beta_jj' F'(u) + (beta_j + beta_jj' x_j')(beta_j' + beta_jj' x_j) F''(u),  u = X_i'beta.
#' beta_jj' is the coefficient of the column of X equal to x_j x_j' (0 when there is none:
#' the cross-partial is then beta_j beta_j' F''(u), not zero).
cross_partial_i <- function(beta, X, j, k, link = c("probit", "logit")) {
  L <- link_fns(link)
  X <- as.matrix(X); j <- col_index(X, j); k <- col_index(X, k)
  jk <- interaction_col(X, j, k)
  b_jk <- if (is.na(jk)) 0 else unname(beta[jk])
  u <- drop(X %*% beta)
  b_jk * L$Fp(u) + (beta[j] + b_jk * X[, k]) * (beta[k] + b_jk * X[, j]) * L$Fpp(u)
}

#' Sample average cross-partial effect (the chapter card's cross_partial(beta, X, j, k, link)).
cross_partial <- function(beta, X, j, k, link = c("probit", "logit")) {
  mean(cross_partial_i(beta, X, j, k, link))
}

#' Inference for a sample-average effect mean_i g_i(theta) whose target is the population average
#' E[g(X_i, theta_0)]. g: function(theta) returning the n per-observation effects; theta: estimate;
#' V: its estimated variance. Two pieces (asymptotically uncorrelated when E[s_i | X_i] = 0):
#'   se_delta   delta method, sqrt(G V G'), G = d mean_i g_i / d theta' (numDeriv): the variance
#'              conditional on the X_i, as marginaleffects reports it;
#'   se_sample  sd of g_i / sqrt(n): the sampling variation of the X_i. With `cluster`, the sd of the
#'              cluster means of g_i - mean(g) scaled for clusters (panel data).
#' se = sqrt(se_delta^2 + se_sample^2) is the standard error for the population APE; ci = est -/+ 1.96 se.
avg_effect_inference <- function(g, theta, V, cluster = NULL) {
  gi <- g(theta)
  est <- mean(gi)
  G <- numDeriv::grad(function(t) mean(g(t)), theta)
  se_delta <- sqrt(drop(t(G) %*% V %*% G))
  if (is.null(cluster)) {
    se_sample <- sd(gi) / sqrt(length(gi))
  } else {
    tot <- rowsum(gi - est, cluster)
    se_sample <- sqrt(sum(tot^2)) / length(gi)
  }
  se <- sqrt(se_delta^2 + se_sample^2)
  list(est = est, se = se, se_delta = se_delta, se_sample = se_sample,
       ci = c(est - qnorm(0.975) * se, est + qnorm(0.975) * se), grad = G)
}

ci95 <- function(est, se) c(est - qnorm(0.975) * se, est + qnorm(0.975) * se)

#' The sampling part of a population-average effect's standard error: sd of the effects divided by sqrt(n)
#' (with `cluster`: the firm-clustered version, sqrt(sum_c (sum_{i in c} (g_i - mean g))^2) / n). This is
#' se_sample of avg_effect_inference, for effects g that are not a function of a parameter vector (the
#' fixed effects logit APE of fe_logit_ape).
sampling_se_clustered <- function(g, cluster = NULL) {
  if (is.null(cluster)) return(sd(g) / sqrt(length(g)))
  tot <- rowsum(g - mean(g), cluster)
  sqrt(sum(tot^2)) / length(g)
}

#' APE of the FE logit MLE fe (alpaca::feglm(Y ~ X | fid, binomial("logit")), one regressor) and of its
#' bias correction bc = alpaca::biasCorr(fe), on the panel with firm index `id` and outcome `y` (rows as in the
#' fit). alpaca::getAPEs()$delta is sum_{kept rows} beta-hat F'(eta-hat) / nobs.full: it already averages over
#' ALL firm-months, the firms that never or always adopt (alpha-hat = -/+ infinity, dropped by feglm) counting
#' as zero effects. So the all-firm APE is delta itself (ape_fe, ape_bc: the chapter's target). It is nearly
#' unbiased in the chapter's designs (-0.7% and -0.2% for the FE MLE at T = 4 and 8): two large errors offset,
#' the kept firms' overstated effects (incidental-parameter bias of the coefficient) and the dropped firms'
#' zero estimated effects. The APE among the firms that move is delta / share_kept for the FE MLE
#' (ape_fe_movers, exact; share_kept = nobs / nobs.full). There is NO movers-only version for the bias-corrected
#' fit: getAPEs(bc)$delta is the mean over all rows minus alpaca's APE bias term, so dividing it by share_kept
#' also rescales that term and is not a consistent movers-only estimator. Standard errors: getAPEs' variance is
#' conditional on the covariates (n.pop = NULL); for the POPULATION APE the sampling variance of the individual
#' effects g_it (zero for dropped firms, beta-hat F'(alpha-hat_i + beta-hat X_it) for the others, from the
#' uncorrected fit) is added, firm-clustered, as avg_effect_inference does in the cross-section:
#' se^2 = se_cond^2 + se_sample^2 (the same se_sample for fe and bc). Returns ape_fe, ape_bc, ape_fe_movers,
#' share_kept, g, se_cond_fe, se_cond_bc, se_sample, se_fe, se_bc, ci_fe, ci_bc.
fe_logit_ape <- function(fe, bc, id, y) {
  kept <- !(ave(y, id) %in% c(0, 1))                                  # firms that move: 0 < sum_t Y_it < T
  stopifnot(sum(kept) == fe$nobs[["nobs"]], length(y) == fe$nobs[["nobs.full"]])
  g <- numeric(length(y)); g[kept] <- unname(coef(fe)) * dlogis(fe$eta)
  a_fe <- alpaca::getAPEs(fe); a_bc <- alpaca::getAPEs(bc)
  share_kept <- fe$nobs[["nobs"]] / fe$nobs[["nobs.full"]]
  se_s <- sampling_se_clustered(g, id)
  se_c_fe <- sqrt(a_fe$vcov[1, 1]); se_c_bc <- sqrt(a_bc$vcov[1, 1])
  se_fe <- sqrt(se_c_fe^2 + se_s^2); se_bc <- sqrt(se_c_bc^2 + se_s^2)
  list(ape_fe = unname(a_fe$delta), ape_bc = unname(a_bc$delta),
       ape_fe_movers = unname(a_fe$delta) / share_kept,
       share_kept = share_kept, g = g, se_cond_fe = se_c_fe, se_cond_bc = se_c_bc, se_sample = se_s,
       se_fe = se_fe, se_bc = se_bc, ci_fe = ci95(unname(a_fe$delta), se_fe), ci_bc = ci95(unname(a_bc$delta), se_bc))
}

## ---- binary MLE helpers -------------------------------------------------------------------------
#' Probit or logit MLE by glm.fit on a design matrix (no formula), with the model-based variance
#' (inverse Hessian; information-matrix equality) and the HC0 sandwich (Wooldridge 2010 eq. 15.18 and
#' the robust form of section 13.5). Returns coef, vcov_model, vcov_sand, loglik.
binary_mle <- function(y, X, link = c("probit", "logit")) {
  link <- match.arg(link)
  X <- as.matrix(X)
  f <- glm.fit(X, y, family = binomial(link), control = glm.control(epsilon = 1e-12, maxit = 100))
  if (!f$converged) stop("binary_mle: glm.fit did not converge")
  b <- f$coefficients
  L <- link_fns(link)
  u <- drop(X %*% b); p <- pmin(pmax(L$F(u), 1e-15), 1 - 1e-15); d <- L$Fp(u)
  w <- d^2 / (p * (1 - p))                                     # expected-information weights
  A <- crossprod(X, X * w)
  S <- X * ((y - p) * d / (p * (1 - p)))                         # scores s_i(beta)
  Ainv <- solve(A)
  list(coef = b, vcov_model = Ainv, vcov_sand = Ainv %*% crossprod(S) %*% Ainv, scores = S,
       loglik = sum(y * log(p) + (1 - y) * log(1 - p)))
}

## ---- Rivers-Vuong control function (Wooldridge 2010 section 15.7.2, eqs. 15.43 to 15.47) ----------
#' One pass of the two steps on the rows `idx`. Step 1: OLS of D_i on (X_i, Z_i) gives pi-hat and the
#' first-stage residual v_i-hat. Step 2: probit of Y_i on (X_i, D_i, v_i-hat) gives the scaled
#' coefficients (beta_e, beta_De, rho_e) = (beta, beta_D, rho) / sigma_e. The APE of D_i comes from the
#' average structural function ASF(x, d) = mean_h Phi(x'beta_e + beta_De d + rho_e v_h-hat) (eq. 15.47):
#'   APE = mean_i mean_h beta_De phi(X_i'beta_e + beta_De D_i + rho_e v_h-hat)   (a double average).
rv_fit <- function(y, X_exog, x_endog, Z, idx = seq_along(y), chunk = 400L) {
  Xe <- X_exog[idx, , drop = FALSE]; d <- x_endog[idx]; yy <- y[idx]
  W <- cbind(Xe, Z[idx, , drop = FALSE])
  pi_hat <- qr.coef(qr(W), d)
  vhat <- d - drop(W %*% pi_hat)
  R <- cbind(Xe, d, vhat)
  f <- glm.fit(R, yy, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
  b <- f$coefficients
  k <- ncol(Xe)
  a <- drop(Xe %*% b[1:k]) + b[k + 1] * d
  cc <- b[k + 2] * vhat
  s <- 0
  for (ch in split(seq_along(a), ceiling(seq_along(a) / chunk))) s <- s + sum(dnorm(outer(a[ch], cc, "+")))
  ape <- unname(b[k + 1]) * s / (length(a) * length(cc))
  list(coef = unname(b), pi = unname(pi_hat), vhat = vhat, ape = ape, glm = f, R = R)
}

#' The chapter card's rivers_vuong(y, X_exog, x_endog, Z, B = 499, seed). X_exog: exogenous regressors
#' X_i including the constant; x_endog: D_i; Z: excluded instruments Z_i (matrix or vector).
#' Both steps are bootstrapped together (pairs bootstrap of whole rows, B draws, indices drawn after
#' set.seed(seed)). Returns
#'   beta     scaled step-2 coefficients (X_i block, beta_D, rho) / sigma_e, named
#'   ape      APE of D_i via the ASF;  ape_se  bootstrap sd;  ape_ci  ape -/+ 1.96 ape_se;
#'   ape_ci_pct   2.5% and 97.5% bootstrap percentiles;  beta_se  bootstrap sd of beta;
#'   pi       first-stage coefficients (pi_1, pi_2);  rho_t  usual probit t statistic of v-hat
#'            (a valid test of exogeneity H0: rho = 0, eq. 15.44);  boot  B x (k + 3) draws (beta, ape).
rivers_vuong <- function(y, X_exog, x_endog, Z, B = 499, seed, cores = 1L) {
  X_exog <- as.matrix(X_exog); Z <- as.matrix(Z)
  n <- length(y)
  full <- rv_fit(y, X_exog, x_endog, Z)
  xn <- colnames(X_exog)
  if (is.null(xn)) xn <- rep("", ncol(X_exog))
  xn[is.na(xn) | !nzchar(xn)] <- paste0("x", which(is.na(xn) | !nzchar(xn)))
  nm <- c(xn, "D", "vhat")
  set.seed(seed)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), n, B)
  one <- function(b) { f <- rv_fit(y, X_exog, x_endog, Z, idx[, b]); c(f$coef, f$ape) }
  draws <- if (cores > 1L) parallel::mclapply(seq_len(B), one, mc.cores = cores) else lapply(seq_len(B), one)
  boot <- do.call(rbind, draws)
  colnames(boot) <- c(nm, "ape")
  ape_se <- sd(boot[, "ape"])
  # usual probit standard error of the v-hat coefficient (valid under H0: rho = 0)
  bm <- binary_mle(y, full$R, "probit")
  rho_t <- full$coef[length(full$coef)] / sqrt(bm$vcov_model[ncol(full$R), ncol(full$R)])
  list(beta = setNames(full$coef, nm), ape = full$ape, ape_se = ape_se,
       ape_ci = ci95(full$ape, ape_se), ape_ci_pct = unname(quantile(boot[, "ape"], c(0.025, 0.975))),
       beta_se = setNames(apply(boot[, nm, drop = FALSE], 2, sd), nm), pi = full$pi, rho_t = unname(rho_t),
       B = B, boot = boot)
}

#' IV probit by maximum likelihood (Wooldridge 2010 eqs. 15.48 to 15.50): the joint density of
#' (Y_i, D_i) given (X_i, Z_i) with (eps_i, v_i) bivariate normal, Var(eps_i) = 1, Var(v_i) = sigma_v^2,
#' Corr = r. Parameters: (beta, beta_D, pi, log sigma_v, atanh r). Returns coef (structural scale),
#' vcov (inverse numerical Hessian) and the APE of D_i, mean_i beta_D phi(X_i'beta + beta_D D_i).
ivprobit_mle <- function(y, X_exog, x_endog, Z, start = NULL) {
  X_exog <- as.matrix(X_exog); Z <- as.matrix(Z)
  W <- cbind(X_exog, Z); k <- ncol(X_exog); q <- ncol(W)
  d <- x_endog
  nll <- function(th) {
    b <- th[1:k]; bD <- th[k + 1]; pi <- th[k + 1 + 1:q]
    sv <- exp(th[k + q + 2]); r <- tanh(th[k + q + 3])
    v <- d - drop(W %*% pi)
    m <- (drop(X_exog %*% b) + bD * d + r * v / sv) / sqrt(1 - r^2)
    -sum(dnorm(v, 0, sv, log = TRUE) + ifelse(y == 1, pnorm(m, log.p = TRUE), pnorm(-m, log.p = TRUE)))
  }
  if (is.null(start)) {
    rv <- rv_fit(y, X_exog, x_endog, Z)
    sv0 <- sd(rv$vhat); re <- rv$coef[k + 2] * sv0; se <- 1 / sqrt(1 + re^2)
    start <- c(rv$coef[1:(k + 1)] * se, rv$pi, log(sv0), atanh(max(min(re * se, 0.95), -0.95)))
  }
  o <- optim(start, nll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-14))
  o <- optim(o$par, nll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-15))
  H <- numDeriv::hessian(nll, o$par)
  V <- solve(H)
  th <- o$par
  ape_fn <- function(t) mean(t[k + 1] * dnorm(drop(X_exog %*% t[1:k]) + t[k + 1] * d))
  G <- numDeriv::grad(ape_fn, th)
  list(coef = th, vcov = V, convergence = o$convergence, loglik = -o$value,
       beta = th[1:(k + 1)], se_beta = sqrt(diag(V))[1:(k + 1)], sigma_v = exp(th[k + q + 2]), r = tanh(th[k + q + 3]),
       ape = ape_fn(th), ape_se_delta = sqrt(drop(t(G) %*% V %*% G)))
}

## ---- panel binary models (Wooldridge 2010 section 15.8) ------------------------------------------
#' The chapter card's cre_probit(y, X, id): pooled correlated random effects probit (Mundlak device,
#' Chamberlain 1980; Wooldridge 2010 eqs. 15.72 to 15.76): pooled probit of Y_it on (1, X_it, Xbar_i),
#' which estimates the scaled coefficients (xi_0, beta, xi) / sqrt(1 + sigma_a^2). The APE of the j-th
#' time-varying regressor averages over all (i, t): mean_it beta_aj phi(xi_a0 + X_it'beta_a + Xbar_i'xi_a).
#' Variance: cluster-robust by unit (the partial-likelihood sandwich of section 15.8.1); the APE standard
#' error adds the clustered sampling variance of the per-observation effects.
#' X: n x k matrix (or vector) of time-varying regressors without a constant. Returns beta (all pooled
#' coefficients, named), ape (one per column of X), ape_se, ape_ci (list), vcov.
cre_probit <- function(y, X, id) {
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  Xbar <- apply(X, 2, function(x) ave(x, id))
  colnames(Xbar) <- paste0(colnames(X), "_bar")
  R <- cbind("(Intercept)" = 1, X, Xbar)
  m <- binary_mle(y, R, "probit")
  b <- setNames(m$coef, colnames(R))
  # cluster-robust variance: A^-1 (sum over units of summed scores) A^-1
  Sc <- rowsum(m$scores, id)
  V <- m$vcov_model %*% crossprod(Sc) %*% m$vcov_model
  apes <- lapply(seq_len(ncol(X)), function(j) {
    jj <- 1 + j
    avg_effect_inference(function(t) dnorm(drop(R %*% t)) * t[jj], m$coef, V, cluster = id)
  })
  names(apes) <- colnames(X)
  list(beta = b, ape = vapply(apes, `[[`, numeric(1), "est"), ape_se = vapply(apes, `[[`, numeric(1), "se"),
       ape_ci = lapply(apes, `[[`, "ci"), vcov = V, se = sqrt(diag(V)))
}

#' Fixed effects logit by conditional MLE (Chamberlain 1980; Wooldridge 2010 eqs. 15.81 to 15.83, where the
#' general-T conditional likelihood is credited to Chamberlain 1984): for unit i with n_i = sum_t Y_it,
#'   l_i(beta) = sum_t Y_it X_it'beta - log sum_{d in S_i} exp(sum_t d_t X_it'beta),
#' S_i the binary sequences of length T with sum n_i (alpha_i drops out). Units with n_i = 0 or T carry no
#' information. Newton iterations; returns coef, vcov_model (inverse Hessian), vcov_robust (sum of unit
#' score outer products in the sandwich: hand-written because survival::clogit refuses robust variances
#' under its exact method), n_units_used.
fe_logit_cmle <- function(y, X, id, tol = 1e-12, maxit = 100) {
  X <- as.matrix(X); k <- ncol(X)
  units <- split(seq_along(y), id)
  prep <- list()
  for (u in units) {
    s <- sum(y[u]); T <- length(u)
    if (s == 0 || s == T) next
    cmb <- utils::combn(T, s)
    Dm <- matrix(0, ncol(cmb), T)
    Dm[cbind(rep(seq_len(ncol(cmb)), each = s), as.vector(cmb))] <- 1
    Sx <- Dm %*% X[u, , drop = FALSE]                 # sum_t d_t X_it for every d in S_i
    prep[[length(prep) + 1L]] <- list(obs = colSums(X[u, , drop = FALSE] * y[u]), Sx = Sx)
  }
  b <- rep(0, k)
  parts <- function(b) {
    g <- rep(0, k); H <- matrix(0, k, k); S <- matrix(0, length(prep), k); ll <- 0
    for (i in seq_along(prep)) {
      e <- drop(prep[[i]]$Sx %*% b); m <- max(e); w <- exp(e - m); w <- w / sum(w)
      mu <- colSums(prep[[i]]$Sx * w)
      si <- prep[[i]]$obs - mu
      C <- crossprod(prep[[i]]$Sx * sqrt(w)) - tcrossprod(mu)
      g <- g + si; H <- H - C; S[i, ] <- si
      ll <- ll + sum(prep[[i]]$obs * b) - (m + log(sum(exp(e - m))))
    }
    list(g = g, H = H, S = S, ll = ll)
  }
  for (it in seq_len(maxit)) {
    p <- parts(b)
    step <- solve(p$H, p$g)
    b <- b - step
    if (max(abs(step)) < tol) break
  }
  p <- parts(b)
  Hinv <- solve(-p$H)
  list(coef = setNames(b, colnames(X)), vcov_model = Hinv, vcov_robust = Hinv %*% crossprod(p$S) %*% Hinv,
       n_units_used = length(prep), iter = it, loglik = p$ll)
}

## ---- multinomial and nested logit (Wooldridge 2010 sections 16.2.1 and 16.2.2) --------------------
#' MNL probabilities P_m(x) for categories 0..M (eq. 16.1); B is M x k (row m = beta_(m)'), W is n x k.
mnl_probs <- function(B, W) {
  E <- cbind(1, exp(W %*% t(B)))
  E / rowSums(E)
}

#' Nested logit probabilities with the base category alone and categories 1..M in one nest with
#' dissimilarity lambda (Wooldridge 2010 eqs. 16.10 to 16.12; utilities V_0 = 0, V_m = W beta_(m)):
#'   P(nest) = exp(lambda IV) / (1 + exp(lambda IV)),  IV = log sum_m exp(V_m / lambda),
#'   P(m) = P(nest) exp(V_m / lambda - IV).  lambda = 1 gives mnl_probs (IIA holds).
nl_probs <- function(B, W, lambda) {
  V <- W %*% t(B)
  mx <- apply(V / lambda, 1, max)
  IV <- mx + log(rowSums(exp(V / lambda - mx)))
  Pn <- plogis(lambda * IV)
  cbind(1 - Pn, Pn * exp(V / lambda - IV))
}

#' One categorical draw per row of the probability matrix P (columns = categories 1..ncol).
draw_choice <- function(P) {
  u <- runif(nrow(P))
  1L + rowSums(u > t(apply(P, 1, cumsum))[, -ncol(P), drop = FALSE])
}

#' APE of regressor j on every category probability under MNL: mean_i P_m(X_i)(beta_(m)j - sum_h P_h beta_(h)j),
#' with beta_(0) = 0 (the derivative of eq. 16.1). Returns a vector over categories 0..M.
mnl_ape <- function(B, W, j) {
  P <- mnl_probs(B, W)
  bj <- c(0, B[, j])
  colMeans(P * (matrix(bj, nrow(P), length(bj), byrow = TRUE) - drop(P %*% bj)))
}

#' Hausman and McFadden (1984) test of IIA with nnet::multinom (tight tolerance; the default
#' reltol = 1e-8 moves the statistic). Full model on all categories; restricted model on the rows that
#' did not choose `drop`. H = (b_r - b_f)'(V_r - V_f)^{-1}(b_r - b_f), df = number of restricted
#' coefficients; V_r - V_f need not be positive semidefinite in a finite sample, so H can be negative
#' (its p-value is then 1). formula: e.g. choice ~ x + s. When `drop` is the base category (the first
#' factor level), the factor is re-levelled to the first remaining category before both fits (the MNL
#' likelihood does not depend on the base): with drop = "Free" the test compares the Enterprise-versus-Pro
#' coefficients beta_Ent - beta_Pro of the full model with the binary logit among the paying firms.
#' The arithmetic is that of mlogit::hmftest (tests/test_fns.R checks it). V_r - V_f is the difference of two
#' nearly equal matrices (dropping one category costs little information about the remaining coefficients), so it
#' is near-singular and often not positive semidefinite in finite samples: that, not the chi-square
#' approximation, is why the test over-rejects under IIA (hm_test_stacked is the contrast that avoids it).
#' Returns stat, df, p, base, min_eigen and max_eigen (smallest and largest eigenvalue of V_r - V_f), full,
#' restricted, and what hm_test_stacked needs to refit nothing: formula, data (re-levelled), resp, drop, nm
#' (names of the compared coefficients, "category:term"), dif (b_r - b_f).
hm_test <- function(formula, data, drop) {
  resp <- all.vars(formula)[1]
  if (identical(drop, levels(data[[resp]])[1])) {
    data[[resp]] <- relevel(data[[resp]], ref = setdiff(levels(data[[resp]]), drop)[1])
  }
  mf <- nnet::multinom(formula, data = data, trace = FALSE, reltol = 1e-12, maxit = 1000, Hess = TRUE)
  ds <- droplevels(data[data[[resp]] != drop, , drop = FALSE])
  mr <- nnet::multinom(formula, data = ds, trace = FALSE, reltol = 1e-12, maxit = 1000, Hess = TRUE)
  keep <- setdiff(levels(data[[resp]])[-1], drop)
  Bf <- coef(mf); Br <- coef(mr)
  if (is.null(dim(Br))) Br <- matrix(Br, nrow = 1, dimnames = list(levels(ds[[resp]])[2], names(Br)))
  nm <- as.vector(t(outer(keep, colnames(Bf), paste, sep = ":")))
  br <- as.vector(t(Br[keep, , drop = FALSE])); bf <- as.vector(t(Bf[keep, , drop = FALSE]))
  Vr <- vcov(mr); Vf <- vcov(mf)[nm, nm]
  dif <- br - bf
  H <- drop(t(dif) %*% solve(Vr - Vf) %*% dif)
  df <- length(dif)
  ev <- eigen(Vr - Vf, symmetric = TRUE, only.values = TRUE)$values
  list(stat = H, df = df, p = pchisq(H, df, lower.tail = FALSE), base = levels(data[[resp]])[1],
       full = mf, restricted = mr, min_eigen = min(ev), max_eigen = max(ev),
       formula = formula, data = data, resp = resp, drop = drop, nm = nm, dif = dif)
}

#' Stacked-score (suest-type) version of the Hausman and McFadden (1984) test, a contrast to hm_test. The full
#' and the restricted multinomial logit are estimated on the same observations (the restricted one with zero
#' score for the rows that chose `drop`); the variance of b_r - b_f is that of the difference of the two
#' influence functions, IF_i = V s_i (V the model-based variance of each fit, s_i its score), summed over
#' i: V_d = sum_i (IF_r,i - IF_f,i)(IF_r,i - IF_f,i)'. V_d is positive semidefinite by construction, so the
#' statistic H = (b_r - b_f)' V_d^{-1} (b_r - b_f) is never negative (hm_test's V_r - V_f is not), and it
#' needs no information-matrix equality for the restricted fit. h: the result of hm_test (it carries the two
#' fits, the re-levelled data and the compared coefficients). Scores of the multinomial logit: s_i,(m) =
#' W_i (1{Y_i = m} - P_m(W_i)) for each non-base category m (the same as sandwich::estfun on a multinom fit;
#' tests/test_fns.R checks it). Returns stat, df, p, min_eigen (of V_d), dif, V_diff, scores (full and restricted).
hm_test_stacked <- function(h) {
  d <- h$data; resp <- h$resp
  W <- model.matrix(delete.response(terms(h$formula)), d)
  scores <- function(fit, rows) {
    ds <- d[rows, , drop = FALSE]; Wd <- W[rows, , drop = FALSE]; yd <- droplevels(ds[[resp]])
    B <- coef(fit)
    if (is.null(dim(B))) B <- matrix(B, nrow = 1, dimnames = list(levels(yd)[2], names(B)))
    P <- mnl_probs(B, Wd)
    S <- do.call(cbind, lapply(seq_len(nrow(B)), function(m) {
      s <- Wd * (as.numeric(yd == rownames(B)[m]) - P[, m + 1L])
      colnames(s) <- paste(rownames(B)[m], colnames(Wd), sep = ":")
      s
    }))
    out <- matrix(0, nrow(d), ncol(S), dimnames = list(NULL, colnames(S)))
    out[rows, ] <- S
    out
  }
  Sf <- scores(h$full, seq_len(nrow(d))); Sr <- scores(h$restricted, which(d[[resp]] != h$drop))
  Vf <- vcov(h$full); Vr <- vcov(h$restricted)
  if (!all(grepl(":", rownames(Vr), fixed = TRUE))) dimnames(Vr) <- list(colnames(Sr), colnames(Sr))   # one restricted category
  stopifnot(identical(rownames(Vf), colnames(Sf)), identical(rownames(Vr), colnames(Sr)), setequal(colnames(Sr), h$nm))
  IFf <- (Sf %*% Vf)[, h$nm, drop = FALSE]; IFr <- (Sr %*% Vr)[, h$nm, drop = FALSE]
  Vd <- crossprod(IFr - IFf)
  H <- drop(t(h$dif) %*% solve(Vd, h$dif))
  list(stat = H, df = length(h$dif), p = pchisq(H, length(h$dif), lower.tail = FALSE),
       min_eigen = min(eigen(Vd, symmetric = TRUE, only.values = TRUE)$values),
       dif = h$dif, V_diff = Vd, scores = list(full = Sf, restricted = Sr))
}

## ---- ordered probit (Wooldridge 2010 section 16.3, eqs. 16.21 to 16.25) ----------------------------
#' Category probabilities of an ordered probit with cut points mu_1 < ... < mu_M and latent standard
#' deviation sigma_i (heteroskedastic: ordinal::clm's scale part): P(Y = m | x) =
#' Phi((mu_m - x'beta)/sigma) - Phi((mu_{m-1} - x'beta)/sigma). index: x_i'beta; returns n x (M + 1).
oprobit_probs <- function(index, cuts, sigma = 1) {
  cdf <- sapply(cuts, function(m) pnorm((m - index) / sigma))
  cdf <- cbind(0, matrix(cdf, nrow = length(index)), 1)
  cdf[, -1, drop = FALSE] - cdf[, -ncol(cdf), drop = FALSE]
}

#' Per-observation effect of regressor x_j on P(Y = top category): beta_j / sigma_i phi((mu_M - x_i'beta)/sigma_i).
oprobit_top_i <- function(beta_j, index, mu_top, sigma = 1) beta_j / sigma * dnorm((mu_top - index) / sigma)

## ---- data generators of the new Helix blocks B2 to B7 (01_dgp.R calls them after set.seed) --------
## p is the block's parameter list (truth12.rds$params); every generator draws in a fixed order, so that
## set.seed(seed) followed by the generator reproduces the saved data exactly (tests/test_dgp.R).

#' B2: two features. x1 = API engagement (z-score), x2 = Helix Deploy pipeline usage (z-score),
#' Corr(x1, x2) = rho_x; probit P(Y = 1 | x) = Phi(b0 + b1 x1 + b2 x2 + b12 x1 x2) (beta_jj' = b12).
gen_B2 <- function(p, n = p$n) {
  x1 <- rnorm(n)
  x2 <- p$rho_x * x1 + sqrt(1 - p$rho_x^2) * rnorm(n)
  b <- p$beta
  p_true <- pnorm(b[1] + b[2] * x1 + b[3] * x2 + b[4] * x1 * x2)
  y <- rbinom(n, 1, p_true)
  data.frame(x1 = x1, x2 = x2, y = y, p_true = p_true)
}

#' B3: endogenous trial engagement D_i, randomised onboarding reminder Z_i. X_i = log(cloud spend /
#' median), spend in thousands of USD; first stage D = pi0 + pix X + piz Z + v, v ~ N(0, sigma_v^2);
#' structural Y = 1{b0 + bx X + beta_D D + eps > 0}, (eps, v) bivariate normal, Var(eps) = 1, Corr = r.
gen_B3 <- function(p, n = p$n) {
  X <- p$sd_logspend * rnorm(n)
  spend_kusd <- p$spend_median_kusd * exp(X)
  Z <- rbinom(n, 1, p$p_reminder)
  v <- p$sigma_v * rnorm(n)
  e <- rnorm(n)
  D <- p$pi[1] + p$pi[2] * X + p$pi[3] * Z + v
  eps <- p$r * v / p$sigma_v + sqrt(1 - p$r^2) * e
  Y <- as.integer(p$beta[1] + p$beta[2] * X + p$beta_D * D + eps > 0)
  data.frame(spend_kusd = spend_kusd, X = X, Z = Z, D = D, Y = Y, v_true = v, eps_true = eps)
}

#' B4: omitted heterogeneity independent of X. Y = 1{b0 + b1 X + a + eps > 0}, a ~ N(0, sigma_a^2)
#' (unobserved budget flexibility), eps ~ N(0, 1), X = engagement z-score, all independent.
gen_B4 <- function(p, n = p$n) {
  X <- rnorm(n)
  a <- p$sigma_a * rnorm(n)
  eps <- rnorm(n)
  Y <- as.integer(p$beta[1] + p$beta[2] * X + a + eps > 0)
  data.frame(X = X, Y = Y, a_true = a)
}

#' B5: balanced monthly panel, N firms, T months. X_it = c_i + w_it (monthly engagement), c_i ~ N(0, sigma_c^2);
#' alpha_i = xi0 + xi Xbar_i + a_i, a_i ~ N(0, sigma_a^2) (Mundlak form, correlated with X);
#' P(Y_it = 1 | X_i, alpha_i) = Lambda(alpha_i + beta X_it) (FE logit is correctly specified).
gen_B5 <- function(p, N = p$N, T = p$T) {
  c_i <- p$sigma_c * rnorm(N)
  X <- rep(c_i, each = T) + rnorm(N * T)
  id <- rep(seq_len(N), each = T)
  Xbar <- ave(X, id)
  a <- p$sigma_a * rnorm(N)
  alpha <- p$xi0 + p$xi * Xbar + rep(a, each = T)
  Y <- rbinom(N * T, 1, plogis(alpha + p$beta * X))
  data.frame(id = id, t = rep(seq_len(T), times = N), X = X, Y = Y, alpha_true = alpha)
}

#' B6: 1 to 5 satisfaction rating. Y* = b1 x1 + b2 bill + delta G + sigma_G eps, eps ~ N(0, 1);
#' x1 = engagement z-score; bill = monthly Helix bill in thousands of USD (lognormal);
#' G = 1 for managed (account-manager) firms, latent sd sigma_G = sigma[G + 1]; rating m when
#' mu_{m-1} < Y* <= mu_m (mu_0 = -Inf, mu_5 = Inf).
gen_B6 <- function(p, n = p$n) {
  x1 <- rnorm(n)
  bill <- exp(log(p$bill_median_kusd) + p$sd_logbill * rnorm(n))
  G <- rbinom(n, 1, p$p_managed)
  ystar <- p$beta[1] * x1 + p$beta[2] * bill + p$delta_G * G + p$sigma[G + 1] * rnorm(n)
  rating <- findInterval(ystar, p$cuts, left.open = TRUE) + 1L
  data.frame(x1 = x1, bill_kusd = bill, G = G, rating = factor(rating, levels = 1:5, ordered = TRUE),
             ystar_true = ystar)
}

#' B7: plan choice Free (base) / Pro / Enterprise. x = engagement z-score, s = standardised log seats;
#' V_Pro = W beta_Pro, V_Ent = W beta_Ent, W = (1, x, s). lambda = 1: MNL (IIA holds); lambda < 1: nested
#' logit with the paid nest {Pro, Enterprise} (IIA fails).
gen_B7 <- function(p, n = p$n, lambda = 1) {
  x <- rnorm(n); s <- rnorm(n)
  W <- cbind(1, x, s)
  P <- if (lambda == 1) mnl_probs(p$coef, W) else nl_probs(p$coef, W, lambda)
  ch <- draw_choice(P)
  data.frame(x = x, s = s, choice = factor(p$levels[ch], levels = p$levels),
             p_free_true = P[, 1], p_pro_true = P[, 2], p_ent_true = P[, 3])
}

#' Gauss-Hermite nodes and weights for a standard normal variable (Golub-Welsch: the eigenvalues and the squared
#' first eigenvector components of the Jacobi matrix of the probabilists' Hermite polynomials); the rule with n
#' nodes integrates polynomials of degree up to 2n - 1 exactly. Returns z (nodes) and w (weights summing to one).
gauss_hermite_normal <- function(n = 120L) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}

#' Exact population truths of block B5 (gen_B5) by two-dimensional Gauss-Hermite quadrature (controller ruling 44).
#' In gen_B5, X_it = c_i + w_it and alpha_i = xi0 + xi Xbar_i + a_i, with Xbar_i the mean of the Tm months' X: the pair
#' (X_it, alpha_i) is bivariate normal with the same law in every month t:
#'   Var X = sigma_c^2 + 1,  Cov(X, alpha) = xi (sigma_c^2 + 1 / Tm),  Var alpha = xi^2 (sigma_c^2 + 1 / Tm) + sigma_a^2,
#' mean (0, xi0), and Y | (X, alpha) ~ Bernoulli(Lambda(alpha + beta X)). So the months-1-to-4 sub-panel of a T = 8
#' panel has exactly the T = 8 values (Tm = 8), and a 4-month DGP is Tm = 4. p: list(sigma_c, sigma_a, xi0, xi,
#' beta); n: nodes per dimension (stable at 60 and 120). Returns the true APE E[beta Lambda'(alpha + beta X)], the rate
#' P(Y = 1), the pooled logit probability limit (pooled_int, pooled_plim: the maximiser of the expected logit
#' log likelihood, Newton iterations on the population score) and the APE of the pooled logit at its limit.
b5_quadrature <- function(p, Tm, n = 120L) {
  vX <- p$sigma_c^2 + 1; cXa <- p$xi * (p$sigma_c^2 + 1 / Tm); va <- p$xi^2 * (p$sigma_c^2 + 1 / Tm) + p$sigma_a^2
  g <- gauss_hermite_normal(n); G <- expand.grid(i = seq_len(n), j = seq_len(n))
  z1 <- g$z[G$i]; z2 <- g$z[G$j]; w <- g$w[G$i] * g$w[G$j]
  X <- sqrt(vX) * z1; al <- p$xi0 + cXa / sqrt(vX) * z1 + sqrt(va - cXa^2 / vX) * z2
  pr <- plogis(al + p$beta * X)
  b <- c(0, 1)
  for (it in 1:100) {
    q <- plogis(b[1] + b[2] * X); sc <- c(sum(w * (pr - q)), sum(w * (pr - q) * X))
    H <- matrix(c(sum(w * q * (1 - q)), sum(w * q * (1 - q) * X), sum(w * q * (1 - q) * X), sum(w * q * (1 - q) * X^2)), 2)
    b <- b + solve(H, sc)
    if (max(abs(sc)) < 1e-15) break
  }
  list(ape = sum(w * p$beta * dlogis(al + p$beta * X)), rate = sum(w * pr), pooled_int = b[1], pooled_plim = b[2],
       pooled_ape_plim = sum(w * dlogis(b[1] + b[2] * X)) * b[2])
}

#' Quadrature nodes for a standard normal variable: an equally spaced grid on [-lim, lim] with weights
#' proportional to phi (normalised to sum to one). With 401 nodes on [-8, 8] the error for smooth
#' integrands is far below 1e-10 (the trapezoid rule is spectrally accurate here).
normal_grid <- function(n = 401L, lim = 8) {
  z <- seq(-lim, lim, length.out = n)
  w <- dnorm(z)
  list(z = z, w = w / sum(w))
}
