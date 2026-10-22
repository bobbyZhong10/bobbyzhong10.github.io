# ------------------------------------------------------------------
# dyn_fns.R -- estimators, bias formulas and data generators for chapter 11 (dynamic panels).
#
# Notation (notation_v3.md 1.4, 1.8 and 1.9): unit i = 1..N (a Nimbus merchant), period t (a month).
# Model
#     Y_it = rho Y_{i,t-1} + X_it' beta + alpha_i + eps_it,
# rho the AR coefficient (bare rho is only the AR coefficient in chapter 11), alpha_i the merchant
# effect (Wooldridge writes c_i), eps_it the idiosyncratic error, e_it an innovation, s a lag order.
# Arguments, everywhere: y is the stacked Y_it; id is the unit index i; t is the period index (integer
# months, consecutive within a unit); X (optional) is the n x k matrix of X_it.
#
# Two period counts, never mixed up:
#   T_obs  months observed per merchant (Y_i1, ..., Y_iT_obs); the legacy code's T (gen_dyn(T = 8))
#          and the instrument count (T_obs - 1)(T_obs - 2)/2 of difference GMM use T_obs;
#   T      the number of REGRESSION periods, T = T_obs - 1 (the first month only supplies a lag).
#          nickell_bias, nickell_first_factor and hk_correct take this T, as in Nickell (1981) and
#          Hahn and Kuersteiner (2001 working paper, model (5): t = 1, ..., T with Y_i0 known).
#
# Card functions:
#   anderson_hsiao(y, id, t, inst)  -> list(est, se, ci, n, inst): rho_hat_AH, IV on the first-differenced
#                                      equation with Y_{i,t-2} ("level") or Delta Y_{i,t-2} ("diff")
#   spj_fe(y, id, t)                -> numeric: rho_hat_SPJ = 2 rho_hat_FE - bar_rho_{1/2}, the halves weighted
#                                      by |S|/T (Dhaene and Jochmans 2015, eqs. 2.2 and 2.4; equal weights for even T)
#   hk_correct(rho_fe, T)           -> numeric: rho_hat_HK = ((T + 1)/T) rho_hat_FE + 1/T
#   n_instruments(fit)              -> integer: r, the number of instrument columns of a plm::pgmm fit
# Helpers: lag_within, fe_ar1, pooled_ar1, ols_diff_ar1, spj_parts, nickell_bias, nickell_first_factor,
# gmm_summary, ci95; data generators gen_D2, gen_D2_ma, gen_D3, gen_D4 (blocks D1, D1p and the Nickell Monte Carlo
# D1T use the inherited legacy generator gen_dyn in 01_dgp.R; block D5 reuses gen_dyn).
# Cluster-robust variances are by unit without a small-sample factor (CR0 = plm's
# vcovHC(method = "arellano", type = "HC0")).
# ------------------------------------------------------------------

## ---- internal helpers ---------------------------------------------------------------------
.unit_index <- function(id) match(id, unique(id))                  # 1..G in order of appearance
.demean <- function(M, g) {                                         # subtract unit means (rows of M)
  M <- as.matrix(M)
  M - (rowsum(M, g) / as.vector(rowsum(rep(1, nrow(M)), g)))[g, , drop = FALSE]
}
# CR0 cluster-robust variance: bread * sum_g (Z_g' res_g)(Z_g' res_g)' * bread, res the residual vector
.cl_vcov <- function(Z, res, g, bread = solve(crossprod(Z))) {
  S <- rowsum(Z * as.vector(res), g)
  V <- bread %*% crossprod(S) %*% bread
  dimnames(V) <- list(colnames(Z), colnames(Z))
  V
}
# numeric key of (i, t) with room for lags and leads of up to `pad` months (no aliasing across units)
.panel_key <- function(id, t, pad = 0) {
  if (any(t != round(t))) stop("t must be integer-valued (months)")
  span <- max(t) - min(t) + pad + 1
  match(id, unique(id)) * span + (t - min(t))
}
.check_panel <- function(y, id, t) {
  if (length(y) != length(id) || length(y) != length(t)) stop("y, id and t must have the same length")
  if (anyNA(id) || anyNA(t)) stop("id and t must not contain NA")
  if (anyDuplicated(.panel_key(id, t))) stop("(id, t) pairs must be unique")
  invisible(TRUE)
}

#' 95% interval est -/+ qnorm(0.975) se.
ci95 <- function(est, se) c(est - qnorm(0.975) * se, est + qnorm(0.975) * se)

#' v at (i, t - k): the k-th lag of v within unit i (k < 0: a lead), NA when month t - k of unit i is not
#' in the data. v, id, t are stacked over (i, t) in any order; t is integer-valued.
lag_within <- function(v, id, t, k = 1L) {
  key <- .panel_key(id, t, pad = abs(k))
  v[match(key - k, key)]
}

## ---- FE, pooled OLS and OLS in differences of the AR(1) panel model -----------------------
#' Within (FE) estimator of Y_it = rho Y_{i,t-1} [+ X_it' beta] + alpha_i + eps_it, i.e. rho_hat_FE
#' (and beta_hat_FE): pooled OLS of the demeaned Y_it on the demeaned Y_{i,t-1} [and X_it], unit means
#' over the regression periods (the months with Y_{i,t-1} observed), as plm(model = "within").
#' Returns est (named: "rho", then the columns of X), se and V (CR0 by unit), n (regression rows),
#' G (units), T (regression periods per unit, when balanced).
#' Arguments: y is $Y_{it}$, id is $i$, t is $t$, X is $X_{it}$; est is ($\hat\rho_{FE}$, $\hat\beta_{FE}$).
fe_ar1 <- function(y, id, t, X = NULL) {
  .check_panel(y, id, t)
  yl <- lag_within(y, id, t, 1L)
  R <- cbind(rho = yl)
  if (!is.null(X)) {
    X <- as.matrix(X); if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
    R <- cbind(R, X)
  }
  ok <- !is.na(y) & stats::complete.cases(R)
  g <- .unit_index(id[ok])
  Yd <- .demean(y[ok], g); Rd <- .demean(R[ok, , drop = FALSE], g); colnames(Rd) <- colnames(R)
  A <- solve(crossprod(Rd)); b <- drop(A %*% crossprod(Rd, Yd)); names(b) <- colnames(R)
  e <- drop(Yd - Rd %*% b)
  V <- .cl_vcov(Rd, e, g, A)
  Tg <- tabulate(g)
  list(est = b, se = sqrt(diag(V)), V = V, n = length(e), G = max(g),
       T = if (length(unique(Tg)) == 1L) Tg[1] else NA_integer_)
}

#' Pooled OLS of Y_it on (1, Y_{i,t-1}) [plus month dummies when time_fe = TRUE]: the level equation of
#' the system without instruments. Returns est (the coefficient on Y_{i,t-1}), se (CR0 by unit), n.
pooled_ar1 <- function(y, id, t, time_fe = FALSE) {
  .check_panel(y, id, t)
  yl <- lag_within(y, id, t, 1L); ok <- !is.na(y) & !is.na(yl)
  Z <- if (time_fe) stats::model.matrix(~ yl + factor(t), data.frame(yl = yl[ok], t = t[ok])) else cbind(1, yl = yl[ok])
  colnames(Z)[2] <- "yl"
  A <- solve(crossprod(Z)); b <- drop(A %*% crossprod(Z, y[ok])); e <- drop(y[ok] - Z %*% b)
  V <- .cl_vcov(Z, e, .unit_index(id[ok]), A)
  list(est = unname(b[2]), se = sqrt(V[2, 2]), n = length(e))
}

#' OLS of Delta Y_it on Delta Y_{i,t-1} (no intercept): the first-differenced equation without
#' instruments, kept as a reference value (block D2 stores it). It is NOT the estimator that
#' proliferation pulls GMM towards (controller Ruling 41, correcting Ruling 15): with the standard
#' one-step weight H = D D' (D the differencing matrix), GLS on the differenced equation equals
#' within-groups (FE) exactly, so difference GMM with many instruments is pulled towards FE; system GMM
#' is pulled mostly towards OLS in levels (a mix with about 0.13 weight on within-groups, D2 sample).
ols_diff_ar1 <- function(y, id, t) {
  .check_panel(y, id, t)
  dy <- y - lag_within(y, id, t, 1L); dyl <- lag_within(dy, id, t, 1L); ok <- !is.na(dy) & !is.na(dyl)
  b <- sum(dy[ok] * dyl[ok]) / sum(dyl[ok]^2)
  e <- dy[ok] - b * dyl[ok]
  s <- rowsum(dyl[ok] * e, .unit_index(id[ok]))
  list(est = b, se = sqrt(sum(s^2)) / sum(dyl[ok]^2), n = sum(ok))
}

## ---- card function: Anderson-Hsiao ---------------------------------------------------------
#' Anderson-Hsiao (1981, 1982) IV estimator rho_hat_AH of the first-differenced AR(1) equation
#'     Delta Y_it = rho Delta Y_{i,t-1} + Delta eps_it,
#' instrumenting Delta Y_{i,t-1} with Y_{i,t-2} (inst = "level", lag order s = 2 of the Arellano-Bond
#' moment E[Y_{i,t-s} Delta eps_it] = 0) or with Delta Y_{i,t-2} (inst = "diff", one more month lost).
#' Just identified, no intercept (no time effects, as pgmm(effect = "individual")); y is Y_it, id is i,
#' t is t. se is CR0 clustered by unit (Delta eps_it is MA(1) within a unit); with inst = "level" it
#' equals the robust SE of pgmm(y ~ lag(y, 1) | lag(y, 2:2), collapse = TRUE, model = "onestep").
#' Arguments: y is $Y_{it}$, id is $i$, t is $t$; inst picks the instrument $Y_{i,t-2}$ or $\Delta Y_{i,t-2}$;
#' est is $\hat\rho_{AH}$.
anderson_hsiao <- function(y, id, t, inst = c("level", "diff")) {
  inst <- match.arg(inst)
  .check_panel(y, id, t)
  dy <- y - lag_within(y, id, t, 1L); dyl <- lag_within(dy, id, t, 1L)
  z <- if (inst == "level") lag_within(y, id, t, 2L) else lag_within(dy, id, t, 2L)
  ok <- !is.na(dy) & !is.na(dyl) & !is.na(z)
  zx <- sum(z[ok] * dyl[ok])
  est <- sum(z[ok] * dy[ok]) / zx
  e <- dy[ok] - est * dyl[ok]
  s <- rowsum(z[ok] * e, .unit_index(id[ok]))
  se <- sqrt(sum(s^2)) / abs(zx)
  list(est = est, se = se, ci = ci95(est, se), n = sum(ok), inst = inst)
}

## ---- card functions: bias formulas and corrections -----------------------------------------
#' Exact probability limit (N -> infinity, T fixed) of rho_hat_FE - rho in the AR(1) panel with fixed
#' effects and a stationary start, Nickell (1981); T = number of regression periods (Y_i0 known, t = 1..T):
#'   -(1 + rho)/(T - 1) * [1 - (1/T)(1 - rho^T)/(1 - rho)]
#'     * [1 - 2 rho/((1 - rho)(T - 1)) * (1 - (1 - rho^T)/(T (1 - rho)))]^{-1}.
#' Checked against large simulated panels in tests/test_fns.R (lit_software §13.3: journal text not read).
#' Arguments: rho is $\rho$, T is $T$ (regression periods); returns plim $\hat\rho_{FE} - \rho$.
nickell_bias <- function(rho, T) {
  if (any(T < 2)) stop("nickell_bias: T (regression periods) must be at least 2")
  if (any(abs(rho) >= 1)) stop("nickell_bias: |rho| must be below 1")
  a <- (1 - rho^T) / (1 - rho)
  -(1 + rho) / (T - 1) * (1 - a / T) / (1 - 2 * rho / ((1 - rho) * (T - 1)) * (1 - a / T))
}

#' The first factor -(1 + rho)/(T - 1) of nickell_bias (T = regression periods). It is not the whole
#' leading-order bias: to order 1/T the bias is -(1 + rho)/T (Hahn and Kuersteiner), which is what the
#' legacy table row "nickell_approx" prints (it uses the number of months T_obs - 1 = T).
nickell_first_factor <- function(rho, T) -(1 + rho) / (T - 1)

#' Hahn and Kuersteiner analytic bias correction of the FE estimator rho_fe = rho_hat_FE:
#'   rho_hat_HK = ((T + 1)/T) rho_hat_FE + 1/T = rho_hat_FE + (1 + rho_hat_FE)/T,
#' Hahn and Kuersteiner (2001 working paper) eq. (6), from their Theorem 2 (journal version: Hahn and
#' Kuersteiner 2002, Econometrica 70(4)); T = number of regression periods. Vectorised.
#' Arguments: rho_fe is $\hat\rho_{FE}$, T is $T$ (regression periods); returns $\hat\rho_{HK}$.
hk_correct <- function(rho_fe, T) ((T + 1) / T) * rho_fe + 1 / T

#' Half-panel jackknife pieces (Dhaene and Jochmans 2015, eqs. 2.2 to 2.4): the FE estimate on all T
#' regression periods and on two consecutive halves in time (each half uses every unit; the first month
#' of a half supplies only the lag, so no regression period is lost).
#' A cut puts the first k regression periods in S1 and the other T - k in S2. Dhaene and Jochmans (2015),
#' eq. (2.2): the pieces of a cut are averaged with weights |S|/T, their shares of the T periods,
#'     bar_rho_cut = (|S1|/T) rho_hat_FE,S1 + (|S2|/T) rho_hat_FE,S2,
#' and eq. (2.4): rho_hat_SPJ = 2 rho_hat_FE - bar_rho_{1/2}, with bar_rho_{1/2} the average of bar_rho_cut over
#' the cuts: one cut (k = T/2) for even T, where the weights are 1/2 and 1/2 and nothing is weighted; two cuts
#' for odd T, (k, T - k) = (floor(T/2), ceiling(T/2)) and (ceiling(T/2), floor(T/2)), where the pieces have
#' (T -/+ 1)/2 periods and the weights are (T -/+ 1)/(2T), not 1/2 (fix round 2: the earlier code gave the
#' four half-panel estimates equal weight 1/4 for odd T, which is not their estimator).
#' Returns list(fe, halves (a 2-column matrix: one row per cut, columns S1 and S2), weights (same shape: |S|/T),
#' half_mean (bar_rho_{1/2}), spj, T).
spj_parts <- function(y, id, t) {
  .check_panel(y, id, t)
  yl <- lag_within(y, id, t, 1L); ok <- !is.na(y) & !is.na(yl)
  periods <- sort(unique(t[ok])); TT <- length(periods)
  if (TT < 2) stop("spj_fe: need at least two regression periods")
  fe_on <- function(keep) {
    g <- .unit_index(id[keep]); Yd <- .demean(y[keep], g); Ld <- .demean(yl[keep], g)
    sum(Ld * Yd) / sum(Ld^2)
  }
  fe <- fe_on(ok)
  cuts <- unique(c(floor(TT / 2), ceiling(TT / 2)))                 # k = |S1|; |S2| = T - k
  halves <- do.call(rbind, lapply(cuts, function(k) {
    s1 <- periods[seq_len(k)]
    c(S1 = fe_on(ok & t %in% s1), S2 = fe_on(ok & !(t %in% s1)))
  }))
  weights <- cbind(S1 = cuts / TT, S2 = (TT - cuts) / TT)             # |S|/T, DJ (2.2)
  half_mean <- mean(rowSums(weights * halves))                        # bar_rho_{1/2}, DJ (2.4)
  list(fe = fe, halves = halves, weights = weights, half_mean = half_mean, spj = 2 * fe - half_mean, T = TT)
}

#' Split-panel (half-panel) jackknife rho_hat_SPJ = 2 rho_hat_FE - bar_rho_{1/2}, the halves in time, Dhaene and
#' Jochmans (2015) eqs. (2.2) and (2.4); bar_rho_{1/2} weights each half by its share |S|/T of the regression
#' periods (equal weights for even T; for odd T the two cuts are averaged, see spj_parts); y is Y_it, id is i,
#' t is t. Returns a numeric scalar.
#' Arguments: y is $Y_{it}$, id is $i$, t is $t$; returns $\hat\rho_{SPJ}$ (halves $\hat\rho_{FE,1}$ and $\hat\rho_{FE,2}$,
#' notation_v3.md 1.9, each weighted by its share |S|/T of the regression periods $T$; equal weights 1/2 for even $T$).
spj_fe <- function(y, id, t) spj_parts(y, id, t)$spj

#' Number of instruments r of a plm::pgmm fit (the columns of each unit's instrument matrix Z_i:
#' GMM-type lags, plus the time dummies when effect = "twoways").
#' Arguments: fit is a pgmm fit; returns $r$, the number of columns of $\mathbf{Z}_i$.
n_instruments <- function(fit) {
  if (!inherits(fit, "pgmm")) stop("n_instruments: fit must be a plm::pgmm object")
  as.integer(ncol(fit$W[[1L]]))
}

#' Summary of a pgmm fit for the coefficient par: est; se = robust SE (vcovHC: the Windmeijer (2005)
#' corrected two-step variance for model = "twosteps", the robust one-step sandwich for "onestep");
#' se_conv = the uncorrected two-step SE (vcov; NA for one-step, whose plm vcov is unscaled); ci from
#' se; r = n_instruments; Hansen J (sargan(fit, "twosteps"): the quadratic form of the moments in the
#' inverse of sum_i Z_i' e_i e_i' Z_i, which summary.pgmm prints as "Sargan test" for a two-step fit and
#' which is the legacy "Sargan p 0.290") with df and p; the m_1 and m_2 serial-correlation tests (mtest
#' with the robust variance, as summary.pgmm).
#' There is no one-step Sargan statistic here: plm::sargan(fit, "onestep") returns the unscaled quadratic
#' form z' A1 z (A1 the one-step weight matrix, without the scale 1/(N sigma_eps^2)), e.g. 12,012.9
#' against a properly scaled 22.79 on the D1 sample, so its p-value is 0 for any sample; it is a plm
#' artifact and not a test, so it is neither computed nor stored (Task 11.1 fix round 1).
#' Arguments: fit is a pgmm fit ($\hat\rho^{(1)}_{AB}$, $\hat\rho^{(2)}_{AB}$ or $\hat\rho_{BB}$), par the coefficient
#' name; r is $r$, J is Hansen's $J$, m1 and m2 are $m_1$ and $m_2$.
gmm_summary <- function(fit, par = "lag(y, 1)") {
  V <- suppressWarnings(plm::vcovHC(fit))
  est <- unname(coef(fit)[par]); se <- unname(sqrt(V[par, par]))
  two <- identical(fit$args$model, "twosteps")
  J <- plm::sargan(fit, "twosteps")
  m1 <- suppressWarnings(plm::mtest(fit, order = 1L, vcov = V))
  m2 <- suppressWarnings(plm::mtest(fit, order = 2L, vcov = V))
  list(est = est, se = se, se_conv = if (two) unname(sqrt(stats::vcov(fit)[par, par])) else NA_real_,
       ci = ci95(est, se), r = n_instruments(fit), N = length(fit$W),
       J = unname(J$statistic), J_df = unname(J$parameter), J_p = unname(J$p.value),
       m1 = unname(m1$statistic), m1_p = unname(m1$p.value),
       m2 = unname(m2$statistic), m2_p = unname(m2$p.value))
}

## ---- data generators (new blocks; each script calls set.seed before a generator) ----------
#' Stack an N x T matrix by unit: data frame (id, t) with t = 1..T and the columns of `cols`.
.stack <- function(N, TT, cols) {
  out <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), times = N))
  for (nm in names(cols)) out[[nm]] <- if (is.matrix(cols[[nm]])) as.vector(t(cols[[nm]])) else rep(cols[[nm]], each = TT)
  out
}

#' Block D2: stationary AR(1) panel of key accounts, Y_it in thousands of USD.
#' m_i ~ N(m_mean, m_sd^2) is the long-run monthly GMV, alpha_i = (1 - rho) m_i, Y_i,start drawn from the
#' stationary distribution m_i + N(0, sigma_e^2/(1 - rho^2)), then Y_it = rho Y_{i,t-1} + alpha_i + eps_it,
#' eps ~ N(0, sigma_e^2), `burn` months discarded. Draw order: m (N), start (N), then eps month by month.
gen_D2 <- function(P) {
  m <- stats::rnorm(P$N, P$m_mean, P$m_sd); alpha <- (1 - P$rho) * m
  Tt <- P$T + P$burn; Y <- matrix(0, P$N, Tt)
  Y[, 1] <- m + stats::rnorm(P$N, 0, P$sigma_e / sqrt(1 - P$rho^2))
  for (s in 2:Tt) Y[, s] <- P$rho * Y[, s - 1] + alpha + stats::rnorm(P$N, 0, P$sigma_e)
  .stack(P$N, P$T, list(y = Y[, (P$burn + 1):Tt], m_true = m, alpha_true = alpha))
}

#' Block D2 misspecified twin (Monte Carlo of 04): the D2 design with MA(1) errors,
#'     eps_it = e_it + ma_coef e_{i,t-1},  e ~ N(0, sigma_e^2) the innovation,
#' so Y_{i,t-2} is correlated with Delta eps_it (through e_{i,t-2}): the lag-2 moment E[Y_{i,t-2} Delta eps_it] = 0
#' and the level moments E[Delta Y_{i,t-1}(alpha_i + eps_it)] = 0 are invalid, and a Hansen J test with power
#' must reject. m_i, alpha_i = (1 - rho) m_i and the burn-in as in gen_D2; the starting deviation has the
#' stationary variance sigma_e^2 (1 + 2 rho ma_coef + ma_coef^2)/(1 - rho^2) of the AR(1) deviation with
#' MA(1) errors. P is the D2 parameter list plus ma_coef. Draw order: m (N), e_{i,0} (N), start (N), then
#' e month by month.
gen_D2_ma <- function(P) {
  m <- stats::rnorm(P$N, P$m_mean, P$m_sd); alpha <- (1 - P$rho) * m
  Tt <- P$T + P$burn; Y <- matrix(0, P$N, Tt)
  e_prev <- stats::rnorm(P$N, 0, P$sigma_e)
  Y[, 1] <- m + stats::rnorm(P$N, 0, P$sigma_e * sqrt((1 + 2 * P$rho * P$ma_coef + P$ma_coef^2) / (1 - P$rho^2)))
  for (s in 2:Tt) {
    e <- stats::rnorm(P$N, 0, P$sigma_e)
    Y[, s] <- P$rho * Y[, s - 1] + alpha + e + P$ma_coef * e_prev
    e_prev <- e
  }
  .stack(P$N, P$T, list(y = Y[, (P$burn + 1):Tt], m_true = m, alpha_true = alpha))
}

#' Block D3: GMV with a lagged dependent variable and a predetermined regressor with feedback.
#' alpha_i ~ N(a_mean, a_sd^2); promotion spend (hundreds of USD)
#'     X_it = x_mean + x_alpha alpha_i + feedback eps_{i,t-1} + eps_x_it,  eps_x ~ N(0, eps_x_sd^2),
#' eps_x_it the promotion shock (a code symbol only: it is not the outcome error eps_it, and u is not an
#' error symbol in notation_v3.md),
#' and GMV (thousands of USD) Y_it = rho Y_{i,t-1} + beta X_it + alpha_i + eps_it, eps ~ N(0, sigma_e^2).
#' X_it is sequentially exogenous (it depends on eps_{i,t-1}, not on eps_it), not strictly exogenous.
#' Start at the long-run mean (alpha_i + beta (x_mean + x_alpha alpha_i))/(1 - rho) + eps, `burn` months
#' discarded. Draw order: alpha (N); month 1: eps, eps_x; then month by month eps, eps_x.
gen_D3 <- function(P) {
  alpha <- stats::rnorm(P$N, P$a_mean, P$a_sd)
  Tt <- P$T + P$burn; Y <- X <- E <- matrix(0, P$N, Tt)
  E[, 1] <- stats::rnorm(P$N, 0, P$sigma_e)
  X[, 1] <- P$x_mean + P$x_alpha * alpha + stats::rnorm(P$N, 0, P$eps_x_sd)
  Y[, 1] <- (alpha + P$beta * (P$x_mean + P$x_alpha * alpha)) / (1 - P$rho) + E[, 1]
  for (s in 2:Tt) {
    E[, s] <- stats::rnorm(P$N, 0, P$sigma_e)
    X[, s] <- P$x_mean + P$x_alpha * alpha + P$feedback * E[, s - 1] + stats::rnorm(P$N, 0, P$eps_x_sd)
    Y[, s] <- P$rho * Y[, s - 1] + P$beta * X[, s] + alpha + E[, s]
  }
  k <- (P$burn + 1):Tt
  .stack(P$N, P$T, list(y = Y[, k], x = X[, k], eps_true = E[, k], alpha_true = alpha))
}

#' Block D4: a cohort of newly onboarded merchants (non-stationary initial condition).
#' m_i ~ N(m_mean, m_sd^2) long-run monthly GMV (thousands of USD), alpha_i = (1 - rho) m_i. Month 1 on
#' Nimbus: Y_i1 = start_share m_i + w_i1, w_i1 ~ N(0, sigma_e^2/(1 - rho^2)) (start_share < 1: the merchant starts
#' below its long-run level, so Y_i1 - m_i is correlated with alpha_i and mean stationarity fails);
#' then Y_it = rho Y_{i,t-1} + alpha_i + eps_it. y_stat_true is the same merchants with the same draws
#' had they joined at their long-run level (start_share = 1): y_stat_true - y = rho^(t - 1) (1 - start_share) m_i.
#' (The start share is not called lambda: lambda_t is the time effect of notation_v3.md 1.8.)
#' Draw order: m (N), w (N), then eps month by month.
gen_D4 <- function(P) {
  m <- stats::rnorm(P$N, P$m_mean, P$m_sd); alpha <- (1 - P$rho) * m
  w <- stats::rnorm(P$N, 0, P$sigma_e / sqrt(1 - P$rho^2))
  Y <- Ys <- matrix(0, P$N, P$T)
  Y[, 1] <- P$start_share * m + w; Ys[, 1] <- m + w
  for (s in 2:P$T) {
    e <- stats::rnorm(P$N, 0, P$sigma_e)
    Y[, s] <- P$rho * Y[, s - 1] + alpha + e; Ys[, s] <- P$rho * Ys[, s - 1] + alpha + e
  }
  .stack(P$N, P$T, list(y = Y, y_stat_true = Ys, m_true = m, alpha_true = alpha))
}
