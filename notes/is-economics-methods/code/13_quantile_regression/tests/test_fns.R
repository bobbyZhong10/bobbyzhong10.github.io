# ------------------------------------------------------------------
# tests/test_fns.R -- unit tests of qr_fns.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# The first block is the chapter card's assertions, verbatim (comments translated);
# the rest checks every hand-written estimator against exact hand computations and,
# where a package implements the same estimator, against that package:
#   quantreg 6.1 (rq), rifreg 1.1.0 (rifreg), qte 2.0.0 (unc_qte, the deprecated CiC, cic).
# ------------------------------------------------------------------
source("qr_fns.R")
suppressPackageStartupMessages(library(quantreg))

## ---- chapter card (verbatim) ---------------------------------------------------
sup <- c(1, 2, 3); Fv <- c(0.2, 0.5, 1)                                             # discrete distribution: the generalized-inverse properties hold exactly
stopifnot(gen_inverse(sup, Fv, 0.5) == 2, gen_inverse(sup, Fv, 0.51) == 3,
          all(sapply(c(0.1, 0.3, 0.7), function(u) Fv[match(gen_inverse(sup, Fv, u), sup)] >= u)))
set.seed(1); y <- rexp(500); ys <- sort(y)                                          # the minimiser lies between the 150th and 151st order statistics
qq <- optimize(function(q) sum(check_loss(y - q, 0.3)), range(y))$minimum
stopifnot(qq >= ys[150] - 1e-4, qq <= ys[151] + 1e-4)

## ---- check function rho_tau(u) = (tau - 1{u < 0}) u ---------------------------------------
stopifnot(check_loss(2, 0.3) == 0.6, check_loss(-2, 0.3) == 1.4, check_loss(0, 0.3) == 0,
          all(check_loss(c(-1, 1), 0.5) == 0.5), all(check_loss(seq(-3, 3, 0.5), 0.9) >= 0))

## ---- generalized inverse: Wacker (2023) Lemma 1 (g)(3) and (i) on random discrete cdfs ------
## F^{-1}(F(y)) <= y for every support point y; F(F^{-1}(u)) >= u for every u in (0, 1]
set.seed(2)
for (r in 1:50) {
  s <- sort(unique(round(rnorm(8), 2))); p <- diff(c(0, sort(runif(length(s) - 1)), 1)); Fs <- cumsum(p)
  Fs[length(Fs)] <- 1
  stopifnot(all(gen_inverse(s, Fs, Fs) <= s))
  u <- runif(40)
  stopifnot(all(Fs[match(gen_inverse(s, Fs, u), s)] >= u))
}
# a flat stretch of F: the generalized inverse picks the smallest y with F(y) >= u
stopifnot(gen_inverse(c(1, 2, 3, 4), c(0.25, 0.5, 0.5, 1), c(0.5, 0.500001)) == c(2, 4))
stopifnot(inherits(try(gen_inverse(c(2, 1), c(0.5, 1), 0.5), silent = TRUE), "try-error"),
          inherits(try(gen_inverse(1:2, c(0.5, 1), 0), silent = TRUE), "try-error"))

## ---- emp_quantile = generalized inverse of the empirical cdf (= quantile(type = 1) in exact arithmetic) ---
## emp_quantile equals gen_inverse(sort(y), (1:n) / n, tau) exactly. R's quantile(type = 1) agrees with it at
## generic tau (the loop below) but not where n * tau rounds just above an integer (the deterministic
## case after the loop), where quantile(type = 1) returns the next order statistic.
set.seed(3)
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
for (n in c(1, 2, 7, 100, 1499, 1500, 2500, 5001)) {
  yy <- rnorm(n); tt <- c(taus, 1 / n, 1, runif(5))
  stopifnot(identical(emp_quantile(yy, tt), unname(quantile(yy, tt, type = 1))),
            identical(emp_quantile(yy, tt), gen_inverse(sort(yy), seq_len(n) / n, tt)))
}
yt <- c(5, 1, 1, 3)                                                                   # ties
stopifnot(identical(emp_quantile(yt, c(0.25, 0.5, 0.51, 0.75, 1)), c(1, 1, 3, 3, 5)))
## deterministic case where emp_quantile and quantile(type = 1) differ: y = 1:800, tau = 0.55. Here
## n * tau = 440.00000000000006 (not 440), so quantile(type = 1) returns the 441st order statistic, while
## F_n(440) = 440 / 800 = 0.55 >= tau gives F_n^{-1}(0.55) = 440 (the definition inf{y : F_n(y) >= tau}).
y800 <- 1:800
stopifnot(800 * 0.55 > 440,
          emp_quantile(y800, 0.55) == 440,
          emp_quantile(y800, 0.55) == gen_inverse(sort(y800), seq_len(800) / 800, 0.55))
## R 4.6.1's quantile(y800, 0.55, type = 1) returns 441 here; reported, not asserted (base R may change this)
message(sprintf("  note: quantile(type = 1) at n = 800, tau = 0.55 gives %g (emp_quantile 440)",
                unname(quantile(y800, 0.55, type = 1))))
## at the levels k / n of the empirical cdf's own jumps, the generalized inverse of y = 1:n is exactly k
## (n = 2500: 145 of these levels are ones where quantile(type = 1) is one order statistic too high)
y25 <- 1:2500; k25 <- seq_len(2500)
stopifnot(identical(emp_quantile(y25, k25 / 2500), as.numeric(k25)),
          identical(emp_quantile(y25, k25 / 2500), gen_inverse(y25, k25 / 2500, k25 / 2500)))
message(sprintf("  note: quantile(type = 1) differs from k at %d of the 2500 levels k / n (R %s)",
                sum(unname(quantile(y25, k25 / 2500, type = 1)) != k25), getRversion()))

## ---- equivariance to a monotone transformation, exact in a saturated design ------------------
## Sample quantiles commute with any increasing transformation (order statistics are preserved), so in a
## saturated design (one binary regressor) the QR fit of log(y) is the log of the QR fit of y exactly.
## Wooldridge (2010) section 12.10.1 (p. 453) states this only for the median; the statement for general
## tau needs another source (Koenker and Bassett 1978 on equivariance; theorem number to be checked).
## Arm sizes 301 and 299 make n tau non-integer in each arm, so every sample quantile is unique.
set.seed(4)
de <- rep(0:1, c(301, 299)); ye <- exp(0.3 * de + rnorm(600, 0, 0.8))
for (t in taus) {
  a <- coef(rq(ye ~ de, tau = t)); b <- coef(rq(log(ye) ~ de, tau = t))
  q0 <- emp_quantile(ye[de == 0], t); q1 <- emp_quantile(ye[de == 1], t)
  stopifnot(max(abs(a - c(q0, q1 - q0))) < 1e-10,                                   # saturated QR = arm quantiles
            max(abs(b - c(log(a[1]), log(a[1] + a[2]) - log(a[1])))) < 1e-10,        # log of the QR fit
            abs(emp_quantile(log(ye), t) - log(emp_quantile(ye, t))) < 1e-14)
}

## ---- qte_marginal: difference of arm quantiles; qte::unc_qte (qte 2.0.0) as cross-check ------
set.seed(5)
dq <- rep(0:1, c(1001, 999)); yq <- ifelse(dq == 1, rnorm(2000, 1, 2), rnorm(2000))
set.seed(51); qm <- qte_marginal(yq, dq, taus, B = 200)
stopifnot(is.data.frame(qm), identical(names(qm)[1:4], c("tau", "est", "lo", "hi")), nrow(qm) == 5,
          max(abs(qm$est - (emp_quantile(yq[dq == 1], taus) - emp_quantile(yq[dq == 0], taus)))) < 1e-14,
          all(qm$lo < qm$est & qm$est < qm$hi), all(abs((qm$hi - qm$lo) / (2 * qnorm(0.975) * qm$se) - 1) < 1e-12))
set.seed(51); qm2 <- qte_marginal(yq, dq, taus, B = 200)
stopifnot(identical(qm, qm2))                                                          # reproducible under set.seed
# qte 2.0.0's unc_qte(xformla = ~1) minimises the check function numerically (BMisc::weighted_quantile calls
# optimize(), tolerance about 1.2e-4), so it agrees with the exact arm quantiles only to that tolerance; with
# odd arm sizes n tau is never an integer and the minimiser is unique (with an even arm, optimize() can return
# any point between two order statistics: differences up to 0.008 on a toy with an arm of 1000)
suppressPackageStartupMessages(library(qte))
uq <- unc_qte(yname = "y", dname = "d", data = data.frame(y = yq, d = dq), xformla = ~1, probs = taus,
              target = "qte", biters = 2, cband = FALSE)
stopifnot(max(abs(uq$qte - qm$est)) < 1e-3)

## ---- rif_quantile: RIF regression = rifreg 1.1.0 exactly ----------------------------------
set.seed(6)
dr <- data.frame(w = rbinom(3001, 1, 0.3), x = rnorm(3001)); dr$y <- 1 + 0.5 * dr$w + 0.6 * dr$x + rnorm(3001, 0, 0.4 + 0.3 * dr$w)
rf <- rifreg::rifreg(y ~ w + x, data = dr, statistic = "quantiles", probs = taus)
mine <- sapply(taus, function(t) coef(lm(rif_quantile(dr$y, t) ~ w + x, data = dr)))
stopifnot(max(abs(unname(rf$estimates) - unname(mine))) < 1e-10)
rf_sj <- rifreg::rifreg(y ~ w + x, data = dr, statistic = "quantiles", probs = 0.9, bw = "SJ")
stopifnot(max(abs(unname(rf_sj$estimates[, 1]) - coef(lm(rif_quantile(dr$y, 0.9, bw = "SJ") ~ w + x, data = dr)))) < 1e-10)
r9 <- rif_quantile(dr$y, 0.9)
q9 <- attr(r9, "quantile"); f9 <- attr(r9, "density")
stopifnot(length(r9) == nrow(dr), abs(mean(r9) - (q9 + (0.9 - mean(dr$y <= q9)) / f9)) < 1e-10,   # sample mean of the RIF
          q9 == quantile(dr$y, 0.9, type = 7, names = FALSE), f9 > 0,
          abs(attr(r9, "bw") - bw.nrd0(dr$y)) < 1e-15)

## ---- ivqr_grid: with D = Z (full compliance) the inverse-QR zero is the QR coefficient ----------
## rq(y - d b ~ x + z) with z = d has z-coefficient beta_QR - b exactly (regression equivariance), so the
## grid estimate is the grid point nearest the ordinary QR coefficient on d. Both arms have an odd size
## (601): with a binary regressor the QR solution is unique only when n_arm tau is not an integer.
set.seed(7)
ni <- 1202; xi <- rnorm(ni); zi <- sample(rep(0:1, c(601, 601))); yi <- 1 + 0.4 * zi + 0.5 * xi + rnorm(ni, 0, 0.7)
b_qr <- coef(rq(yi ~ zi + xi, tau = 0.5))[["zi"]]
step <- 0.002; gr <- seq(0, 0.8, by = step)
iv <- ivqr_grid(yi, zi, zi, cbind(x = xi), 0.5, gr)
stopifnot(abs(iv$est - b_qr) <= step / 2 + 1e-12, iv$ci[1] < b_qr, iv$ci[2] > b_qr, iv$connected, !iv$at_boundary,
          max(abs(iv$coef_z - (b_qr - gr))) < 1e-8)
# the same with X = NULL (intercept only)
iv0 <- ivqr_grid(yi, zi, zi, NULL, 0.5, gr)
stopifnot(abs(iv0$est - coef(rq(yi ~ zi, tau = 0.5))[["zi"]]) <= step / 2 + 1e-12)

## ---- lqte_iv: with D = Z the complier distributions are the arm distributions ----------------
lq <- lqte_iv(yq, dq, dq, taus)
stopifnot(max(abs(lq$lqte - qm$est)) < 1e-12)

## ---- dist_reg: one binary-response GLM per threshold -------------------------------------
Xd <- cbind(w = dr$w, x = dr$x); th <- quantile(dr$y, c(0.2, 0.5, 0.8), names = FALSE)
g <- dist_reg(dr$y, Xd, th, link = "probit")
stopifnot(is.matrix(g), all(dim(g) == c(3, 3)), identical(colnames(g), c("(Intercept)", "w", "x")))
for (j in seq_along(th))
  stopifnot(max(abs(g[j, ] - coef(glm(I(dr$y <= th[j]) ~ w + x, data = dr, family = binomial("probit"))))) < 1e-8)
stopifnot(all(diff(g[, 1]) > 0))                                                      # intercept rises with the threshold
Fd <- dr_cdf(g, Xd, "probit")
stopifnot(length(Fd) == 3, all(diff(Fd) > 0), max(abs(Fd - c(0.2, 0.5, 0.8))) < 0.02)
gl <- dist_reg(dr$y, Xd, th)                                                          # default link: logit
stopifnot(max(abs(gl[2, ] - coef(glm(I(dr$y <= th[2]) ~ w + x, data = dr, family = binomial("logit"))))) < 1e-8)

## ---- canay_qr and fe_dummy_qr on a small panel ---------------------------------------------
set.seed(8)
Np <- 20; Tp <- 5; idp <- rep(seq_len(Np), each = Tp); ap <- rnorm(Np)[idp]
xp <- 0.5 * ap + rnorm(Np * Tp); yp <- ap + 0.7 * xp + rnorm(Np * Tp, 0, 0.5)
cq <- canay_qr(yp, cbind(x = xp), idp, 0.333)
bw <- coef(lm(yp ~ xp + factor(idp)))[["xp"]]                                         # within (FE) slope
ah <- ave(yp - bw * xp, idp)
stopifnot(abs(cq$beta_mean - bw) < 1e-10, max(abs(cq$alpha_hat - ah)) < 1e-10,
          max(abs(cq$beta - coef(rq(I(yp - ah) ~ xp, tau = 0.333)))) < 1e-10, identical(names(cq$beta), c("(Intercept)", "x")))
fd <- fe_dummy_qr(yp, cbind(x = xp), idp, 0.333)
dense <- coef(rq(yp ~ xp + factor(idp) - 1, tau = 0.333))[["xp"]]
stopifnot(abs(fd$beta[["x"]] - dense) < 1e-5)

## ---- cic_qtt and qdid: a hand-computed toy ---------------------------------------------
## Q_{d=0}(y) = F^{-1}_{Y01}(F_{Y00}(y)): Q(2.5) = F01^{-1}(0.5) = 20, Q(3.5) = 30, Q(4.5) = 40
h00 <- c(1, 2, 3, 4); h01 <- c(10, 20, 30, 40); h10 <- c(2.5, 3.5, 4.5); h11 <- c(25, 36, 50)
ch <- cic_qtt(h00, h01, h10, c(0.2, 0.5, 0.9), h11)
stopifnot(is.data.frame(ch), identical(names(ch)[1:2], c("tau", "qtt")),
          identical(ch$qtt, c(5, 6, 10)), identical(ch$cf, c(20, 30, 40)))
stopifnot(identical(qdid(h00, h01, h10, h11, 0.5)$qdid, 36 - 3.5 - 20 + 2))
# the off-by-one case fixed in chapter 20's .fz_qq: equal cell sizes, y_d00 = 1:25, y_d01 = 101:125
stopifnot(identical(qq_map(7, 1:25, 101:125), 107L), identical(qq_map(c(0, 25, 26), 1:25, 101:125), c(101L, 125L, 125L)))

## ---- cic_qtt against qte 2.0.0 -------------------------------------------------------
## Toy in the spirit of lit_software 13.5 item 8: Y = exp(0.5 U + 0.3 1{t = 2}) + 1{treated, t = 2}.
## Cell sizes 1499, 1500, 1501, 1502 (n00 and n01 coprime, so n01 F00(y) is never an integer and the
## deprecated CiC's floating-point quantile(type = 1) cannot hit the off-by-one case): the deprecated CiC must
## equal cic_qtt to 1e-10; the new cic() uses another quantile convention and is only checked to 0.02.
set.seed(9)
nc <- c(1499, 1500, 1501, 1502)
uu <- list(runif(nc[1]), runif(nc[2]), runif(nc[3])^0.7, runif(nc[4])^0.7)
yc <- c(exp(0.5 * uu[[1]]), exp(0.5 * uu[[2]] + 0.3), exp(0.5 * uu[[3]]), exp(0.5 * uu[[4]] + 0.3) + 1)
dc <- data.frame(y = yc, treat = rep(c(0, 0, 1, 1), nc), t = rep(c(1, 2, 1, 2), nc))
dc$g <- ifelse(dc$treat == 1, 2, 0)
cell <- function(s, t) dc$y[dc$treat == s & dc$t == t]
mine_c <- cic_qtt(cell(0, 1), cell(0, 2), cell(1, 1), taus, cell(1, 2))
old <- suppressWarnings(CiC(y ~ treat, t = 2, tmin1 = 1, tname = "t", data = dc, panel = FALSE, se = FALSE, probs = taus))
new <- suppressMessages(cic(yname = "y", gname = "g", tname = "t", data = dc, panel = FALSE, probs = taus,
                            gt_type = "qtt", biters = 2, cband = FALSE))
stopifnot(max(abs(unname(old$qte) - mine_c$qtt)) < 1e-10, max(abs(new$overall$qtt - mine_c$qtt)) < 0.02)

## ---- cic_lqte: sharp design and the counterfactual numerator -----------------------------------
## With D = S T (no always-takers, every treated-group unit switches) the switchers' LQTE is the CIC QTT:
## F_CIC,1 = F_Y11 and F_CIC,0 = the cdf of Q_{d=0}(Y10).
ds <- data.frame(y = yc, g = rep(c(0, 0, 1, 1), nc), t = rep(c(0, 1, 0, 1), nc)); ds$d <- ds$g * ds$t
for (t in taus) {
  lt <- cic_lqte(ds, "y", "d", "g", "t", t, B = 0)
  stopifnot(abs(lt$est - mine_c$qtt[mine_c$tau == t]) < 1e-10, is.null(lt$ci))
}
num_direct <- mean(cell(1, 2)) - mean(qq_map(cell(1, 1), cell(0, 1), cell(0, 2)))   # E(Y11) - E(Q_{D10}(Y10))
stopifnot(abs(lt$num - num_direct) < 1e-12, abs(lt$first_stage - 1) < 1e-15, abs(lt$wald - num_direct) < 1e-12)
# fuzzy toy: bootstrap interval brackets the estimate and is reproducible under set.seed
set.seed(10)
nf <- 800; df <- data.frame(g = rep(c(0, 0, 1, 1), each = nf), t = rep(c(0, 1, 0, 1), each = nf))
df$v <- runif(nrow(df)); df$z <- rnorm(nrow(df)) + 0.5 * df$g
df$d <- as.integer((df$g == 0 & df$v > 0.8) | (df$g == 1 & (df$v > 0.85 | (df$t == 1 & df$v > 0.4))))
df$y <- exp(0.5 * df$z + 0.1 * df$t) * (1 + 0.2 * df$d)
set.seed(11); f1 <- cic_lqte(df, "y", "d", "g", "t", 0.5, B = 99)
set.seed(11); f2 <- cic_lqte(df, "y", "d", "g", "t", 0.5, B = 99)
stopifnot(identical(f1, f2), f1$ci[1] < f1$est, f1$est < f1$ci[2], abs(f1$first_stage - 0.45) < 0.06)
