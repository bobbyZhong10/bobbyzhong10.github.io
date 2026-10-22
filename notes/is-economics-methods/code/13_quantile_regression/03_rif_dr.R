# ------------------------------------------------------------------
# 03_rif_dr.R -- block Q2: conditional QR against the unconditional quantile effect (chapter 13)
#   cqr   QR of log revenue on FBM (W) and log catalogue size (X); coefficient on W with nid standard
#         errors (Hendricks-Koenker sandwich, Hall-Sheather bandwidth): the effect on the conditional
#         quantile Quant_tau(log Y | W, X); truth 0.5 at every tau
#   rif   OLS of RIF(log Y; tau) on (1, W, X) (Firpo, Fortin and Lemieux 2009): the coefficient on W
#         estimates UQPE(tau), the effect of a larger FBM share on the unconditional tau-quantile.
#         Standard error: bootstrap of the whole procedure (quantile, kernel density, OLS), 999 draws
#         that resample merchants, seed 13201; se_ols is the OLS sandwich (HC0) that treats the RIF as
#         known (what summary(rifreg) reports, with its warning), shown for comparison
#   dr    distribution regression (probit, regressors 1, W, X, W X: correctly specified here) at 199
#         thresholds (sample quantiles 0.005 to 0.995); the counterfactual distribution with the FBM
#         share raised by dp = 0.1, holding F(Y | W, X) and F(X) fixed; quantiles by linear
#         interpolation of the cdf after a monotone (cummax) fix; uqe_discrete = (q_cf - q_fit) / dp. Point estimates only.
# Output: res_extra.rds element Q2 (save_extra).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("qr_fns.R")
suppressPackageStartupMessages(library(quantreg))
tr <- readRDS("truth13.rds"); P <- tr$params; P2 <- P$Q2; taus <- P$taus; taus_fine <- P$taus_fine
Q2 <- readRDS("marigold_Q2.rds")
z975 <- qnorm(0.975)
n <- nrow(Q2)
M <- cbind(1, fbm = Q2$fbm, x = Q2$x)
qrM <- qr(M)

## ---- RIF regression: point estimates and analytic (RIF-as-known) standard errors on the fine grid ----
rif_coef <- function(ly, Mq) {
  R <- vapply(taus_fine, function(t) rif_quantile(ly, t), numeric(length(ly)))
  qr.coef(Mq, R)[2, ]                                                    # coefficient on W, one per tau
}
rif_est <- rif_coef(Q2$logy, qrM)
rif_se_ols <- vapply(taus_fine, function(t) {
  fit <- lm(rif_quantile(Q2$logy, t) ~ fbm + x, data = Q2)
  sqrt(sandwich::vcovHC(fit, type = "HC0")["fbm", "fbm"])
}, numeric(1))
set.seed(P$seed$Q2_boot)
boot <- vapply(seq_len(999), function(b) {
  i <- sample.int(n, replace = TRUE)
  rif_coef(Q2$logy[i], qr(M[i, ]))
}, numeric(length(taus_fine)))
rif_se <- apply(boot, 1, sd)

## ---- conditional QR -------------------------------------------------------------------------
cqr <- t(vapply(taus_fine, function(t) {
  fit <- .muffle_nonunique(rq(logy ~ fbm + x, tau = t, data = Q2))
  summary(fit, se = "nid")$coefficients["fbm", 1:2]
}, numeric(2)))

grid <- data.frame(tau = taus_fine, rif_est = rif_est, rif_se = rif_se, rif_lo = rif_est - z975 * rif_se,
                   rif_hi = rif_est + z975 * rif_se, rif_se_ols = rif_se_ols,
                   cqr_est = cqr[, 1], cqr_se = cqr[, 2], cqr_lo = cqr[, 1] - z975 * cqr[, 2], cqr_hi = cqr[, 1] + z975 * cqr[, 2],
                   uqpe_true = tr$Q2$fine$uqpe, cqr_true = tr$Q2$fine$cqr)
k <- which(taus_fine == P2$tau_uqe)

## ---- distribution regression and the discrete counterfactual -------------------------------------
th <- quantile(Q2$logy, seq(0.005, 0.995, by = 0.005), type = 1, names = FALSE)
XD <- cbind(fbm = Q2$fbm, x = Q2$x, fbm_x = Q2$fbm * Q2$x)
gam <- dist_reg(Q2$logy, XD, th, link = "probit")
X1 <- cbind(fbm = 1, x = Q2$x, fbm_x = Q2$x); X0 <- cbind(fbm = 0, x = Q2$x, fbm_x = 0)
F1 <- dr_cdf(gam, X1, "probit"); F0 <- dr_cdf(gam, X0, "probit")
p_hat <- mean(Q2$fbm)
dr_quant <- function(Fv, t) {
  # monotone (cummax) fix of the estimated cdf: running maximum over the ascending thresholds. It is not the
  # Chernozhukov, Fernandez-Val and Galichon (2010) rearrangement, which sorts the values; the two differ in
  # general, but F_fit and F_cf below are already non-decreasing on this sample (asserted after they are
  # built), so cummax, sorting and no fix at all give the same cdf here.
  Fm <- cummax(Fv)
  approx(Fm, th, xout = t, ties = "ordered")$y
}
F_fit <- p_hat * F1 + (1 - p_hat) * F0; F_cf <- (p_hat + P2$dp) * F1 + (1 - p_hat - P2$dp) * F0
stopifnot(identical(cummax(F_fit), sort(F_fit)), identical(cummax(F_cf), sort(F_cf)))   # cummax fix = sorting here
dr <- data.frame(tau = taus, q_fit = vapply(taus, dr_quant, numeric(1), Fv = F_fit),
                 q_cf = vapply(taus, dr_quant, numeric(1), Fv = F_cf))
dr$uqe_discrete <- (dr$q_cf - dr$q_fit) / P2$dp
dr$uqe_discrete_true <- tr$Q2$uqe_discrete
dr$q_sample <- emp_quantile(Q2$logy, taus)

Q2res <- list(tau = P2$tau_uqe,
              rif = list(est = rif_est[k], se = rif_se[k], ci = c(rif_est[k] - z975 * rif_se[k], rif_est[k] + z975 * rif_se[k]),
                         se_ols = rif_se_ols[k]),
              cqr = list(est = cqr[k, 1], se = cqr[k, 2], ci = c(cqr[k, 1] - z975 * cqr[k, 2], cqr[k, 1] + z975 * cqr[k, 2])),
              grid = grid, dr = dr, dr_thresholds = length(th), dr_n_extreme = attr(gam, "n_extreme"), B = ncol(boot),
              sample = list(n = n, n_fbm = sum(Q2$fbm), share_fbm = p_hat, median_y = emp_quantile(Q2$y, 0.5),
                            p90_y = emp_quantile(Q2$y, 0.9), p10_y = emp_quantile(Q2$y, 0.1)))
save_extra("Q2", Q2res)
cat(sprintf("Saved res_extra.rds: Q2 (%d of %d DR threshold fits had fitted probabilities numerically 0 or 1)\n", attr(gam, "n_extreme"), length(th)))
print(grid[grid$tau %in% taus, ], digits = 4)
print(dr, digits = 4)
