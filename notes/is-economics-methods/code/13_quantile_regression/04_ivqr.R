# ------------------------------------------------------------------
# 04_ivqr.R -- block Q3: endogenous AdBoost adoption, IV-QR and the complier QTE (chapter 13)
#   first_stage  E(D | Z = 1) - E(D | Z = 0)
#   qr     naive QR of log revenue on D and X (nid standard errors): biased, adopters have higher ranks
#   ivqr   Chernozhukov-Hansen inverse QR (ivqr_grid): for each b on the grid -0.30, -0.298, ..., 0.60,
#          QR of log Y - D b on (1, X, Z); estimate where the Z coefficient's t statistic (nid) is closest
#          to zero; 95% confidence set = grid points with |t| <= 1.96 (dual inference). The target is the
#          structural effect beta_D(tau) of the whole population under rank similarity (iv_tab$struct_true),
#          not the complier effect; the LQTE (iv_tab$lqte_true) equals beta_D only at tau = 0.5 (both 0.10 by
#          the symmetric design), so the IV-QR interval is compared with the LQTE at the median only
#   lqte   complier QTE from the Wald-type complier cdfs (lqte_iv, no covariates: Z is randomised),
#          bootstrap that resamples within each Z arm, 999 draws, seed 13301; the same draws give the
#          Wald (LATE) estimate's standard error
# Output: res_extra.rds element Q3 (save_extra).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("qr_fns.R")
suppressPackageStartupMessages(library(quantreg))
tr <- readRDS("truth13.rds"); P <- tr$params; P3 <- P$Q3; taus <- P$taus
Q3 <- readRDS("marigold_Q3.rds")
z975 <- qnorm(0.975)
grid <- round(seq(-0.30, 0.60, by = P3$grid_step), 3)

fs <- mean(Q3$D[Q3$Z == 1]) - mean(Q3$D[Q3$Z == 0])

## ---- naive QR and IV-QR at every tau ---------------------------------------------------------
qr_tab <- do.call(rbind, lapply(taus, function(t) {
  fit <- .muffle_nonunique(rq(logy ~ D + x, tau = t, data = Q3))
  cf <- .muffle_nonunique(summary(fit, se = "nid"))$coefficients
  data.frame(tau = t, est = cf["D", 1], se = cf["D", 2], lo = cf["D", 1] - z975 * cf["D", 2], hi = cf["D", 1] + z975 * cf["D", 2])
}))
iv_list <- lapply(taus, function(t) ivqr_grid(Q3$logy, Q3$D, Q3$Z, cbind(x = Q3$x), t, grid))
iv_tab <- data.frame(tau = taus, est = vapply(iv_list, `[[`, numeric(1), "est"),
                     lo = vapply(iv_list, function(r) r$ci[1], numeric(1)), hi = vapply(iv_list, function(r) r$ci[2], numeric(1)),
                     connected = vapply(iv_list, `[[`, logical(1), "connected"),
                     at_boundary = vapply(iv_list, `[[`, logical(1), "at_boundary"),
                     struct_true = tr$Q3$struct, lqte_true = tr$Q3$lqte)
if (any(iv_tab$at_boundary)) stop("04_ivqr.R: an IV-QR confidence set reaches the grid edge; widen the grid")
# the accepted set can have small holes where the t statistic of the Z coefficient (continuous but not
# monotone in b: the coefficient is piecewise linear and its SE moves with b) crosses 1.96 more than once;
# the reported interval is then its convex hull (range of accepted points)
if (!all(iv_tab$connected)) message("04_ivqr.R: accepted IV-QR set not an interval at tau = ",
                                    paste(iv_tab$tau[!iv_tab$connected], collapse = ", "), " (convex hull reported)")
iv_tab$n_accepted <- vapply(iv_list, function(r) sum(abs(r$tstat) <= z975), numeric(1))
iv_tab$n_hull <- round((iv_tab$hi - iv_tab$lo) / P3$grid_step) + 1

## ---- complier QTE and the Wald estimate, bootstrap -------------------------------------------
wald <- function(ly, d, z) (mean(ly[z == 1]) - mean(ly[z == 0])) / (mean(d[z == 1]) - mean(d[z == 0]))
lq <- lqte_iv(Q3$logy, Q3$D, Q3$Z, taus)
w_est <- wald(Q3$logy, Q3$D, Q3$Z)
i1 <- which(Q3$Z == 1); i0 <- which(Q3$Z == 0)
set.seed(P$seed$Q3_boot)
boot <- vapply(seq_len(999), function(b) {
  i <- c(i1[sample.int(length(i1), replace = TRUE)], i0[sample.int(length(i0), replace = TRUE)])
  c(lqte_iv(Q3$logy[i], Q3$D[i], Q3$Z[i], taus)$lqte, wald(Q3$logy[i], Q3$D[i], Q3$Z[i]))
}, numeric(length(taus) + 1))
se_b <- apply(boot, 1, sd)
lqte_tab <- data.frame(tau = taus, est = lq$lqte, se = se_b[seq_along(taus)], lo = lq$lqte - z975 * se_b[seq_along(taus)],
                       hi = lq$lqte + z975 * se_b[seq_along(taus)], q1c = lq$q1c, q0c = lq$q0c, truth = tr$Q3$lqte)

k <- which(taus == P3$tau)
Q3res <- list(tau = P3$tau, first_stage = fs,
              qr = list(est = qr_tab$est[k], se = qr_tab$se[k], ci = c(qr_tab$lo[k], qr_tab$hi[k])),
              ivqr = list(est = iv_tab$est[k], ci = c(iv_tab$lo[k], iv_tab$hi[k]), grid_step = P3$grid_step,
                          grid_range = range(grid), tstat_at_truth = iv_list[[k]]$tstat[which.min(abs(grid - tr$Q3$lqte_median))]),
              lqte = list(est = lqte_tab$est[k], se = lqte_tab$se[k], ci = c(lqte_tab$lo[k], lqte_tab$hi[k])),
              wald = list(est = w_est, se = se_b[length(se_b)], ci = w_est + c(-1, 1) * z975 * se_b[length(se_b)]),
              qr_tab = qr_tab, iv_tab = iv_tab, lqte_tab = lqte_tab, B = ncol(boot),
              sample = list(n = nrow(Q3), n_z1 = sum(Q3$Z), share_d_z1 = mean(Q3$D[Q3$Z == 1]), share_d_z0 = mean(Q3$D[Q3$Z == 0]),
                            share_complier_true = mean(Q3$type_true == "complier"),
                            median_y = emp_quantile(Q3$y, 0.5)))
save_extra("Q3", Q3res)
cat(sprintf("Saved res_extra.rds: Q3 (first stage %.3f)\n", fs))
print(qr_tab, digits = 4); print(iv_tab, digits = 4); print(lqte_tab, digits = 4)
cat(sprintf("Wald (LATE, log points) %.4f se %.4f\n", w_est, se_b[length(se_b)]))
