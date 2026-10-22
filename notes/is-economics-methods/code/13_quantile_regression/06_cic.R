# ------------------------------------------------------------------
# 06_cic.R -- blocks Q5 and Q5f: changes-in-changes (chapter 13)
#   cic, qdid  Q5 (sharp): CIC QTT (cic_qtt, Athey and Imbens 2006) and quantile DiD (qdid) on the grid
#              0.05, ..., 0.95; bootstrap that resamples within each region-quarter cell, 999 draws,
#              seed 13501, the same draws for both estimators; intervals est -/+ 1.96 se.
#              x$Q5$cic$ci and x$Q5$qdid$ci: lists of c(lo, hi), one per tau in taus (the card's test)
#   lqte       Q5f (fuzzy): the switchers' LQTE at tau = 0.5 (cic_lqte, dCDH 2018 Theorem 2.3), bootstrap
#              within cells, 999 draws, seed 13502; point estimates at every tau (no bootstrap); the
#              Wald-CIC (switchers' LATE) = num / first_stage from the same call
#   cross      toy data sets (columns y, d, g, t) and this chapter's Wald-CIC numerators
#              E(Y11) - E(Q_{D10}(Y10)) for the cross-check against chapter 20's wald_cic, which only
#              tests/test_conclusions.R reads (the card: ../20_did_beyond_binary/ is read only in tests)
# Output: res_extra.rds elements Q5 and cross (save_extra).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("qr_fns.R")
tr <- readRDS("truth13.rds"); P <- tr$params; taus <- P$taus; taus_fine <- P$taus_fine
Q5 <- readRDS("marigold_Q5.rds"); Q5f <- readRDS("marigold_Q5f.rds")
z975 <- qnorm(0.975)

## ---- Q5: sharp CIC and QDID ------------------------------------------------------------------
cells <- lapply(list(c00 = c(0, 0), c01 = c(0, 1), c10 = c(1, 0), c11 = c(1, 1)), function(st) Q5$y[Q5$S == st[1] & Q5$T == st[2]])
both <- function(cl) c(cic_qtt(cl$c00, cl$c01, cl$c10, taus_fine, cl$c11)$qtt, qdid(cl$c00, cl$c01, cl$c10, cl$c11, taus_fine)$qdid)
est <- both(cells)
set.seed(P$seed$Q5_boot)
boot <- vapply(seq_len(999), function(b) both(lapply(cells, function(v) v[sample.int(length(v), replace = TRUE)])),
               numeric(2 * length(taus_fine)))
se <- apply(boot, 1, sd)
nf <- length(taus_fine)
grid <- data.frame(tau = taus_fine, cic_est = est[1:nf], cic_se = se[1:nf], qdid_est = est[nf + 1:nf], qdid_se = se[nf + 1:nf],
                   qtt_true = tr$Q5$fine$qtt, qdid_true = tr$Q5$fine$qdid)
grid$cic_lo <- grid$cic_est - z975 * grid$cic_se; grid$cic_hi <- grid$cic_est + z975 * grid$cic_se
grid$qdid_lo <- grid$qdid_est - z975 * grid$qdid_se; grid$qdid_hi <- grid$qdid_est + z975 * grid$qdid_se
g5 <- grid[grid$tau %in% taus, ]
cic_q <- cic_qtt(cells$c00, cells$c01, cells$c10, taus, cells$c11)
stopifnot(max(abs(cic_q$qtt - g5$cic_est)) < 1e-12)

## ---- Q5f: switchers' LQTE --------------------------------------------------------------------
pf <- P$Q5f
set.seed(P$seed$Q5f_boot)
lq <- cic_lqte(Q5f, "y", "D", "S", "T", pf$tau, B = 999)
lq_grid <- vapply(taus, function(t) cic_lqte(Q5f, "y", "D", "S", "T", t, B = 0)$est, numeric(1))
drates <- with(Q5f, tapply(D, list(S, T), mean))

Q5res <- list(taus = taus, B = ncol(boot),
              cic = list(est = g5$cic_est, se = g5$cic_se, ci = lapply(seq_len(nrow(g5)), function(i) c(g5$cic_lo[i], g5$cic_hi[i])),
                         cf = cic_q$cf, q11 = cic_q$q11),
              qdid = list(est = g5$qdid_est, se = g5$qdid_se, ci = lapply(seq_len(nrow(g5)), function(i) c(g5$qdid_lo[i], g5$qdid_hi[i]))),
              grid = grid,
              lqte = list(tau = pf$tau, est = lq$est, se = lq$se, ci = lq$ci, num = lq$num, first_stage = lq$first_stage,
                          wald = lq$wald, se_wald = lq$se_wald, ci_wald = lq$ci_wald, B = lq$B),
              lqte_grid = data.frame(tau = taus, est = lq_grid, truth = tr$Q5$lqte_grid),
              sample = list(nc = P$Q5$nc, means = vapply(cells, mean, numeric(1)), medians = vapply(cells, emp_quantile, numeric(1), tau = 0.5),
                            nc_f = pf$nc, d_rate = c(d00 = drates["0", "0"], d01 = drates["0", "1"], d10 = drates["1", "0"], d11 = drates["1", "1"])))
save_extra("Q5", Q5res)
cat("Saved res_extra.rds: Q5\n")
print(g5, digits = 4)
cat(sprintf("Q5f LQTE_sw(%.2f) = %.4f [%.4f, %.4f] (truth %.4f); Wald-CIC %.4f [%.4f, %.4f] (truth LATE_sw %.4f); first stage %.4f\n",
            pf$tau, lq$est, lq$ci[1], lq$ci[2], tr$Q5$lqte, lq$wald, lq$ci_wald[1], lq$ci_wald[2], tr$Q5$late_sw, lq$first_stage))
print(cbind(tau = taus, est = lq_grid, truth = tr$Q5$lqte_grid), digits = 4)

## ---- cross-chapter toys: Wald-CIC numerator E(Y11) - E(Q_{D10}(Y10)) --------------------------
## toy_int: integer outcomes with ties and always-takers in both groups (hand-sized);
## toy_equal: equal cell sizes, y_d00 = 1:25 and y_d01 = 101:125 among the untreated (the floating-point
##   off-by-one case of quantile(type = 1)); toy_random: unequal cells, continuous outcomes, seed 13601;
## q5f: the chapter's own fuzzy data set.
toy_int <- data.frame(
  g = rep(c(0, 0, 1, 1), c(10, 10, 6, 7)), t = rep(c(0, 1, 0, 1), c(10, 10, 6, 7)),
  y = c(1:10, (1:10) * 2 + 1, c(2, 4, 4, 7, 9, 11), c(5, 8, 10, 10, 15, 20, 22)),
  d = c(rep(0:1, c(7, 3)), rep(0:1, c(7, 3)), c(0, 0, 0, 0, 1, 1), c(0, 0, 1, 1, 1, 1, 1)))
toy_equal <- data.frame(
  g = rep(c(0, 0, 1, 1), each = 30), t = rep(c(0, 1, 0, 1), each = 30),
  y = c(c(1:25, 40:44), c(101:125, 140:144), c(3:24, 41:48), c(110:131, 145:152)),
  d = c(rep(0:1, c(25, 5)), rep(0:1, c(25, 5)), rep(0:1, c(22, 8)), rep(0:1, c(10, 20))))
set.seed(P$seed$cross_toy)
nr <- c(37, 41, 43, 47)
toy_random <- data.frame(g = rep(c(0, 0, 1, 1), nr), t = rep(c(0, 1, 0, 1), nr))
toy_random$v <- runif(nrow(toy_random))
toy_random$d <- as.integer((toy_random$g == 0 & toy_random$v > 0.75) |
                             (toy_random$g == 1 & (toy_random$v > 0.8 | (toy_random$t == 1 & toy_random$v > 0.35))))
toy_random$y <- exp(0.4 * rnorm(nrow(toy_random)) + 0.2 * toy_random$t + 0.3 * toy_random$g) * (1 + 0.3 * toy_random$d)
toy_random$v <- NULL
q5f_toy <- data.frame(y = Q5f$y, d = Q5f$D, g = Q5f$S, t = Q5f$T)
toys <- list(toy_int = toy_int, toy_equal = toy_equal, toy_random = toy_random, q5f = q5f_toy)
num13 <- vapply(toys, function(d) cic_lqte(d, "y", "d", "g", "t", 0.5, B = 0)$num, numeric(1))
save_extra("cross", list(toys = toys, num13 = num13))
cat("Saved res_extra.rds: cross (Wald-CIC numerators", paste(signif(num13, 6), collapse = ", "), ")\n")
