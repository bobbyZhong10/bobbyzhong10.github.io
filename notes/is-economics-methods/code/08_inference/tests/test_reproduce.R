# ------------------------------------------------------------------
# tests/test_reproduce.R -- the inherited Halcyon code (v2, seeds 626, 6261, 313, 99, 101, 202,
# 11, 12, 414, 818) reproduces every number chapter 8 prints, at the printed rounding.
# Field names follow the res_*.rds files (03 stores s40/s10 as named vectors iid, cr1_z, cr1_t).
# One deliberate change from v2: the wild cluster bootstrap p is 0.0742 (76 of the 1024 sign
# patterns of G = 10 clusters, all enumerated, ties counted) and prints as 0.074, not the v2 text's 0.073
# (the frozen v2 run counted 145 of 1999 random draws): the all +1 and all -1 patterns reproduce |t_obs|
# exactly, and floating point decided whether such a tie counted (numbers_def.R). The conclusion (not
# significant at 5%) is the same.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
H <- readRDS("halcyon.rds"); dl <- readRDS("res_delta.rds"); cl <- readRDS("res_cluster.rds")
hc <- readRDS("res_hac.rds"); te <- readRDS("res_test.rds")
# x printed with d decimals as v
at <- function(x, v, d) isTRUE(abs(as.numeric(x) - v) <= 0.5 * 10^-d + 1e-12)

## DGP (01): truths of section 6.1 and the Moulton factor of section 3.4
tr <- H$truth
stopifnot(at(tr$WTP, 0.20, 2), at(tr$b_time, -0.08, 2), at(tr$b_fee, -0.40, 2), at(tr$tau, 0.25, 2),
          at(tr$rho, 0.083, 3), at(tr$moulton, 4.18, 2), nrow(H$clus_main) == 8000, tr$G == 40, tr$n_c == 200)

## Delta against bootstrap (02): section 6.2 and the table of 6.7
m <- dl$main; w <- dl$weak
stopifnot(at(m$wtp_hat, 0.207, 3), at(m$se_delta, 0.0105, 4), at(m$se_boot, 0.0106, 4),
          at(m$ci_delta[1], 0.186, 3), at(m$ci_delta[2], 0.227, 3),
          at(m$ci_perc[1], 0.188, 3), at(m$ci_perc[2], 0.229, 3),
          at(m$ci_pt[1], 0.188, 3), at(m$ci_pt[2], 0.228, 3))
stopifnot(at(dl$base_w["bf"], -0.246, 3), at(dl$base_w["se_bf"], 0.105, 3), at(w$se_delta, 0.138, 3),
          at(w$skew, 8.0, 1), at(w$ci_delta[1], 0.05, 2), at(w$ci_delta[2], 0.59, 2),
          at(w$ci_perc[1], 0.17, 2), at(w$ci_perc[2], 1.51, 2), at(w$ci_pt[1], 0.24, 2), at(w$ci_pt[2], 0.82, 2))

## Cluster (03): sections 1, 3.4, 4.2, 6.3, 6.4
stopifnot(at(cl$tau_hat, 0.318, 3), at(cl$se_iid, 0.0231, 4), at(cl$se_hc1, 0.0231, 4),
          at(cl$se_cr1, 0.0908, 4), at(cl$se_cr2, 0.0920, 4),
          at(cl$tau_hat / cl$se_iid, 13.76, 2), at(cl$tau_hat / cl$se_hc1, 13.76, 2),
          at(cl$tau_hat / cl$se_cr1, 3.50, 2), at(cl$se_cr1 / cl$se_iid, 3.93, 2))
p_iid <- 2 * pnorm(-abs(cl$tau_hat / cl$se_iid))
stopifnot(abs(p_iid - 4e-43) <= 0.5e-43)                          # printed 4 x 10^-43
stopifnot(at(cl$s40[["iid"]], 0.645, 3), at(cl$s40[["cr1_z"]], 0.070, 3), at(cl$s40[["cr1_t"]], 0.064, 3),
          at(cl$s10[["iid"]], 0.636, 3), at(cl$s10[["cr1_z"]], 0.103, 3), at(cl$s10[["cr1_t"]], 0.067, 3))
stopifnot(at(cl$t10, 2.21, 2), at(cl$p10_z, 0.027, 3), at(cl$p10_t, 0.054, 3))
stopifnot(abs(cl$p10_wcb - 76 / 1024) < 1e-12, at(cl$p10_wcb, 0.074, 3), cl$p10_wcb > 0.05)

## HAC (04): sections 3.5 and 6.5
stopifnot(at(hc$b, 0.189, 3), at(hc$se_iid, 0.0484, 4), at(hc$se_nw, 0.0613, 4), hc$L == 5,
          at(hc$b / hc$se_iid, 3.91, 2), at(hc$b / hc$se_nw, 3.08, 2), at(hc$se_nw / hc$se_iid, 1.27, 2),
          at(hc$cov_iid, 0.849, 3), at(hc$cov_nw, 0.920, 3))

## Testing (05): sections 4.1 and 6.6
stopifnot(at(te$W, 35.58, 2), at(te$LR, 36.93, 2), at(te$LM, 36.47, 2),
          at(te$W_direct, 10.68, 2), at(te$W_repar, 6.79, 2), at(te$pWd, 0.0011, 4), at(te$pWr, 0.0092, 4))
