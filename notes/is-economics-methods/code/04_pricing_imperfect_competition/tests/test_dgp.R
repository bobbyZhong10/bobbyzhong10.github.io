# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Cirrus data and truth04.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("pricing_fns.R")
tr <- readRDS("truth04.rds")
C1 <- readRDS("cirrus_C1.rds"); C2 <- readRDS("cirrus_C2.rds"); C3 <- readRDS("cirrus_C3.rds")
P <- tr$params
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))

## ---- truth04 fields named on the chapter card ---------------------------------------------------
stopifnot(fin(tr$C1$alpha), tr$C1$alpha > 0, fin(tr$C1$delta, 2L), identical(names(tr$C1$delta), c("Cirrus", "Garnet")),
          fin(tr$C2$wtp, P$C2$N), identical(tr$C2$wtp, C2$wtp), identical(tr$C3$J, c(1, 2, 5, 20, 100)))
stopifnot(identical(tr$C1$alpha, P$C1$alpha), identical(tr$C1$delta, P$C1$delta))
# C4 is the C1 duopoly at its baseline (no shocks)
stopifnot(identical(P$C4$delta, P$C1$delta), identical(P$C4$alpha, P$C1$alpha), identical(P$C4$c, P$C1$c_base))
# seeds: integers, all different, one per simulated block (C4 and C0 use no random numbers)
sd <- unlist(P$seed)
stopifnot(identical(names(sd), c("C1", "C2", "C3", "C2_mc")), !anyDuplicated(sd), all(sd == round(sd)))

## ---- the generators reproduce the saved data from the stored seeds ------------------------------------
set.seed(P$seed$C1); stopifnot(identical(gen_C1(P$C1), C1))
set.seed(P$seed$C2); stopifnot(identical(gen_C2(P$C2), C2))
set.seed(P$seed$C3); stopifnot(identical(gen_C3(P$C3), C3))

## ---- C0: hand block, no random numbers -----------------------------------------------------------------
stopifnot(is.list(tr$C0), fin(tr$C0$pm), fin(tr$C0$profit_m), tr$C0$pm == 7)

## ---- C1: logit duopoly markets ---------------------------------------------------------------------------
M <- P$C1$M
stopifnot(M >= 500, nrow(C1) == M, identical(C1$market, seq_len(M)),
          all(c("c1", "c2", "xi1_true", "xi2_true", "a_hom_true", "p1_bertrand", "p2_bertrand", "s1_bertrand", "s2_bertrand",
                "p1_joint", "p2_joint", "s1_joint", "s2_joint", "p_cournot_hom", "q1_cournot_hom", "q2_cournot_hom",
                "p_bertrand_hom") %in% names(C1)))
hw <- P$C1$cost_halfwidth
stopifnot(all(abs(C1$c1 - P$C1$c_base[1]) <= hw), all(abs(C1$c2 - P$C1$c_base[2]) <= hw),
          abs(cor(C1$c1, C1$c2)) < 0.15, abs(cor(C1$c1, C1$xi1_true)) < 0.15, abs(cor(C1$c2, C1$xi1_true)) < 0.15,
          abs(sd(C1$xi1_true) / P$C1$xi_sd - 1) < 0.15)
# every market's prices satisfy the first-order conditions of their conduct, and the shares are the logit shares
for (m in seq_len(M)) {
  dm <- P$C1$delta + c(C1$xi1_true[m], C1$xi2_true[m]); cm <- c(C1$c1[m], C1$c2[m])
  pb <- c(C1$p1_bertrand[m], C1$p2_bertrand[m]); pj <- c(C1$p1_joint[m], C1$p2_joint[m])
  stopifnot(max(abs(logit_foc(pb, dm, P$C1$alpha, cm, 0))) < 1e-10, max(abs(logit_foc(pj, dm, P$C1$alpha, cm, 1))) < 1e-10,
            max(abs(logit_share(pb, dm, P$C1$alpha) - c(C1$s1_bertrand[m], C1$s2_bertrand[m]))) < 1e-14,
            max(abs(logit_share(pj, dm, P$C1$alpha) - c(C1$s1_joint[m], C1$s2_joint[m]))) < 1e-14)
}
stopifnot(all(C1$p1_joint > C1$p1_bertrand), all(C1$p2_joint > C1$p2_bertrand), all(C1$s1_bertrand + C1$s2_bertrand < 1))
# homogeneous-product variant: p = a_m - Q; Cournot p = (a + c1 + c2) / 3; Bertrand: the low-cost firm prices at the
# rival's cost (below its own monopoly price), so p = max(c1, c2)
stopifnot(max(abs(C1$p_cournot_hom - (C1$a_hom_true + C1$c1 + C1$c2) / 3)) < 1e-12,
          max(abs(C1$q1_cournot_hom - (C1$a_hom_true - 2 * C1$c1 + C1$c2) / 3)) < 1e-12,
          max(abs(C1$q2_cournot_hom - (C1$a_hom_true - 2 * C1$c2 + C1$c1) / 3)) < 1e-12,
          all(C1$q1_cournot_hom > 0), all(C1$q2_cournot_hom > 0),
          identical(C1$p_bertrand_hom, pmax(C1$c1, C1$c2)), all((C1$a_hom_true + pmin(C1$c1, C1$c2)) / 2 > C1$p_bertrand_hom))
# truth: pass-through matrices at the baseline and averaged over the markets
for (k in c("pt_bertrand_base", "pt_joint_base", "pt_bertrand_mean", "pt_joint_mean"))
  stopifnot(is.matrix(tr$C1[[k]]), identical(dim(tr$C1[[k]]), c(2L, 2L)), all(is.finite(tr$C1[[k]])))
stopifnot(max(abs(tr$C1$pt_bertrand_base - passthrough_logit(P$C1$delta, P$C1$alpha, P$C1$c_base, 0))) < 1e-12)

## ---- C2: merchants and Cirrus Pro ----------------------------------------------------------------------------
N <- P$C2$N
stopifnot(N >= 20000, nrow(C2) == N, identical(C2$customer, seq_len(N)), sum(C2$train) == N / 2,
          identical(names(C2), c("customer", "sales_z", "tenure_z", "category", "wtp", "train", "index_true", "eps_true")),
          identical(levels(C2$category), P$C2$cat_levels))
b <- P$C2$beta
idx <- b[["const"]] + b[["sales"]] * C2$sales_z + b[["tenure"]] * C2$tenure_z +
  b[["electronics"]] * (C2$category == "electronics") + b[["home"]] * (C2$category == "home")
stopifnot(max(abs(C2$index_true - idx)) < 1e-12, identical(C2$wtp, pmax(1, round(exp(C2$index_true + C2$eps_true)))),
          all(C2$wtp >= 1), all(C2$wtp == round(C2$wtp)))
stopifnot(max(abs(as.numeric(table(C2$category)) / N - P$C2$cat_probs)) < 0.02, abs(sd(C2$eps_true) / P$C2$sigma - 1) < 0.03,
          abs(cor(C2$sales_z, C2$tenure_z)) < 0.03, abs(cor(C2$sales_z, C2$eps_true)) < 0.03)
# the random split balances the features
stopifnot(abs(mean(C2$sales_z[C2$train]) - mean(C2$sales_z[!C2$train])) < 0.05,
          abs(mean(C2$wtp[C2$train]) / mean(C2$wtp[!C2$train]) - 1) < 0.05)
# willingness to pay: median near exp(const), features explain about half of log WTP
stopifnot(abs(median(C2$wtp) / exp(b[["const"]]) - 1) < 0.1,
          var(C2$index_true) / var(C2$index_true + C2$eps_true) > 0.4, var(C2$index_true) / var(C2$index_true + C2$eps_true) < 0.6)

## ---- C3: independent U[0, 1] values, bundle sums -----------------------------------------------------------
stopifnot(P$C3$N >= 50000, is.matrix(C3), identical(dim(C3), c(as.integer(P$C3$N), 5L)),
          identical(colnames(C3), paste0("J", P$C3$J)))
for (k in seq_along(P$C3$J)) {
  J <- P$C3$J[k]; x <- C3[, k]
  stopifnot(all(x >= 0), all(x <= J), abs(mean(x) - J / 2) < 4 * sqrt(J / 12 / P$C3$N), abs(var(x) / (J / 12) - 1) < 0.05)
  if (k > 1) stopifnot(all(x >= C3[, k - 1]))                     # nested bundles: the first J goods
}
