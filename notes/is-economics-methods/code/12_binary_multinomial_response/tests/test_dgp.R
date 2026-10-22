# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Helix data and truth12.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("bin_fns.R")
tr <- readRDS("truth12.rds"); P <- tr$params
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))
md5_lines <- function(lines) { f <- tempfile(); on.exit(unlink(f)); writeLines(lines, f); unname(tools::md5sum(f)) }

## ---- truth12 fields named on the chapter card ---------------------------------------------------
stopifnot(fin(tr$B1$beta, 2L), fin(tr$B1$ame), fin(tr$B2$avg_cross), fin(tr$B3$ape), fin(tr$B4$beta), fin(tr$B4$ape),
          fin(tr$B5$beta), fin(tr$B6$ape_top), is.matrix(tr$B7$coef), all(dim(tr$B7$coef) == c(2, 3)), all(is.finite(tr$B7$coef)))
# every block has its own seed, all different; the inherited block keeps 512, 9001 and 778
sd <- unlist(P$seed)
stopifnot(length(sd) >= 22, !anyDuplicated(sd), all(sd == round(sd)),
          P$seed$B1 == 512, P$seed$B1_oracle == 9001, P$seed$B1_consistency == 778)
# Monte Carlo seed ranges (base + r) do not overlap each other or any block seed
mc <- list(B5_T4 = P$seed$B5_mc_T4 + seq_len(P$B5$R_mc), B5_T8 = P$seed$B5_mc_T8 + seq_len(P$B5$R_mc),
           B7 = P$seed$B7_mc + seq_len(P$B7$R_mc), B6 = P$seed$B6_mc + seq_len(P$B6$R_mc))
stopifnot(!anyDuplicated(c(unlist(mc), sd[!grepl("_mc", names(sd))])))

## ---- inherited block B1: legacy code byte-identical and first; data regenerate exactly ------------------
src <- readLines("01_dgp.R")
stopifnot(md5_lines(src[1:77]) == "d20a281baafc8f4eb82a437d451859fc",          # = the legacy 01_dgp.R
          md5_lines(readLines("02_lpm_naive.R")) == "fcca5906181707afefb11e8e49a5a283",
          md5_lines(readLines("03_mle_logit.R")) == "a9d62a7df1bbf4063caf818d047daa13",
          grep("set.seed", src)[1] == 30, grep("set.seed", src)[2] == 53, grep("^gen_helix <- function", src) == 39,
          grep("set.seed\\(778\\)", readLines("03_mle_logit.R")) == 61)
e1 <- new.env(); eval(parse(text = src[33:49]), e1)                          # legacy parameters and gen_helix
set.seed(512); h1 <- e1$gen_helix(e1$n_main)
hd <- readRDS("helix_data.rds"); or <- readRDS("oracle.rds")
stopifnot(identical(h1, hd), nrow(hd) == 4000, identical(names(hd), c("x", "D", "pD", "y", "mu")),
          identical(c(or$a0, or$a1), tr$B1$beta), tr$B1$ame == or$ame_logit, tr$B1$base_rate == or$base_rate)
# B1 truth: quadrature agrees with the legacy 5,000,000-draw oracle and with stats::integrate
ig <- function(f) integrate(f, -Inf, Inf, rel.tol = 1e-12)$value
b <- tr$B1$beta
stopifnot(abs(tr$B1$ame_exact - ig(function(x) b[2] * dlogis(b[1] + b[2] * x) * dnorm(x))) < 1e-10,
          abs(tr$B1$base_rate_exact - ig(function(x) plogis(b[1] + b[2] * x) * dnorm(x))) < 1e-10,
          abs(tr$B1$ame_exact - tr$B1$ame) < 1e-3, abs(tr$B1$base_rate_exact - tr$B1$base_rate) < 1e-3,
          abs(tr$B1$lpm_slope_pop - tr$B1$ame_exact) < 1e-10,                  # Stein's lemma (X normal)
          abs(tr$B1$mem - b[2] * dlogis(b[1])) < 1e-15)
# probit pseudo-true value: first-order conditions of the expected probit log likelihood hold
pq <- tr$B1$probit_pseudo
foc <- sapply(1:2, function(k) integrate(function(x) { u <- pq[1] + pq[2] * x; pl <- plogis(b[1] + b[2] * x)
  (pl - pnorm(u)) * dnorm(u) / (pnorm(u) * pnorm(-u)) * (if (k == 1) 1 else x) * dnorm(x) }, -10, 10, rel.tol = 1e-12)$value)
stopifnot(max(abs(foc)) < 1e-7, abs(tr$B1$logit_over_probit_pseudo - b[2] / pq[2]) < 1e-12,
          tr$B1$logit_over_probit_pseudo > 1.5, tr$B1$logit_over_probit_pseudo < 1.8,
          abs(tr$B1$probit_ape_pseudo - tr$B1$ame_exact) < 0.01)            # APEs agree across links, coefficients do not

## ---- the generators reproduce the saved data from the stored seeds --------------------------------
D <- list(B2 = readRDS("helix_B2.rds"), B3 = readRDS("helix_B3.rds"), B4 = readRDS("helix_B4.rds"),
          B5 = readRDS("helix_B5.rds"), B6 = readRDS("helix_B6.rds"), B7 = readRDS("helix_B7.rds"),
          B7n = readRDS("helix_B7n.rds"))
for (bk in c("B2", "B3", "B4", "B5", "B6", "B7")) {
  set.seed(P$seed[[bk]])
  stopifnot(identical(get(paste0("gen_", bk))(P[[bk]]), D[[bk]]))
}
set.seed(P$seed$B7_nested)
stopifnot(identical(gen_B7(P$B7, lambda = P$B7$lambda_nested), D$B7n))

## ---- B2: two features, interaction --------------------------------------------------------------
B2 <- D$B2; P2 <- P$B2
stopifnot(nrow(B2) == P2$n, setequal(unique(B2$y), 0:1), abs(cor(B2$x1, B2$x2) - P2$rho_x) < 0.05,
          P2$beta[4] > 0, tr$B2$avg_cross < 0,                                     # the Ai-Norton sign reversal in the truth
          tr$B2$share_pos > 0.2, tr$B2$share_pos < 0.5,                              # the sign varies across firms
          tr$B2$rate > 0.65, tr$B2$rate < 0.8, abs(mean(B2$y) - tr$B2$rate) < 0.03,
          tr$B2$naive_int_me > 0,
          tr$B2$large$N >= 1e6, abs(tr$B2$large$avg_cross - tr$B2$avg_cross) < 5e-4,
          abs(tr$B2$large$share_pos - tr$B2$share_pos) < 5e-3, abs(tr$B2$large$rate - tr$B2$rate) < 2e-3)

## ---- B3: endogenous engagement, randomised reminder ---------------------------------------------------
B3 <- D$B3; P3 <- P$B3
stopifnot(nrow(B3) == P3$n, setequal(unique(B3$Z), 0:1), abs(mean(B3$Z) - 0.5) < 0.03,
          max(abs(B3$spend_kusd - P3$spend_median_kusd * exp(B3$X))) < 1e-12,
          abs(median(B3$spend_kusd) - 5) < 0.5,                                   # median monthly cloud spend about USD 5,000
          abs(cor(B3$eps_true, B3$v_true) - P3$r) < 0.05, abs(var(B3$eps_true) - 1) < 0.08,
          max(abs(B3$D - (P3$pi[1] + P3$pi[2] * B3$X + P3$pi[3] * B3$Z + B3$v_true))) < 1e-12,
          identical(B3$Y, as.integer(P3$beta[1] + P3$beta[2] * B3$X + P3$beta_D * B3$D + B3$eps_true > 0)))
fs <- summary(lm(D ~ X + Z, data = B3))$coefficients
stopifnot(fs["Z", "t value"]^2 > 100)                                              # a strong instrument
stopifnot(abs(tr$B3$large$ape - tr$B3$ape) < 1e-3, abs(tr$B3$large$rate - tr$B3$rate) < 2e-3,
          tr$B3$large$naive_ape > 1.4 * tr$B3$ape,                                 # endogeneity overstates the APE
          abs(tr$B3$scaled[3] - P3$beta_D / sqrt(1 - P3$r^2)) < 1e-12, abs(tr$B3$rho - P3$r / P3$sigma_v) < 1e-12)

## ---- B4: omitted heterogeneity independent of X ---------------------------------------------------
B4 <- D$B4; P4 <- P$B4
stopifnot(nrow(B4) == P4$n, abs(cor(B4$X, B4$a_true)) < 0.05, abs(sd(B4$a_true) - P4$sigma_a) < 0.05,
          abs(tr$B4$beta_att - P4$beta[2] / sqrt(1 + P4$sigma_a^2)) < 1e-12, tr$B4$beta_att < 0.9 * tr$B4$beta,
          tr$B4$large$N >= 1e6, abs(tr$B4$large$ape - tr$B4$ape) < 1e-3, abs(tr$B4$large$rate - tr$B4$rate) < 2e-3,
          max(abs(tr$B4$large$short_coef - tr$B4$beta_att_vec)) < 5e-3)              # eq. 15.35 in a large sample
# the APE in closed form = the average of the attenuated-probit derivative by quadrature
s4 <- tr$B4$scale
ng <- normal_grid()
stopifnot(abs(tr$B4$ape - sum(ng$w * (P4$beta[2] / s4) * dnorm((P4$beta[1] + P4$beta[2] * ng$z) / s4))) < 1e-10)

## ---- B5: monthly adoption panel ---------------------------------------------------------------------
B5 <- D$B5; P5 <- P$B5
stopifnot(nrow(B5) == P5$N * P5$T, all(table(B5$id) == P5$T), all(B5$t == rep(1:P5$T, P5$N)),
          setequal(unique(B5$Y), 0:1))
ab <- tapply(B5$alpha_true, B5$id, function(z) z[1]); xb <- tapply(B5$X, B5$id, mean)
stopifnot(all(tapply(B5$alpha_true, B5$id, function(z) diff(range(z))) == 0),   # alpha_i is constant over months
          cor(ab, xb) > 0.4, abs(coef(lm(ab ~ xb))[2] - P5$xi) < 0.15)             # correlated with Xbar_i (Mundlak)
stopifnot(tr$B5$T8$pooled_plim > P5$beta + 0.1, tr$B5$T4$pooled_plim > P5$beta + 0.1,    # pooled logit is biased
          tr$B5$T8$N == P5$N_large, tr$B5$ape > 0.1, tr$B5$ape < 0.15, abs(tr$B5$ape - tr$B5$ape_T4) < 0.01,
          abs(mean(B5$Y) - tr$B5$T8$rate) < 0.03, P5$R_mc >= 1000)
mv <- tapply(B5$Y, B5$id, function(z) sum(z) > 0 && sum(z) < P5$T)
stopifnot(abs(mean(mv) - tr$B5$T8$share_movers) < 0.05, tr$B5$T4$share_movers < tr$B5$T8$share_movers)

## ---- B6: satisfaction rating -------------------------------------------------------------------------
B6 <- D$B6; P6 <- P$B6
stopifnot(nrow(B6) == P6$n, is.ordered(B6$rating), identical(levels(B6$rating), as.character(1:5)),
          all(table(B6$rating) > 0.05 * P6$n), abs(median(B6$bill_kusd) - P6$bill_median_kusd) < 0.1,
          identical(as.integer(B6$rating), findInterval(B6$ystar_true, P6$cuts, left.open = TRUE) + 1L),
          all(diff(P6$cuts) > 0))
stopifnot(abs(sum(tr$B6$shares) - 1) < 1e-12, max(abs(rowSums(tr$B6$shares_G) - 1)) < 1e-12,
          max(abs(tr$B6$large$shares - tr$B6$shares)) < 2e-3, abs(tr$B6$large$ape_top - tr$B6$ape_top) < 5e-4,
          max(abs(as.vector(table(B6$rating)) / P6$n - tr$B6$shares)) < 0.03,
          # managed accounts: more 1s and more 5s (the larger latent spread), a higher latent mean
          tr$B6$shares_G["managed", "r1"] > tr$B6$shares_G["self_serve", "r1"],
          tr$B6$shares_G["managed", "r5"] > tr$B6$shares_G["self_serve", "r5"],
          tr$B6$bl_latent_diff > 0, tr$B6$bl_log_ratio_h > 0)                       # Bond and Lang: the ranking flips under h

## ---- B7: plan choice ---------------------------------------------------------------------------------
B7 <- D$B7; B7n <- D$B7n; P7 <- P$B7
stopifnot(nrow(B7) == P7$n, nrow(B7n) == P7$n, identical(levels(B7$choice), P7$levels),
          identical(rownames(tr$B7$coef), c("Pro", "Enterprise")), identical(tr$B7$coef, P7$coef),
          max(abs(as.vector(table(B7$choice)) / P7$n - tr$B7$shares)) < 0.03,
          max(abs(as.vector(table(B7n$choice)) / P7$n - tr$B7$shares_nested)) < 0.03,
          max(abs(tr$B7$large$shares - tr$B7$shares)) < 2e-3, max(abs(tr$B7$large$ape_x - tr$B7$ape_x)) < 1e-3,
          abs(sum(tr$B7$ape_x)) < 1e-12, all(tr$B7$shares > 0.2),
          P7$lambda_nested < 1, P7$drop == P7$levels[1], P7$R_mc >= 1000)
# the nested variant changes substitution, not the utilities: same x, s distribution, other shares
stopifnot(max(abs(c(mean(B7n$x), mean(B7n$s)))) < 0.05, abs(tr$B7$shares_nested[3] - tr$B7$shares[3]) > 0.03)
# the coefficient of engagement on Enterprise is positive, its APE on P(Enterprise) is negative
stopifnot(P7$coef["Enterprise", "x"] > 0, tr$B7$ape_x["Enterprise"] < 0)

## ---- fix round 1 and 2: two B5 designs at T = 4 (controller rulings 42 and 44) ------------------------------
## "T4" = a 4-month DGP (alpha_i built on the 4-month mean; used by the Monte Carlo and dgp_B5_T4_*);
## "T8first4" = months 1 to 4 of the T = 8 panel (alpha_i built on the 8-month mean; the one-sample T = 4 panel).
## Ruling 44: (X_it, alpha_i) has the same bivariate normal law in every month, so months 1 to 4 of the T = 8 design have
## exactly the T = 8 pooled-logit limit, true APE, pooled-APE limit and rate. The truths are the 2-D Gauss-Hermite
## quadrature values (b5_quadrature); the large draws are kept as checks ($draw), within their Monte Carlo error.
tf <- tr$B5$T8first4; t8 <- tr$B5$T8; t4 <- tr$B5$T4
stopifnot(!is.null(tf), tf$T == 4, tf$T_dgp == P5$T, tf$N == P5$N_large, t8$N == P5$N_large, t4$N == P5$N_large,
          is.finite(tf$ape), is.finite(tf$pooled_plim), is.finite(tf$pooled_ape_plim), length(tf$cre_probit_plim) == 3)
q8 <- b5_quadrature(P5, 8, 120); q4 <- b5_quadrature(P5, 4, 120)
stopifnot(abs(t8$ape - q8$ape) < 1e-12, abs(t8$pooled_plim - q8$pooled_plim) < 1e-12, abs(t8$pooled_ape_plim - q8$pooled_ape_plim) < 1e-12,
          abs(t8$rate - q8$rate) < 1e-12, abs(tr$B5$ape - q8$ape) < 1e-12,
          abs(t4$ape - q4$ape) < 1e-12, abs(t4$pooled_plim - q4$pooled_plim) < 1e-12, abs(t4$pooled_ape_plim - q4$pooled_ape_plim) < 1e-12,
          abs(tr$B5$ape_T4 - q4$ape) < 1e-12,
          # months 1 to 4 of the T = 8 panel share the T = 8 values
          tf$ape == t8$ape, tf$pooled_plim == t8$pooled_plim, tf$pooled_ape_plim == t8$pooled_ape_plim, tf$rate == t8$rate,
          # the values of the controller's re-review (quad.R)
          abs(t8$pooled_plim - 0.97690) < 5e-6, abs(t8$ape - 0.124527) < 5e-7, abs(t8$pooled_ape_plim - 0.18877) < 5e-6,
          abs(t4$pooled_plim - 1.04498) < 5e-6, abs(t4$ape - 0.121431) < 5e-7)
# the two designs differ: the pooled-logit limit is farther from beta in the 4-month DGP, the true APE is lower
stopifnot(abs(tf$pooled_plim - t4$pooled_plim) > 0.02, abs(tf$ape - tr$B5$ape_T4) > 0.002, P$seed$B5_large_first4 == 1235)
# stable in the number of nodes
q8s <- b5_quadrature(P5, 8, 60)
stopifnot(max(abs(unlist(q8) - unlist(q8s))) < 1e-9)
# the large draws (checks only): within their Monte Carlo error of the quadrature values (sd about 0.003 for the pooled limit,
# 0.00013 for the APE at N_large = 200,000; the bounds are about five sd) and the stored draws are labelled
for (L in list(t8, t4, tf)) {
  stopifnot(!is.null(L$draw), abs(L$draw$pooled_plim - L$pooled_plim) < 0.015, abs(L$draw$ape - L$ape) < 7e-4,
            abs(L$draw$rate - L$rate) < 0.005, abs(L$draw$pooled_ape_plim - L$pooled_ape_plim) < 0.01)
}
stopifnot(abs(tf$draw$pooled_plim - t8$pooled_plim) < 0.015)                       # draw of months 1 to 4 vs the T = 8 limit
# the stored one-sample panel's first four months are the first four months of the T = 8 panel; the sample truth is close to tf
mm4 <- B5[B5$t <= 4, ]
stopifnot(nrow(mm4) == P5$N * 4, abs(mean(P5$beta * dlogis(mm4$alpha_true + P5$beta * mm4$X)) - tf$ape) < 0.01)
stopifnot(tf$share_movers < t8$share_movers)
