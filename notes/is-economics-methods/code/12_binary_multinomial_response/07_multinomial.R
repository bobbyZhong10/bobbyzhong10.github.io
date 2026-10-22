# ------------------------------------------------------------------
# 07_multinomial.R -- blocks B6 (ordered rating) and B7 (plan choice) -> res_mult.rds
#   B6   ordinal::clm ordered probit with a scale part for managed accounts (correct specification):
#        APE of engagement on P(rating = 5) with a standard error (delta method plus the sampling
#        variance of the X_i); the homoskedastic ordered probit (clm, = MASS::polr) for comparison;
#        the Bond and Lang (2019) latent-scale comparison from the estimated location and scale; a Monte
#        Carlo of the APE (R = 1,000 samples, seeds 1260000 + r; ruling 29, the one-sample interval misses)
#   B7   MNL by nnet::multinom on the IIA sample (coefficients beta_Pro, beta_Ent, base Free; APEs of
#        engagement and seats on each plan's probability); the Hausman and McFadden test that drops Free
#        (hm_test; = mlogit::hmftest, checked below on both samples) on the IIA sample and on the nested
#        sample; Monte Carlo of the test's rejection rate under IIA (R = 4,000 samples of n = 5,000,
#        seeds 1270000 + r), with the share of negative statistics and of draws where V_r - V_f is not
#        positive semidefinite. The over-rejection under IIA (7.0% at 5%, about 3.9% at 1%) comes from that
#        near-singular, often non-PSD V_r - V_f, not from a coding error: the stacked-score (suest-type)
#        statistic hm_test_stacked, PSD by construction, is computed on the same Monte Carlo draws as a contrast
#        (rejection rates at 1, 5 and 10%, share of negative statistics = 0). The main test stays hm_test
#        (= mlogit::hmftest).
# Reads helix_B6.rds, helix_B7.rds, helix_B7n.rds, truth12.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("bin_fns.R")
suppressPackageStartupMessages({ library(ordinal); library(nnet); library(mlogit); library(marginaleffects) })
tr <- readRDS("truth12.rds"); seed <- tr$params$seed; P6 <- tr$params$B6; P7 <- tr$params$B7
cores <- max(1L, min(8L, parallel::detectCores() - 2L))
covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])

## =================================================================
## B6: ordered probit
## =================================================================
d6 <- readRDS("helix_B6.rds")
f6 <- clm(rating ~ x1 + bill_kusd + G, scale = ~ G, data = d6, link = "probit")
th <- coef(f6); V6 <- vcov(f6)
stopifnot(identical(names(th), c("1|2", "2|3", "3|4", "4|5", "x1", "bill_kusd", "G", "G")))
Xo <- cbind(d6$x1, d6$bill_kusd, d6$G)
top_i <- function(t) oprobit_top_i(t[5], drop(Xo %*% t[5:7]), t[4], exp(t[8] * d6$G))
top_inf <- avg_effect_inference(top_i, unname(th), V6)
ms6 <- avg_slopes(f6, variables = "x1")
stopifnot(abs(ms6$estimate[ms6$group == "5"] - top_inf$est) < 1e-5)
f6h <- clm(rating ~ x1 + bill_kusd + G, data = d6, link = "probit")
thh <- coef(f6h); Vh <- vcov(f6h)
toph_inf <- avg_effect_inference(function(t) oprobit_top_i(t[5], drop(Xo %*% t[5:7]), t[4], 1), unname(thh), Vh)
p6 <- MASS::polr(rating ~ x1 + bill_kusd + G, data = d6, method = "probit", Hess = TRUE)
stopifnot(max(abs(coef(p6) - thh[5:7])) < 1e-4)
lr_scale <- 2 * (as.numeric(logLik(f6)) - as.numeric(logLik(f6h)))
sig1 <- exp(unname(th[8])); dG <- unname(th[7])
B6r <- list(coef = setNames(unname(th), c("mu1", "mu2", "mu3", "mu4", "x1", "bill_kusd", "G", "log_sigma_G")),
            se = setNames(sqrt(diag(V6)), c("mu1", "mu2", "mu3", "mu4", "x1", "bill_kusd", "G", "log_sigma_G")),
            sigma_managed = sig1,
            ape_top = top_inf$est, ape_top_se = top_inf$se, ape_top_ci = top_inf$ci, ape_top_se_delta = top_inf$se_delta,
            ape_by_rating = setNames(ms6$estimate, paste0("r", ms6$group)),
            homo = list(coef = setNames(unname(thh), c("mu1", "mu2", "mu3", "mu4", "x1", "bill_kusd", "G")),
                        se = setNames(sqrt(diag(Vh)), c("mu1", "mu2", "mu3", "mu4", "x1", "bill_kusd", "G")),
                        ape_top = toph_inf$est, ape_top_ci = toph_inf$ci),
            lr_scale = lr_scale, p_lr_scale = pchisq(lr_scale, 1, lower.tail = FALSE),
            # Bond and Lang (2019): for the same x, latent mean difference (managed minus self-serve) on the Y*
            # scale, and the log ratio E[h | x1, bill, G = 1] / E[h | x1, bill, G = 0] under h(Y*) = -exp(-Y*), = -delta + (sigma_1^2 - 1) / 2
            # (CONDITIONAL on x1 and bill: the unconditional E[h | G] is minus infinity, see 01_dgp.R)
            bl_latent_diff = dG, bl_log_ratio_h = -dG + (sig1^2 - 1) / 2,
            shares = as.vector(table(d6$rating)) / nrow(d6),
            shares_G = rbind(self_serve = as.vector(table(d6$rating[d6$G == 0])) / sum(d6$G == 0),
                             managed = as.vector(table(d6$rating[d6$G == 1])) / sum(d6$G == 1)),
            mean_rating_G = c(self_serve = mean(as.integer(d6$rating[d6$G == 0])), managed = mean(as.integer(d6$rating[d6$G == 1]))),
            n = nrow(d6), share_managed = mean(d6$G), mean_bill_kusd = mean(d6$bill_kusd), median_bill_kusd = median(d6$bill_kusd))
B6r$covers_truth <- covers_num(B6r$ape_top_ci, tr$B6$ape_top)
# Monte Carlo of the B6 design (ruling 29: the one-sample interval above misses its truth): R = 1,000 samples,
# seeds 1260000 + r; mean of the APE estimate and coverage of its nominal 95% interval
fit_B6 <- function(d) {
  f <- clm(rating ~ x1 + bill_kusd + G, scale = ~ G, data = d, link = "probit")
  X <- cbind(d$x1, d$bill_kusd, d$G)
  a <- avg_effect_inference(function(t) oprobit_top_i(t[5], drop(X %*% t[5:7]), t[4], exp(t[8] * d$G)), unname(coef(f)), vcov(f))
  c(est = a$est, se = a$se, lo = a$ci[1], hi = a$ci[2], b1 = unname(coef(f)[5]))
}
t0 <- proc.time()[["elapsed"]]
M6 <- do.call(rbind, parallel::mclapply(seq_len(P6$R_mc), function(r) { set.seed(seed$B6_mc + r); fit_B6(gen_B6(P6)) },
                                        mc.cores = cores))
B6r$mc <- list(R = nrow(M6), mean_est = mean(M6[, "est"]), sd_est = sd(M6[, "est"]), mean_se = mean(M6[, "se"]),
               cover = mean(M6[, "lo"] <= tr$B6$ape_top & tr$B6$ape_top <= M6[, "hi"]),
               bias = mean(M6[, "est"]) - tr$B6$ape_top, mean_b1 = mean(M6[, "b1"]),
               share_est_above_sample = mean(M6[, "est"] >= B6r$ape_top))
cat(sprintf("B6 Monte Carlo (R = %d, %.1f s): mean APE %.5f (truth %.5f), sd %.5f, mean SE %.5f, coverage %.3f\n",
            B6r$mc$R, proc.time()[["elapsed"]] - t0, B6r$mc$mean_est, tr$B6$ape_top, B6r$mc$sd_est, B6r$mc$mean_se, B6r$mc$cover))
cat(sprintf("B6: APE of engagement on P(5) %.5f [%.5f, %.5f] (truth %.5f); homoskedastic %.5f; delta_G %.3f, sigma_managed %.3f, h log ratio %.3f\n",
            B6r$ape_top, B6r$ape_top_ci[1], B6r$ape_top_ci[2], tr$B6$ape_top, toph_inf$est, dG, sig1, B6r$bl_log_ratio_h))

## =================================================================
## B7: plan choice
## =================================================================
d7 <- readRDS("helix_B7.rds"); d7n <- readRDS("helix_B7n.rds")
mnl <- multinom(choice ~ x + s, data = d7, trace = FALSE, reltol = 1e-12, maxit = 1000, Hess = TRUE)
B <- coef(mnl); V7 <- vcov(mnl)
seB <- matrix(sqrt(diag(V7)), 2, 3, byrow = TRUE, dimnames = dimnames(B))
W7 <- cbind(1, d7$x, d7$s)
ape7 <- function(j) {
  f <- function(t) { Bt <- matrix(t, 2, 3, byrow = TRUE); mnl_ape(Bt, W7, j) }
  est <- f(as.vector(t(B)))
  J <- numDeriv::jacobian(f, as.vector(t(B)))
  P <- mnl_probs(B, W7); bj <- c(0, B[, j])
  gi <- P * (matrix(bj, nrow(P), 3, byrow = TRUE) - drop(P %*% bj))
  se <- sqrt(diag(J %*% V7 %*% t(J)) + apply(gi, 2, var) / nrow(P))
  list(est = setNames(est, P7$levels), se = setNames(se, P7$levels))
}
ax <- ape7(2); as_ <- ape7(3)
ms7 <- avg_slopes(mnl, variables = "x")
stopifnot(max(abs(ms7$estimate[match(P7$levels, ms7$group)] - ax$est)) < 1e-5)
hm_iia <- hm_test(choice ~ x + s, d7, P7$drop)
hm_nest <- hm_test(choice ~ x + s, d7n, P7$drop)
hm_iia_ent <- hm_test(choice ~ x + s, d7, "Enterprise")
hm_nest_ent <- hm_test(choice ~ x + s, d7n, "Enterprise")
# stacked-score (suest-type, PSD) contrast on the same two samples (drop Free)
hs_iia <- hm_test_stacked(hm_iia); hs_nest <- hm_test_stacked(hm_nest)
# mlogit::hmftest on the same samples (long format), as a check of hm_test
to_long <- function(d) {
  long <- do.call(rbind, lapply(P7$levels, function(a) data.frame(id = seq_len(nrow(d)), alt = a, chosen = d$choice == a, x = d$x, s = d$s)))
  dfidx::dfidx(long[order(long$id, match(long$alt, P7$levels)), ], idx = c("id", "alt"), choice = "chosen")
}
l7 <- to_long(d7); l7n <- to_long(d7n)
ml_iia <- mlogit(chosen ~ 0 | x + s, data = l7, reflevel = "Pro")
ml_nest <- mlogit(chosen ~ 0 | x + s, data = l7n, reflevel = "Pro")
hmf_iia <- hmftest(ml_iia, mlogit(chosen ~ 0 | x + s, data = l7, reflevel = "Pro", alt.subset = c("Pro", "Enterprise")))
hmf_nest <- hmftest(ml_nest, mlogit(chosen ~ 0 | x + s, data = l7n, reflevel = "Pro", alt.subset = c("Pro", "Enterprise")))
stopifnot(abs(hm_iia$stat - unname(hmf_iia$statistic)) < 1e-3 * max(1, abs(hm_iia$stat)),
          abs(hm_nest$stat - unname(hmf_nest$statistic)) < 1e-3 * max(1, abs(hm_nest$stat)))
mnl_n <- multinom(choice ~ x + s, data = d7n, trace = FALSE, reltol = 1e-12, maxit = 1000, Hess = TRUE)
paid <- droplevels(d7n[d7n$choice != "Free", ])
paid$choice <- relevel(paid$choice, ref = "Pro")
gpaid <- glm(I(choice == "Enterprise") ~ x + s, data = paid, family = binomial(), control = glm.control(epsilon = 1e-12))
log10p <- function(h) pchisq(h$stat, h$df, lower.tail = FALSE, log.p = TRUE) / log(10)

t0 <- proc.time()[["elapsed"]]
mc_one <- function(r) {
  set.seed(seed$B7_mc + r)
  d <- gen_B7(P7)
  h <- hm_test(choice ~ x + s, d, P7$drop)
  hs <- hm_test_stacked(h)                         # same draw, same two fits: the stacked-score contrast
  c(stat = h$stat, p = h$p, min_eigen = h$min_eigen, max_eigen = h$max_eigen,
    s_stat = hs$stat, s_p = hs$p, s_min_eigen = hs$min_eigen)
}
M7 <- do.call(rbind, parallel::mclapply(seq_len(P7$R_mc), mc_one, mc.cores = cores))
cat(sprintf("B7 Monte Carlo (R = %d) done in %.1f s\n", nrow(M7), proc.time()[["elapsed"]] - t0))
iia <- list(coef = B, se = seB, ci_lb = B - qnorm(0.975) * seB, ci_ub = B + qnorm(0.975) * seB,
            covers = (B - qnorm(0.975) * seB <= tr$B7$coef) & (tr$B7$coef <= B + qnorm(0.975) * seB),
            ape_x = ax$est, ape_x_se = ax$se, ape_s = as_$est, ape_s_se = as_$se,
            hm_stat = hm_iia$stat, hm_df = hm_iia$df, hm_p = hm_iia$p, hm_stat_mlogit = unname(hmf_iia$statistic),
            hm_ent_stat = hm_iia_ent$stat, hm_ent_p = hm_iia_ent$p,
            shares = setNames(as.vector(table(d7$choice)) / nrow(d7), P7$levels),
            mc_R = nrow(M7), mc_reject = mean(M7[, "p"] < 0.05), mc_reject01 = mean(M7[, "p"] < 0.01),
            mc_reject10 = mean(M7[, "p"] < 0.10), mc_neg_share = mean(M7[, "stat"] < 0),
            mc_mean_stat = mean(M7[, "stat"]), mc_median_stat = median(M7[, "stat"]),
            mc_mean_stat_pos = mean(M7[M7[, "stat"] >= 0, "stat"]), mc_stats = M7[, "stat"],
            mc_mcse = sqrt(0.05 * 0.95 / nrow(M7)),
            # why the classic test over-rejects: V_r - V_f is near-singular (smallest over largest eigenvalue, median
            # over draws) and not positive semidefinite in a share of the draws
            mc_nonpsd_share = mean(M7[, "min_eigen"] < 0), mc_median_eigen_ratio = median(M7[, "min_eigen"] / M7[, "max_eigen"]),
            # contrast: stacked-score (suest-type) statistic on the same draws; PSD by construction, never negative
            suest = list(hm_stat = hs_iia$stat, hm_df = hs_iia$df, hm_p = hs_iia$p,
                         mc_R = nrow(M7), mc_reject01 = mean(M7[, "s_p"] < 0.01), mc_reject = mean(M7[, "s_p"] < 0.05),
                         mc_reject10 = mean(M7[, "s_p"] < 0.10), mc_neg_share = mean(M7[, "s_stat"] < 0),
                         mc_mean_stat = mean(M7[, "s_stat"]), mc_median_stat = median(M7[, "s_stat"]),
                         mc_min_eigen_min = min(M7[, "s_min_eigen"]), mc_stats = M7[, "s_stat"]))
nested <- list(coef = coef(mnl_n), hm_stat = hm_nest$stat, hm_df = hm_nest$df, hm_p = hm_nest$p, hm_log10p = log10p(hm_nest),
               hm_stat_mlogit = unname(hmf_nest$statistic), hm_ent_stat = hm_nest_ent$stat, hm_ent_p = hm_nest_ent$p,
               hm_ent_log10p = log10p(hm_nest_ent),
               suest = list(hm_stat = hs_nest$stat, hm_df = hs_nest$df, hm_p = hs_nest$p,
                            hm_log10p = pchisq(hs_nest$stat, hs_nest$df, lower.tail = FALSE, log.p = TRUE) / log(10)),
               full_ent_minus_pro = coef(mnl_n)["Enterprise", ] - coef(mnl_n)["Pro", ],
               paid_logit = coef(gpaid), shares = setNames(as.vector(table(d7n$choice)) / nrow(d7n), P7$levels))
iia$full_ent_minus_pro <- B["Enterprise", ] - B["Pro", ]
cat(sprintf("B7 IIA sample: HM (drop Free) %.3f, p = %.3f; MC rejection at 5%% %.4f (1%% %.4f, 10%% %.4f), negative %.3f, V_r - V_f not PSD in %.3f of the draws\n",
            iia$hm_stat, iia$hm_p, iia$mc_reject, iia$mc_reject01, iia$mc_reject10, iia$mc_neg_share, iia$mc_nonpsd_share))
cat(sprintf("B7 IIA, stacked-score contrast: statistic %.3f, p = %.3f; MC rejection at 1%% %.4f, 5%% %.4f, 10%% %.4f; negative %.3f; mean statistic %.3f (chi2_%d: %d)\n",
            iia$suest$hm_stat, iia$suest$hm_p, iia$suest$mc_reject01, iia$suest$mc_reject, iia$suest$mc_reject10,
            iia$suest$mc_neg_share, iia$suest$mc_mean_stat, iia$suest$hm_df, iia$suest$hm_df))
cat(sprintf("B7 nested sample: HM (drop Free) %.2f, p = %.3g; (drop Enterprise) %.2f, p = %.3g; stacked-score %.2f, log10 p = %.2f\n",
            nested$hm_stat, nested$hm_p, nested$hm_ent_stat, nested$hm_ent_p, nested$suest$hm_stat, nested$suest$hm_log10p))
cat(sprintf("    engagement APE on P(Free, Pro, Ent): %s (truth %s)\n", paste(sprintf("%.4f", ax$est), collapse = ", "),
            paste(sprintf("%.4f", tr$B7$ape_x), collapse = ", ")))

saveRDS(list(B6 = B6r, B7 = list(iia = iia, nested = nested)), "res_mult.rds")
cat("Saved res_mult.rds\n")
