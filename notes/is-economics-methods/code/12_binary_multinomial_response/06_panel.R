# ------------------------------------------------------------------
# 06_panel.R -- block B5 (monthly adoption panel) -> res_panel.rds
# One sample (N = 1,000 firms, T = 8 months) and, from the same panel, its first four months:
#   pooled   pooled logit of Y_it on X_it, cluster-robust by firm (a partial MLE that ignores alpha_i)
#   clogit   FE logit by conditional MLE: survival::clogit (exact) and the hand-written CMLE with the
#            firm-clustered sandwich (clogit refuses robust variances under its exact method)
#   fe, bc   FE logit with firm dummies (alpaca::feglm, the incidental-parameters MLE) and its analytical
#            bias correction (alpaca::biasCorr: the correction of Fernandez-Val and Weidner 2016, with
#            bandwidth L = 0 for strictly exogenous X; with firm effects only, the one-way case, it coincides
#            with the analytical correction of Hahn and Newey 2004); APEs by alpaca::getAPEs
#            (marginaleffects does not support feglm objects)
#   cre      pooled CRE probit (cre_probit: Mundlak device) and the CRE random-effects probit by MLE
#            (lme4::glmer, nAGQ = 8) on (X_it, Xbar_i)
# TWO DESIGNS AT T = 4 (controller ruling 42), never mixed:
#   one sample, "T = 4" (keys panel_T4_*)   months 1 to 4 of the T = 8 panel in helix_B5.rds: alpha_i was built
#            on the 8-month mean of X_it. Its truths: dgp_B5_T8first4_*; the law of (X_it, alpha_i) is the same in
#            every month, so they EQUAL the T = 8 truths exactly (2-D Gauss-Hermite quadrature, controller ruling 44:
#            pooled-logit limit 0.97690, true APE 0.124527, pooled-APE limit 0.18877); every one-sample T = 4
#            estimate is compared with them. The large draw of 01_dgp.R only checks them.
#   Monte Carlo, T = 4 (keys panel_mc_T4_*)   a 4-month DGP, gen_B5(T = 4): alpha_i built on the 4-month mean.
#            Its truths (also exact quadrature values): dgp_B5_T4_*, dgp_B5_ape_T4 (pooled limit 1.04498, APE 0.121431).
# Monte Carlo (R = 1,000 panels per T in {4, 8}, N = 1,000, seeds 1254000 + r and 1258000 + r): the bias and
# sd of the FE MLE, the bias-corrected FE and the CMLE, interval coverage, and the APE bias, sd, standard error
# and coverage.
# APE convention: the chapter's target is the APE over ALL firm-months, E[beta Lambda'(alpha_i + beta X_it)].
# alpaca::getAPEs()$delta already averages over all firm-months (it divides by nobs.full, which includes the
# firms that never or always adopt: alpha-hat = -/+ infinity, dropped by feglm, a zero estimated effect). So
# ape_fe = delta and ape_bc = delta of the bias-corrected fit ARE the all-firm APEs (tests/test_fns.R checks
# delta = sum(beta-hat F'(eta-hat)) / nobs.full); the APEs among the firms that move are delta / share_kept
# (ape_fe_kept for the FE MLE only; "kept" = the firms feglm keeps), compared with the average true effect over the
# movers' firm-months of the same panel (ape_true_movers). There is no movers-only bias-corrected APE:
# getAPEs(bc)$delta is the mean over all rows minus alpaca's APE bias term, so dividing by share_kept also rescales
# that term (controller ruling 44). The all-firm FE APE is nearly unbiased (-0.7% at T = 4, -0.2% at T = 8) because
# two large errors offset, which the Monte Carlo stores as ape_fe_kept_part (the kept firms' overstated effects,
# positive: sum over kept rows of (estimated - true effect) / all rows) and ape_fe_drop_part (the dropped firms'
# zero estimated effect, negative: minus the true effects of their rows / all rows); with ape_fe_sampling_part (the
# mean of the per-panel true APE minus the population truth) they add up to ape_fe_bias exactly.
# APE standard errors (fe_logit_ape in bin_fns.R): getAPEs' variance is conditional on the covariates (n.pop =
# NULL, no finite-population term); the standard errors stored here are for the POPULATION APE, adding the
# firm-clustered sampling variance of the individual effects (se^2 = se_cond^2 + se_sample^2), the pattern of
# avg_effect_inference in the cross-section blocks. ape_*_se_cond is the conditional part alone.
# Reads helix_B5.rds, truth12.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("bin_fns.R")
suppressPackageStartupMessages({ library(survival); library(alpaca); library(sandwich) })
tr <- readRDS("truth12.rds"); P5 <- tr$params$B5; seed <- tr$params$seed
cores <- max(1L, min(8L, parallel::detectCores() - 2L))
covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])
beta <- P5$beta

## ---- estimators on one panel ---------------------------------------------------------------------
fit_panel <- function(d, full = TRUE) {
  d$fid <- factor(d$id)
  fe <- alpaca::feglm(Y ~ X | fid, data = d, family = binomial("logit"))
  bc <- alpaca::biasCorr(fe)
  mov <- !(ave(d$Y, d$id) %in% c(0, 1))                       # firm-months of firms that move
  me_true <- beta * dlogis(d$alpha_true + beta * d$X)          # true marginal effects (simulation truth)
  ap <- fe_logit_ape(fe, bc, d$id, d$Y)                       # all-firm APEs = getAPEs' delta; movers-only = delta / share_kept
  cl <- survival::clogit(Y ~ X + strata(id), data = d)
  pl <- glm(Y ~ X, data = d, family = binomial("logit"), control = glm.control(epsilon = 1e-12, maxit = 100))
  out <- list(fe = list(est = unname(coef(fe)), se = sqrt(vcov(fe)[1, 1])),
              bc = list(est = unname(coef(bc)), se = sqrt(vcov(bc)[1, 1])),
              clogit = list(est = unname(coef(cl)), se = sqrt(vcov(cl)[1, 1])),
              pooled = list(est = unname(coef(pl)[2])),
              share_kept = ap$share_kept, n_firms_kept = length(unique(d$id[mov])),
              # all-firm APEs (the chapter's target) and their population-APE standard errors
              ape_fe = ap$ape_fe, ape_bc = ap$ape_bc, ape_fe_se = ap$se_fe, ape_bc_se = ap$se_bc,
              ape_fe_se_cond = ap$se_cond_fe, ape_bc_se_cond = ap$se_cond_bc, ape_se_sample = ap$se_sample,
              # APEs among the firms that move (delta / share_kept), against ape_true_movers
              ape_fe_kept = ap$ape_fe_movers,
              # error of the all-firm FE APE against this panel's true APE = kept part + dropped part (exact)
              ape_fe_kept_part = sum(ap$g - me_true * mov) / nrow(d), ape_fe_drop_part = -sum(me_true[!mov]) / nrow(d),
              ape_true_sample = mean(me_true), ape_true_movers = mean(me_true[mov]))
  for (k in c("fe", "bc", "clogit")) out[[k]]$ci <- ci95(out[[k]]$est, out[[k]]$se)
  out$fe$ape_ci <- ap$ci_fe; out$bc$ape_ci <- ap$ci_bc
  if (full) {
    Vpl <- sandwich::vcovCL(pl, cluster = d$id, type = "HC0", cadjust = FALSE)
    Xp <- cbind(1, d$X)
    pa <- avg_effect_inference(function(t) dlogis(drop(Xp %*% t)) * t[2], unname(coef(pl)), Vpl, cluster = d$id)
    out$pooled <- list(est = unname(coef(pl)[2]), se = sqrt(Vpl[2, 2]), ci = ci95(coef(pl)[2], sqrt(Vpl[2, 2])),
                       se_model = sqrt(vcov(pl)[2, 2]), ape = pa$est, ape_se = pa$se, ape_ci = pa$ci)
    cm <- fe_logit_cmle(d$Y, cbind(X = d$X), d$id)
    out$cmle <- list(est = unname(cm$coef), se_model = sqrt(cm$vcov_model[1, 1]), se_robust = sqrt(cm$vcov_robust[1, 1]),
                     ci_robust = ci95(cm$coef, sqrt(cm$vcov_robust[1, 1])), n_units_used = cm$n_units_used,
                     max_diff_clogit = abs(unname(cm$coef) - unname(coef(cl))))
  }
  out
}

## ---- one-sample estimates against the truth of THEIR design (T8: dgp_B5_*; months 1 to 4 of T = 8: dgp_B5_T8first4_*) ----
vs_truth <- function(s, truth, design) {
  s$truth <- list(design = design, ape = truth$ape, pooled_plim = truth$pooled_plim, pooled_ape_plim = truth$pooled_ape_plim)
  s$vs <- list(pooled_minus_plim = s$pooled$est - truth$pooled_plim,
               pooled_covers_plim = covers_num(s$pooled$ci, truth$pooled_plim),
               pooled_ape_minus_truth = s$pooled$ape - truth$ape, pooled_ape_covers_truth = covers_num(s$pooled$ape_ci, truth$ape),
               ape_fe_minus_truth = s$ape_fe - truth$ape, ape_bc_minus_truth = s$ape_bc - truth$ape,
               ape_fe_covers_truth = covers_num(s$fe$ape_ci, truth$ape), ape_bc_covers_truth = covers_num(s$bc$ape_ci, truth$ape))
  s
}

## ---- one sample: T = 8 and T = 4 (months 1 to 4) ------------------------------------------------------
d5 <- readRDS("helix_B5.rds")
s8 <- vs_truth(fit_panel(d5), tr$B5$T8, "T = 8 panel")
d54 <- d5[d5$t <= P5$T_short, ]                                 # months 1 to 4 of the T = 8 panel (not a 4-month DGP)
s4 <- vs_truth(fit_panel(d54), tr$B5$T8first4, "months 1 to 4 of the T = 8 panel")
cre <- cre_probit(d5$Y, cbind(X = d5$X), d5$id)
d5$Xbar <- ave(d5$X, d5$id)
gm <- lme4::glmer(Y ~ X + Xbar + (1 | id), data = d5, family = binomial("probit"), nAGQ = 8)
sig_a <- as.numeric(attr(lme4::VarCorr(gm)$id, "stddev"))
bgm <- lme4::fixef(gm)
cre_re <- list(coef = bgm, se = sqrt(diag(as.matrix(vcov(gm)))), sigma_a = sig_a,
               scaled_coef = bgm / sqrt(1 + sig_a^2),
               # APE of the RE probit: E over a_i ~ N(0, sigma_a^2) of beta phi(.) = scaled beta phi(scaled index)
               ape = mean(dnorm(drop(cbind(1, d5$X, d5$Xbar) %*% bgm) / sqrt(1 + sig_a^2))) * unname(bgm[2]) / sqrt(1 + sig_a^2))
cat(sprintf("B5 (T = 8): pooled %.4f [%.4f, %.4f]; clogit %.4f [%.4f, %.4f]; FE MLE %.4f; BC %.4f; truth %.2f\n",
            s8$pooled$est, s8$pooled$ci[1], s8$pooled$ci[2], s8$clogit$est, s8$clogit$ci[1], s8$clogit$ci[2],
            s8$fe$est, s8$bc$est, beta))
cat(sprintf("B5 (months 1 to 4 of the T = 8 panel): clogit %.4f; FE MLE %.4f; BC %.4f; pooled %.4f (its limit %.4f)\n",
            s4$clogit$est, s4$fe$est, s4$bc$est, s4$pooled$est, tr$B5$T8first4$pooled_plim))
cat(sprintf("B5 APE over all firm-months (T = 8; truth %.4f): pooled %.4f, CRE probit %.4f, CRE RE probit %.4f, FE %.4f, BC FE %.4f [%.4f, %.4f]\n",
            tr$B5$ape, s8$pooled$ape, cre$ape, cre_re$ape, s8$ape_fe, s8$ape_bc, s8$bc$ape_ci[1], s8$bc$ape_ci[2]))
cat(sprintf("B5 APE, months 1 to 4 (truth %.4f): pooled %.4f, FE %.4f, BC FE %.4f [%.4f, %.4f]; movers only (%.1f%% of firm-months kept): FE %.4f (average true effect of the movers %.4f)\n",
            tr$B5$T8first4$ape, s4$pooled$ape, s4$ape_fe, s4$ape_bc, s4$bc$ape_ci[1], s4$bc$ape_ci[2], 100 * s4$share_kept,
            s4$ape_fe_kept, s4$ape_true_movers))

## ---- Monte Carlo ---------------------------------------------------------------------------------------
mc_one <- function(r, T) {
  set.seed((if (T == P5$T_short) seed$B5_mc_T4 else seed$B5_mc_T8) + r)
  d <- gen_B5(P5, N = P5$N, T = T)
  f <- fit_panel(d, full = FALSE)
  c(fe = f$fe$est, fe_se = f$fe$se, bc = f$bc$est, bc_se = f$bc$se, cl = f$clogit$est, cl_se = f$clogit$se,
    pooled = f$pooled$est, ape_fe = f$ape_fe, ape_bc = f$ape_bc, ape_fe_se = f$ape_fe_se, ape_bc_se = f$ape_bc_se,
    share_kept = f$share_kept, ape_fe_kept = f$ape_fe_kept,
    ape_fe_kept_part = f$ape_fe_kept_part, ape_fe_drop_part = f$ape_fe_drop_part,
    ape_true_sample = f$ape_true_sample, ape_true_movers = f$ape_true_movers)
}
t0 <- proc.time()[["elapsed"]]
mc <- list()
for (T in c(P5$T_short, P5$T)) {
  M <- do.call(rbind, parallel::mclapply(seq_len(P5$R_mc), mc_one, T = T, mc.cores = cores))
  ape_true <- if (T == P5$T) tr$B5$ape else tr$B5$ape_T4
  cov_to <- function(est, se, v) mean(abs(est - v) <= qnorm(0.975) * se)
  cov <- function(est, se) cov_to(est, se, beta)
  mc[[paste0("T", T)]] <- list(T = T, R = nrow(M), draws = M,
    fe_mean = mean(M[, "fe"]), fe_sd = sd(M[, "fe"]), fe_bias = mean(M[, "fe"]) - beta,
    bc_mean = mean(M[, "bc"]), bc_sd = sd(M[, "bc"]), bc_bias = mean(M[, "bc"]) - beta,
    cl_mean = mean(M[, "cl"]), cl_sd = sd(M[, "cl"]), cl_bias = mean(M[, "cl"]) - beta,
    pooled_mean = mean(M[, "pooled"]), pooled_bias = mean(M[, "pooled"]) - beta,
    cover_fe = cov(M[, "fe"], M[, "fe_se"]), cover_bc = cov(M[, "bc"], M[, "bc_se"]), cover_cl = cov(M[, "cl"], M[, "cl_se"]),
    mean_se_cl = mean(M[, "cl_se"]), mean_se_bc = mean(M[, "bc_se"]),
    # APE over all firm-months (getAPEs' delta) against the DGP truth: bias, sd, mean population-APE standard error, coverage
    ape_fe_bias = mean(M[, "ape_fe"]) - ape_true, ape_bc_bias = mean(M[, "ape_bc"]) - ape_true,
    ape_fe_rel_bias = (mean(M[, "ape_fe"]) - ape_true) / ape_true, ape_bc_rel_bias = (mean(M[, "ape_bc"]) - ape_true) / ape_true,
    ape_fe_sd = sd(M[, "ape_fe"]), ape_bc_sd = sd(M[, "ape_bc"]),
    ape_fe_mean_se = mean(M[, "ape_fe_se"]), ape_bc_mean_se = mean(M[, "ape_bc_se"]),
    ape_fe_cover = cov_to(M[, "ape_fe"], M[, "ape_fe_se"], ape_true), ape_bc_cover = cov_to(M[, "ape_bc"], M[, "ape_bc_se"], ape_true),
    ape_true_sample_mean = mean(M[, "ape_true_sample"]),
    fe_rel_bias = (mean(M[, "fe"]) - beta) / beta, bc_rel_bias = (mean(M[, "bc"]) - beta) / beta,
    # APE among the firms that move (the sample the FE estimator uses; delta / share_kept), against its own true value
    ape_movers_true = mean(M[, "ape_true_movers"]),
    ape_fe_kept_bias = mean(M[, "ape_fe_kept"] - M[, "ape_true_movers"]),
    ape_fe_kept_rel_bias = mean(M[, "ape_fe_kept"] - M[, "ape_true_movers"]) / mean(M[, "ape_true_movers"]),
    # why the all-firm FE APE is nearly unbiased: ape_fe_bias = kept part + dropped part + sampling part (exact)
    ape_fe_kept_part = mean(M[, "ape_fe_kept_part"]), ape_fe_drop_part = mean(M[, "ape_fe_drop_part"]),
    ape_fe_sampling_part = mean(M[, "ape_true_sample"]) - ape_true,
    share_kept = mean(M[, "share_kept"]), ape_true = ape_true)
}
cat(sprintf("B5 Monte Carlo done in %.1f s\n", proc.time()[["elapsed"]] - t0))
for (k in names(mc)) with(mc[[k]], cat(sprintf(
  "  %s: FE bias %+.4f (sd %.4f, cover %.3f); BC bias %+.4f (sd %.4f, cover %.3f); CMLE bias %+.4f (sd %.4f, cover %.3f)\n      APE over all firm-months (truth %.4f): FE %+.4f (%+.1f%%, sd %.4f, mean SE %.4f, cover %.3f), BC %+.4f (%+.1f%%, sd %.4f, mean SE %.4f, cover %.3f)\n      APE among the movers (truth %.4f): FE %+.4f (%+.1f%%); error of the all-firm FE APE: kept firms %+.4f, dropped firms %+.4f, sampling %+.4f\n",
  k, fe_bias, fe_sd, cover_fe, bc_bias, bc_sd, cover_bc, cl_bias, cl_sd, cover_cl,
  ape_true, ape_fe_bias, 100 * ape_fe_rel_bias, ape_fe_sd, ape_fe_mean_se, ape_fe_cover,
  ape_bc_bias, 100 * ape_bc_rel_bias, ape_bc_sd, ape_bc_mean_se, ape_bc_cover,
  ape_movers_true, ape_fe_kept_bias, 100 * ape_fe_kept_rel_bias, ape_fe_kept_part, ape_fe_drop_part, ape_fe_sampling_part)))

res <- list(T8 = s8, T4 = s4, cre = cre, cre_re = cre_re, mc = mc,
            sample = list(N = P5$N, T = P5$T, rate = mean(d5$Y),
                          share_movers = mean(tapply(d5$Y, d5$id, function(z) sum(z) > 0 && sum(z) < P5$T)),
                          share_movers_T4 = mean(tapply(d54$Y, d54$id, function(z) sum(z) > 0 && sum(z) < P5$T_short))))
res$T8$clogit$covers_truth <- covers_num(s8$clogit$ci, beta); res$T8$pooled$covers_truth <- covers_num(s8$pooled$ci, beta)
res$T4$clogit$covers_truth <- covers_num(s4$clogit$ci, beta)
saveRDS(res, "res_panel.rds")
cat("Saved res_panel.rds\n")
