# ------------------------------------------------------------------
# 05_cross_section.R -- blocks B2, B3 and B4 (cross sections) -> res_cs.rds
#   B2   probit with an interaction: beta_jj' and its tests, the average cross-partial effect
#        (Ai and Norton 2003) with a standard error, the APEs, the share of firms with a positive
#        cross-partial, the "marginal effect of the interaction term" that is not the cross-partial
#   B3   naive probit of Y on (X, D); Rivers-Vuong control function with both steps bootstrapped
#        together (B = 499, seed B3_boot) and the APE of D from the average structural function;
#        the exogeneity t test on v-hat; IV probit by maximum likelihood; linear 2SLS for comparison
#   B4   short probit on X (attenuated coefficient, APE), the long probit with a_i (infeasible),
#        the LPM slope
# APE standard errors: delta method plus the sampling variance of the X_i (avg_effect_inference in
# bin_fns.R), so that the intervals are for the population APE.
# Reads helix_B2.rds, helix_B3.rds, helix_B4.rds, truth12.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("bin_fns.R")
suppressPackageStartupMessages({ library(sandwich); library(marginaleffects) })
tr <- readRDS("truth12.rds"); seed <- tr$params$seed
cores <- max(1L, min(8L, parallel::detectCores() - 2L))
covers_num <- function(ci, v) as.numeric(ci[1] <= v && v <= ci[2])

## =================================================================
## B2: interaction
## =================================================================
d2 <- readRDS("helix_B2.rds")
X2 <- cbind("(Intercept)" = 1, x1 = d2$x1, x2 = d2$x2, "x1:x2" = d2$x1 * d2$x2)
m2 <- binary_mle(d2$y, X2, "probit")
g2 <- glm(y ~ x1 * x2, data = d2, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
stopifnot(max(abs(m2$coef - coef(g2))) < 1e-8)
b2 <- m2$coef; V2 <- m2$vcov_model
se_int <- sqrt(V2[4, 4]); se_int_rob <- sqrt(m2$vcov_sand[4, 4])
g2r <- glm(y ~ x1 + x2, data = d2, family = binomial("probit"), control = glm.control(epsilon = 1e-12, maxit = 100))
lr_int <- 2 * (as.numeric(logLik(g2)) - as.numeric(logLik(g2r)))
cross_inf <- avg_effect_inference(function(t) cross_partial_i(t, X2, 2, 3, "probit"), b2, V2)
ape1_inf <- avg_effect_inference(function(t) dnorm(drop(X2 %*% t)) * (t[2] + t[4] * d2$x2), b2, V2)
ape2_inf <- avg_effect_inference(function(t) dnorm(drop(X2 %*% t)) * (t[3] + t[4] * d2$x1), b2, V2)
ms2 <- avg_slopes(g2, variables = c("x1", "x2"))
stopifnot(abs(ms2$estimate[ms2$term == "x1"] - ape1_inf$est) < 1e-6, abs(ms2$estimate[ms2$term == "x2"] - ape2_inf$est) < 1e-6)
ci_i <- cross_partial_i(b2, X2, 2, 3, "probit")
naive_inf <- avg_effect_inference(function(t) dnorm(drop(X2 %*% t)) * t[4], b2, V2)
# the cross-partial along the latent index: averages within deciles of X_i'beta-hat (for the figure)
u2 <- drop(X2 %*% b2)
B2r <- list(coef = setNames(b2, colnames(X2)), se = setNames(sqrt(diag(V2)), colnames(X2)),
            coef_interaction = unname(b2[4]), se_interaction = se_int, se_interaction_robust = se_int_rob,
            z_interaction = unname(b2[4]) / se_int, p_interaction = 2 * pnorm(-abs(unname(b2[4]) / se_int)),
            ci_interaction = ci95(unname(b2[4]), se_int), lr_interaction = lr_int,
            p_lr_interaction = pchisq(lr_int, 1, lower.tail = FALSE),
            avg_cross = cross_inf$est, avg_cross_se = cross_inf$se, avg_cross_se_delta = cross_inf$se_delta,
            avg_cross_ci = cross_inf$ci, avg_cross_z = cross_inf$est / cross_inf$se,
            share_pos = mean(ci_i > 0), cross_min = min(ci_i), cross_max = max(ci_i),
            ape1 = ape1_inf$est, ape1_se = ape1_inf$se, ape1_ci = ape1_inf$ci,
            ape2 = ape2_inf$est, ape2_se = ape2_inf$se, ape2_ci = ape2_inf$ci,
            ape1_se_me = ms2$std.error[ms2$term == "x1"],
            naive_int_me = naive_inf$est, naive_int_me_ci = naive_inf$ci,
            rate = mean(d2$y), n = nrow(d2), index = u2, cross_i = ci_i)
B2r$covers_truth <- covers_num(B2r$avg_cross_ci, tr$B2$avg_cross)
cat(sprintf("B2: beta_jj' = %.4f (SE %.4f), average cross-partial %.5f (SE %.5f, truth %.5f), share positive %.3f\n",
            B2r$coef_interaction, se_int, B2r$avg_cross, B2r$avg_cross_se, tr$B2$avg_cross, B2r$share_pos))

## =================================================================
## B3: endogenous engagement, control function
## =================================================================
d3 <- readRDS("helix_B3.rds")
Xe <- cbind("(Intercept)" = 1, X = d3$X)
Xn <- cbind(Xe, D = d3$D)
mn <- binary_mle(d3$Y, Xn, "probit")
naive_inf <- avg_effect_inference(function(t) dnorm(drop(Xn %*% t)) * t[3], mn$coef, mn$vcov_model)
naive <- list(coef = setNames(mn$coef, colnames(Xn)), se = setNames(sqrt(diag(mn$vcov_model)), colnames(Xn)),
              ape = naive_inf$est, ape_se = naive_inf$se, ape_ci = naive_inf$ci)
fs <- lm(D ~ X + Z, data = d3)
fs_t <- coef(fs)["Z"] / sqrt(vcovHC(fs, type = "HC0")["Z", "Z"])
t0 <- proc.time()[["elapsed"]]
cf <- rivers_vuong(d3$Y, Xe, d3$D, cbind(Z = d3$Z), B = 499, seed = seed$B3_boot, cores = cores)
cat(sprintf("B3: bootstrap (B = %d) done in %.1f s\n", cf$B, proc.time()[["elapsed"]] - t0))
cf$boot_ape <- cf$boot[, "ape"]; cf$boot <- NULL
cf$coef_scaled_truth <- tr$B3$scaled
cf$beta_ci <- lapply(seq_along(cf$beta), function(j) ci95(cf$beta[j], cf$beta_se[j])); names(cf$beta_ci) <- names(cf$beta)
cf$rho_p <- 2 * pnorm(-abs(cf$rho_t))
cf$first_stage_F <- unname(fs_t^2)
cf$first_stage_coef <- unname(coef(fs))
iv <- ivprobit_mle(d3$Y, Xe, d3$D, cbind(Z = d3$Z))
iv_inf <- avg_effect_inference(function(t) t[3] * dnorm(drop(Xe %*% t[1:2]) + t[3] * d3$D), iv$coef, iv$vcov)
ivm <- list(beta = iv$beta, se_beta = iv$se_beta, sigma_v = iv$sigma_v, r = iv$r, ape = iv_inf$est, ape_se = iv_inf$se,
            ape_ci = iv_inf$ci, convergence = iv$convergence, loglik = iv$loglik,
            ci_beta_D = ci95(iv$beta[3], iv$se_beta[3]))
# linear 2SLS of Y on (1, X, D) with Z as the instrument (a linear probability IV), HC0
Dhat <- fitted(fs)
Xh <- cbind(Xe, D = Dhat)
b_iv <- solve(crossprod(Xh, Xn), crossprod(Xh, d3$Y))
e_iv <- drop(d3$Y - Xn %*% b_iv)
A_iv <- solve(crossprod(Xh, Xn))
V_iv <- A_iv %*% crossprod(Xh * e_iv) %*% t(A_iv)
tsls <- list(est = b_iv[3], se = sqrt(V_iv[3, 3]), ci = ci95(b_iv[3], sqrt(V_iv[3, 3])))
lpm_ols <- lm(Y ~ X + D, data = d3)
ols_lpm <- list(est = unname(coef(lpm_ols)["D"]), se = sqrt(vcovHC(lpm_ols, type = "HC0")["D", "D"]))
B3r <- list(naive = naive, cf = cf, ivmle = ivm, tsls = tsls, ols_lpm = ols_lpm, rate = mean(d3$Y), n = nrow(d3),
            mean_spend_kusd = mean(d3$spend_kusd), median_spend_kusd = median(d3$spend_kusd))
B3r$naive$covers_truth <- covers_num(naive$ape_ci, tr$B3$ape); B3r$cf$covers_truth <- covers_num(cf$ape_ci, tr$B3$ape)
B3r$cf$covers_truth_pct <- covers_num(cf$ape_ci_pct, tr$B3$ape); B3r$ivmle$covers_truth <- covers_num(ivm$ape_ci, tr$B3$ape)
B3r$tsls$covers_ape <- covers_num(tsls$ci, tr$B3$ape)
cat(sprintf("B3: naive APE %.4f [%.4f, %.4f]; CF APE %.4f (boot SE %.4f) [%.4f, %.4f]; IV-MLE APE %.4f; 2SLS %.4f; truth %.4f\n",
            naive$ape, naive$ape_ci[1], naive$ape_ci[2], cf$ape, cf$ape_se, cf$ape_ci[1], cf$ape_ci[2], ivm$ape,
            tsls$est, tr$B3$ape))
cat(sprintf("    rho t = %.2f, first-stage F = %.1f\n", cf$rho_t, cf$first_stage_F))

## =================================================================
## B4: omitted heterogeneity independent of X
## =================================================================
d4 <- readRDS("helix_B4.rds")
Xs <- cbind("(Intercept)" = 1, X = d4$X)
ms <- binary_mle(d4$Y, Xs, "probit")
short_inf <- avg_effect_inference(function(t) dnorm(drop(Xs %*% t)) * t[2], ms$coef, ms$vcov_model)
Xl <- cbind(Xs, a = d4$a_true)
ml <- binary_mle(d4$Y, Xl, "probit")
long_inf <- avg_effect_inference(function(t) dnorm(drop(Xl %*% t)) * t[2], ml$coef, ml$vcov_model)
lpm4 <- lm(Y ~ X, data = d4)
B4r <- list(beta = unname(ms$coef[2]), beta_se = sqrt(ms$vcov_model[2, 2]), beta_ci = ci95(ms$coef[2], sqrt(ms$vcov_model[2, 2])),
            coef = setNames(ms$coef, colnames(Xs)),
            ape = short_inf$est, ape_se = short_inf$se, ape_ci = short_inf$ci,
            long = list(coef = setNames(ml$coef, colnames(Xl)), se = setNames(sqrt(diag(ml$vcov_model)), colnames(Xl)),
                        ape = long_inf$est, ape_se = long_inf$se, ape_ci = long_inf$ci),
            ratio_short_long = unname(ms$coef[2] / ml$coef[2]),
            lpm_slope = unname(coef(lpm4)[2]), lpm_se = sqrt(vcovHC(lpm4, type = "HC0")[2, 2]),
            rate = mean(d4$Y), n = nrow(d4))
B4r$lpm_ci <- ci95(B4r$lpm_slope, B4r$lpm_se)
B4r$covers_att <- covers_num(B4r$beta_ci, tr$B4$beta_att); B4r$covers_ape <- covers_num(B4r$ape_ci, tr$B4$ape)
B4r$long$covers_ape <- covers_num(B4r$long$ape_ci, tr$B4$ape)
cat(sprintf("B4: short probit beta %.4f (truth beta %.2f, attenuated %.4f); APE %.4f [%.4f, %.4f] (truth %.4f); long beta %.4f, APE %.4f\n",
            B4r$beta, tr$B4$beta, tr$B4$beta_att, B4r$ape, B4r$ape_ci[1], B4r$ape_ci[2], tr$B4$ape, ml$coef[2], long_inf$est))

saveRDS(list(B2 = B2r, B3 = B3r, B4 = B4r), "res_cs.rds")
cat("Saved res_cs.rds\n")
