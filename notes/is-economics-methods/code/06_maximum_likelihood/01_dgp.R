# ------------------------------------------------------------------
# 01_dgp.R -- the Sorrel simulation world for chapter 6 (maximum likelihood)
#
# Sorrel is a (fictional) meal-kit subscription platform: a subscriber picks a plan with a
# monthly fee of 45, 60, 75 or 90 USD (two to five meals a week), the first box ships on the
# first day of the signup month, and the subscription runs until the subscriber cancels.
# Sorrel opened 36 months before the data pull and has grown about 4% a month in signups.
# The outcome is tenure in months; a subscriber still active at the pull is right-censored
# after c = 36 - (signup month) + 1 months, independently of how long he or she would stay.
# Money is in US dollars.
# Four blocks (all new in v3; chapter 6 inherits no code from legacy chapter 5):
#   C1   Weibull tenure, right-censored at the data pull: naive OLS on observed tenure
#        against the censored-data MLE; information matrix equality; W, LR, LM.
#   C2   a second market whose tenure is a mixture of two log-logistic segments (deal-seekers
#        and regulars), fitted by a Weibull: QMLE, pseudo-true value, Hessian against OPG
#        against sandwich, Monte Carlo coverage, Vuong test.
#   C3   C1 regrouped into monthly person-period rows (first 24 months of tenure): pooled
#        complementary log-log = grouped-duration MLE.
#   C4   two-step: a usage-intensity index estimated by a Poisson CMLE on first-week in-app
#        check-ins enters a Weibull hazard as a generated regressor.
# Every number in chapter 6 comes from these data and from truth06.rds.
#
# Outputs (in this directory):
#   sorrel_C1.rds   data.frame: id signup_month fee_usd referred c_months t_months delta t_star_true
#   sorrel_C2.rds   data.frame: id signup_month fee_usd referred c_months t_months delta deal_true t_star_true
#   sorrel_C3.rds   data.frame: id m y fee_usd referred (person-period rows built from C1)
#   sorrel_C4.rds   data.frame: id signup_month fee_usd hh_size quiz_std usage_wk1 c_months t_months
#                   delta u_true t_star_true
#   truth06.rds     list(C1, C2, C3, C4, params): truth and parameters
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 1.5): t_i* (t_star_true), c_i (c_months), t_i = min(t_i*, c_i) (t_months),
# delta_i = 1{t_i* <= c_i} (delta); h(t | x) = exp(x'beta) kappa t^(kappa - 1), S(t | x) = exp(-exp(x'beta) t^kappa);
# grouped intervals m = 1, ..., M with cuts a_m = m months, y_im, eta_m; theta_* pseudo-true value;
# gamma first-step parameter (C4). Field names follow the chapter card: truth06$C1$beta (the fee
# coefficient), $C1$shape (kappa), $C2$pseudo (the pseudo-true shape kappa_*), $C4$beta (the
# coefficient on the usage index).
#
# Seeds (every random block its own set.seed, nothing inherited): C1 6101, C2 6102, C2 large-sample
# pseudo-true value 6112, C4 6104; the Monte Carlo seeds used by 03 and 05 are stored here too
# (C2 6122, C4 6124). C3 has no random draws.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("ml_fns.R")

# Seeds were fixed before any result was seen and are not tuned afterwards (Phase 3 ruling 25). With these
# draws the C2 sample's sandwich interval for kappa_* misses the pseudo-true value (kappa_hat is 2.14
# sandwich SEs below it), an ordinary 5% event. Ruling 29 keeps seed 6102 and replaces the card assertion
# covers(c2$qmle$ci, tr$C2$pseudo) by its Monte Carlo counterpart (03_qmle.R carries the coverage claim).
seed <- list(C1 = 6101, C2 = 6102, C2_large = 6112, C4 = 6104, C2_mc = 6122, C4_mc = 6124)

fees <- c(45, 60, 75, 90)                  # monthly fee, USD (2, 3, 4, 5 meals a week)
p_fee <- c(0.30, 0.35, 0.20, 0.15)
months <- 36
growth <- 1.04
p_c <- growth^(months - seq_len(months) + 1)   # P(c), c = 1, ..., 36 (signup month k = 37 - c)
p_c <- p_c / sum(p_c)

## =================================================================
## C1: Weibull tenure, right-censored at the data pull
## =================================================================
## 5,000 subscribers. x = (1, fee_usd, referred): 30% joined through a friend's referral code.
## h(t | x) = exp(-3.44 + 0.02 fee - 0.40 referred) 0.9 t^(-0.1): each extra dollar of monthly fee
## raises the cancellation hazard by 2% (a 15-dollar plan step by 35%), a referral lowers it by 33%,
## and kappa = 0.9 < 1: the hazard falls slowly with tenure. Median tenure on the 60-dollar plan
## without referral: 8.0 months. About 42% of the spells are still running at the pull (population share).
## truth06$C1$beta = 0.02 (the fee coefficient, PH scale); beta_aft = -beta / kappa is the same effect on
## log tenure (what OLS of log t* on x estimates on complete data).
P1 <- list(n = 5000, months = months, growth = growth, fees = fees, p_fee = p_fee, p_ref = 0.30,
           beta = c(-3.44, 0.02, -0.40), kappa = 0.9)
set.seed(seed$C1)
C1 <- gen_C1(P1)

eta60 <- P1$beta[1] + 60 * P1$beta[2]
cells1 <- expand.grid(c = seq_len(months), fee = fees, ref = 0:1)
cells1$w <- p_c[cells1$c] * p_fee[match(cells1$fee, fees)] * ifelse(cells1$ref == 1, P1$p_ref, 1 - P1$p_ref)
S1 <- exp(-exp(P1$beta[1] + P1$beta[2] * cells1$fee + P1$beta[3] * cells1$ref) * cells1$c^P1$kappa)
truth_C1 <- list(beta = P1$beta[2], beta_vec = P1$beta, shape = P1$kappa, beta_aft = -P1$beta[2] / P1$kappa,
                 median_fee60 = (log(2) / exp(eta60))^(1 / P1$kappa), S12_fee60 = exp(-exp(eta60) * 12^P1$kappa),
                 share_cens_pop = sum(cells1$w * S1))

## =================================================================
## C2: log-logistic mixture fitted by a Weibull (QMLE)
## =================================================================
## 8,000 subscribers in Sorrel's second market, same plans, referral share and censoring design.
## Two latent segments: deal-seekers (joined on a first-box discount; 37.8% of non-referred and 18.2%
## of referred subscribers, logit -0.5 - 1.0 referred) and regulars. Tenure is log-logistic in each
## segment, S(t) = 1 / (1 + (t / a)^p): deal-seekers have median 1.5 months and shape 0.5 (very
## dispersed: many cancel after the first box, some linger for years; decreasing hazard); regulars
## have median 12 months and shape 5 (hazard rising to a peak at 15.8 months, then falling). The fee
## rescales time in both segments, a = median exp(-0.02 (fee - 60)). The analyst fits the Weibull PH
## model with x = (1, fee_usd, referred). No Weibull is the truth: the QMLE converges to the
## pseudo-true value theta_* that minimises the expected Kullback-Leibler divergence for this
## censoring design (White 1982). truth06$C2$pseudo is kappa_* from the Weibull MLE on 2,000,000
## draws (the card's definition); the quadrature of the population objective (c2_cells,
## pseudo_true_weibull) is the check. Monte Carlo (03_qmle.R): R_mc = 2000 replications.
P2 <- list(n = 8000, months = months, growth = growth, fees = fees, p_fee = p_fee, p_ref = 0.30,
           g0 = -0.5, g_ref = -1.0, med_deal = 1.5, p_deal = 0.5, med_reg = 12, p_reg = 5,
           theta_fee = 0.02, fee_ref = 60, N_large = 2e6, R_mc = 2000)
set.seed(seed$C2)
C2 <- gen_C2(P2)

set.seed(seed$C2_large)
L2 <- gen_C2(P2, n = P2$N_large)
X2L <- cbind(1, L2$fee_usd, L2$referred)
fl <- fit_weibull(L2$t_months, L2$delta, X2L, engine = "newton")
stopifnot(fl$converged)
share_cens_large <- 1 - mean(L2$delta)
rm(L2, X2L)
qd <- pseudo_true_weibull(P2)
stopifnot(qd$converged)
cells2 <- expand.grid(c = seq_len(months), fee = fees, ref = 0:1)
cells2$w <- p_c[cells2$c] * p_fee[match(cells2$fee, fees)] * ifelse(cells2$ref == 1, P2$p_ref, 1 - P2$p_ref)
med_ref2 <- uniroot(function(s) c2_surv(s, 60, 0, P2) - 0.5, c(0.01, 100), tol = 1e-12)$root
hg <- c2_hazard(seq(0.25, 36, by = 0.25), 60, 0, P2)
truth_C2 <- list(pseudo = exp(fl$par[4]), pseudo_vec = fl$par, pseudo_fee = fl$par[2],
                 large = list(N = P2$N_large, par = fl$par, se_sand = sqrt(diag(fl$V_sand)), converged = fl$converged,
                              share_cens = share_cens_large),
                 quad = list(par = qd$par, kappa = qd$kappa, value = qd$value, converged = qd$converged),
                 deal_share_nonref = plogis(P2$g0), deal_share_ref = plogis(P2$g0 + P2$g_ref),
                 median_fee60 = med_ref2,
                 share_cens_pop = sum(cells2$w * c2_surv(cells2$c, cells2$fee, cells2$ref, P2)),
                 hazard_fee60 = data.frame(t = seq(0.25, 36, by = 0.25), h = hg))

## =================================================================
## C3: C1 as monthly person-period rows (grouped durations)
## =================================================================
## Intervals m = 1, ..., 24 with a_m = m months: the first two years of tenure. A spell censored at
## c_i months (a whole number: censoring happens at the cut points, as Wooldridge 2010 section 22.4
## assumes) contributes min(c_i, 24) rows; an exit in month m contributes m rows with y = 1 in the
## last; a subscriber still active after 24 months is censored at 24. Because C1's hazard is
## proportional, the grouped model with one free eta_m per month is exactly right:
##   eta_m = beta_0 + log(m^kappa - (m - 1)^kappa)   (the log of the integrated baseline hazard),
## and beta (fee, referral) is C1's.
P3 <- list(cuts = 0:24, M = 24L)
C3 <- person_period(C1$t_months, C1$delta, cbind(fee_usd = C1$fee_usd, referred = C1$referred), P3$cuts)
truth_C3 <- list(eta = P1$beta[1] + log(P3$cuts[-1]^P1$kappa - P3$cuts[-(P3$M + 1)]^P1$kappa), beta = P1$beta[2:3],
                 cuts = P3$cuts)

## =================================================================
## C4: two-step, a generated usage index in the hazard
## =================================================================
## 5,000 subscribers who joined after the Sorrel app launched; same plans and censoring design.
## Onboarding gives z = (1, hh_size, quiz_std): household size (1 to 5; shares 15, 35, 20, 20, 10%)
## and a standardised cooking-quiz score. The usage-intensity index u = -1.0 + 0.25 hh_size + 0.35
## quiz_std is the log of expected in-app recipe check-ins in the first week; the platform sees one
## week of check-ins, v ~ Poisson(exp(u)) (mean 0.82: a sparse signal). The hazard depends on the
## index, not on the noisy count: h(t | x) = exp(-3.70 + 0.02 fee - 0.80 u) 0.9 t^(-0.1), so a
## subscriber with e-times the expected usage has a 55% lower cancellation hazard. Given z, v and
## t* are independent (Wooldridge 2010 eq. 13.66 holds), but the second-step objective does not
## involve v, so the section 13.10.2 efficiency gain does not arise: Gamma_0 != 0 adds variance.
## truth06$C4$beta = -0.80. Monte Carlo (05_twostep.R): R_mc = 1000 replications.
P4 <- list(n = 5000, months = months, growth = growth, fees = fees, p_fee = p_fee,
           hh_sizes = 1:5, p_hh = c(0.15, 0.35, 0.20, 0.20, 0.10), gamma = c(-1.0, 0.25, 0.35),
           beta = c(-3.70, 0.02, -0.80), kappa = 0.9, R_mc = 1000)
set.seed(seed$C4)
C4 <- gen_C4(P4)
mean_usage <- exp(P4$gamma[1]) * sum(P4$p_hh * exp(P4$gamma[2] * P4$hh_sizes)) * exp(P4$gamma[3]^2 / 2)
truth_C4 <- list(beta = P4$beta[3], beta_vec = P4$beta, gamma = P4$gamma, shape = P4$kappa, mean_usage = mean_usage,
                 var_u = P4$gamma[2]^2 * (sum(P4$p_hh * P4$hh_sizes^2) - sum(P4$p_hh * P4$hh_sizes)^2) + P4$gamma[3]^2)

## =================================================================
## Save
## =================================================================
saveRDS(C1, "sorrel_C1.rds"); saveRDS(C2, "sorrel_C2.rds"); saveRDS(C3, "sorrel_C3.rds"); saveRDS(C4, "sorrel_C4.rds")
tr <- list(C1 = truth_C1, C2 = truth_C2, C3 = truth_C3, C4 = truth_C4,
           params = list(seed = seed, C1 = P1, C2 = P2, C3 = P3, C4 = P4))
saveRDS(tr, "truth06.rds")

cat(sprintf("C1: %d subscribers; censored %.3f (population %.3f); median tenure (60 USD, no referral) %.2f months\n",
            nrow(C1), 1 - mean(C1$delta), truth_C1$share_cens_pop, truth_C1$median_fee60))
cat(sprintf("C2: %d subscribers; censored %.3f; deal-seekers %.3f; pseudo-true kappa_* %.4f (2e6 draws) %.4f (quadrature)\n",
            nrow(C2), 1 - mean(C2$delta), mean(C2$deal_true), truth_C2$pseudo, truth_C2$quad$kappa))
cat(sprintf("    pseudo-true theta_* (large sample): %s; quadrature: %s\n",
            paste(sprintf("%.5f", fl$par), collapse = " "), paste(sprintf("%.5f", qd$par), collapse = " ")))
cat(sprintf("C3: %d person-period rows, %d exits in months 1 to 24 (fewest in a month: %d)\n",
            nrow(C3), sum(C3$y), min(tapply(C3$y, C3$m, sum))))
cat(sprintf("C4: %d subscribers; mean first-week check-ins %.3f (population %.3f); censored %.3f\n",
            nrow(C4), mean(C4$usage_wk1), mean_usage, 1 - mean(C4$delta)))
cat("Saved sorrel_C1.rds ... sorrel_C4.rds, truth06.rds\n")
