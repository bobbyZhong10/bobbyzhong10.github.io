# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Sorrel data and truth06.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("ml_fns.R")
tr <- readRDS("truth06.rds")
C1 <- readRDS("sorrel_C1.rds"); C2 <- readRDS("sorrel_C2.rds"); C3 <- readRDS("sorrel_C3.rds"); C4 <- readRDS("sorrel_C4.rds")
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))

## ---- truth06 fields named on the chapter card ---------------------------------------------------
stopifnot(fin(tr$C1$beta), fin(tr$C1$shape), fin(tr$C2$pseudo), fin(tr$C4$beta))
# every random block has its own seed, all different (C3 has none: it is C1 regrouped)
sd <- unlist(tr$params$seed)
stopifnot(length(sd) >= 6, !anyDuplicated(sd), all(sd == round(sd)), is.null(tr$params$seed$C3))

## ---- the generators reproduce the saved data from the stored seeds ---------------------------------
for (b in c("C1", "C2", "C4")) {
  set.seed(tr$params$seed[[b]])
  regen <- get(paste0("gen_", b))(tr$params[[b]])
  stopifnot(isTRUE(all.equal(regen, get(b), tolerance = 0, check.attributes = TRUE)))
}

## ---- C1: Weibull subscription tenure, right-censored at the data pull ---------------------------------
P1 <- tr$params$C1
stopifnot(is.data.frame(C1), nrow(C1) == P1$n,
          all(c("id", "signup_month", "fee_usd", "referred", "c_months", "t_months", "delta", "t_star_true") %in% names(C1)))
stopifnot(all(C1$c_months == P1$months - C1$signup_month + 1), all(C1$c_months %in% seq_len(P1$months)),
          identical(C1$t_months, pmin(C1$t_star_true, C1$c_months)),
          identical(C1$delta, as.numeric(C1$t_star_true <= C1$c_months)),
          all(C1$fee_usd %in% P1$fees), setequal(unique(C1$referred), 0:1))
cs1 <- 1 - mean(C1$delta)
stopifnot(cs1 > 0.30, cs1 < 0.55)                                                    # a growing platform: many recent signups
stopifnot(tr$C1$beta == P1$beta[2], tr$C1$shape == P1$kappa, identical(tr$C1$beta_vec, P1$beta),
          abs(tr$C1$beta_aft - (-P1$beta[2] / P1$kappa)) < 1e-15, P1$beta[2] > 0, P1$kappa < 1)
# latent tenure is Weibull: S(12 | fee 60, not referred) = exp(-exp(x'beta) 12^kappa)
sel <- C1$fee_usd == 60 & C1$referred == 0
S12 <- exp(-exp(sum(P1$beta * c(1, 60, 0))) * 12^P1$kappa)
stopifnot(abs(mean(C1$t_star_true[sel] > 12) - S12) < 0.04, abs(tr$C1$S12_fee60 - S12) < 1e-15,
          tr$C1$median_fee60 > 5, tr$C1$median_fee60 < 12)                             # median tenure in months
# independent censoring (Wooldridge 2010 eq. 22.22): the signup month carries no information on t*
stopifnot(abs(cor(C1$signup_month, C1$t_star_true, method = "spearman")) < 0.05)

## ---- C2: log-logistic mixture (deal-seekers and regulars), fitted by a Weibull --------------------------
P2 <- tr$params$C2
stopifnot(is.data.frame(C2), nrow(C2) == P2$n,
          all(c("id", "signup_month", "fee_usd", "referred", "c_months", "t_months", "delta", "deal_true", "t_star_true") %in% names(C2)))
stopifnot(all(C2$c_months == P2$months - C2$signup_month + 1), identical(C2$t_months, pmin(C2$t_star_true, C2$c_months)),
          identical(C2$delta, as.numeric(C2$t_star_true <= C2$c_months)), setequal(unique(C2$deal_true), 0:1))
stopifnot(abs(mean(C2$deal_true[C2$referred == 0]) - plogis(P2$g0)) < 0.03,
          abs(mean(C2$deal_true[C2$referred == 1]) - plogis(P2$g0 + P2$g_ref)) < 0.03)
# the mixture's latent survival function is right (share of t* > 6 months among fee 60, not referred)
sel2 <- C2$fee_usd == 60 & C2$referred == 0
stopifnot(abs(mean(C2$t_star_true[sel2] > 6) - c2_surv(6, 60, 0, P2)) < 0.04)
# pseudo-true value: large sample (the card's definition) against quadrature of the population objective
stopifnot(tr$C2$large$N >= 1e6, tr$C2$large$converged, fin(tr$C2$pseudo_vec, 4L), fin(tr$C2$quad$par, 4L),
          abs(tr$C2$pseudo - exp(tr$C2$pseudo_vec[4])) < 1e-12, abs(tr$C2$pseudo_fee - tr$C2$pseudo_vec[2]) < 1e-15,
          all(abs(tr$C2$pseudo_vec - tr$C2$quad$par) < 4 * tr$C2$large$se_sand))
q20 <- pseudo_true_weibull(P2, n_gl = 20)                                            # quadrature is stable
stopifnot(q20$converged, max(abs(q20$par - tr$C2$quad$par)) < 1e-7)
# misspecification shows: the pseudo-true shape says "hazard falls with tenure", while the reference
# subscriber's true hazard falls, rises and falls again (not a Weibull)
stopifnot(tr$C2$pseudo < 0.9)
hg <- c2_hazard(seq(0.25, 36, by = 0.25), 60, 0, P2)
stopifnot(any(diff(hg) > 0), any(diff(hg) < 0))
stopifnot(P2$R_mc >= 1000)

## ---- C3: C1 regrouped into monthly person-period rows (first 24 months of tenure) ------------------------
P3 <- tr$params$C3
stopifnot(identical(P3$cuts, 0:P3$M))
pp <- person_period(C1$t_months, C1$delta, cbind(fee_usd = C1$fee_usd, referred = C1$referred), P3$cuts)
stopifnot(identical(pp, C3), identical(names(C3), c("id", "m", "y", "fee_usd", "referred")))
ex <- tapply(C3$y, C3$m, sum)
stopifnot(length(ex) == P3$M, all(ex >= 5))                                          # every interval has exits: the MLE exists
stopifnot(sum(C3$y) == sum(C1$delta == 1 & C1$t_months <= P3$M),
          nrow(C3) == sum(pmin(ifelse(C1$delta == 1, ceiling(C1$t_months), floor(C1$t_months)), P3$M)))
eta_true <- P1$beta[1] + log(P3$cuts[-1]^P1$kappa - P3$cuts[-(P3$M + 1)]^P1$kappa)   # eta_m = log of the integrated baseline hazard
stopifnot(max(abs(tr$C3$eta - eta_true)) < 1e-12, identical(tr$C3$beta, P1$beta[2:3]))

## ---- C4: first-week usage (Poisson) and the usage-intensity index in the hazard ---------------------------
P4 <- tr$params$C4
stopifnot(is.data.frame(C4), nrow(C4) == P4$n,
          all(c("id", "signup_month", "fee_usd", "hh_size", "quiz_std", "usage_wk1", "c_months", "t_months", "delta",
                "u_true", "t_star_true") %in% names(C4)))
stopifnot(max(abs(C4$u_true - (P4$gamma[1] + P4$gamma[2] * C4$hh_size + P4$gamma[3] * C4$quiz_std))) < 1e-12,
          all(C4$usage_wk1 >= 0), all(C4$usage_wk1 == round(C4$usage_wk1)),
          mean(C4$usage_wk1) > 0.5, mean(C4$usage_wk1) < 1.5,                         # sparse first-week signal
          identical(C4$t_months, pmin(C4$t_star_true, C4$c_months)),
          identical(C4$delta, as.numeric(C4$t_star_true <= C4$c_months)))
stopifnot(tr$C4$beta == P4$beta[3], identical(tr$C4$beta_vec, P4$beta), identical(tr$C4$gamma, P4$gamma),
          P4$beta[3] != 0, P4$R_mc >= 1000)
