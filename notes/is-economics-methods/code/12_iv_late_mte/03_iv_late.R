# ============================================================================
# 03_iv_late.R  --  first stage, ITT (reduced form), Wald = 2SLS = complier LATE
# ----------------------------------------------------------------------------
# The randomized encouragement Z is the instrument.  We show:
#   (1) first stage : Z shifts adoption D (strong), with F and effective F
#   (2) ITT / reduced form : Z's intention-to-treat effect on Y (clean but
#       DILUTED by never/always-takers)
#   (3) Wald = ITT / first-stage = 2SLS(just-identified) = complier LATE
#   (4) compare to TRUE complier-LATE, and describe who the compliers are.
# ============================================================================

suppressMessages(library(fixest))

d <- readRDS("meridian_data.rds")
dat <- d$dat; truth <- d$truth; par <- d$params

## ---- (1) first stage: D ~ Z -----------------------------------------------
fs <- feols(D ~ Z, data = dat)
fs_coef <- coef(fs)["Z"]
# classical (IID) first-stage F and heteroskedasticity-robust F (effective F).
# For a single endogenous regressor with a single instrument the robust
# first-stage Wald equals the Montiel-Olea-Pflueger effective F.
F_iid <- fitstat(fs, "wald", vcov = "iid")$wald$stat
F_hc  <- fitstat(fs, "wald", vcov = "hetero")$wald$stat
eff_F <- F_hc

## ---- (2) ITT / reduced form: Y ~ Z ----------------------------------------
rf <- feols(Y ~ Z, data = dat)
itt      <- coef(rf)["Z"]
itt_se   <- se(rf)["Z"]

## ---- (3) Wald estimator = ITT / first stage -------------------------------
wald <- itt / fs_coef

## 2SLS via feols IV syntax:  Y ~ 1 | D ~ Z
iv <- feols(Y ~ 1 | D ~ Z, data = dat)
b_2sls  <- coef(iv)["fit_D"]
se_2sls <- se(iv)["fit_D"]

## delta-method SE for the Wald ratio (should match 2SLS SE closely)
itt_var <- itt_se^2
fs_se   <- se(fs)["Z"]
wald_se <- sqrt(itt_var / fs_coef^2)   # leading term (fs very precise here)

## ---- (4) complier share and complier characteristics ----------------------
# complier share = first stage (Abadie): P(complier) = E[D|Z=1]-E[D|Z=0]
complier_share_hat <- fs_coef
# Abadie kappa weights to describe compliers on an observed covariate W:
#   E[g(X)|complier] = E[kappa * g(X)] / E[kappa]
pz_hat <- mean(dat$Z)                       # P(Z=1), randomized ~ 0.5
kappa <- 1 - dat$D * (1 - dat$Z) / (1 - pz_hat) - (1 - dat$D) * dat$Z / pz_hat
w_mean_all      <- mean(dat$W_size)
w_mean_complier <- sum(kappa * dat$W_size) / sum(kappa)
# oracle check: true complier covariate mean
w_mean_complier_oracle <- mean(dat$W_size[dat$is_complier])

cat("=== First stage ===\n")
cat(sprintf("E[D|Z=1]-E[D|Z=0]      = %.4f  (SE %.4f)\n", fs_coef, fs_se))
cat(sprintf("First-stage F (IID)    = %.2f\n", F_iid))
cat(sprintf("Effective F (robust)   = %.2f\n", eff_F))
cat("\n=== ITT (reduced form) ===\n")
cat(sprintf("ITT  b(Z on Y)         = %.4f  (SE %.4f)\n", itt, itt_se))
cat("\n=== Wald / 2SLS = complier LATE ===\n")
cat(sprintf("Wald = ITT/firststage  = %.4f  (SE %.4f)\n", wald, wald_se))
cat(sprintf("2SLS (feols IV)        = %.4f  (SE %.4f)\n", b_2sls, se_2sls))
cat(sprintf("TRUE complier-LATE     = %.4f\n", truth$LATE))
cat(sprintf("(for contrast) TRUE ATE= %.4f , TRUE ATT = %.4f\n", truth$ATE, truth$ATT))
cat("\n=== Compliers ===\n")
cat(sprintf("Complier share (=1st stage)      = %.4f\n", complier_share_hat))
cat(sprintf("E[W | complier] (Abadie kappa)   = %.4f\n", w_mean_complier))
cat(sprintf("E[W | complier] (oracle)         = %.4f\n", w_mean_complier_oracle))
cat(sprintf("E[W | all firms]                 = %.4f\n", w_mean_all))

saveRDS(list(fs_coef = fs_coef, fs_se = fs_se, F_iid = F_iid, eff_F = eff_F,
             itt = itt, itt_se = itt_se, wald = wald, wald_se = wald_se,
             b_2sls = b_2sls, se_2sls = se_2sls,
             complier_share = complier_share_hat,
             w_complier = w_mean_complier, w_complier_oracle = w_mean_complier_oracle,
             w_all = w_mean_all),
        "res_iv.rds")
