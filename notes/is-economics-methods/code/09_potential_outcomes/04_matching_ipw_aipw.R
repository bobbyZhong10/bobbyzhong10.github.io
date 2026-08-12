# ============================================================================
# 04_matching_ipw_aipw.R  —  Study B (observational, CIA arm), ATT estimators
#   - propensity score (logit)
#   - nearest-neighbour matching (MatchIt, 1:1 on PS, ATT) point estimate
#     with the Abadie-Imbens SE (via Matching::Match on the same PS)
#   - IPW: Horvitz-Thompson (un-normalized) and Hajek (normalized), ATT
#   - AIPW / doubly-robust ATT
#   Point estimate + SE for each.
# ============================================================================
suppressMessages({
  library(MatchIt); library(Matching); library(marginaleffects)
})

d   <- readRDS("northwind.rds")
obs <- d$obs
obs$D <- obs$D_cia
obs$Y <- obs$Y_cia
Xrhs <- ~ size + dig + eff + ind

# ---- propensity score (logit) ---------------------------------------------
ps_fit <- glm(update(Xrhs, D ~ .), data = obs, family = binomial)
obs$e  <- predict(ps_fit, type = "response")

# ============================================================================
# (1) Nearest-neighbour matching on the PS (ATT)
# ============================================================================
# MatchIt: 1:1 NN on the logit of the PS, with replacement (ATT target)
mi <- matchit(update(Xrhs, D ~ .), data = obs, method = "nearest",
              distance = "glm", link = "logit", estimand = "ATT",
              replace = TRUE)
mdat <- match.data(mi)
fit_m <- lm(Y ~ D, data = mdat, weights = weights)
match_pt <- coef(fit_m)["D"]

# genuine Abadie-Imbens variance via Matching::Match on the same linear PS
lps <- predict(ps_fit)            # linear predictor (logit PS)
mout <- Match(Y = obs$Y, Tr = obs$D, X = lps, estimand = "ATT",
              M = 1, replace = TRUE, ties = TRUE, BiasAdjust = FALSE)
match_ai_pt <- mout$est[1]
match_ai_se <- mout$se[1]         # Abadie-Imbens SE

# a naive (wrong) SE for contrast: OLS SE on the matched sample
match_naive_se <- summary(fit_m)$coefficients["D","Std. Error"]

# ============================================================================
# (2) IPW for the ATT  (weights: treated 1 ; control e/(1-e))
# ============================================================================
ipw_att <- function(Y, D, e) {
  w <- e / (1 - e)
  m1 <- mean(Y[D == 1])                                   # treated mean
  # Horvitz-Thompson (un-normalized control weights, divide by n1)
  ht  <- m1 - sum((1 - D) * w * Y) / sum(D)
  # Hajek (normalized control weights)
  haj <- m1 - sum((1 - D) * w * Y) / sum((1 - D) * w)
  c(ht = ht, hajek = haj)
}
ipw_pt <- ipw_att(obs$Y, obs$D, obs$e)

# ============================================================================
# (3) AIPW / doubly-robust ATT
#   tau = mean{ D(Y-m0) - (1-D) e/(1-e) (Y-m0) } / mean(D)
#   m0 = E[Y | X, D=0]
# ============================================================================
aipw_att <- function(dat, e) {
  m0_fit <- lm(update(Xrhs, Y ~ .), data = dat[dat$D == 0, ])
  m0 <- predict(m0_fit, newdata = dat)
  p  <- mean(dat$D)
  w  <- e / (1 - e)
  psi <- (dat$D * (dat$Y - m0) - (1 - dat$D) * w * (dat$Y - m0)) / p
  mean(psi)
}
aipw_pt <- aipw_att(obs, obs$e)

# ============================================================================
# SEs for IPW / AIPW via nonparametric bootstrap (PS + outcome refit each draw)
#   (bootstrap is valid for IPW/AIPW; note: NOT valid for NN matching)
# ============================================================================
set.seed(202)
B <- 1000L
boot <- t(replicate(B, {
  ix <- sample(nrow(obs), replace = TRUE)
  db <- obs[ix, ]
  psb <- glm(update(Xrhs, D ~ .), data = db, family = binomial)
  eb  <- predict(psb, type = "response")
  ip  <- ipw_att(db$Y, db$D, eb)
  ai  <- aipw_att(db, eb)
  c(ip, aipw = ai)
}))
se_boot <- apply(boot, 2, sd)

cat("==== Study B : matching / IPW / AIPW (CIA arm) ====\n")
cat(sprintf("propensity range        : [%.4f, %.4f]\n", min(obs$e), max(obs$e)))
cat(sprintf("NN matching (MatchIt)    : %.4f\n", match_pt))
cat(sprintf("NN matching (Matching)   : %.4f  AI-SE %.4f  (naive OLS-SE %.4f)\n",
            match_ai_pt, match_ai_se, match_naive_se))
cat(sprintf("IPW Horvitz-Thompson     : %.4f  SE %.4f\n", ipw_pt["ht"],    se_boot["ht"]))
cat(sprintf("IPW Hajek                : %.4f  SE %.4f\n", ipw_pt["hajek"], se_boot["hajek"]))
cat(sprintf("AIPW (doubly robust)     : %.4f  SE %.4f\n", aipw_pt,         se_boot["aipw"]))
cat(sprintf("true ATT (CIA)           : %.4f\n", d$truths$true_ATT_cia))

saveRDS(list(
  match_pt = match_pt, match_ai_pt = match_ai_pt,
  match_ai_se = match_ai_se, match_naive_se = match_naive_se,
  ipw_ht = ipw_pt["ht"], se_ht = se_boot["ht"],
  ipw_hajek = ipw_pt["hajek"], se_hajek = se_boot["hajek"],
  aipw = aipw_pt, se_aipw = se_boot["aipw"],
  e_range = range(obs$e)),
  "res_matching_ipw_aipw.rds")
