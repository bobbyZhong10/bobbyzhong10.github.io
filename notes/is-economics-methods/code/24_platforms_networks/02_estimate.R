# ------------------------------------------------------------------
# 02_estimate.R -- estimating state dependence (switching cost) in
# the Cadence panel, three ways, to expose spurious state dependence.
#
#   (1) naive pooled logit: choiceA ~ dfee + Lprev
#         no control for theta_i -> delta badly overstated.
#   (2) random-effects logit WITHOUT initial-conditions fix:
#         choiceA ~ dfee + Lprev + (1|merchant)
#         still biased because the month-1 choice is correlated with
#         theta_i (Heckman/Wooldridge initial-conditions problem).
#   (3) Wooldridge (2005) correlated-random-effects logit: add the
#         initial choice initA as an auxiliary regressor for theta_i:
#         choiceA ~ dfee + Lprev + initA + (1|merchant)
#         -> recovers delta.
# Reported switching cost in dollars = delta_hat / |beta_hat|.
# ------------------------------------------------------------------
library(data.table)
library(lme4)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dt    <- readRDS(file.path(code_dir, "cadence_panel.rds"))
truth <- readRDS(file.path(code_dir, "truth24.rds"))
est   <- dt[month >= 2]

sc <- function(delta, beta) delta / abs(beta)

## ---- (1) naive pooled logit -------------------------------------
m_naive <- glm(choiceA ~ dfee + Lprev, family = binomial, data = est)
b_naive <- coef(m_naive)
se_naive <- sqrt(diag(vcov(m_naive)))

## ---- (2) RE logit, no initial-conditions correction -------------
m_re <- glmer(choiceA ~ dfee + Lprev + (1 | merchant),
              family = binomial, data = est,
              control = glmerControl(optimizer = "bobyqa"),
              nAGQ = 1)
b_re  <- fixef(m_re)
se_re <- sqrt(diag(as.matrix(vcov(m_re))))

## ---- (3) Wooldridge CRE logit (initial choice as aux regressor) -
m_cre <- glmer(choiceA ~ dfee + Lprev + initA + (1 | merchant),
               family = binomial, data = est,
               control = glmerControl(optimizer = "bobyqa"),
               nAGQ = 1)
b_cre  <- fixef(m_cre)
se_cre <- sqrt(diag(as.matrix(vcov(m_cre))))

## ---- collect ----------------------------------------------------
tab <- data.table(
  model = c("True (DGP)", "Naive pooled logit",
            "RE logit (no IC fix)", "Wooldridge CRE logit"),
  beta  = c(truth$beta, b_naive["dfee"], b_re["dfee"], b_cre["dfee"]),
  delta = c(truth$delta, b_naive["Lprev"], b_re["Lprev"], b_cre["Lprev"]),
  delta_se = c(NA, se_naive["Lprev"], se_re["Lprev"], se_cre["Lprev"]),
  sigma_theta = c(truth$sigma_theta, NA,
                  sqrt(unlist(VarCorr(m_re))[1]),
                  sqrt(unlist(VarCorr(m_cre))[1]))
)
tab[, sc_dollars := sc(delta, beta)]

cat("== State dependence estimates (delta_true = 0.90) ==\n")
print(tab[, .(model,
              beta = round(beta, 3),
              delta = round(delta, 3),
              delta_se = round(delta_se, 3),
              sigma_theta = round(sigma_theta, 3),
              sc_dollars = round(sc_dollars, 2))])

cat(sprintf("\nNaive delta overstates truth by %.0f%% (%.3f vs %.3f)\n",
            100 * (b_naive["Lprev"] / truth$delta - 1),
            b_naive["Lprev"], truth$delta))
cat(sprintf("Naive implied switching cost $%.2f/mo vs true $%.2f/mo\n",
            sc(b_naive["Lprev"], b_naive["dfee"]), truth$sc_dollars))
cat(sprintf("Wooldridge CRE delta = %.3f (SE %.3f), switching cost $%.2f/mo\n",
            b_cre["Lprev"], se_cre["Lprev"], sc(b_cre["Lprev"], b_cre["dfee"])))

saveRDS(list(tab = tab,
             naive = list(b = b_naive, se = se_naive),
             re = list(b = b_re, se = se_re),
             cre = list(b = b_cre, se = se_cre)),
        file.path(code_dir, "res24_estimate.rds"))
cat("Saved res24_estimate.rds\n")
