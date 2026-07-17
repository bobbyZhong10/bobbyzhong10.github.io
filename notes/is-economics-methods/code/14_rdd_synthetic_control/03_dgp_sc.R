# =====================================================================
# 14_rdd_synthetic_control / 03_dgp_sc.R
# Aster marketplace -- synthetic control design line (single-metro fee
# reform). One treated metro adopts a seller fee reform at month T0;
# ~20 donor metros do not. We want the reform's effect on metro-level
# log GMV.
#
# The panel is generated with an INTERACTIVE FIXED-EFFECTS / FACTOR
# structure:  Y_it = delta_i + lambda_i . f_t + tau_it * treat_it + eps
#   - 2 common factors f_t (a rising adoption trend + a cyclical term);
#   - metro-specific loadings lambda_i on those factors.
# Because loadings differ across metros, metros have DIFFERENT trends,
# so DiD's parallel-trends assumption FAILS. But the treated metro's
# loadings are a convex combination of a few donors' loadings, so a
# convex weight vector reproduces its factor exposure -> synthetic
# control fits the pre-period tightly and recovers the effect.
# =====================================================================

set.seed(12021)

n_donors <- 20
n_units  <- n_donors + 1        # unit 1 = treated ("Aster-Metro")
Tt       <- 36                   # months
T0       <- 25                   # reform starts at month 25 (12 post months)

# ---- two common factors --------------------------------------------
t_idx <- 1:Tt
# factor 1: a rising platform-adoption trend (concave growth)
f1 <- 2.0 * (1 - exp(-t_idx / 12)) + 0.05 * t_idx
# factor 2: a cyclical / seasonal component
f2 <- sin(2 * pi * t_idx / 12) + 0.02 * t_idx
F  <- cbind(f1, f2)

# ---- donor loadings and levels -------------------------------------
delta_d  <- runif(n_donors, 4.5, 6.0)          # metro base log GMV
lam1_d   <- runif(n_donors, 0.20, 1.60)         # loading on trend factor
lam2_d   <- runif(n_donors, -0.60, 0.60)        # loading on cyclical factor

# ---- treated loadings = convex combo of a few "anchor" donors ------
# guarantees the treated metro sits INSIDE the donor convex hull, so a
# convex synthetic control can reproduce its factor exposure.
anchors  <- c(3, 7, 12, 18)
w_true   <- c(0.35, 0.30, 0.20, 0.15)          # true generating weights
delta_tr <- sum(w_true * delta_d[anchors])
lam1_tr  <- sum(w_true * lam1_d[anchors])
lam2_tr  <- sum(w_true * lam2_d[anchors])

delta <- c(delta_tr, delta_d)
lam1  <- c(lam1_tr,  lam1_d)
lam2  <- c(lam2_tr,  lam2_d)

# ---- true treatment effect path (fee reform lifts GMV) -------------
# zero pre-reform, then a smooth ramp to a plateau of ~0.12 log points.
tau_path <- numeric(Tt)
post <- t_idx >= T0
tau_path[post] <- 0.12 * (1 - exp(-(t_idx[post] - T0 + 1) / 3))

# ---- assemble the panel --------------------------------------------
sigma_e <- 0.02
Y <- matrix(NA_real_, n_units, Tt)
for (i in 1:n_units) {
  common <- delta[i] + lam1[i] * f1 + lam2[i] * f2
  Y[i, ] <- common + rnorm(Tt, 0, sigma_e)
}
# add the reform effect to the treated unit's post periods
Y[1, ] <- Y[1, ] + tau_path

# true effect (measured on the noise-free counterfactual)
true_effect_path <- tau_path
true_post_avg    <- mean(tau_path[post])

# ---- long data frame ----------------------------------------------
unit_names <- c("Aster-Metro", sprintf("Donor-%02d", 1:n_donors))
panel <- data.frame(
  metro = factor(rep(unit_names, times = Tt), levels = unit_names),
  unit  = rep(1:n_units, times = Tt),
  month = rep(t_idx, each = n_units),
  log_gmv = as.vector(Y)
)
panel$treated  <- as.integer(panel$unit == 1)
panel$post     <- as.integer(panel$month >= T0)
panel$treat_it <- panel$treated * panel$post

# ---- what a naive DiD would say (to show it FAILS) -----------------
# 2x2 DiD: treated vs simple donor average, pre vs post.
tr_pre  <- mean(Y[1, !post]);          tr_post  <- mean(Y[1, post])
do_pre  <- mean(colMeans(Y[-1, !post])); do_post <- mean(colMeans(Y[-1, post]))
did_est <- (tr_post - tr_pre) - (do_post - do_pre)

# pre-period trend gap (parallel-trends check): treated vs donor avg
tr_pre_slope <- coef(lm(Y[1, !post] ~ t_idx[!post]))[2]
do_pre_slope <- coef(lm(colMeans(Y[-1, !post]) ~ t_idx[!post]))[2]

meta_sc <- list(
  n_donors = n_donors, Tt = Tt, T0 = T0,
  anchors = anchors, w_true = w_true,
  true_post_avg = true_post_avg,
  true_effect_path = true_effect_path,
  did_est = did_est, true_did_target = true_post_avg,
  tr_pre_slope = as.numeric(tr_pre_slope),
  do_pre_slope = as.numeric(do_pre_slope),
  sigma_e = sigma_e
)

out_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(out_dir) == 0 || out_dir == "") out_dir <- "."
saveRDS(list(panel = panel, Y = Y, meta = meta_sc,
             unit_names = unit_names),
        file.path(out_dir, "data_sc.rds"))

cat(sprintf("SC DGP: %d donors + 1 treated, T=%d, reform at T0=%d\n",
            n_donors, Tt, T0))
cat(sprintf("  true post-period avg effect     = %.3f\n", true_post_avg))
cat(sprintf("  true effect at last month       = %.3f\n", tau_path[Tt]))
cat(sprintf("  naive 2x2 DiD estimate          = %.3f  (target %.3f)\n",
            did_est, true_post_avg))
cat(sprintf("  DiD bias                        = %.3f\n",
            did_est - true_post_avg))
cat(sprintf("  pre-trend slope treated/donors  = %.4f / %.4f  (parallel? no)\n",
            meta_sc$tr_pre_slope, meta_sc$do_pre_slope))
cat(sprintf("  true generating weights (donors %s) = %s\n",
            paste(anchors, collapse=","), paste(w_true, collapse=",")))
