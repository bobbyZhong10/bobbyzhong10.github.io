# =====================================================================
# 22_synthetic_control_factor_models / 01_dgp.R
# The simulated Aster marketplace for chapter 22 (synthetic control and
# panel factor models). One script, four data blocks, one truth file:
#
#   S1  aster_S1.rds  single-metro seller fee reform, treated metro INSIDE
#                     the donor hull (the legacy block from the old
#                     03_dgp_sc.R; seed seeds22[["S1"]])
#   S2  aster_S2.rds  the same reform launched in Aster's highest-level metro
#                     (Aster-Capital), which lies OUTSIDE the donor hull
#                     (seed seeds22[["S2"]])
#   S3  aster_S3.rds  staggered rollout of a seller relief program in
#                     three waves; the program goes to the metros hit
#                     hardest by a demand downturn (largest cycle
#                     loading), with wave order nearly random among them,
#                     so parallel trends between adopters and never-
#                     treated metros fail (gen_S3() in dgp_s3.R; seed seeds22[["S3"]],
#                     chosen by a seed screen documented there)
#   S4  aster_S4.rds  pilot design: pick 5 of 40 metros in advance for a
#                     6-month commission pilot, with replicate pilot-period
#                     futures for RMSE comparisons (seed seeds22[["S4"]])
#   truth22.rds       list(S1 = list(att_post, tau_path), S2 = list(att_post,
#                     tau_path), S3 = list(att, att_gt), S4 = list(tau), seeds)
#
# Every block is a factor model for untreated log GMV (notation_v3.md 2.9(b)):
#   Y_it(0) = alpha_i + lambda_t + mu_i' f_t + u_it,
# alpha_i the metro level, lambda_t a common time effect, f_t the common
# factors, mu_i the metro's loadings, u_it idiosyncratic noise. The legacy
# S1 code calls the level delta and the loadings lam1, lam2; S2 keeps those
# names so the two single-unit blocks read alike.
#
# Random-number discipline: the S1 block is the legacy code with two
# changes, its set.seed() argument read from seeds22 (value 12021, the
# legacy literal) and its output file renamed; nothing draws random
# numbers before it, so its stream is the legacy stream draw for draw.
# Each new block starts with its own set.seed() (S2 and S4 inside local(),
# S3 inside gen_S3()), so S2 to S4 never move S1 and do not move one
# another. The four seeds are defined once, in seeds22 below, and reused
# by the set.seed() calls, by gen_S3() and in truth22$seeds.
#
# Time convention: S1 and S2 keep the legacy meta$T0 = 25, the FIRST reform
# month (the legacy comment "reform starts at month T0"). The chapter's T_0
# (notation 2.9(a)) is the LAST untreated month, 24; S2 stores it as
# meta$n_pre, and the functions in sc_fns.R take T0 = 24.
#
# Run from this directory:  Rscript 01_dgp.R   (a few seconds)
# =====================================================================

# The block seeds, in one place. S1 = 12021 is the inherited legacy seed:
# changing it changes every inherited chapter 22 number.
seeds22 <- c(S1 = 12021, S2 = 22002, S3 = 21, S4 = 22004)

# =====================================================================
# Block S1 (legacy 03_dgp_sc.R, verbatim except two lines: the seed is
# written set.seed(seeds22[["S1"]]), still 12021, and the output file is
# aster_S1.rds instead of data_sc.rds). Its original header:
#
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

set.seed(seeds22[["S1"]])

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
        file.path(out_dir, "aster_S1.rds"))

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
rm(F)   # the legacy F <- cbind(f1, f2) masks FALSE; no random numbers involved

# =====================================================================
# After S1: its truth, and the hull LP from sc_fns.R (sourcing it draws
# no random numbers). Ruling S1: tau_path holds the 12 post-period true
# effects (months 25..36) and att_post is their mean.
# =====================================================================
truth_S1 <- list(att_post = true_post_avg, tau_path = tau_path[post])
source(file.path(out_dir, "sc_fns.R"))

# =====================================================================
# Block S2: the fee reform in Aster's highest-level metro (outside the hull)
# ---------------------------------------------------------------------
# Same world as S1: two factors (f1 a platform-adoption trend, f2 a
# seasonal cycle), 20 donor metros whose levels and loadings come from
# the S1 distributions, 36 months, reform from month 25 (12 post months).
# This time the reform launches in "Aster-Capital", Aster's flagship
# metro. Its base log GMV is 6.20 while donor levels are U(4.5, 6.0), so
# no donor and no convex combination of donors reaches its level: the
# treated metro is outside the donor convex hull. Its loadings are
# typical (lam1 = 0.90 and lam2 = 0.00 are the centres of the donor
# ranges): the problem is the level, as it is for the largest unit in
# most real panels. The effect ramps to a plateau of 0.10 log points.
#
# meta$treated_outside_hull: hull membership is judged in the space of
# the NOISE-FREE generating parameters (delta_i, lam1_i, lam2_i), the
# metro level plus its two factor loadings, by LP feasibility:
# in_hull_lp() in sc_fns.R minimises sum |X0 w - x1| over the simplex
# with lpSolve, and the metro is inside iff that minimum is 0. (The
# observed 24-month pre-period outcome vector is almost never inside the
# hull of 20 noisy donor vectors, so testing it would make the flag
# vacuous.) tests/test_dgp.R cross-checks the flag with
# pensynth::in_convex_hull.
#
# Built for the Task 10 conclusions on S2. Task 9 feasibility run at the
# shipped level 6.20 (seeds 1 to 100, otherwise this design; SC on the 24
# pre-period outcomes with V = I, augsynth ridge ASC, synthdid): SC
# pre-period RMSPE above 3 x 0.0165 in 100/100 seeds (mean 0.165, minimum
# 0.086); synthetic DiD closer to the truth than SC in 100/100 and ridge
# ASC in 99/100 (mean absolute errors 0.020 and 0.018 vs 0.219 for SC);
# a negative ridge-ASC weight in 100/100. Penalized SC used no more donors
# than SC in 98%, 90%, 87% and 97% of seeds at pen = 0.001, 0.01, 0.1, 1.
# (A first run at level 6.25 gave the same picture; 6.20 was kept as the
# smaller gap above the donor range.) The Task 9 report gives the numbers
# for the stored S2 seed, including its penalty path.
# =====================================================================
truth_S2 <- local({
  set.seed(seeds22[["S2"]])

  n_donors <- 20
  n_units  <- n_donors + 1        # unit 1 = treated ("Aster-Capital")
  Tt       <- 36                  # months
  T0       <- 25                  # FIRST reform month (legacy S1 convention)
  n_pre    <- T0 - 1              # the chapter's T_0: last untreated month

  t_idx <- 1:Tt
  f1 <- 2.0 * (1 - exp(-t_idx / 12)) + 0.05 * t_idx
  f2 <- sin(2 * pi * t_idx / 12) + 0.02 * t_idx

  delta_d <- runif(n_donors, 4.5, 6.0)          # donor base log GMV
  lam1_d  <- runif(n_donors, 0.20, 1.60)         # loading on trend factor
  lam2_d  <- runif(n_donors, -0.60, 0.60)        # loading on cyclical factor
  delta_tr <- 6.20                               # above every possible donor level
  lam1_tr  <- 0.90
  lam2_tr  <- 0.00
  delta <- c(delta_tr, delta_d)
  lam1  <- c(lam1_tr,  lam1_d)
  lam2  <- c(lam2_tr,  lam2_d)

  post <- t_idx >= T0
  tau_path <- numeric(Tt)
  tau_path[post] <- 0.10 * (1 - exp(-(t_idx[post] - T0 + 1) / 3))

  sigma_e <- 0.02
  Y <- matrix(NA_real_, n_units, Tt)
  for (i in 1:n_units) {
    Y[i, ] <- delta[i] + lam1[i] * f1 + lam2[i] * f2 + rnorm(Tt, 0, sigma_e)
  }
  Y[1, ] <- Y[1, ] + tau_path

  # hull membership of the treated metro in (level, loading1, loading2) space
  params <- rbind(level = delta, loading1 = lam1, loading2 = lam2)
  hull   <- in_hull_lp(params[, 1], params[, -1])

  unit_names <- c("Aster-Capital", sprintf("Donor-%02d", 1:n_donors))
  colnames(params) <- unit_names
  panel <- data.frame(
    metro = factor(rep(unit_names, times = Tt), levels = unit_names),
    unit  = rep(1:n_units, times = Tt),
    month = rep(t_idx, each = n_units),
    log_gmv = as.vector(Y)
  )
  panel$treated  <- as.integer(panel$unit == 1)
  panel$post     <- as.integer(panel$month >= T0)
  panel$treat_it <- panel$treated * panel$post

  # naive 2x2 DiD and pre-trend slopes, as in S1
  tr_pre  <- mean(Y[1, !post]);            tr_post <- mean(Y[1, post])
  do_pre  <- mean(colMeans(Y[-1, !post])); do_post <- mean(colMeans(Y[-1, post]))
  did_est <- (tr_post - tr_pre) - (do_post - do_pre)
  tr_pre_slope <- unname(coef(lm(Y[1, !post] ~ t_idx[!post]))[2])
  do_pre_slope <- unname(coef(lm(colMeans(Y[-1, !post]) ~ t_idx[!post]))[2])

  meta <- list(
    n_donors = n_donors, Tt = Tt, T0 = T0, n_pre = n_pre,
    true_post_avg = mean(tau_path[post]), true_effect_path = tau_path,
    did_est = did_est, tr_pre_slope = tr_pre_slope, do_pre_slope = do_pre_slope,
    sigma_e = sigma_e, f1 = f1, f2 = f2,
    params = params,                                   # noise-free level and loadings, by metro
    level_gap = delta_tr - max(delta_d),               # treated level minus the highest donor level
    treated_outside_hull = !as.vector(hull),
    hull_l1_gap = attr(hull, "l1_gap"),
    hull_method = paste("LP feasibility (in_hull_lp in sc_fns.R, lpSolve) in the space of the",
                        "noise-free generating parameters (level, loading1, loading2);",
                        "outside iff min sum|X0 w - x1| over the simplex > 0")
  )
  rownames(Y) <- unit_names
  saveRDS(list(panel = panel, Y = Y, meta = meta, unit_names = unit_names),
          file.path(out_dir, "aster_S2.rds"))

  cat(sprintf("\nS2 DGP: Aster-Capital + %d donors, T=%d, reform from month %d (n_pre = %d)\n",
              n_donors, Tt, T0, n_pre))
  cat(sprintf("  treated level %.2f vs highest donor level %.3f (gap %.3f)\n",
              delta_tr, max(delta_d), meta$level_gap))
  cat(sprintf("  treated outside donor hull (LP, parameter space) = %s  (L1 gap %.4f)\n",
              meta$treated_outside_hull, meta$hull_l1_gap))
  cat(sprintf("  true post-period avg effect     = %.4f\n", meta$true_post_avg))
  cat(sprintf("  naive 2x2 DiD estimate          = %.4f\n", did_est))

  list(att_post = mean(tau_path[post]), tau_path = tau_path[post])
})

# =====================================================================
# Block S3: staggered rollout with selection on factor exposure
# ---------------------------------------------------------------------
# The generator gen_S3(seed) lives in dgp_s3.R, with the full design
# notes (the cycle loading drives which metros get the program; wave
# order is nearly random among them) and the seed screen behind seed 21,
# so that Task 10 can redraw the block without regenerating S1, S2, S4.
# =====================================================================
source(file.path(out_dir, "dgp_s3.R"))
truth_S3 <- local({
  s3 <- gen_S3(seeds22[["S3"]])
  saveRDS(s3$data, file.path(out_dir, "aster_S3.rds"))

  m  <- s3$data$meta
  ad <- is.finite(m$cohort)
  cat(sprintf("\nS3 DGP: %d metros, T=%d, waves at months %s (%d metros each), %d contract metros\n",
              m$N, m$Tt, paste(m$waves, collapse = ", "), m$n_wave, m$n_contract))
  cat("  mean cycle loading by cohort (29, 33, 37, never):",
      paste(sprintf("%.3f", m$loading1_by_cohort), collapse = " / "), "\n")
  cat(sprintf("  adopters vs never-treated mean loading1 = %.3f vs %.3f\n",
              mean(m$loading1[ad]), mean(m$loading1[!ad])))
  cat(sprintf("  cor(loading1, adopt_rank)       = %.3f  (within adopters, cor(loading1, cohort) = %.3f)\n",
              cor(m$loading1, m$adopt_rank), cor(m$loading1[ad], m$cohort[ad])))
  cat(sprintf("  true ATT (all treated cells)    = %.4f\n", s3$truth$att))

  s3$truth
})

# =====================================================================
# Block S4: designing a pilot (synthetic experimental design)
# ---------------------------------------------------------------------
# 40 eligible mid-sized metros. When the pilot team must name 5 pilot
# metros, 24 months of history are available (months 1 to 24); the pilot,
# a lower seller commission, runs in months 25 to 30. Untreated outcome =
# log monthly GMV in USD millions:
#   Y_it(0) = alpha_i + lambda_t + mu_i1 f1_t + mu_i2 f2_t + mu_i3 g_t + e_it
#   alpha_i ~ U(-0.5, 1.0); lambda_t = 0.3 (1 - exp(-t / 18)) + 0.005 t
#   f1_t = the S1 adoption trend, centred over the 30 months; mu_i1 ~ U(0.4, 0.8)
#   f2_t = sin(2 pi t / 12);                                  mu_i2 ~ U(-0.15, 0.15)
#   g_t  = AR(1) common demand shock, rho 0.7, innovation sd 0.05,
#          stationary start;                                  mu_i3 ~ U(0, 1)
#   e_it ~ N(0, 0.02^2)
# The pilot raises log GMV by tau = 0.04 in every pilot metro-month: the
# analyst observes Y0 + tau * D_i in months 25 to 30, where D_i = 1 for
# the five metros the design picks (Y0 here is the untreated outcome;
# the design decides who gets tau).
# Replicates: Y0_pilot_reps[r, , ] (r = 1..500, metros x pilot months) are
# independent draws of the pilot window given the same history: new g_t
# innovations starting from g_24, new e_it. Replicate 1 is the pilot
# window stored in Y0 and panel. RMSEs over replicates compare designs.
# Feasibility (Task 9, seeds 1 to 5, a crude design: the best of 3000
# random 5-metro sets, with simplex weights on the chosen metros and on the
# rest each fitting the 40-metro average history): RMSE of the synthetic
# design 0.005 to 0.009, against 0.06 to 0.085 for random assignment of 5
# metros with a difference-in-differences estimator (0.21 to 0.23 with a
# pilot-period difference in means).
# =====================================================================
truth_S4 <- local({
  set.seed(seeds22[["S4"]])

  N <- 40; T_pre <- 24; T_pilot <- 6; Tt <- T_pre + T_pilot; t_idx <- 1:Tt
  m_pilot <- 5
  R <- 500
  tau <- 0.04
  sigma_e <- 0.02
  rho_g <- 0.7; sd_g <- 0.05

  alpha <- runif(N, -0.5, 1.0)
  mu1   <- runif(N, 0.4, 0.8)
  mu2   <- runif(N, -0.15, 0.15)
  mu3   <- runif(N, 0, 1)
  lambda_t <- 0.3 * (1 - exp(-t_idx / 18)) + 0.005 * t_idx
  f1 <- 2.0 * (1 - exp(-t_idx / 12)) + 0.05 * t_idx
  f1 <- f1 - mean(f1)
  f2 <- sin(2 * pi * t_idx / 12)

  g <- numeric(Tt)
  g[1] <- rnorm(1, 0, sd_g / sqrt(1 - rho_g^2))
  for (t in 2:T_pre) g[t] <- rho_g * g[t - 1] + rnorm(1, 0, sd_g)
  base <- outer(alpha, rep(1, Tt)) + outer(rep(1, N), lambda_t) + outer(mu1, f1) + outer(mu2, f2)
  Y0_pre <- base[, 1:T_pre] + outer(mu3, g[1:T_pre]) + matrix(rnorm(N * T_pre, 0, sigma_e), N, T_pre)

  pilot <- T_pre + 1:T_pilot
  g_reps <- matrix(NA_real_, R, T_pilot)
  Y0_pilot_reps <- array(NA_real_, c(R, N, T_pilot))
  for (r in 1:R) {
    prev <- g[T_pre]
    for (k in 1:T_pilot) {
      prev <- rho_g * prev + rnorm(1, 0, sd_g)
      g_reps[r, k] <- prev
    }
    Y0_pilot_reps[r, , ] <- base[, pilot] + outer(mu3, g_reps[r, ]) +
      matrix(rnorm(N * T_pilot, 0, sigma_e), N, T_pilot)
  }
  g[pilot] <- g_reps[1, ]

  unit_names <- sprintf("Metro-%02d", 1:N)
  Y0 <- cbind(Y0_pre, Y0_pilot_reps[1, , ])
  dimnames(Y0) <- list(unit_names, NULL)
  dimnames(Y0_pilot_reps) <- list(NULL, unit_names, NULL)
  panel <- data.frame(
    metro  = factor(rep(unit_names, times = Tt), levels = unit_names),
    unit   = rep(1:N, times = Tt),
    month  = rep(t_idx, each = N),
    log_gmv0 = as.vector(Y0),                     # untreated outcome Y_it(0)
    phase  = rep(ifelse(t_idx <= T_pre, "history", "pilot"), each = N)
  )
  meta <- list(
    N = N, T_pre = T_pre, T_pilot = T_pilot, Tt = Tt, m_pilot = m_pilot, R = R, tau = tau,
    alpha = alpha, loading1 = mu1, loading2 = mu2, loading3 = mu3,
    lambda_t = lambda_t, f1 = f1, f2 = f2, g = g, g_pilot_reps = g_reps,
    rho_g = rho_g, sd_g = sd_g, sigma_e = sigma_e,
    observed_rule = "pilot-window outcome = Y0 + tau * D_i, D_i = 1 for the m_pilot metros the design picks"
  )
  saveRDS(list(panel = panel, Y0 = Y0, Y0_pilot_reps = Y0_pilot_reps, meta = meta, unit_names = unit_names),
          file.path(out_dir, "aster_S4.rds"))

  cat(sprintf("\nS4 DGP: %d metros, %d history months + %d pilot months, pick %d pilot metros; %d replicate pilot windows\n",
              N, T_pre, T_pilot, m_pilot, R))
  cat(sprintf("  true pilot effect tau           = %.3f\n", tau))

  list(tau = tau)
})

# =====================================================================
# Truth file
# =====================================================================
truth22 <- list(S1 = truth_S1, S2 = truth_S2, S3 = truth_S3, S4 = truth_S4,
                seeds = seeds22)
saveRDS(truth22, file.path(out_dir, "truth22.rds"))
cat(sprintf("\ntruth22.rds: S1 att_post %.4f | S2 att_post %.4f | S3 att %.4f (%d ATT(g,t) cells) | S4 tau %.3f\n",
            truth22$S1$att_post, truth22$S2$att_post, truth22$S3$att, nrow(truth22$S3$att_gt), truth22$S4$tau))
