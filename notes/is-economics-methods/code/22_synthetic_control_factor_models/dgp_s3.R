# =====================================================================
# 22_synthetic_control_factor_models / dgp_s3.R
# gen_S3(seed): the generator of block S3, sourced by 01_dgp.R (which
# calls gen_S3(seeds22[["S3"]]) and writes aster_S3.rds) and by
# 05b_stag_screen.R, which redraws the block for its seed screen without
# regenerating S1, S2 or S4. Pure function apart from the RNG: it calls
# set.seed(seed) (so it resets the session's random stream), writes no
# file and prints nothing. Returns
#   list(data  = list(panel, Y, meta, unit_names),   # the aster_S3.rds object
#        truth = list(att, att_gt))                  # truth22$S3
#
# Block S3: staggered rollout with selection on factor exposure
# ---------------------------------------------------------------------
# 80 mid-sized metros, 48 months; outcome = log monthly GMV in USD
# millions (alpha_i ~ U(-0.5, 1.0): about 0.6 to 2.7 million USD a month).
# Untreated outcomes:
#   Y_it(0) = alpha_i + lambda_t + mu_i1 f1_t + mu_i2 f2_t + u_it
#   lambda_t = 0.6 (1 - exp(-t / 18)) + 0.01 t     common platform growth
#   f1_t = cos(2 pi (t - 1) / 48)                  demand cycle: boom in
#          month 1, downturn to a trough in month 25, recovery to a new
#          boom by month 48
#   f2_t = sin(2 pi t / 12)                        seasonality
#   mu_i1 ~ U(0, 0.25)  cyclical exposure (meta$loading1)
#   mu_i2 ~ U(0, 0.15)  seasonal amplitude (meta$loading2)
#   u_it = e_it + a_it: e_it ~ N(0, 0.02^2) iid, a_it a stationary AR(1)
#          with rho = 0.8 and sd 0.04 (persistent metro-level shocks)
# Rollout: Aster's seller relief program (fee rebates and promotion
# credits) starts in the recovery, in three waves of 5 metros at months
# 29, 33 and 37 (4, 8 and 12 months after the trough). 20 metros on a
# legacy partner contract are not eligible (drawn at random). Aster ranks
# the other 60 by how hard the downturn hit them, measured with noise:
# priority_i = mu_i1 + N(0, 0.05^2). The top 15 by priority get the
# program, in three waves of 5 in priority order; the other 45, and the
# 20 contract metros, never get it (cohort Inf).
#
# What the loadings drive: SELECTION INTO THE PROGRAM (who adopts), not
# the order of the waves. In the stored draw (seed seeds22[["S3"]] in 01_dgp.R, 21 when this was written) the 15 adopters
# have mean cycle loading 0.208 against 0.113 for the 65 never-treated
# metros, so the metros that fell furthest in the downturn, and rebound
# furthest in the recovery, are the ones that get the program: selection
# on factor exposure (not Ashenfelter's dip, which chapters 19 and 22 keep
# for selection on transitory shocks), and parallel trends between
# adopters and never-treated metros fail. Among the adopters, wave order
# is nearly random: the top-15 cut leaves a narrow band of high loadings
# (sd 0.032) and the priority noise (sd 0.05) swamps the differences
# within it, so the wave means are 0.212, 0.210, 0.203 and, within
# adopters, cor(loading1, cohort) = -0.113. The brief's check
# cor(loading1, adopt_rank) = -0.486 is mostly the treated / never-treated
# contrast (cor(loading1, treated indicator) = 0.487), because the 65
# never-treated metros tie at the last rank. The DiD-family bias in S3
# therefore comes from comparing adopters with never-treated metros
# (selection on factor exposure), not from a timing gradient.
# Effect: constant, tau_it = 0.06 in every treated metro-month, so every
# aggregation (TWFE, Callaway-Sant'Anna "simple", fect ATT.avg, the
# multisynth average) targets the same number, and a miss is bias, not a
# different weighting.
#
# Seed screen behind the stored seed, 21 (Task 9 feasibility run; the seed numbers in this
# record, including "seed 21" below, are the values at the time of writing: if seeds22[["S3"]]
# is retuned, edit them by hand; Ruling S6 moves the
# reproducible screen to 05b_stag_screen.R, whose rerun reproduced every
# tally below and the covering seeds exactly; cached in res_stag_screen.rds). Candidates: seeds 1, 2, ..., 24, data
# = gen_S3(seed) (identical to the exploration generator for every
# candidate). True ATT 0.06; covered = 95% interval contains 0.06.
# Recipes, run on the panel with `unit` as the unit index:
#   TWFE  fixest::feols(log_gmv ~ treat_it | unit + month, cluster = ~unit),
#         confint() at 95%.
#   CS    set.seed(1); did::att_gt("log_gmv", "month", "unit", gname =
#         cohort with Inf -> 0, control_group = "nevertreated", est_method =
#         "reg", bstrap = TRUE, cband = FALSE, biters = 999); aggte(type =
#         "simple", bstrap = TRUE, biters = 999, cband = FALSE);
#         interval = estimate -/+ qnorm(0.975) * se.
#   IFE   fect::fect(log_gmv ~ treat_it, index = c("unit", "month"),
#         method = "ife", force = "two-way", CV = TRUE, r = c(0, 4),
#         se = TRUE, nboots = 100, parallel = FALSE, seed = 1);
#         interval = est.avg CI.lower, CI.upper.
#   MC default  the same call with method = "mc" (CV over lambda, no r).
#   MC refit    fect(..., method = "mc", force = "two-way", CV = FALSE,
#         lambda = <lambda.cv of the MC default run>, se = TRUE, nboots =
#         100, parallel = FALSE, seed = 1, tol = 1e-7, max.iteration = 50000).
#   multisynth  augsynth::multisynth(log_gmv ~ treat_it, unit, month, panel)
#         with defaults; summary(); row Level == "Average", Time NA;
#         lower_bound, upper_bound.
# Tallies over the 24 candidates: TWFE covered 0/24 and CS 0/24; IFE
# 24/24 (CV picked r = 2 in 23, r = 3 at seed 4); multisynth 24/24; MC
# default 6/24 (seeds 1, 7, 12, 21, 23, 24); MC refit 14/24 (seeds 1, 3,
# 6, 7, 9, 10, 12, 13, 16, 18, 19, 21, 23, 24). In 24/24 seeds both MC
# estimates lay strictly between the truth and the TWFE estimate. Seed 21
# is one of the six candidates where every Task 10 S3 conclusion held
# under both MC recipes, and the one whose default MC interval covered
# with the widest margin when 21, 24, 23 and 7 were rerun with fect seeds
# 11, 22, 33 and nboots = 500. MC coverage is therefore a property of
# this draw, not of the design (Ruling S6). Choices that help MC: the
# post-adoption cycle values stay inside the range each treated metro
# shows before adoption (the month-1 boom); the random contract
# exclusions leave high-exposure controls; persistent AR(1) noise widens
# all intervals alike; and levels near zero (USD millions) tighten fect's
# relative convergence tolerance on the factor part (fect 2.4.5 runs its
# CV, and its final MC fit after CV, with tol of at least 1e-3).
# =====================================================================
gen_S3 <- function(seed) {
  set.seed(seed)

  N  <- 80
  Tt <- 48
  t_idx <- 1:Tt
  waves      <- c(29, 33, 37)     # first treated month of each wave (cohort G_i)
  n_wave     <- 5
  n_contract <- 20
  tau     <- 0.06
  sigma_e <- 0.02
  ar_rho  <- 0.8
  ar_sd   <- 0.04
  prio_sd <- 0.05

  lambda_t <- 0.6 * (1 - exp(-t_idx / 18)) + 0.01 * t_idx
  f1 <- cos(2 * pi * (t_idx - 1) / 48)
  f2 <- sin(2 * pi * t_idx / 12)

  alpha <- runif(N, -0.5, 1.0)
  mu1   <- runif(N, 0, 0.25)
  mu2   <- runif(N, 0, 0.15)
  priority <- mu1 + rnorm(N, 0, prio_sd)
  rank_index <- priority
  contract <- rep(FALSE, N)
  excl <- sample.int(N, n_contract)
  contract[excl] <- TRUE
  rank_index[excl] <- -Inf
  ord <- order(rank_index, decreasing = TRUE)
  cohort <- rep(Inf, N)
  for (k in seq_along(waves)) cohort[ord[((k - 1) * n_wave + 1):(k * n_wave)]] <- waves[k]

  Y0 <- outer(alpha, rep(1, Tt)) + outer(rep(1, N), lambda_t) + outer(mu1, f1) + outer(mu2, f2) +
    matrix(rnorm(N * Tt, 0, sigma_e), N, Tt)
  a <- matrix(0, N, Tt)
  a[, 1] <- rnorm(N, 0, ar_sd)
  for (t in 2:Tt) a[, t] <- ar_rho * a[, t - 1] + rnorm(N, 0, ar_sd * sqrt(1 - ar_rho^2))
  Y0 <- Y0 + a

  D <- outer(cohort, t_idx, function(g, t) as.integer(t >= g))
  tau_it <- tau * D
  Y <- Y0 + tau_it

  unit_names <- sprintf("Metro-%02d", 1:N)
  rownames(Y) <- unit_names
  panel <- data.frame(
    metro    = factor(rep(unit_names, times = Tt), levels = unit_names),
    unit     = rep(1:N, times = Tt),
    month    = rep(t_idx, each = N),
    log_gmv  = as.vector(Y),
    cohort   = rep(cohort, times = Tt),          # first treated month; Inf = never treated
    treat_it = as.vector(D),
    contract = rep(contract, times = Tt),
    tau_true = as.vector(tau_it)                 # simulation truth, not observable
  )
  adopt_rank <- rank(cohort, ties.method = "average")

  att_gt <- expand.grid(t = t_idx, g = waves)
  att_gt <- att_gt[att_gt$t >= att_gt$g, c("g", "t")]
  att_gt$e   <- att_gt$t - att_gt$g
  att_gt$att <- mapply(function(g, t) mean(tau_it[cohort == g, t]), att_gt$g, att_gt$t)
  rownames(att_gt) <- NULL
  att <- mean(tau_it[D == 1])

  meta <- list(
    N = N, Tt = Tt, waves = waves, n_wave = n_wave, n_contract = n_contract, tau = tau,
    loading1 = mu1, loading2 = mu2, alpha = alpha,       # by unit, in the order of unit_names
    cohort = cohort, adopt_rank = adopt_rank, contract = contract, priority = priority,
    lambda_t = lambda_t, f1 = f1, f2 = f2,
    sigma_e = sigma_e, ar_rho = ar_rho, ar_sd = ar_sd, prio_sd = prio_sd,
    loading1_by_cohort = tapply(mu1, cohort, mean)
  )
  list(data  = list(panel = panel, Y = Y, meta = meta, unit_names = unit_names),
       truth = list(att = att, att_gt = att_gt))
}
