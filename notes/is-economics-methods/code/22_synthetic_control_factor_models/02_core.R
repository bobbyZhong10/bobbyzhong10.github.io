# =====================================================================
# 22_synthetic_control_factor_models / 02_core.R
# The inherited chapter 22 analysis of block S1 (the single-metro fee
# reform, aster_S1.rds), as the legacy 04_sc_estimate.R ran it with
# tidysynth 0.2.1 (Ruling S2: this script replaces it):
#   - the naive 2x2 DiD (treated metro vs the simple donor average);
#   - donor weights (convex, sum to 1) from a predictor-based SC fit
#     with tidysynth's own V optimisation (predictors: the pre-period
#     mean and log GMV in months 6, 12, 18, 24);
#   - pre-period fit (RMSPE), treated-vs-synthetic gap = effect path;
#   - in-space placebos and the permutation p-value from the post/pre
#     MSPE ratio (tidysynth's grab_significance, which counts the
#     first reform month 25 as pre: time_unit <= i_time);
#   - in-time placebo (fake reform 8 months early).
# The tidysynth code is the legacy code line for line: its conventions
# are what give the inherited numbers (MSPE ratio 32.777, pre-RMSPE
# 0.0165), which tests/test_reproduce.R locks. The hand-written SC of
# sc_fns.R (V = I on outcomes) is used from 03_s1_modern.R on.
# Writes res_sc_core.rds. Run from this directory:  Rscript 02_core.R
# =====================================================================

suppressWarnings(suppressMessages({
  library(tidysynth)
  library(dplyr)
}))

in_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(in_dir) == 0 || in_dir == "") in_dir <- "."
obj   <- readRDS(file.path(in_dir, "aster_S1.rds"))
panel <- obj$panel
meta  <- obj$meta
T0    <- meta$T0          # legacy convention: the FIRST reform month (25)
Y     <- obj$Y

# ---- naive 2x2 DiD, the legacy algorithm (as in 01_dgp.R) ------------
post    <- seq_len(meta$Tt) >= T0
tr_pre  <- mean(Y[1, !post]);            tr_post <- mean(Y[1, post])
do_pre  <- mean(colMeans(Y[-1, !post])); do_post <- mean(colMeans(Y[-1, post]))
naive_did <- (tr_post - tr_pre) - (do_post - do_pre)
stopifnot(abs(naive_did - meta$did_est) < 1e-12)      # two codings of one formula

# ---- build the synthetic control (legacy code) ---------------------
sc <- panel %>%
  synthetic_control(outcome = log_gmv, unit = metro, time = month,
                    i_unit = "Aster-Metro", i_time = T0,
                    generate_placebos = TRUE) %>%
  generate_predictor(time_window = 1:(T0 - 1), mean_pre = mean(log_gmv)) %>%
  generate_predictor(time_window = 6,          gmv_m06  = log_gmv) %>%
  generate_predictor(time_window = 12,         gmv_m12  = log_gmv) %>%
  generate_predictor(time_window = 18,         gmv_m18  = log_gmv) %>%
  generate_predictor(time_window = 24,         gmv_m24  = log_gmv) %>%
  generate_weights(optimization_window = 1:(T0 - 1)) %>%
  generate_control()

# ---- donor weights -------------------------------------------------
wts <- sc %>% grab_unit_weights() %>% arrange(desc(weight))
top_w <- head(wts, 6)
w_all <- setNames(wts$weight, wts$unit)
w_all <- w_all[sprintf("Donor-%02d", 1:meta$n_donors)]     # all 20 donors, in donor order

# predictor weights (tidysynth's optimised V, the chapter's Omega)
v_pred <- sc %>% grab_predictor_weights()
v_pred <- setNames(v_pred$weight, v_pred$variable)

# ---- observed vs synthetic path & the gap --------------------------
traj <- sc %>% grab_synthetic_control(placebo = FALSE)
# columns: time_unit, real_y, synth_y
traj$gap <- traj$real_y - traj$synth_y

# full placebo trajectories (all units) for the figure's gap panel
traj_all <- sc %>% grab_synthetic_control(placebo = TRUE)

pre_idx  <- traj$time_unit <  T0
post_idx <- traj$time_unit >= T0
pre_rmspe  <- sqrt(mean(traj$gap[pre_idx]^2))
post_avg   <- mean(traj$gap[post_idx])
post_last  <- traj$gap[traj$time_unit == max(traj$time_unit)]

# ---- permutation inference (in-space placebos) ---------------------
sig <- sc %>% grab_significance()          # per-unit post/pre MSPE ratio + p
treated_row <- sig %>% filter(unit_name == "Aster-Metro")
perm_p    <- treated_row$fishers_exact_pvalue
rmspe_rat <- treated_row$mspe_ratio
n_perm    <- nrow(sig)
rank_perm <- treated_row$rank
# how many placebos the figure keeps (pre-RMSPE <= 5 x the treated unit's, 07_figures.R)
pre_fit_all <- traj_all %>% filter(time_unit < T0) %>%
  group_by(.id) %>% summarise(r = sqrt(mean((real_y - synth_y)^2)), .groups = "drop")
n_fig_kept <- sum(pre_fit_all$r <= 5 * pre_fit_all$r[pre_fit_all$.id == "Aster-Metro"]) - 1L

# ---- in-time placebo: fake reform 8 months before the real one -----
T_fake <- T0 - 8
sc_time <- panel %>%
  synthetic_control(outcome = log_gmv, unit = metro, time = month,
                    i_unit = "Aster-Metro", i_time = T_fake,
                    generate_placebos = FALSE) %>%
  generate_predictor(time_window = 1:(T_fake - 1), mean_pre = mean(log_gmv)) %>%
  generate_predictor(time_window = 6,  gmv_m06 = log_gmv) %>%
  generate_predictor(time_window = 12, gmv_m12 = log_gmv) %>%
  generate_weights(optimization_window = 1:(T_fake - 1)) %>%
  generate_control()
traj_t <- sc_time %>% grab_synthetic_control(placebo = FALSE)
traj_t$gap <- traj_t$real_y - traj_t$synth_y
# "effect" in the window between the fake date and the real date (still
# pre-reform, so a valid SC should show ~0 here)
intime_idx <- traj_t$time_unit >= T_fake & traj_t$time_unit < T0
intime_avg <- mean(traj_t$gap[intime_idx])

res <- list(
  # the fields tests/test_reproduce.R reads
  naive_did     = naive_did,
  att_post      = post_avg,          # SC estimate: average post-period gap (months 25-36)
  gap_last      = post_last,         # SC gap in month 36
  pre_rmspe     = pre_rmspe,         # months 1-24
  mspe_ratio    = rmspe_rat,         # tidysynth: months 1-25 count as pre
  p_perm        = perm_p,
  intime_gap    = intime_avg,
  w             = w_all,             # SC weights of all 20 donors, named
  # the rest of the inherited numbers and the figure inputs
  rank_perm     = rank_perm,
  n_perm        = n_perm,
  n_fig_placebos_kept = n_fig_kept,
  v_pred        = v_pred,
  T_fake        = T_fake,
  true_post_avg = meta$true_post_avg,
  true_last     = meta$true_effect_path[meta$Tt],
  did_bias      = naive_did - meta$true_post_avg,
  did_ratio     = naive_did / meta$true_post_avg,
  tr_pre_slope  = meta$tr_pre_slope,
  do_pre_slope  = meta$do_pre_slope,
  w_true        = meta$w_true,
  anchors       = meta$anchors,
  sigma_e       = meta$sigma_e,
  n_donors      = meta$n_donors,
  Tt            = meta$Tt,
  top_weights   = top_w,
  sig           = sig,
  traj          = traj,
  traj_all      = traj_all,
  T0            = T0
)
saveRDS(res, file.path(in_dir, "res_sc_core.rds"))

# ---- report --------------------------------------------------------
cat("================ Synthetic control (S1, inherited) ================\n")
cat(sprintf("true post-period avg effect   : %.3f\n", res$true_post_avg))
cat(sprintf("naive DiD (fails)             : %.3f  (bias %.3f, %.1f x the truth)\n",
            res$naive_did, res$did_bias, res$did_ratio))
cat(sprintf("pre-trend slope treated/donors: %.4f / %.4f\n", res$tr_pre_slope, res$do_pre_slope))
cat(sprintf("SC post-period avg gap        : %.3f  (<- estimate)\n", res$att_post))
cat(sprintf("SC gap at last month          : %.3f  (true %.3f)\n",
            res$gap_last, res$true_last))
cat(sprintf("pre-period RMSPE              : %.4f  (tight fit)\n", res$pre_rmspe))
cat(sprintf("in-time placebo avg gap       : %.4f  (should be ~0)\n", res$intime_gap))
cat(sprintf("permutation p-value           : %.3f  (rank %d of %d units)\n",
            res$p_perm, res$rank_perm, res$n_perm))
cat(sprintf("post/pre MSPE ratio (treated) : %.3f\n", res$mspe_ratio))
cat("---- top donor weights (true anchors: donors 3,7,12,18 = .35/.30/.20/.15) ----\n")
print(top_w)
