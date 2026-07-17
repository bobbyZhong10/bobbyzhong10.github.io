# =====================================================================
# 14_rdd_synthetic_control / 04_sc_estimate.R
# Synthetic control for the Aster fee reform via tidysynth:
#   - donor weights (convex, sum to 1);
#   - pre-period fit (RMSPE);
#   - treated-vs-synthetic gap = the estimated effect path;
#   - in-space placebo (run SC pretending each donor was treated);
#   - in-time placebo (fake reform date in the pre-period);
#   - permutation p-value from the post/pre RMSPE-ratio distribution.
# =====================================================================

suppressWarnings(suppressMessages({
  library(tidysynth)
  library(dplyr)
}))

in_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(in_dir) == 0 || in_dir == "") in_dir <- "."
obj   <- readRDS(file.path(in_dir, "data_sc.rds"))
panel <- obj$panel
meta  <- obj$meta
T0    <- meta$T0

# ---- build the synthetic control -----------------------------------
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
sig <- sc %>% grab_significance()          # per-unit pre/post RMSPE ratio + p
treated_row <- sig %>% filter(unit_name == "Aster-Metro")
perm_p    <- treated_row$fishers_exact_pvalue
rmspe_rat <- treated_row$mspe_ratio
n_perm    <- nrow(sig)

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
  true_post_avg = meta$true_post_avg,
  true_last     = meta$true_effect_path[meta$Tt],
  did_est       = meta$did_est,
  sc_post_avg   = post_avg,
  sc_post_last  = post_last,
  pre_rmspe     = pre_rmspe,
  perm_p        = perm_p,
  rmspe_ratio   = rmspe_rat,
  n_perm        = n_perm,
  intime_avg    = intime_avg,
  top_weights   = top_w,
  w_true        = meta$w_true,
  anchors       = meta$anchors,
  traj          = traj,
  traj_all      = traj_all,
  T0            = T0
)
saveRDS(res, file.path(in_dir, "res_sc.rds"))

# ---- report --------------------------------------------------------
cat("================ Synthetic control ================\n")
cat(sprintf("true post-period avg effect   : %.3f\n", res$true_post_avg))
cat(sprintf("naive DiD (fails)             : %.3f\n", res$did_est))
cat(sprintf("SC post-period avg gap        : %.3f  (<- estimate)\n", res$sc_post_avg))
cat(sprintf("SC gap at last month          : %.3f  (true %.3f)\n",
            res$sc_post_last, res$true_last))
cat(sprintf("pre-period RMSPE              : %.4f  (tight fit)\n", res$pre_rmspe))
cat(sprintf("in-time placebo avg gap       : %.4f  (should be ~0)\n", res$intime_avg))
cat(sprintf("permutation p-value           : %.3f  (of %d units)\n",
            res$perm_p, res$n_perm))
cat(sprintf("post/pre RMSPE ratio (treated): %.2f\n", res$rmspe_ratio))
cat("---- top donor weights (true anchors: donors 3,7,12,18 = .35/.30/.20/.15) ----\n")
print(top_w)
