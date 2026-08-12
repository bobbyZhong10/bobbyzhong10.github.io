# =====================================================================
# 14_rdd_synthetic_control / 01_dgp_rdd.R
# Aster marketplace -- RDD design line (Preferred badge).
#
# Sellers have a continuous quality score X (the running variable),
# cutoff c = 80. Sellers with score >= 80 get a "Preferred" badge.
# We simulate the effect of the badge on log monthly GMV.
#
# Y(0) is a SMOOTH NONLINEAR function of the score (higher-scoring
# sellers sell more anyway). The badge adds a JUMP at the cutoff.
#   - sharp variant:  everyone above 80 gets the badge (D = 1[X>=80]).
#   - fuzzy variant:  badge probability JUMPS at 80 but not 0 -> 1
#                     (exceptions below, non-adopters above).
#
# The DGP is tuned so the NAIVE above-vs-below difference is MUCH
# larger than the true cutoff jump, because Y(0) rises steeply in X.
# =====================================================================

set.seed(101)

N <- 12000
c_cut <- 80

# ---- running variable: continuous quality score in [45,100] --------
# Drawn from a smooth truncated-normal so the density is continuous at
# the cutoff (no manipulation) -> McCrary/rddensity should NOT reject.
raw   <- rnorm(N, mean = 74, sd = 11)
score <- pmin(pmax(raw, 45), 100)

# ---- smooth nonlinear baseline potential outcome Y(0) --------------
# centered, scaled running variable
s <- (score - c_cut) / 10

# mu0 is increasing and genuinely WIGGLY (continuous at the cutoff).
# The steep linear slope makes the naive above-vs-below gap huge (high
# scorers sell far more anyway). The curvature/wiggle is smooth AT the
# cutoff, so a local-linear fit in a narrow window stays ~unbiased,
# but a high-order GLOBAL polynomial has to bend across the whole
# support and its cutoff extrapolation swings with the order chosen --
# exactly the Gelman & Imbens (2019) overfitting warning.
mu0 <- 5.60 + 0.95 * s + 0.20 * s^2 + 0.50 * sin(4 * s)

# ---- true cutoff effect (the sharp RDD estimand) -------------------
tau_true <- 0.20              # true jump in log GMV at c = 80
sigma_e  <- 0.30             # idiosyncratic noise (log GMV)

eps <- rnorm(N, 0, sigma_e)

# ---- SHARP variant -------------------------------------------------
D_sharp <- as.integer(score >= c_cut)
y_sharp <- mu0 + tau_true * D_sharp + eps

# ---- FUZZY variant -------------------------------------------------
# Badge probability jumps at the cutoff but not from 0 to 1:
#   below 80: p = 0.10 (exceptions granted)
#   at/above 80: p = 0.85 (not every eligible seller activates)
# => first-stage compliance jump ~ 0.75. Structural badge effect is
# the SAME tau_true, so the fuzzy-RD LATE should recover ~ tau_true.
p_badge <- ifelse(score >= c_cut, 0.85, 0.10)
D_fuzzy <- rbinom(N, 1, p_badge)
eps2    <- rnorm(N, 0, sigma_e)
y_fuzzy <- mu0 + tau_true * D_fuzzy + eps2

first_stage_jump_true <- 0.85 - 0.10   # = 0.75

dat <- data.frame(
  score   = score,
  s       = s,
  mu0     = mu0,
  D_sharp = D_sharp,
  y_sharp = y_sharp,
  D_fuzzy = D_fuzzy,
  y_fuzzy = y_fuzzy
)

# ---- naive comparison (the "trap") ---------------------------------
naive_diff <- mean(dat$y_sharp[dat$D_sharp == 1]) -
              mean(dat$y_sharp[dat$D_sharp == 0])

meta <- list(
  N = N, c_cut = c_cut, tau_true = tau_true, sigma_e = sigma_e,
  first_stage_jump_true = first_stage_jump_true,
  naive_diff = naive_diff,
  n_above = sum(dat$D_sharp == 1), n_below = sum(dat$D_sharp == 0),
  frac_badge_below = mean(dat$D_fuzzy[dat$score < c_cut]),
  frac_badge_above = mean(dat$D_fuzzy[dat$score >= c_cut])
)

out_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(out_dir) == 0 || out_dir == "") out_dir <- "."
saveRDS(list(dat = dat, meta = meta), file.path(out_dir, "data_rdd.rds"))

cat(sprintf("RDD DGP: N=%d, cutoff=%d, tau_true=%.3f\n", N, c_cut, tau_true))
cat(sprintf("  naive above-vs-below diff = %.3f (should be >> tau_true)\n", naive_diff))
cat(sprintf("  n_below=%d, n_above=%d\n", meta$n_below, meta$n_above))
cat(sprintf("  fuzzy badge frac: below=%.3f, above=%.3f (jump=%.3f)\n",
            meta$frac_badge_below, meta$frac_badge_above,
            meta$frac_badge_above - meta$frac_badge_below))
