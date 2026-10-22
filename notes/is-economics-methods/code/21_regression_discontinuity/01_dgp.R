# =====================================================================
# 21_regression_discontinuity / 01_dgp.R
# Aster marketplace: every simulated data block of chapter 21.
# Replaces 01_dgp_rdd.R. Run from this directory:  Rscript 01_dgp.R  (about a second)
#
# Block R1 (inherited, FIRST, verbatim from 01_dgp_rdd.R, set.seed(101)): the Preferred
#   badge at quality score 80; sharp and fuzzy variants; true jump 0.20. Its random stream
#   is untouched: the code below is the legacy code line for line, only the output file
#   is renamed (data_rdd.rds -> aster_R1.rds, legacy list(dat, meta) layout). Every
#   inherited number of the chapter (1.891, 0.200, [0.133, 0.267], ...) comes from it.
#   Reading note on the legacy comment "smooth truncated-normal" in the score draw below: the
#   code draws a normal and clamps it to [45, 100] with pmin/pmax, so the two ends carry small
#   point masses far from the cutoff and the density is smooth at 80. The comment keeps its
#   legacy wording so that the block stays verbatim; the chapter text explains the comment.
#
# New blocks, each from its generator in rd_fns.R with its own seed (ruling R4), saved as
# a plain list of equal-length column vectors plus meta (ruling R2). y is log monthly GMV
# (USD) unless stated; columns ending in _true are simulation truth, not observable.
#
#   R1c (seed 2101, N 12,000) same world as R1 (score, Y(0) mean, 0.20 jump at 80) plus two
#       predetermined covariates: gmv_prev (log monthly GMV in the year before the badge
#       launched, USD) and tenure (years on Aster). They explain about 37% of the outcome
#       variance left after the score and do not jump at 80.
#       Task 7 needs: with covs, CI length / without <= 0.9 and |est change| < 0.05
#       (covariates buy precision, not a new estimand).
#   R1h (seed 2102, N 12,000) the R1 world a year later, score reported with one decimal; sellers
#       just below 80 push their reported score to exactly 80.0 (about 3.5% of sellers),
#       their outcome still follows the true score. True jump 0.20.
#       Task 7 needs: density test p < 0.05 at 80, donut estimate closer to 0.20 than the
#       naive estimate (the heapers drag the right-hand limit down).
#   R2  (seed 2103, N 8,000) discrete running variable: the average of a seller's 10 most
#       recent review stars (support 3.0, 3.1, ..., 5.0, 21 points, at least 42 sellers
#       each); the "4.5 stars & up" search filter starts at 4.25. True jump 0.15
#       (truth21$R2$jump); the filter gain rises just right of the cutoff (first results
#       page), so local-linear extrapolation from the support points is biased while the
#       two points 4.2 and 4.3 differ by the jump plus about 0.01. meta$M_bound (7.79)
#       bounds |second derivative| on each side. tenure is a predetermined covariate that
#       balances only in the window [4.2, 4.3] (rdwinselect's choice).
#       Task 7 needs: Monte Carlo (>= 500 draws of gen_R2) with honest coverage >= 0.93
#       and rdrobust conventional coverage <= 0.88; honest and local-randomization CIs
#       covering 0.15 in this draw.
#   R3  (seed 2104, N 60,000) regression kink: the commission rate b(x) is 12% up to USD 1,500 of
#       last quarter's average monthly GMV x, then falls 0.6 pp per USD 100 to a 6% floor
#       from USD 2,500 (meta$b_fun, meta$cutoff = 1500). Kink effect -0.05 per percentage
#       point (truth21$R3$kink_effect), no level jump.
#       Task 7 needs: kink CI covering -0.05, level jump at 1,500 within 2 se of zero.
#   R4  (seed 2105, N 24,000 = 12,800 / 8,000 / 3,200) multi-cutoff: badge cutoffs 70
#       (vintage), 80 (home), 90 (electronics), effects 0.40, 0.25, 0.10
#       (truth21$R4$tau_c); parallel Y(0) across categories, so
#       the vintage effect at the interior score 75 is 0.325 (truth21$R4$tau_extrap).
#       Lemma 1 weights about 0.57, 0.31, 0.12, pooled target about 0.32, simple mean 0.25.
#       Task 7 needs: each cutoff CI and the extrapolation CI covering their truths; the
#       pooled estimate more than 0.02 away from mean(tau_c) (0.25).
#   R4g (seed 2106, N 40,000) geographic RD: same-day delivery zone = within 25 km of the hub;
#       coordinates x1, x2 (km east, north); six boundary points with effects 0.30,
#       0.25, 0.15, 0.10, 0.15, 0.25 (truth21$R4g$tau_b).
#       Task 7 needs: boundary-point CIs covering at least 80% of their truths.
#   Power (Task 7): rdpower on R1 for a 0.20 jump (Task 6 prototype: power 0.9995).
#
# Outputs (this directory): aster_R1.rds, aster_R1c.rds, aster_R1h.rds, aster_R2.rds,
# aster_R3.rds, aster_R4.rds, aster_R4g.rds, truth21.rds (true parameters of every block).
# =====================================================================

# ---------------------------------------------------------------------
# Block R1: copied verbatim from 01_dgp_rdd.R (lines 19 to 90), then saved as aster_R1.rds
# ---------------------------------------------------------------------
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
saveRDS(list(dat = dat, meta = meta), file.path(out_dir, "aster_R1.rds"))

cat(sprintf("RDD DGP: N=%d, cutoff=%d, tau_true=%.3f\n", N, c_cut, tau_true))
cat(sprintf("  naive above-vs-below diff = %.3f (should be >> tau_true)\n", naive_diff))
cat(sprintf("  n_below=%d, n_above=%d\n", meta$n_below, meta$n_above))
cat(sprintf("  fuzzy badge frac: below=%.3f, above=%.3f (jump=%.3f)\n",
            meta$frac_badge_below, meta$frac_badge_above,
            meta$frac_badge_above - meta$frac_badge_below))

# ---------------------------------------------------------------------
# New blocks (generators in rd_fns.R; each sets its own seed, so they do not depend on
# the R1 stream above or on each other)
# ---------------------------------------------------------------------
source(file.path(out_dir, "rd_fns.R"))

seeds <- c(R1c = 2101, R1h = 2102, R2 = 2103, R3 = 2104, R4 = 2105, R4g = 2106)

# each gen_*() starts with set.seed(seed), so every block has its own stream
R1c <- gen_R1c(n = 12000, seed = seeds[["R1c"]])
R1h <- gen_R1h(n = 12000, seed = seeds[["R1h"]])
R2  <- gen_R2(n = 8000,   seed = seeds[["R2"]])
R3  <- gen_R3(n = 60000,  seed = seeds[["R3"]])
R4  <- gen_R4(n = 24000,  seed = seeds[["R4"]])
R4g <- gen_R4g(n = 40000, seed = seeds[["R4g"]])

blocks <- list(R1c = R1c, R1h = R1h, R2 = R2, R3 = R3, R4 = R4, R4g = R4g)
for (nm in names(blocks)) saveRDS(blocks[[nm]], file.path(out_dir, paste0("aster_", nm, ".rds")))

truth21 <- list(
  R1  = list(jump = meta$tau_true, first_stage_jump = meta$first_stage_jump_true, cutoff = meta$c_cut),
  R1c = list(jump = R1c$meta$jump, cutoff = R1c$meta$cutoff),
  R1h = list(jump = R1h$meta$jump, cutoff = R1h$meta$cutoff, share_heaped = R1h$meta$share_heaped),
  R2  = list(jump = R2$meta$jump, cutoff = R2$meta$cutoff, M_bound = R2$meta$M_bound, cef = R2$meta$cef),
  R3  = list(kink_effect = R3$meta$kink_effect, level_jump = R3$meta$level_jump, cutoff = R3$meta$cutoff,
             dslope_b = R3$meta$slope_right - R3$meta$slope_left),
  R4  = list(tau_c = R4$meta$tau_c, tau_extrap = R4$meta$tau_extrap, x_extrap = R4$meta$x_extrap,
             extrap_cutoff = R4$meta$extrap_cutoff, ref_cutoff = R4$meta$ref_cutoff,
             omega = R4$meta$omega, tau_pooled = sum(R4$meta$omega * R4$meta$tau_c)),
  R4g = list(tau_b = R4g$meta$tau_b, boundary = R4g$meta$boundary, radius = R4g$meta$radius),
  seeds = c(R1 = 101, seeds)
)
saveRDS(truth21, file.path(out_dir, "truth21.rds"))

f3 <- function(v) paste(sprintf("%.3f", v), collapse = ", ")
cat("\nNew blocks (seeds ", paste(names(seeds), seeds, sep = " = ", collapse = ", "), ")\n", sep = "")
cat(sprintf("  R1c: N=%d, jump %.3f, covariates gmv_prev, tenure\n", R1c$meta$N, truth21$R1c$jump))
cat(sprintf("  R1h: N=%d, jump %.3f, heaped at 80.0: %d sellers (%.1f%%), at 79.0: %d\n", R1h$meta$N,
            truth21$R1h$jump, sum(R1h$x == 80), 100 * R1h$meta$share_heaped, sum(R1h$x == 79)))
cat(sprintf("  R2:  N=%d, %d support points, min mass %d, jump %.3f, M_bound %.3f\n", R2$meta$N,
            length(unique(R2$x)), min(table(R2$x)), truth21$R2$jump, truth21$R2$M_bound))
cat(sprintf("  R3:  N=%d, kink at USD %d, policy slope change %.3f pp/USD, kink effect %.3f\n", R3$meta$N,
            R3$meta$cutoff, truth21$R3$dslope_b, truth21$R3$kink_effect))
cat(sprintf("  R4:  N=%d, cutoffs %s, tau_c %s, Lemma 1 weights %s, pooled %.3f vs mean %.3f, tau(75) %.3f\n",
            R4$meta$N, paste(R4$meta$cutoffs, collapse = "/"), f3(truth21$R4$tau_c), f3(truth21$R4$omega),
            truth21$R4$tau_pooled, mean(truth21$R4$tau_c), truth21$R4$tau_extrap))
cat(sprintf("  R4g: N=%d, radius %d km, inside %d, tau_b %s\n", R4g$meta$N, R4g$meta$radius, sum(R4g$d),
            f3(truth21$R4g$tau_b)))
