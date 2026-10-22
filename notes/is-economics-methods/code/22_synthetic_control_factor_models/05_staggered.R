# =====================================================================
# 22_synthetic_control_factor_models / 05_staggered.R
# Block S3 (aster_S3.rds, seed seeds22[["S3"]]): a seller relief program rolled out
# in three waves of 5 metros (months 29, 33, 37) to the metros the demand
# downturn hit hardest. The cycle loading drives WHO adopts (selection on
# factor exposure; wave order is nearly random among adopters), so the
# adopters' untreated paths are not parallel to the never-treated metros'.
# True effect 0.06 in every treated metro-month.
#
# Estimators (recipes in stag_fns.R; every bootstrap seeded, so the fits
# are deterministic):
#   twfe         fixest two-way FE, clustered by metro
#   cs           Callaway-Sant'Anna, never-treated controls, "simple"
#   cs_nyt       the same with not-yet-treated controls (extra)
#   ife          fect interactive FE, CV over r in 0..4
#   mc           fect matrix completion: CV over lambda (fect defaults),
#                then a refit at the CV lambda with CV = FALSE, tol 1e-7,
#                max.iteration 50000 (Ruling S6a: fect's default fit after
#                CV stops at tol >= 1e-3, Task 9 report 7.1)
#   mc_default   the fect-default MC fit (CV, default tolerance)
#   multisynth   augsynth::multisynth defaults (Ben-Michael, Feller,
#                Rothstein 2022)
#   sdid_cohort  SDID per cohort, aggregated by treated metro-months
#                (AAHIW 2021, Section 8), unit-bootstrap interval
# fect bootstraps use nboots = 500 here (the screen uses 100, as Task 9
# recorded); every fit uses est_seeds22[["s3_fits"]] (est_seeds.R, value 1):
# fect seed, set.seed before CS (biters 999) and before the multisynth
# summary, SDID bootstrap (B = 500).
#
# The 24-seed screen (Ruling S6c) is computed by 05b_stag_screen.R and
# cached in res_stag_screen.rds; when the cache exists this script copies
# it into res_stag.rds$screen (05b does the same after computing it).
#
# READ BEFORE QUOTING THE MC NUMBERS (Task 10 review, Important 1). The
# stored draw (seed seeds22[["S3"]], 21 when this was written) is the most
# favourable of the 24 screened draws for BOTH matrix-completion recipes: its
# MC estimates are the lowest of the 24, so the share of the TWFE bias they
# remove ranks 1st of 24 and lies far above the mean over the draws. MC
# overestimates the true effect in every one of the 24 draws under both
# recipes (it always lies between the truth and TWFE), and its interval covers
# the truth in only some of them. The other estimators' stored estimates are
# not extreme. Quote the values from numbers.json, not from this header (keys
# s3_screen_stored_rank_*, s3_screen_n_est_above_truth_*,
# s3_screen_*_share_twfe_bias_removed_*, s3_screen_covered_*); seed 21 was
# chosen by Task 9's screen (dgp_s3.R).
# Writes res_stag.rds. Run from this directory:  Rscript 05_staggered.R
# =====================================================================

suppressWarnings(suppressMessages({
  library(fixest)
  library(did)
  library(fect)
  library(augsynth)
  library(synthdid)
}))
source("stag_fns.R")
source("est_seeds.R")

S3 <- readRDS("aster_S3.rds"); tr <- readRDS("truth22.rds")
p <- S3$panel
nb <- 500L
cores <- max(1L, min(8L, parallel::detectCores() - 1L))
tm <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  out <- expr
  cat(sprintf("  %-12s %6.1f s\n", label, proc.time()[["elapsed"]] - t0))
  out
}

cat(sprintf("fitting S3 (seed %d):\n", tr$seeds[["S3"]]))
twfe       <- tm("twfe", fit_twfe(p))
sd3 <- est_seeds22[["s3_fits"]]
cs         <- tm("cs", fit_cs(p, "nevertreated", biters = 999, seed = sd3))
cs_nyt     <- tm("cs_nyt", fit_cs(p, "notyettreated", biters = 999, seed = sd3))
ife        <- tm("ife", fit_fect(p, "ife", nboots = nb, seed = sd3))
mc_default <- tm("mc_default", fit_fect(p, "mc", nboots = nb, seed = sd3))
mc         <- tm("mc_refit", fit_mc_refit(p, lambda = mc_default$lambda_cv, nboots = nb, seed = sd3))
multisynth <- tm("multisynth", fit_multisynth(p, seed = sd3))
sdid_cohort <- tm("sdid_cohort", fit_sdid_cohort(p, B = 500, seed = sd3, cores = cores))

m <- S3$meta; ad <- is.finite(m$cohort)
design <- list(n_units = m$N, n_months = m$Tt, n_adopters = sum(ad), n_never = sum(!ad), n_contract = m$n_contract,
               waves = m$waves, n_wave = m$n_wave, n_treated_cells = sum(p$treat_it),
               loading1_adopters = mean(m$loading1[ad]), loading1_never = mean(m$loading1[!ad]),
               loading1_by_cohort = m$loading1_by_cohort,
               cor_loading1_adopt_rank = cor(m$loading1, m$adopt_rank),
               cor_loading1_cohort_adopters = cor(m$loading1[ad], m$cohort[ad]),
               cor_loading1_treated = cor(m$loading1, as.numeric(ad)))

res <- list(twfe = twfe, cs = cs, cs_nyt = cs_nyt, ife = ife, mc = mc, mc_default = mc_default,
            multisynth = multisynth, sdid_cohort = sdid_cohort, truth = tr$S3$att, design = design,
            settings = list(fect_nboots = nb, fect_seed = sd3, cs_biters = 999, cs_seed = sd3, multisynth_seed = sd3,
                            sdid_B = 500, sdid_seed = sd3),
            screen = if (file.exists("res_stag_screen.rds")) readRDS("res_stag_screen.rds") else NULL)
saveRDS(res, "res_stag.rds")

f4 <- function(v) sprintf("%.4f", v)
cat("================ S3: staggered rollout, selection on factor exposure (truth", f4(tr$S3$att), ") ================\n")
for (k in c("twfe", "cs", "cs_nyt", "ife", "mc", "mc_default", "multisynth", "sdid_cohort")) {
  o <- res[[k]]
  cat(sprintf("  %-12s %s  [%s, %s]  covers truth: %s\n", k, f4(o$est), f4(o$ci[1]), f4(o$ci[2]),
              o$ci[1] <= tr$S3$att & tr$S3$att <= o$ci[2]))
}
cat("  IFE CV r =", ife$r_cv, "| MC lambda.cv =", format(mc_default$lambda_cv, digits = 4),
    "(lambda.norm", format(mc_default$lambda_norm, digits = 4), ") | refit iterations", mc$niter, "\n")
cat("  SDID by cohort:", paste(names(sdid_cohort$by_cohort), f4(sdid_cohort$by_cohort), sep = ": ", collapse = " | "), "\n")
cat(if (is.null(res$screen)) "  screen cache res_stag_screen.rds missing: run 05b_stag_screen.R\n"
    else "  screen: copied from res_stag_screen.rds\n")
