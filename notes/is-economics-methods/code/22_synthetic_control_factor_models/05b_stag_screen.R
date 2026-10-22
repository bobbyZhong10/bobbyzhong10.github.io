# =====================================================================
# 22_synthetic_control_factor_models / 05b_stag_screen.R
# The S3 seed screen (Ruling S6c): rerun the 24 candidate draws that
# Task 9 screened before it fixed seed 21, with the recipes recorded in
# dgp_s3.R, and tally how often each estimator's 95% interval covers the
# true effect 0.06. Data = gen_S3(seed) for seed = 1, ..., 24.
#
#   twfe        fit_twfe (stag_fns.R)
#   cs          fit_cs, never-treated controls, biters 999, set.seed(1)
#   ife         fit_fect "ife", CV r in 0..4, nboots 100, fect seed 1
#   mc_default  fit_fect "mc", CV over lambda, nboots 100, fect seed 1
#   mc          fit_mc_refit at the mc_default lambda.cv: CV = FALSE,
#               tol 1e-7, max.iteration 50000, nboots 100, fect seed 1
#   multisynth  fit_multisynth, defaults; set.seed(1) before summary()
#               (Task 9's record did not seed this bootstrap; seeding it
#               makes the screen reproducible)
# Also counted: the draws in which each MC estimate lies strictly between
# the truth and the TWFE estimate.
#
# SLOW (several minutes on 8 cores). The result is cached in
# res_stag_screen.rds (kept with the code) and recomputed only when the
# cache is missing or the environment variable RUN_SLOW is "1":
#     RUN_SLOW=1 Rscript 05b_stag_screen.R
# Either way the script then checks the screen against res_stag.rds (the
# stored seed's six point estimates must match within 1e-8) and against
# the installed versions of fect, did, fixest and augsynth, and stops
# with a RUN_SLOW=1 hint when they disagree (a stale cache); it then
# copies the screen into res_stag.rds$screen (when res_stag.rds exists)
# and prints the tallies next to Task 9's. Seeds: est_seeds22[["s3_screen"]]
# (est_seeds.R; value 1, Task 9's recorded seed).
# Run from this directory:  Rscript 05b_stag_screen.R
# =====================================================================

cache <- "res_stag_screen.rds"
recompute <- !file.exists(cache) || identical(Sys.getenv("RUN_SLOW"), "1")

# Task 9's record (dgp_s3.R header), for comparison only
task9 <- list(covered = list(twfe = 0, cs = 0, ife = 24, mc_default = 6, mc = 14, multisynth = 24),
              mc_default_seeds = c(1, 7, 12, 21, 23, 24),
              mc_seeds = c(1, 3, 6, 7, 9, 10, 12, 13, 16, 18, 19, 21, 23, 24),
              mc_between = list(mc_default = 24, mc = 24))

if (recompute) {
  suppressWarnings(suppressMessages({
    library(fixest)
    library(did)
    library(fect)
    library(augsynth)
  }))
  source("stag_fns.R")
  source("est_seeds.R")
  source("dgp_s3.R")
  seeds <- 1:24
  truth <- 0.06
  nb <- 100L
  cores <- max(1L, min(8L, parallel::detectCores() - 1L))
  one_seed <- function(s) {
    g <- gen_S3(s)
    p <- g$data$panel
    stopifnot(abs(g$truth$att - truth) < 1e-12)
    tw  <- fit_twfe(p)
    sds <- est_seeds22[["s3_screen"]]
    cs  <- fit_cs(p, "nevertreated", biters = 999, seed = sds)
    ife <- fit_fect(p, "ife", nboots = nb, seed = sds)
    mcd <- fit_fect(p, "mc", nboots = nb, seed = sds)
    mcr <- fit_mc_refit(p, lambda = mcd$lambda_cv, nboots = nb, seed = sds)
    ms  <- fit_multisynth(p, seed = sds)
    fits <- list(twfe = tw, cs = cs, ife = ife, mc_default = mcd, mc = mcr, multisynth = ms)
    row <- data.frame(seed = s)
    for (k in names(fits)) {
      row[[paste0(k, "_est")]] <- fits[[k]]$est
      row[[paste0(k, "_lo")]]  <- fits[[k]]$ci[1]
      row[[paste0(k, "_hi")]]  <- fits[[k]]$ci[2]
      row[[paste0(k, "_cover")]] <- fits[[k]]$ci[1] <= truth & truth <= fits[[k]]$ci[2]
    }
    row$ife_r_cv <- ife$r_cv
    row$mc_lambda_cv <- mcd$lambda_cv
    row$mc_lambda_norm <- mcd$lambda_norm
    row
  }
  t0 <- proc.time()[["elapsed"]]
  rows <- parallel::mclapply(seeds, one_seed, mc.cores = cores)
  bad <- vapply(rows, inherits, logical(1), "try-error")
  if (any(bad)) stop("05b_stag_screen.R: seeds ", paste(seeds[bad], collapse = ", "), " failed: ",
                     paste(unique(vapply(rows[bad], as.character, character(1))), collapse = " | "), call. = FALSE)
  per_seed <- do.call(rbind, rows)
  wall <- proc.time()[["elapsed"]] - t0
  est_keys <- c("twfe", "cs", "ife", "mc_default", "mc", "multisynth")
  between <- function(k) {
    lo <- pmin(truth, per_seed$twfe_est); hi <- pmax(truth, per_seed$twfe_est)
    per_seed[[paste0(k, "_est")]] > lo & per_seed[[paste0(k, "_est")]] < hi
  }
  screen <- list(
    seeds = seeds, n_seeds = length(seeds), truth = truth, per_seed = per_seed,
    covered = setNames(lapply(est_keys, function(k) sum(per_seed[[paste0(k, "_cover")]])), est_keys),
    covered_seeds = setNames(lapply(est_keys, function(k) seeds[per_seed[[paste0(k, "_cover")]]]), est_keys),
    mc_between = list(mc_default = sum(between("mc_default")), mc = sum(between("mc"))),
    ife_r_cv_counts = table(per_seed$ife_r_cv),
    settings = list(fect_nboots = nb, fect_seed = est_seeds22[["s3_screen"]], cs_biters = 999,
                    cs_seed = est_seeds22[["s3_screen"]], multisynth_seed = est_seeds22[["s3_screen"]],
                    mc_refit_tol = 1e-7, mc_refit_max_iter = 50000, cores = cores),
    wall_seconds = wall, computed = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    versions = vapply(c("fect", "did", "fixest", "augsynth"), function(k) as.character(utils::packageVersion(k)), ""),
    task9_record = task9)
  saveRDS(screen, cache)
  cat(sprintf("screen computed: %d seeds in %.0f s on %d cores; cached in %s\n", length(seeds), wall, cores, cache))
} else {
  screen <- readRDS(cache)
  cat(sprintf("screen read from the cache %s (computed %s in %.0f s; RUN_SLOW=1 recomputes it)\n",
              cache, screen$computed, screen$wall_seconds))
}

if (file.exists("res_stag.rds")) {
  st <- readRDS("res_stag.rds")
  # staleness guard (stag_fns.R): the cached screen must reproduce the stored draw's fits in
  # res_stag.rds and must have been computed with the installed package versions
  if (!exists("screen_staleness")) source("stag_fns.R")
  probs <- screen_staleness(screen, st, readRDS("truth22.rds")$seeds[["S3"]])
  if (length(probs))
    stop("05b_stag_screen.R: the screen in ", cache, " does not match res_stag.rds or the installed packages:\n  ",
         paste(probs, collapse = "\n  "),
         if (recompute) "\nRerun 05_staggered.R, then this script." else
           "\nThe cache is stale: recompute it with  RUN_SLOW=1 Rscript 05b_stag_screen.R  (or RUN_SLOW=1 Rscript 00_run_all.R)",
         call. = FALSE)
  st$screen <- screen
  saveRDS(st, "res_stag.rds")
  cat("cache checked against res_stag.rds and the installed packages; copied into res_stag.rds$screen\n")
}

cat(sprintf("%-11s %8s %8s\n", "estimator", "covered", "Task 9"))
for (k in names(screen$covered))
  cat(sprintf("%-11s %5d/24 %5d/24\n", k, screen$covered[[k]], screen$task9_record$covered[[k]]))
cat("MC default covers at seeds:", paste(screen$covered_seeds$mc_default, collapse = ", "), "\n")
cat("MC refit covers at seeds:  ", paste(screen$covered_seeds$mc, collapse = ", "), "\n")
cat(sprintf("MC strictly between truth and TWFE: default %d/24, refit %d/24 (Task 9: %d, %d)\n",
            screen$mc_between$mc_default, screen$mc_between$mc,
            screen$task9_record$mc_between$mc_default, screen$task9_record$mc_between$mc))
