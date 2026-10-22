# tests/test_conclusions.R
# The teaching conclusions of chapter 22 that rest on the four data blocks, checked against
# the res_*.rds files written by 03_s1_modern.R (S1), 04_hull.R (S2), 05_staggered.R and
# 05b_stag_screen.R (S3) and 06_design.R (S4). Run from the code directory after them:
#     Rscript tests/test_conclusions.R
#
# Part 1 is the Task 10 brief's test, code verbatim, comments translated to English (the
# English-comments rule; the direction of a check never changes with the results). Part 2
# adds the Ruling S6 check on matrix completion, the Ruling S8 lock on the moving-block
# conformal p-value, the fix-round-1 checks (MC refit overestimates in every screened draw;
# the cached screen is not stale) and guards that keep Part 1 from passing vacuously (the
# interval objects have the promised shape, the two pre-period RMSPEs come from the same
# recipe, the S3 seed screen is stored). Part 3 and Part 4 lock seed-dependent claims that
# the chapter text states (final review, Task 14). Part 5 checks the staggered scpi run of
# 05c_stag_scpi.R (Ruling S12) and locks the coverage count the text states.
source("../../tools/R/chapter_kit.R")

## ---- Part 1: the brief ------------------------------------------------------------
tr <- readRDS("truth22.rds"); s1 <- readRDS("res_s1_ext.rds"); s2 <- readRDS("res_hull.rds")
s3 <- readRDS("res_stag.rds"); s4 <- readRDS("res_design.rds")
stopifnot(abs(s1$did_uniform - s1$naive_did) < 1e-8)                     # SDID with uniform weights is DiD (numerical evidence for ledger item 19->22)
stopifnot(covers(s1$sdid$ci, tr$S1$att_post), s1$conformal$p <= 0.10)
stopifnot(mean(s1$scpi$lo <= tr$S1$tau_path & tr$S1$tau_path <= s1$scpi$hi) >= 10/12)
stopifnot(s2$sc$pre_rmspe > 3 * s1$sc_pre_rmspe)                         # outside the hull: the SC fit is clearly worse
err <- function(e) abs(e - tr$S2$att_post)
stopifnot(err(s2$asc$est) < err(s2$sc$est), err(s2$sdid$est) < err(s2$sc$est), any(s2$asc$w < 0))
stopifnot(all(s2$pen$w >= -1e-10), abs(sum(s2$pen$w) - 1) < 1e-8, s2$pen$n_nonzero <= s2$sc$n_nonzero,
          s2$pen$max_diff_across_starts < 1e-6)                          # penalized SC: on the simplex, sparser, unique
stopifnot(!covers(s3$twfe$ci, tr$S3$att), !covers(s3$cs$ci, tr$S3$att))  # DiD family fails when factor exposure drives who adopts
stopifnot(covers(s3$ife$ci, tr$S3$att), covers(s3$mc$ci, tr$S3$att), covers(s3$multisynth$ci, tr$S3$att))
stopifnot(s4$rmse_synthetic_design < s4$rmse_random_design)

## ---- Part 2: Ruling S6 and guards ---------------------------------------------------
# Ruling S6(b): the reported matrix completion (CV, then refit at the CV lambda to tight
# convergence) removes part of the TWFE bias. Coverage of the stored draw is checked above;
# across draws it is not a property of the design (s3$screen holds the 24-seed tallies).
stopifnot(abs(s3$mc$est - tr$S3$att) < abs(s3$twfe$est - tr$S3$att))
# every interval is a finite c(lo, hi) with lo < hi, and every estimate a finite number
ci_ok <- function(o) is.list(o) && is.numeric(o$est) && length(o$est) == 1L && is.finite(o$est) &&
  is.numeric(o$ci) && length(o$ci) == 2L && all(is.finite(o$ci)) && o$ci[1] < o$ci[2]
stopifnot(ci_ok(s1$sdid), all(vapply(s3[c("twfe", "cs", "ife", "mc", "mc_default", "multisynth", "sdid_cohort")], ci_ok, logical(1))))
# scpi: one prediction interval per post month, aligned with the 12-month truth (Ruling S1)
stopifnot(length(tr$S1$tau_path) == 12L, length(s1$scpi$lo) == 12L, length(s1$scpi$hi) == 12L, all(s1$scpi$lo < s1$scpi$hi))
# the reported conformal p-value uses iid permutations (identity plus n_perm random ones, so
# it lives on the grid 1/(n_perm + 1), ..., 1); the moving-block p-value is stored next to
# it, on the grid of the T = 36 cyclic shifts (1/36, ..., 1)
np1 <- s1$conformal$n_perm + 1
stopifnot(identical(s1$conformal$scheme, "iid"), s1$conformal$n_perm >= 999,
          s1$conformal$p >= 1/np1, abs(s1$conformal$p * np1 - round(s1$conformal$p * np1)) < 1e-6,
          s1$conformal$p_moving_block >= 1/36,
          abs(s1$conformal$p_moving_block * 36 - round(s1$conformal$p_moving_block * 36)) < 1e-9)
# Ruling S8: the conformal conclusion depends on the permutation scheme. With iid permutations
# (valid under CWZ Assumption 2.1, which S1 meets) the test rejects at 10% (Part 1); with the 36
# moving blocks it does not. Both results are teaching conclusions and both are locked.
stopifnot(s1$conformal$p_moving_block > 0.10)
# the S1 and S2 pre-period RMSPEs compared above use the same recipe: sc_weights on the 24
# pre-period outcomes with V = I (recomputed here from the data blocks)
source("sc_fns.R")
pre_rmspe_v1 <- function(Y) {
  X1 <- Y[1, 1:24]; X0 <- t(Y[-1, 1:24]); w <- sc_weights(X1, X0, diag(24))
  sqrt(mean((X1 - drop(X0 %*% w))^2))
}
stopifnot(abs(s1$sc_pre_rmspe - pre_rmspe_v1(readRDS("aster_S1.rds")$Y)) < 1e-10,
          abs(s2$sc$pre_rmspe - pre_rmspe_v1(readRDS("aster_S2.rds")$Y)) < 1e-10)
# penalized SC: weights named by donor, a stated penalty rule, at least 10 donor orderings
stopifnot(length(s2$pen$w) == 20L, !is.null(names(s2$pen$w)), is.character(s2$pen$rule), s2$pen$pen > 0,
          s2$pen$n_orderings >= 10L, length(s2$asc$w) == 20L)
# the S3 seed screen (Ruling S6c) is stored with the 24 candidates and every tally
scr <- s3$screen
stopifnot(is.list(scr), scr$n_seeds == 24L, identical(as.integer(scr$seeds), 1:24),
          all(c("twfe", "cs", "ife", "mc_default", "mc", "multisynth") %in% names(scr$covered)),
          all(vapply(scr$covered, function(x) x >= 0 && x <= 24, logical(1))))
# Task 10 review, Important 1: across the 24 screened draws the refit MC overestimates the
# true effect in every draw (it is biased upward in this design; the stored draw is the most
# favourable of the 24 for both MC recipes, see 05_staggered.R and the s3_screen_stored_* keys)
stopifnot(sum(scr$per_seed$mc_est > tr$S3$att) == scr$n_seeds)
# Task 10 review, Important 3: the cached screen must still describe the current recipes and
# packages. Its row for the stored seed must reproduce res_stag.rds's six point estimates, the
# package versions it recorded must be the installed ones, and res_stag.rds must carry the cache.
source("stag_fns.R")
stale <- screen_staleness(scr, s3, tr$seeds[["S3"]])
if (length(stale))
  stop("the cached S3 screen is stale:\n  ", paste(stale, collapse = "\n  "),
       "\nRecompute it:  RUN_SLOW=1 Rscript 05b_stag_screen.R  (or RUN_SLOW=1 Rscript 00_run_all.R)", call. = FALSE)
stopifnot(identical(scr, readRDS("res_stag_screen.rds")))
# the guard itself must not pass vacuously: a screen without recorded package versions, or without
# a point-estimate column, is reported as a problem (Task 10 minor, Task 12a)
scr_nov <- scr; scr_nov$versions <- NULL
scr_noc <- scr; scr_noc$per_seed$mc_est <- NULL
stopifnot(length(screen_staleness(scr_nov, s3, tr$seeds[["S3"]])) >= 1L,
          length(screen_staleness(scr_noc, s3, tr$seeds[["S3"]])) >= 1L,
          length(screen_staleness(scr, s3, tr$seeds[["S3"]])) == 0L)
# S4: both designs are evaluated on the same 500 replicate pilot windows
stopifnot(s4$reps == 500L, s4$m_pilot == 5L, length(s4$synthetic$treated) >= 1L, length(s4$synthetic$treated) <= 5L)

## ---- Part 3: final-review locks on seed-dependent teaching claims ---------------------
# S4: the comparison the chapter teaches is the synthetic design against random assignment
# plus DiD (the 8.8x RMSE ratio), not only against the difference in means of Part 1
stopifnot(s4$rmse_synthetic_design < s4$rmse_random_did)
# S3: CS with not-yet-treated controls misses the truth too (that pool is mostly never-treated metros)
stopifnot(!covers(s3$cs_nyt$ci, tr$S3$att))
# S3: multisynth's bootstrap interval is wide enough to cover zero
stopifnot(covers(s3$multisynth$ci, 0))
# S2: the SDID interval (placebo standard error) covers the true effect
stopifnot(covers(s2$sdid$ci, tr$S2$att_post))
# S1: the full scpi intervals cover the true effect in all 12 post months, the in-sample-only
# intervals in fewer (9): the out-of-sample uncertainty cannot be dropped
stopifnot(sum(s1$scpi$covered) == 12L, sum(s1$scpi$covered_in) < sum(s1$scpi$covered))
# S2: the ASC conformal p-value (augsynth's default inference) does not reject at 10%
stopifnot(s2$asc$p_conformal > 0.10)

## ---- Part 4: Task 14 final-review locks on claims the text states --------------------
# S3: fect's cross-validation picks r = 2, the number of factors in the DGP
stopifnot(s3$ife$r_cv == 2L)
# S3 screen: the tallies the text quotes. Intervals covering 0.06 in 0 (TWFE, CS), 24 (IFE,
# multisynth), 6 (MC default) and 14 (MC refit) of the 24 draws; r = 2 chosen in 23 of them;
# TWFE and CS above the truth in every draw
cov <- vapply(c("twfe", "cs", "ife", "mc_default", "mc", "multisynth"), function(k) as.numeric(scr$covered[[k]]), 0)
stopifnot(identical(unname(cov), c(0, 0, 24, 6, 14, 24)), sum(scr$per_seed$ife_r_cv == 2) == 23L,
          all(scr$per_seed$twfe_est > tr$S3$att), all(scr$per_seed$cs_est > tr$S3$att))
# S1: SDID is farther from the truth than the hand-written SC (V = I)
stopifnot(abs(s1$sdid$est - tr$S1$att_post) > abs(s1$sc_hand$est - tr$S1$att_post))
# S2: outside the hull, SC and penalized SC under the holdout rule both get the sign wrong
stopifnot(s2$sc$est < 0, s2$pen$holdout$est < 0)

## ---- Part 5: Ruling S12, scpi on the staggered rollout (05c_stag_scpi.R) --------------
# res_stag_scpi.rds exists and holds one prediction interval per event month 0..11 for the 15
# adopters against the 65 never-treated donors; every interval contains its point estimate,
# and the reported mean lies inside the range of the bounds
stopifnot(file.exists("res_stag_scpi.rds"))
sc3 <- readRDS("res_stag_scpi.rds"); be <- sc3$by_event
# the text (6.5, third reading point): the separate fits leave the average path's pre-period gap
# (Theorem 7's q^pool) as a common bias, here negative before and after adoption
stopifnot(length(sc3$pool_gap) == 12L, sc3$pre_gap_mean < 0, sc3$est < tr$S3$att, sc3$q_pool12 > 0)
stopifnot(is.data.frame(be), nrow(be) == 12L, identical(as.integer(be$e), 0:11), sc3$n_months == 12L,
          sc3$n_treated == 15L, sc3$n_donors == 65L, sc3$n_cells == 180L,
          all(is.finite(c(be$est, be$lo, be$hi))), all(be$lo < be$est & be$est < be$hi),
          abs(sc3$est - mean(be$est)) < 1e-12, sc3$lo_min <= sc3$est, sc3$est <= sc3$hi_max,
          sc3$settings$sims == 200L, sc3$settings$failed_sims == 0)
# the text states that 2 of the 12 intervals cover 0.06
stopifnot(sum(vapply(seq_len(nrow(be)), function(i) covers(c(be$lo[i], be$hi[i]), tr$S3$att), logical(1))) == 2L,
          sc3$n_covered == 2L)

cat("test_conclusions: OK\n")
