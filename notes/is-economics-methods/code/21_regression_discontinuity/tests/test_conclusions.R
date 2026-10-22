# tests/test_conclusions.R
# The teaching conclusions of chapter 21 that rest on the new blocks, checked against the
# res_*.rds files written by 03_discrete_honest.R, 04_kink.R, 05_multi.R and
# 06_covs_donut_power.R. Run from the code directory:  Rscript tests/test_conclusions.R
#
# Part 1 is the Task 7 brief's test, code verbatim, comments translated to English (the
# direction of a check never changes with the results; block parameters were fixed in
# Task 6), except one line: the brief's identity check on the pooled estimator (true by
# construction in 05_multi.R) is replaced by a coverage check (Task 14 final review). Part 2
# adds the ruling R3 check and guards that keep Part 1 from passing vacuously. Part 3 locks
# the directions that the chapter text asserts.
source("../../tools/R/chapter_kit.R")

## ---- Part 1: the brief ------------------------------------------------------------
tr <- readRDS("truth21.rds"); d <- readRDS("res_discrete.rds"); k <- readRDS("res_kink.rds")
m <- readRDS("res_multi.rds"); x <- readRDS("res_extra.rds")
stopifnot(d$mc$cover_honest >= 0.93, d$mc$cover_conv <= 0.88)           # discrete running variable: honest covers, conventional undercovers
stopifnot(covers(d$honest$ci, tr$R2$jump), covers(d$locrand$ci, tr$R2$jump))
stopifnot(covers(k$rk$ci, tr$R3$kink_effect), abs(k$level_jump$est) < 2 * k$level_jump$se)   # no level jump, the slope change is identified
stopifnot(covers(m$pooled$ci, sum(m$weights * tr$R4$tau_c)))             # the Lemma 1 weighted estimator covers its target, the sample-weighted truth
stopifnot(all(mapply(covers, m$by_cutoff$ci, tr$R4$tau_c)), covers(m$extrap$ci, tr$R4$tau_extrap))
stopifnot(abs(m$pooled$est - mean(tr$R4$tau_c)) > 0.02)                 # pooled is not the simple mean: who sets the weights must be explained
stopifnot(mean(mapply(covers, m$geo$ci, tr$R4g$tau_b)) >= 0.8)
stopifnot(x$covs$ci_len / x$nocovs$ci_len <= 0.9, abs(x$covs$est - x$nocovs$est) < 0.05)   # covariates buy precision, not a new estimand
stopifnot(x$heap$density_p < 0.05, abs(x$donut$est - tr$R1h$jump) < abs(x$heap_naive$est - tr$R1h$jump))
stopifnot(x$power$power > 0.8)

## ---- Part 2: ruling R3 ------------------------------------------------------------
# The pooled RD on the normalized score X - c_i (what applied work runs) targets the
# Lemma 1 weighted average of the cutoff effects, with the weights of res_multi$weights.
stopifnot(covers(m$pooled_norm$ci, sum(m$weights * tr$R4$tau_c)))
# guards: the weights are the ruling R3 weights (named by cutoff, nonnegative, sum to one) and
# line up with the by-cutoff rows and the truth, so the coverage checks above compare like with like
stopifnot(identical(names(m$weights), c("70", "80", "90")), all(m$weights >= 0), abs(sum(m$weights) - 1) < 1e-12,
          identical(as.character(m$by_cutoff$cutoff), names(m$weights)), identical(names(tr$R4$tau_c), names(m$weights)),
          length(m$geo$ci) == length(tr$R4g$tau_b), length(m$geo$ci) == 6L, d$mc$reps >= 500)

## ---- Part 3: directions asserted in the text (Task 14 final review) ------------------
# R2 Monte Carlo (Section 3.4): clustering by support point covers less often than the
# conventional interval, bias correction does not restore nominal coverage, and plugging in
# the data-based lower confidence limit for M undercovers. "Below nominal" means more than
# three Monte Carlo standard errors below 0.95.
mc_floor <- 0.95 - 3 * sqrt(0.95 * 0.05 / d$mc$reps)
stopifnot(d$mc$cover_crv < d$mc$cover_conv, d$mc$cover_robust < mc_floor, d$mc$cover_honest_lb < mc_floor)
# R3 kink (Section 3.6): local linear misses the true kink effect by more than local quadratic
stopifnot(abs(k$rk_p1$est - tr$R3$kink_effect) > abs(k$rk$est - tr$R3$kink_effect))

cat("test_conclusions: OK\n")
