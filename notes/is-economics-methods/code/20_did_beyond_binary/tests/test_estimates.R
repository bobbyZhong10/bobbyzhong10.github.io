# tests/test_estimates.R
# Locks the teaching conclusions of chapter 20 against the stored estimates.
# Run from the code directory AFTER 02_fuzzy.R, 03_dose.R and 04_switching.R:
#     Rscript tests/test_estimates.R
#
# Part 1 is the Task 11 brief's test, verbatim (only the field access follows
# what the scripts store). Part 2 locks further conclusions the chapter states;
# they come from the DGP theory and from Task 10's rehearsal numbers, and were
# written before 02-04 were run.
f  <- readRDS("res_fuzzy.rds"); dz <- readRDS("res_dose.rds"); sw <- readRDS("res_switch.rds")
tr <- readRDS("truth20.rds")
covers <- function(x, v) x$lo <= v && v <= x$hi

## ---- Part 1: the brief --------------------------------------------------------------
stopifnot(covers(f$tc, tr$A$late), !covers(f$did, tr$A$late))          # Wald-TC covers the truth, Wald-DID does not
stopifnot(covers(dz$cont$overall, tr$B$att_overall))                    # continuous treatment: overall ATT covered
stopifnot(abs(sw$twfe$est - tr$C$avg_effect) > 3 * sw$twfe$se)           # TWFE is off under switching treatments
stopifnot(covers(sw$dyn$l1, tr$C$delta[1]))                              # DID_1 covers the true effect at horizon 1
stopifnot(sw$weights$share_negative > 0)

## ---- Part 2, block A: fuzzy take-up ---------------------------------------------------
stopifnot(covers(f$cic, tr$A$late))                                      # Wald-CIC covers too (Assumptions 2, 7, 8 hold)
stopifnot(f$did$est < tr$A$late, f$did$hi < tr$A$late)                  # Wald-DID is biased NEGATIVELY (pA1 < pA0)
stopifnot(f$did$lo <= f$did$est, f$did$est <= f$did$hi, f$tc$lo < f$tc$est, f$tc$est < f$tc$hi,
          f$cic$lo < f$cic$est, f$cic$est < f$cic$hi)
# variant A': rising control take-up violates Assumption 2; Wald-TC and Wald-CIC overshoot
stopifnot(f$rising$tc$est - tr$A_rising$late > 2, f$rising$cic$est - tr$A_rising$late > 2,
          abs(f$rising$did$est - tr$A_rising$late) > 2)
stopifnot(!covers(f$rising$tc, tr$A_rising$late), !covers(f$rising$cic, tr$A_rising$late))

## ---- Part 2, block B: continuous dose --------------------------------------------------
# AS / WAS (did_multiplegt_stat) and the HAD WAS cover their own targets
stopifnot(covers(dz$stat$aoss, tr$B$as), covers(dz$stat$waoss, tr$B$was), covers(dz$had, tr$B$had_was))
# contdid's overall ACRT recovers the parallel-trends-only slope, not the true ACRT^o
a <- dz$cont$acrt_overall$est
stopifnot(abs(a - tr$B$acrt_overall_pt_pop) < abs(a - tr$B$acrt_overall))
# the estimated ATT(d|d) curve is stored on a dose grid inside the support of the true curve
stopifnot(is.data.frame(dz$cont$curve), all(c("dose", "est", "lo", "hi") %in% names(dz$cont$curve)),
          all(dz$cont$curve$dose >= min(tr$B$att_curve$dose)), all(dz$cont$curve$dose <= max(tr$B$att_curve$dose)))

## ---- Part 2, block C: switching treatment ---------------------------------------------
stopifnot(sw$twfe$est < tr$C$avg_effect)                                 # TWFE is biased downward
stopifnot(sw$weights$n_negative > 0)
for (l in 1:8) stopifnot(covers(sw$dyn[[paste0("l", l)]], tr$C$delta[l]))   # every DID_l covers its truth
stopifnot(all(vapply(1:8, function(l) covers(sw$dyn$norm[[paste0("l", l)]], tr$C$delta_norm[l]), TRUE)))
stopifnot(covers(sw$dyn$avg_total, tr$C$avg_total_effect))
stopifnot(sw$dyn$p_jointplacebo > 0.05)                                  # joint placebo test does not reject

cat("test_estimates: OK\n")
