# tests/test_dgp.R
# Checks on the simulated Aster world written by 01_dgp.R (blocks S1 to S4 and
# truth22.rds). Run from the code directory AFTER 01_dgp.R:
#     Rscript tests/test_dgp.R
#
# Part 1 is the Task 9 brief's test, verbatim except that its Chinese comments are
# translated to English (Ruling S7; assertions unchanged). Part 2 checks the truth
# file, the single seed definition, that S1 is the legacy block, the structure each
# new block promises, the S3 selection mechanism (Ruling S5) and that gen_S3() in
# dgp_s3.R reproduces aster_S3.rds (Ruling S6). Part 2 was written after the block
# designs were fixed, to lock them.

## ---- Part 1: the brief ------------------------------------------------------------
S2 <- readRDS("aster_S2.rds"); S3 <- readRDS("aster_S3.rds")
stopifnot(isTRUE(S2$meta$treated_outside_hull))                          # 01_dgp.R decides this by LP feasibility
stopifnot(length(unique(S3$panel$cohort[is.finite(S3$panel$cohort)])) >= 3,
          abs(cor(S3$meta$loading1, S3$meta$adopt_rank)) > 0.3)          # loadings related to adoption (brief: "loadings drive adoption timing"; in S3 mostly who adopts, see Part 2)

## ---- Part 2: truth22.rds ------------------------------------------------------------
source("sc_fns.R")
tr <- readRDS("truth22.rds"); S1 <- readRDS("aster_S1.rds"); S4 <- readRDS("aster_S4.rds")
fin1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
stopifnot(fin1(tr$S1$att_post), fin1(tr$S2$att_post), fin1(tr$S3$att), fin1(tr$S4$tau))
stopifnot(length(tr$S1$tau_path) == 12, abs(tr$S1$att_post - mean(tr$S1$tau_path)) < 1e-15,      # Ruling S1
          identical(tr$S1$tau_path, S1$meta$true_effect_path[25:36]), tr$S1$att_post == S1$meta$true_post_avg)
stopifnot(is.data.frame(tr$S3$att_gt), all(c("g", "t", "att") %in% names(tr$S3$att_gt)))
# the seeds are defined once, as seeds22 at the top of 01_dgp.R; truth22 must carry them
seed_line <- grep("^seeds22 <- c\\(", readLines("01_dgp.R"), value = TRUE)
stopifnot(length(seed_line) == 1)
seeds22 <- eval(parse(text = seed_line))
stopifnot(identical(tr$seeds, seeds22), seeds22[["S1"]] == 12021)            # S1 keeps the legacy seed

## ---- Part 2: S1 is the legacy block (layout and printed numbers) -------------------------
stopifnot(identical(names(S1), c("panel", "Y", "meta", "unit_names")), all(dim(S1$Y) == c(21, 36)), S1$meta$T0 == 25,
          identical(names(S1$panel), c("metro", "unit", "month", "log_gmv", "treated", "post", "treat_it")))
stopifnot(abs(S1$meta$true_post_avg - 0.095) < 5e-4, abs(S1$meta$true_effect_path[36] - 0.118) < 5e-4,
          abs(S1$meta$did_est - 0.244) < 5e-4, abs(S1$meta$tr_pre_slope - 0.1306) < 5e-5,
          abs(S1$meta$do_pre_slope - 0.1165) < 5e-5)

## ---- Part 2: S2, the highest-level metro, outside the donor hull --------------------------------
par2 <- S2$meta$params
stopifnot(identical(names(S2), names(S1)), identical(names(S2$panel), names(S1$panel)), all(dim(S2$Y) == c(21, 36)),
          S2$unit_names[1] == "Aster-Capital", S2$meta$T0 == 25, S2$meta$n_pre == 24)
stopifnot(par2["level", 1] > max(par2["level", -1]), S2$meta$hull_l1_gap > 0.1,
          isFALSE(as.vector(in_hull_lp(par2[, 1], par2[, -1]))),
          isTRUE(as.vector(in_hull_lp(par2[, 2], par2[, -1]))))                  # a donor is a vertex
stopifnot(isFALSE(pensynth::in_convex_hull(par2[, 1], par2[, -1], verbose = FALSE,       # independent solver (clarabel)
                                           control = clarabel::clarabel_control(verbose = FALSE))))
stopifnot(all(S2$panel$treat_it == as.integer(S2$panel$unit == 1 & S2$panel$month >= 25)),
          identical(unname(S2$Y[cbind(S2$panel$unit, S2$panel$month)]), S2$panel$log_gmv))
stopifnot(length(tr$S2$tau_path) == 12, tr$S2$att_post == S2$meta$true_post_avg,
          identical(tr$S2$tau_path, S2$meta$true_effect_path[25:36]))

## ---- Part 2: S3, staggered waves; the cycle loading drives who adopts (Ruling S5) ---------
m3 <- S3$meta; p3 <- S3$panel
coh <- as.vector(tapply(p3$cohort, p3$unit, function(x) x[1]))
stopifnot(nrow(p3) == m3$N * m3$Tt, all(table(p3$unit) == m3$Tt),
          all(tapply(p3$cohort, p3$unit, function(x) length(unique(x))) == 1), identical(coh, m3$cohort))
stopifnot(all(p3$treat_it == as.integer(p3$month >= p3$cohort)),
          all(table(m3$cohort[is.finite(m3$cohort)]) == m3$n_wave), sum(is.infinite(m3$cohort)) == m3$N - 3 * m3$n_wave)
stopifnot(length(m3$loading1) == m3$N, identical(m3$adopt_rank, rank(m3$cohort, ties.method = "average")),
          cor(m3$loading1, m3$adopt_rank) < -0.3,                                 # more exposure -> lower adopt_rank (mostly: adopters rank before the tied never-treated)
          all(tapply(m3$loading1, m3$cohort, mean)[1:3] > mean(m3$loading1[is.infinite(m3$cohort)])),
          all(is.infinite(m3$cohort[m3$contract])))                              # contract metros never treated
stopifnot(abs(tr$S3$att - mean(p3$tau_true[p3$treat_it == 1])) < 1e-15, all(p3$tau_true[p3$treat_it == 0] == 0))
gt <- tr$S3$att_gt
stopifnot(nrow(gt) == sum(m3$Tt - m3$waves + 1), all(gt$t >= gt$g), all(abs(gt$att - tr$S3$att) < 1e-15))
# What the cycle loading drives (Ruling S5): selection into the program, not the wave
# order. Stored draw: adopters' mean loading1 0.208 vs 0.113 for the never-treated (a gap
# of 0.095, about 1.3 sd of U(0, 0.25)); the test asks for a gap above 0.06 (0.8 sd).
# Within the 15 adopters cor(loading1, cohort) = -0.113 (wave means 0.212, 0.210, 0.203):
# wave order is nearly random, and the test asks for |cor| < 0.3.
ad3 <- is.finite(m3$cohort)
stopifnot(mean(m3$loading1[ad3]) - mean(m3$loading1[!ad3]) > 0.06,
          abs(cor(m3$loading1[ad3], m3$cohort[ad3])) < 0.3)
# gen_S3() in dgp_s3.R redraws exactly the stored block and its truth (Ruling S6)
source("dgp_s3.R")
g3 <- gen_S3(seeds22[["S3"]])
stopifnot(identical(g3$data, S3), identical(g3$truth, tr$S3))

## ---- Part 2: S4, the pilot design ---------------------------------------------------------
m4 <- S4$meta
stopifnot(all(dim(S4$Y0) == c(m4$N, m4$T_pre + m4$T_pilot)), all(dim(S4$Y0_pilot_reps) == c(m4$R, m4$N, m4$T_pilot)),
          identical(unname(S4$Y0[, m4$T_pre + 1:m4$T_pilot]), unname(S4$Y0_pilot_reps[1, , ])),
          identical(unname(S4$Y0[cbind(S4$panel$unit, S4$panel$month)]), S4$panel$log_gmv0))
stopifnot(m4$m_pilot == 5, tr$S4$tau == m4$tau, sd(S4$Y0_pilot_reps[, 1, 1]) > 0.01)     # replicates differ

cat("test_dgp.R: all checks passed\n")
