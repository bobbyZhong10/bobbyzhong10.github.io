# ------------------------------------------------------------------
# tests/test_reproduce.R -- the inherited blocks D1 (dyn_main, gen_dyn(seed = 808)), D1p (dyn_persist,
# seed 7372) and D1T (the Nickell Monte Carlo, set.seed(303)), run by the legacy 01_dgp.R and
# 03_dynamic_nickell.R unchanged, reproduce the legacy numbers. The first stopifnot is the chapter
# card, verbatim (the legacy rds stores pool, fe, ab, bb as named one-element vectors).
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
d <- readRDS("res_dynamic.rds")
stopifnot(abs(d$pool - 0.9461) < 5e-5, abs(d$fe - 0.4246) < 5e-5, abs(d$ab - 0.7045) < 5e-5, abs(d$bb - 0.7005) < 5e-5)
# Arellano-Bond specification tests (legacy text: Sargan p = 0.29, AR(1) p = 0.000, AR(2) p = 0.33) and
# the rho -> 1 block D1p (difference GMM 0.809, SE 0.119; system GMM 0.911, SE 0.017)
stopifnot(abs(d$sargan_p - 0.290) < 5e-4, abs(d$ar2_p - 0.334) < 5e-4, d$ar1_p < 5e-4,
          abs(d$ab9 - 0.8094) < 5e-5, abs(d$se_ab9 - 0.1185) < 5e-5, abs(d$bb9 - 0.9113) < 5e-5,
          abs(d$se_bb9 - 0.0169) < 5e-5)
# the legacy Nickell table D1T (rows fe, ab, nickell_approx; columns T_obs = 4, 8, 16, 32)
nk <- rbind(fe = c(0.0821, 0.4169, 0.5703, 0.6400), ab = c(0.7050, 0.6864, 0.6995, 0.7014),
            nickell_approx = c(0.1333, 0.4571, 0.5867, 0.6452))
stopifnot(identical(dimnames(d$nick), list(c("fe", "ab", "nickell_approx"), c("4", "8", "16", "32"))),
          max(abs(d$nick - nk)) < 5e-5)
