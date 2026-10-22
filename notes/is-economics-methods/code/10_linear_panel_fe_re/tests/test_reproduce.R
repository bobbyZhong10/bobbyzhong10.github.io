# ------------------------------------------------------------------
# tests/test_reproduce.R -- the inherited static block P1 (legacy Nimbus code, set.seed(737))
# reproduces the legacy numbers (chapter card, verbatim except for the Hausman field name:
# the legacy 02_static_fe_re.R stores the statistic as haus_stat, not hausman).
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
r <- readRDS("res_static.rds")
stopifnot(abs(r$pool - 0.9511) < 5e-5, abs(r$re - 0.7390) < 5e-5, abs(r$fe - 0.5063) < 5e-5, abs(r$haus_stat - 3362.02) < 0.01)
