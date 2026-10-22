# tests/test_reproduce.R
# v2 verified numbers; they must reproduce under R 4.6.1.
# Run from the code directory:  Rscript tests/test_reproduce.R
s <- readRDS("res_twfe.rds"); m <- readRDS("res_modern.rds"); b <- readRDS("res_bacon.rds")
stopifnot(abs(s$static$est - 0.0227) < 5e-5, abs(s$static$se - 0.0087) < 5e-5)
# res_modern.rds stores the BJS imputation estimator under `imp` (brief's `bjs`)
stopifnot(abs(m$cs_simple$est - 0.0818) < 5e-5,
          abs(m$sa_att$est   - 0.0806) < 5e-5,
          abs(m$imp$est      - 0.0777) < 5e-5)
# Goodman-Bacon weights by comparison type (v2: 0.3042 / 0.3061 / 0.3898)
w <- sort(round(b$summary$weight, 4))
stopifnot(isTRUE(all.equal(w, sort(c(0.3042, 0.3061, 0.3898)), tolerance = 1e-9)))
cat("test_reproduce: OK\n")
