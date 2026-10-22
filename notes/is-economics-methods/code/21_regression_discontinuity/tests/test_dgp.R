# tests/test_dgp.R
# Checks on the simulated Aster blocks written by 01_dgp.R.
# Run from the code directory AFTER 01_dgp.R:  Rscript tests/test_dgp.R   (a few seconds)
#
# Part 1 is the Task 6 brief's test, verbatim. Part 2 adds the non-vacuous guards of
# ruling R2 (every vector a check reads exists and is non-empty; every block has the
# R2 layout) and the design properties each block promises in its 01_dgp.R header.
source("rd_fns.R"); source("tests/helpers.R")      # helpers.R: is_block()

R1 <- readRDS("aster_R1.rds"); R1c <- readRDS("aster_R1c.rds"); tr <- readRDS("truth21.rds")
R4g <- readRDS("aster_R4g.rds")

## ---- Part 0: guards (ruling R2) ---------------------------------------------------------
R2 <- readRDS("aster_R2.rds"); R3 <- readRDS("aster_R3.rds"); R4 <- readRDS("aster_R4.rds"); H <- readRDS("aster_R1h.rds")
nonempty <- function(v) !is.null(v) && length(v) > 0
stopifnot(nonempty(R2$x), nonempty(R4$cutoff), nonempty(H$x), is.function(R3$meta$b_fun),
          is.numeric(R3$meta$cutoff), length(R3$meta$cutoff) == 1, is.finite(R3$meta$cutoff),
          nonempty(tr$R4$tau_c))
stopifnot(any(H$x == 79), any(H$x == 80))          # the heaping comparison below is not against an empty cell
stopifnot(length(unique(R2$x)) >= 10)              # R2 is discrete but not degenerate

## ---- Part 1: the brief ------------------------------------------------------------------
stopifnot(length(unique(R2$x)) <= 60, min(table(R2$x)) >= 10)                 # finite support, mass at every point
bf <- R3$meta$b_fun; cc <- R3$meta$cutoff                                     # policy continuous at cutoff, slope changes
stopifnot(abs(bf(cc + 1e-9) - bf(cc - 1e-9)) < 1e-6, abs((bf(cc + 1) - bf(cc)) - (bf(cc) - bf(cc - 1))) > 1e-6)
stopifnot(length(unique(R4$cutoff)) >= 3, diff(range(readRDS("truth21.rds")$R4$tau_c)) > 0.1)
stopifnot(mean(H$x == 80) > 3 * mean(H$x == 79))                                 # heaping is there

## ---- Part 2: layout (ruling R2) ------------------------------------------------------------
blocks <- list(R1c = R1c, R1h = H, R2 = R2, R3 = R3, R4 = R4, R4g = R4g)
for (nm in names(blocks)) {
  B <- blocks[[nm]]
  stopifnot(is_block(B), all(c("x", "y") %in% names(B)) || nm == "R4g", nonempty(B$y), all(is.finite(B$y)),
            B$meta$N == length(B$y), is.numeric(B$meta$seed))
  # the stored block is what its generator gives for the stored n and seed: doubles (and
  # meta) at tolerance 1e-9, since libm may drift across machines (as for ruling R1);
  # integer and character columns exactly
  G <- get(paste0("gen_", nm))(n = B$meta$N, seed = B$meta$seed)
  stopifnot(identical(names(G), names(B)), isTRUE(all.equal(G$meta, B$meta, tolerance = 1e-9)))
  for (v in setdiff(names(B), "meta")) {
    if (is.double(B[[v]])) stopifnot(is.double(G[[v]]), isTRUE(all.equal(G[[v]], B[[v]], tolerance = 1e-9)))
    else stopifnot(identical(G[[v]], B[[v]]))
  }
}
# the seven new and inherited seeds are distinct
seeds <- c(101, vapply(blocks, function(B) B$meta$seed, 1))
stopifnot(!anyDuplicated(seeds))

## ---- Part 2: truth21 fields -----------------------------------------------------------------
fin <- function(v, n = 1) is.numeric(v) && length(v) >= n && all(is.finite(v))
stopifnot(fin(tr$R1$jump), fin(tr$R1h$jump), fin(tr$R2$jump), fin(tr$R3$kink_effect),
          fin(tr$R4$tau_c, 3), fin(tr$R4$tau_extrap), fin(tr$R4g$tau_b, 3))
stopifnot(abs(tr$R1$jump - 0.20) < 1e-12, abs(tr$R1h$jump - 0.20) < 1e-12)

## ---- Part 2: R1, the inherited block (legacy layout and printed numbers) ---------------------
stopifnot(identical(names(R1), c("dat", "meta")), is.data.frame(R1$dat), nrow(R1$dat) == 12000,
          identical(names(R1$dat), c("score", "s", "mu0", "D_sharp", "y_sharp", "D_fuzzy", "y_fuzzy")),
          R1$meta$c_cut == 80, abs(R1$meta$tau_true - 0.20) < 1e-12,
          abs(R1$meta$naive_diff - 1.891) < 5e-4, R1$meta$n_below == 8436, R1$meta$n_above == 3564)
# the fuzzy half of the R1 stream (drawn after the sharp half): badge shares at the printed
# precision of 01_dgp_rdd.R (0.098, 0.849) and column checksums (integer exact, doubles at
# relative tolerance), so a shifted stream fails here before Task 7's res_core.rds exists
stopifnot(abs(R1$meta$frac_badge_below - 0.098) < 5e-4, abs(R1$meta$frac_badge_above - 0.849) < 5e-4,
          sum(R1$dat$D_fuzzy) == 3848L, sum(R1$dat$D_sharp) == 3564L,
          isTRUE(all.equal(sum(R1$dat$y_fuzzy), 64967.19029935, tolerance = 1e-9)),
          isTRUE(all.equal(sum(R1$dat$y_sharp), 64866.83104977, tolerance = 1e-9)),
          isTRUE(all.equal(sum(R1$dat$score), 888923.35478813, tolerance = 1e-9)))

## ---- Part 2: R1c, same world plus predetermined covariates -------------------------------------
stopifnot(R1c$meta$cutoff == 80, identical(R1c$d, as.integer(R1c$x >= 80)),
          all(c("gmv_prev", "tenure") %in% names(R1c)), min(R1c$x) >= 45, max(R1c$x) <= 100)
# the covariates genuinely predict the outcome: near the cutoff (|x - 80| <= 5, where a cubic
# follows the R1 wiggle), the residual of y net of d and x correlates with that of gmv_prev
# (design value about 0.56), and the two covariates explain at least 20% of the residual
# variance left after d and x (partial R^2; design value about 0.37)
w5 <- abs(R1c$x - 80) <= 5
D5 <- as.data.frame(R1c[c("x", "d", "y", "gmv_prev", "tenure")])[w5, ]
r_y <- resid(lm(y ~ d * poly(x, 3), data = D5))
r_g <- resid(lm(gmv_prev ~ poly(x, 3), data = D5))
stopifnot(nrow(D5) > 2000, cor(r_y, r_g) > 0.4)
rss_covs <- sum(resid(lm(y ~ d * poly(x, 3) + gmv_prev + tenure, data = D5))^2)
stopifnot(1 - rss_covs / sum(r_y^2) > 0.2)
# and are predetermined: no jump at the cutoff. gmv_prev inherits the wiggly R1 mean mu0(x),
# whose curvature biases a local-linear fit at h = 5, so its check nets out mu0 first
stopifnot(abs(ll_jump(R1c$gmv_prev - .r1_mu0(R1c$x), R1c$x, 80, 5)$est) < 0.08,
          abs(ll_jump(R1c$tenure, R1c$x, 80, 5)$est) < 0.5)

## ---- Part 2: R1h, heaping at 80 ---------------------------------------------------------------
stopifnot(all(abs(H$x * 10 - round(H$x * 10)) < 1e-9),           # scores on a one-decimal grid
          identical(H$d, as.integer(H$x >= 80)), H$meta$cutoff == 80,
          all(H$x[H$heaped_true == 1] == 80), all(H$x_true[H$heaped_true == 1] < 80),
          all(H$x[H$heaped_true == 0] == H$x_true[H$heaped_true == 0]), mean(H$heaped_true) > 0.01)

## ---- Part 2: R2, discrete running variable ------------------------------------------------------
sx <- sort(unique(R2$x))
stopifnot(all(abs(sx * 10 - round(sx * 10)) < 1e-9), min(sx) >= 3, max(sx) <= 5,   # averages of 10 integer ratings
          R2$meta$cutoff == 4.25, !any(abs(sx - 4.25) < 1e-9),
          identical(R2$d, as.integer(R2$x >= 4.25)), sum(sx < 4.25) >= 5, sum(sx > 4.25) >= 5,
          fin(R2$meta$M_bound), R2$meta$M_bound > 0)

## ---- Part 2: R3, kinked commission schedule --------------------------------------------------------
stopifnot(all(abs(R3$b - bf(R3$x)) < 1e-12), sum(R3$x < cc) > 1000, sum(R3$x >= cc) > 1000,
          abs((bf(cc + 100) - bf(cc)) / 100 - R3$meta$slope_right) < 1e-12,
          abs((bf(cc) - bf(cc - 100)) / 100 - R3$meta$slope_left) < 1e-12,
          tr$R3$kink_effect < 0)                                    # a higher commission lowers GMV
# serialised policy function is small (no captured data in its environment)
stopifnot(length(serialize(bf, NULL)) < 20000)

## ---- Part 2: R4, category-specific cutoffs ------------------------------------------------------------
cuts <- sort(unique(R4$cutoff))
stopifnot(identical(names(tr$R4$tau_c), as.character(cuts)), identical(R4$d, as.integer(R4$x >= R4$cutoff)),
          all(tapply(R4$category, R4$cutoff, function(v) length(unique(v))) == 1),
          tr$R4$x_extrap > min(cuts), tr$R4$x_extrap < max(cuts), !(tr$R4$x_extrap %in% cuts),
          tr$R4$extrap_cutoff %in% cuts, tr$R4$ref_cutoff %in% cuts,
          tr$R4$extrap_cutoff < tr$R4$x_extrap, tr$R4$x_extrap < tr$R4$ref_cutoff)
# the true pooled target (CKTV 2016 Lemma 1 weights) is not the simple mean
stopifnot(fin(tr$R4$tau_pooled), abs(tr$R4$tau_pooled - mean(tr$R4$tau_c)) > 0.04,
          abs(sum(tr$R4$omega) - 1) < 1e-12, identical(names(tr$R4$omega), names(tr$R4$tau_c)),
          abs(tr$R4$tau_pooled - sum(tr$R4$omega * tr$R4$tau_c)) < 1e-12)
# omega matches the sampled DGP: scores are clamped (not truncated), so the density at an
# interior cutoff is dnorm; and the sample shares within h = 1 of the own cutoff (the
# pooled_cutoff_weights analog) are close to omega (sampling sd about 0.01)
Cg <- R4$meta$categories
om <- dnorm(Cg$cutoff, Cg$mu, Cg$sd) * Cg$share
stopifnot(all(abs(tr$R4$omega - om / sum(om)) < 1e-12),
          all(abs(pooled_cutoff_weights(R4$cutoff, R4$x, h = 1) - tr$R4$omega) < 0.03))

## ---- Part 2: R4g, delivery-zone boundary ---------------------------------------------------------------
bp <- R4g$meta$boundary
stopifnot(all(c("x1", "x2", "d", "y") %in% names(R4g)), is.data.frame(bp), nrow(bp) == length(tr$R4g$tau_b),
          nrow(bp) >= 3, identical(R4g$d, as.integer(sqrt((R4g$x1 - R4g$meta$hub[1])^2 + (R4g$x2 - R4g$meta$hub[2])^2) <= R4g$meta$radius)),
          all(abs(sqrt((bp$b1 - R4g$meta$hub[1])^2 + (bp$b2 - R4g$meta$hub[2])^2) - R4g$meta$radius) < 1e-9),
          diff(range(tr$R4g$tau_b)) > 0.05)
# enough sellers near every boundary point on both sides
for (k in seq_len(nrow(bp))) {
  dk <- sqrt((R4g$x1 - bp$b1[k])^2 + (R4g$x2 - bp$b2[k])^2)
  stopifnot(sum(dk <= 4 & R4g$d == 1) >= 100, sum(dk <= 4 & R4g$d == 0) >= 100)
}

cat("test_dgp: OK\n")
