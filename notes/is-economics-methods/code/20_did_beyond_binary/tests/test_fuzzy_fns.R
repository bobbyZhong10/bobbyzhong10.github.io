# tests/test_fuzzy_fns.R
# Unit tests for the hand-coded fuzzy DiD estimators in fuzzy_did_fns.R
# (Wald-DID, Wald-TC, Wald-CIC; de Chaisemartin & D'Haultfoeuille 2018).
# Run from the code directory:  Rscript tests/test_fuzzy_fns.R
# (the Step 5 cross-check at the end reads tessera_A.rds and tessera_A_rising.rds,
# written by 01_dgp.R, and runs only if Rfuzzydid is installed; otherwise it is
# skipped with a message)
#
# Hand calculations for the toy data (cells: group first index, period second;
# Y_dgt = outcome among units with D = d in cell (g, t)):
#   E(Y00) = 11, E(Y01) = 13.8, E(Y10) = 11.5, E(Y11) = 20.3
#   E(D00) = 0.1, E(D01) = 0.2, E(D10) = 0.1, E(D11) = 0.6
#   Wald-DID = [(20.3 - 11.5) - (13.8 - 11)] / [(0.6 - 0.1) - (0.2 - 0.1)] = 6.0 / 0.4 = 15.0
#   Wald-TC:  delta_0 = 12 - 10 = 2, delta_1 = 21 - 20 = 1
#             E(Y10 + delta_D10) = 11.5 + 0.9 * 2 + 0.1 * 1 = 13.4
#             Wald-TC = (20.3 - 13.4) / (0.6 - 0.1) = 13.8
#   Wald-CIC: Q_0(10.5) = F^{-1}_{Y001}(F_{Y000}(10.5)) = F^{-1}_{Y001}(1) = 12
#             Q_1(20.5) = F^{-1}_{Y101}(F_{Y100}(20.5)) = F^{-1}_{Y101}(1) = 21
#             E(Q_D10(Y10)) = 0.9 * 12 + 0.1 * 21 = 12.9
#             Wald-CIC = (20.3 - 12.9) / 0.5 = 14.8
#   (F^{-1}(q) = inf{y : F(y) >= q}, the generalized inverse used in the paper.)
source("fuzzy_did_fns.R")

## ---- toy data: 10 units per cell -------------------------------------------
toy <- rbind(
  data.frame(g = 0, t = 0, d = c(rep(0, 9), 1),       y = c(rep(10, 9), 20)),
  data.frame(g = 0, t = 1, d = c(rep(0, 8), 1, 1),    y = c(rep(12, 8), 21, 21)),
  data.frame(g = 1, t = 0, d = c(rep(0, 9), 1),       y = c(rep(10.5, 9), 20.5)),
  data.frame(g = 1, t = 1, d = c(rep(0, 4), rep(1, 6)), y = c(rep(12.5, 4), rep(25.5, 6)))
)
wd <- wald_did(toy, "y", "d", "g", "t")
stopifnot(abs(wd$estimate - 15.0) < 1e-10, abs(wd$did_y - 6.0) < 1e-10, abs(wd$did_d - 0.4) < 1e-10)
stopifnot(abs(wald_tc(toy, "y", "d", "g", "t")$estimate - 13.8) < 1e-10)
stopifnot(abs(wald_cic(toy, "y", "d", "g", "t")$estimate - 14.8) < 1e-10)
# Wald-TC also reports delta_0 and delta_1
tc <- wald_tc(toy, "y", "d", "g", "t")
stopifnot(abs(tc$delta[["0"]] - 2) < 1e-10, abs(tc$delta[["1"]] - 1) < 1e-10)

# column names are passed as strings; other names and data.table input work too
toy2 <- setNames(toy, c("grp", "per", "takeup", "orders"))
stopifnot(abs(wald_tc(toy2, "orders", "takeup", "grp", "per")$estimate - 13.8) < 1e-10)
if (requireNamespace("data.table", quietly = TRUE))
  stopifnot(abs(wald_did(data.table::as.data.table(toy), "y", "d", "g", "t")$estimate - 15.0) < 1e-10)

## ---- pure-shift data: CIC and TC coincide ------------------------------------
# control T1 = control T0 + delta_d observation by observation (stable control
# treatment rate); treatment-group T0 values are drawn from the control-group
# T0 values with the same D. Then Q_d(y) = y + delta_d on the support, so the
# Wald-CIC equals the Wald-TC exactly.
set.seed(20)
y_c0_d0 <- round(rnorm(40, 50, 10), 3); y_c0_d1 <- round(rnorm(10, 60, 10), 3)
del0 <- 3.7; del1 <- 5.2
shift <- rbind(
  data.frame(g = 0, t = 0, d = 0, y = y_c0_d0),
  data.frame(g = 0, t = 0, d = 1, y = y_c0_d1),
  data.frame(g = 0, t = 1, d = 0, y = y_c0_d0 + del0),
  data.frame(g = 0, t = 1, d = 1, y = y_c0_d1 + del1),
  data.frame(g = 1, t = 0, d = 0, y = sample(y_c0_d0, 45, replace = TRUE)),
  data.frame(g = 1, t = 0, d = 1, y = sample(y_c0_d1, 5, replace = TRUE)),
  data.frame(g = 1, t = 1, d = rep(0:1, c(25, 25)), y = round(rnorm(50, 62, 10), 3))
)
stopifnot(abs(wald_cic(shift, "y", "d", "g", "t")$estimate -
              wald_tc(shift, "y", "d", "g", "t")$estimate) < 1e-10)
stopifnot(abs(wald_tc(shift, "y", "d", "g", "t")$delta[["0"]] - del0) < 1e-10)

## ---- generalized inverse in the quantile map (fix round 1) -------------------
# Q_d(y) = F^{-1}_{Y_d01}(F_{Y_d00}(y)) with F^{-1}(q) = inf{y : F(y) >= q}.
# Regression: F_{Y00}(7) = 7/25 and F^{-1}_{Y01}(7/25) = 107. The first version
# used stats::quantile(type = 1) at probs = k/n, where n * (k/n) can land just
# above k in floating point and return the (k+1)-th order statistic (108).
stopifnot(.fz_qq(7, 1:25, 101:125) == 107)
for (n in 2:120) stopifnot(identical(.fz_qq(1:n, 1:n, 1000 + 1:n), 1000 + 1:n))
# ties, unequal cell sizes, values below min(y00) (q = 0 maps to min(y01)) and above max(y00)
y00 <- c(1, 1, 2, 2, 2, 3); y01 <- c(10, 10, 10, 20, 30, 30, 40)
stopifnot(identical(.fz_qq(c(0, 1, 1.5, 2, 2.5, 3, 99), y00, y01), c(10, 10, 10, 30, 30, 40, 40)))
# exact reference by integer cross-multiplication: smallest v in y01 with
# #{y01 <= v} * n0 >= #{y00 <= y} * n1; tie-heavy draws with unequal sizes
gi_ref <- function(yeval, y00, y01) {
  n0 <- length(y00); n1 <- length(y01); u <- sort(unique(y01))
  cnt <- vapply(u, function(v) sum(y01 <= v), numeric(1))
  vapply(yeval, function(y) {
    k0 <- sum(y00 <= y)
    if (k0 == 0) return(min(y01))
    u[which(cnt * n0 >= k0 * n1)[1]]
  }, numeric(1))
}
set.seed(7)
for (r in 1:400) {
  n0 <- sample(1:70, 1); n1 <- sample(1:70, 1)
  y00 <- sample(1:15, n0, replace = TRUE) / 3
  y01 <- sample(1:15, n1, replace = TRUE) / 7
  ye <- c(sort(unique(y00)), min(y00) - 1, max(y00) + 1, runif(5, 0, 6))
  stopifnot(identical(.fz_qq(ye, y00, y01), gi_ref(ye, y00, y01)))
}

## ---- equal change in treatment rates: NA with a warning ---------------------
# control 0.1 -> 0.3, treatment group 0.2 -> 0.4: the first-stage DiD is zero
flat <- rbind(
  data.frame(g = 0, t = 0, d = rep(0:1, c(9, 1)), y = rep(c(10, 20), c(9, 1))),
  data.frame(g = 0, t = 1, d = rep(0:1, c(7, 3)), y = rep(c(12, 22), c(7, 3))),
  data.frame(g = 1, t = 0, d = rep(0:1, c(8, 2)), y = rep(c(11, 21), c(8, 2))),
  data.frame(g = 1, t = 1, d = rep(0:1, c(6, 4)), y = rep(c(13, 24), c(6, 4)))
)
w <- tryCatch(wald_did(flat, "y", "d", "g", "t"), warning = function(w) w)
stopifnot(inherits(w, "warning"), grepl("near\\) zero", conditionMessage(w)))
quiet <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
for (f in list(wald_did, wald_tc, wald_cic)) {
  w <- tryCatch(f(flat, "y", "d", "g", "t"), warning = function(w) w)
  stopifnot(inherits(w, "warning"), grepl("near\\) zero", conditionMessage(w)))
  r <- quiet(f(flat, "y", "d", "g", "t"))
  stopifnot(is.list(r), identical(r$estimate, NA_real_))
}

## ---- delta_d undefined: NA with a warning ------------------------------------
# treatment-group units with D = 1 at T = 0 but no control-group unit with D = 1
nod1 <- toy
nod1$d[nod1$g == 0] <- 0
nod1$d[nod1$g == 1 & nod1$t == 1] <- rep(0:1, c(2, 8))
w <- tryCatch(wald_tc(nod1, "y", "d", "g", "t"), warning = function(w) w)
stopifnot(inherits(w, "warning"), grepl("undefined", conditionMessage(w)))
stopifnot(identical(quiet(wald_tc(nod1, "y", "d", "g", "t"))$estimate, NA_real_),
          identical(quiet(wald_cic(nod1, "y", "d", "g", "t"))$estimate, NA_real_))

## ---- input validation ---------------------------------------------------------
fails <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
bad_d <- toy; bad_d$d[1] <- 2
bad_g <- toy; bad_g$g[1] <- 2
bad_t <- toy; bad_t$t <- bad_t$t + 1
bad_na <- toy; bad_na$y[3] <- NA
for (f in list(wald_did, wald_tc, wald_cic)) {
  stopifnot(fails(f(bad_d, "y", "d", "g", "t")), fails(f(bad_g, "y", "d", "g", "t")),
            fails(f(bad_t, "y", "d", "g", "t")), fails(f(bad_na, "y", "d", "g", "t")),
            fails(f(toy, "nope", "d", "g", "t")))
}

## ---- Step 5: cross-check against Rfuzzydid ----------------------------------------
# Rfuzzydid (CRAN third-party R port of the Stata fuzzydid command; version
# 1.3.0 was already installed, this task installs nothing). Point estimates
# only (nose = TRUE). The hand-coded functions stay authoritative: they follow
# dCDH (2018) Sections 2.2-2.4 and 4 directly.
# Finding: Rfuzzydid 1.3.0 evaluates F^{-1}(k/n) like stats::quantile(type = 1)
# and so has the floating-point off-by-one fixed in .fz_qq above. Its Wald-CIC
# equals ours computed with that convention; the two differ only when
# n1 * k0 / n0 is an integer, i.e. mostly when a control cell has the same size
# in both periods (tessera_A: stable control group). Wald-DID and Wald-TC agree.
if (!requireNamespace("Rfuzzydid", quietly = TRUE)) {
  message("Rfuzzydid not installed: cross-check skipped")
} else {
  rfz <- function(X) {
    r <- suppressMessages(Rfuzzydid::fuzzydid(as.data.frame(X), y ~ d, group = "g", time = "t",
                                              did = TRUE, tc = TRUE, cic = TRUE, nose = TRUE))$late
    setNames(r$estimate, r$estimator)[c("W_DID", "W_TC", "W_CIC")]
  }
  ours <- function(X) c(W_DID = wald_did(X, "y", "d", "g", "t")$estimate,
                        W_TC  = wald_tc(X, "y", "d", "g", "t")$estimate,
                        W_CIC = wald_cic(X, "y", "d", "g", "t")$estimate)
  type1_qq <- function(yeval, y_d00, y_d01)      # the convention of the first version
    stats::quantile(y_d01, probs = stats::ecdf(y_d00)(yeval), type = 1, names = FALSE)
  cic_type1 <- function(X) {
    fixed <- .fz_qq; .fz_qq <<- type1_qq; on.exit(.fz_qq <<- fixed)
    wald_cic(X, "y", "d", "g", "t")$estimate
  }
  stopifnot(max(abs(rfz(toy) - c(15.0, 13.8, 14.8))) < 1e-6)
  # a cell design where the off-by-one bites: Q_0(7) = 107 (paper) vs 108 (type 1)
  obo <- rbind(data.frame(g = 0, t = 0, d = 0, y = 1:25),        data.frame(g = 0, t = 0, d = 1, y = 51:55),
               data.frame(g = 0, t = 1, d = 0, y = 101:125),     data.frame(g = 0, t = 1, d = 1, y = 161:165),
               data.frame(g = 1, t = 0, d = 0, y = rep(7, 20)),  data.frame(g = 1, t = 0, d = 1, y = c(53, 53)),
               data.frame(g = 1, t = 1, d = rep(0:1, c(10, 12)), y = c(rep(110, 10), rep(170, 12))))
  # by hand: E(Y11) = 3140/22, E(Q_D10(Y10)) = (20 * 107 + 2 * 163)/22, E(D11) - E(D10) = 10/22
  stopifnot(abs(wald_cic(obo, "y", "d", "g", "t")$estimate - (3140 - 20 * 107 - 2 * 163) / 10) < 1e-10)   # 67.4
  r_obo <- rfz(obo)
  stopifnot(max(abs(r_obo[1:2] - ours(obo)[1:2])) < 1e-6, abs(r_obo[["W_CIC"]] - cic_type1(obo)) < 1e-6)  # Rfuzzydid: 65.4
  for (f in c("tessera_A.rds", "tessera_A_rising.rds")) {
    X <- readRDS(f); r <- rfz(X); o <- ours(X)
    stopifnot(max(abs(r[1:2] - o[1:2])) < 1e-6,                    # Wald-DID, Wald-TC
              abs(r[["W_CIC"]] - cic_type1(X)) < 1e-6,              # Wald-CIC, Rfuzzydid's convention
              abs(r[["W_CIC"]] - o[["W_CIC"]]) < 0.01)              # convention gap is small
  }
}

cat("test_fuzzy_fns: OK\n")
