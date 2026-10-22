# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Nimbus panels and truth10.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("panel_fns.R")
tr <- readRDS("truth10.rds"); P <- tr$params
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))
md5_lines <- function(lines) { f <- tempfile(); on.exit(unlink(f)); writeLines(lines, f); unname(tools::md5sum(f)) }

## ---- truth10 fields named on the chapter card ------------------------------------------
stopifnot(fin(tr$P1$beta), fin(tr$P2$beta), fin(tr$P3$beta), fin(tr$P4$gamma), fin(tr$P5$avg_slope),
          fin(tr$P5$fe_plim), fin(tr$P6$beta))
# every block has its own seed, all different; the inherited block keeps set.seed(737)
sd <- unlist(P$seed)
stopifnot(length(sd) >= 14, !anyDuplicated(sd), all(sd == round(sd)), P$seed$P1 == 737)

## ---- inherited block P1: legacy code byte-identical and first; data regenerate exactly ------
src <- readLines("01_dgp.R")
stopifnot(md5_lines(src[1:37]) == "c2322ed9a70d45e9f4f32ff4d60f8ee3",          # = the legacy 01_dgp.R
          md5_lines(readLines("02_static_fe_re.R")) == "e0fb652d40f925c1ba92a8ad80ffcca3",
          grep("set.seed", src)[1] == 18, grep("^gen_static <- function", src) == 21)
ex <- parse(text = src[21:28])                                                   # the legacy gen_static definition
e1 <- new.env(); eval(ex, e1)
set.seed(737); st <- e1$gen_static()
nb <- readRDS("nimbus.rds")
stopifnot(identical(st, nb$static), identical(nb$truth, list(beta = 0.5, ca = 0.7)),
          tr$P1$beta == 0.5, nrow(nb$static) == 600 * 6)
# P1 probability limits (closed form) agree with a large sample (seed P1_large)
L1 <- tr$P1$large
stopifnot(L1$N >= 1e5, abs(L1$pols - tr$P1$plim_pols) < 0.01, abs(L1$be - tr$P1$plim_be) < 0.01,
          abs(L1$fe - tr$P1$beta) < 0.01, abs(L1$re - tr$P1$plim_re) < 0.01, abs(L1$varphi - tr$P1$varphi_plim) < 0.01,
          abs(tr$P1$plim_pols - (0.5 + 0.7 / (0.7^2 + 1))) < 1e-12,
          tr$P1$plim_re > tr$P1$beta, tr$P1$plim_re < tr$P1$plim_pols)

## ---- the generators reproduce the saved data from the stored seeds -----------------------
blocks <- c(P2 = "nimbus_P2.rds", P3 = "nimbus_P3.rds", P4 = "nimbus_P4.rds", P5 = "nimbus_P5.rds",
            P5b = "nimbus_P5b.rds", P6 = "nimbus_P6.rds", P7 = "nimbus_P7.rds")
D <- lapply(blocks, readRDS)
for (b in names(blocks)) {
  set.seed(P$seed[[b]])
  gen <- get(paste0("gen_", if (b == "P5b") "P5" else b))
  stopifnot(isTRUE(all.equal(gen(P[[b]]), D[[b]], tolerance = 0, check.attributes = TRUE)))
}
for (b in names(blocks)) stopifnot(is.data.frame(D[[b]]), all(c("id", "t", "y") %in% names(D[[b]])) || b == "P2")

## ---- P2: same promotions, iid errors versus random-walk errors -------------------------------
d <- D$P2; P2 <- P$P2
stopifnot(nrow(d) == P2$N * P2$T, all(c("x", "y_iid", "y_rw", "alpha_true", "eps_iid_true", "eps_rw_true") %in% names(d)),
          tr$P2$beta == P2$beta, P$P2$R_mc >= 1000)
stopifnot(max(abs(d$y_iid - (P2$mu0 + P2$beta * d$x + d$alpha_true + d$eps_iid_true))) < 1e-12,
          max(abs(d$y_rw - (P2$mu0 + P2$beta * d$x + d$alpha_true + d$eps_rw_true))) < 1e-12)
de <- d$eps_rw_true - ave(d$eps_rw_true, d$id, FUN = function(v) c(0, v[-length(v)]))   # innovations e_it
stopifnot(abs(sd(de) - P2$sigma_e) < 0.05,
          var(d$eps_rw_true[d$t == P2$T]) / var(d$eps_rw_true[d$t == 1]) > P2$T - 2,      # variance grows with t
          abs(cor(d$eps_iid_true[d$t > 1], d$eps_iid_true[d$t < P2$T])) < 0.06,
          cor(d$x, d$alpha_true) > 0.2, min(d$x) > 0, min(d$y_iid) > 0, min(d$y_rw) > 0)

## ---- P3: feedback, promotions respond to last month's shock --------------------------------
d <- D$P3; P3 <- P$P3
stopifnot(nrow(d) == P3$N * P3$T, tr$P3$beta == P3$beta, P3$feedback < 0,
          max(abs(d$y - (P3$mu0 + P3$beta * d$x + d$alpha_true + d$eps_true))) < 1e-12)
o <- d$t < P3$T                                                      # Cov(X_{i,t+1}, eps_it) = feedback * sigma_e^2
xl <- d$x[which(o) + 1]
stopifnot(abs(cov(xl, d$eps_true[o]) - P3$feedback * P3$sigma_e^2) < 0.1, abs(cor(d$x, d$eps_true)) < 0.05)
# closed-form probability limits under feedback: FE bias shrinks like 1/T, FD bias does not
vw <- P3$feedback^2 * P3$sigma_e^2 + P3$x_sd^2
stopifnot(abs(tr$P3$fe_plim - (P3$beta - P3$feedback * P3$sigma_e^2 / (P3$T * vw))) < 1e-12,
          abs(tr$P3$fd_plim - (P3$beta - P3$feedback * P3$sigma_e^2 / (2 * vw))) < 1e-12,
          abs(tr$P3$large$fe - tr$P3$fe_plim) < 0.01, abs(tr$P3$large$fd - tr$P3$fd_plim) < 0.01)

## ---- P4: certified merchants (time-invariant W_i) with an exogenous time-varying instrument ----
d <- D$P4; P4 <- P$P4
stopifnot(nrow(d) == P4$N * P4$T, tr$P4$gamma == P4$gamma, setequal(unique(d$w), 0:1),
          all(tapply(d$w, d$id, function(v) length(unique(v))) == 1),                   # W_i does not vary over t
          all(tapply(d$x1, d$id, sd) > 0),
          max(abs(d$y - (P4$mu0 + P4$beta1 * d$x1 + P4$beta2 * d$x2 + P4$gamma * d$w + d$alpha_true + d$eps_true))) < 1e-12)
u <- d[d$t == 1, ]; x1b <- tapply(d$x1, d$id, mean)
stopifnot(mean(u$w) > 0.25, mean(u$w) < 0.6, cor(u$w, u$alpha_true) > 0.2,              # certification correlates with alpha_i
          abs(cor(d$x1, d$alpha_true)) < 0.05,                                          # X_1 is uncorrelated with alpha_i
          cor(u$w, x1b[as.character(u$id)]) > 0.3, cor(d$x2, d$alpha_true) > 0.2, min(d$y) > 0)

## ---- P5 and P5b: random slopes ---------------------------------------------------------------
P5 <- P$P5; P5b <- P$P5b
Es2 <- sum(P5$s_probs * P5$s_types^2); Vs2 <- sum(P5$s_probs * P5$s_types^4) - Es2^2
stopifnot(P5$design == "within_sd", P5b$design == "level", P5$avg_slope == P5b$avg_slope,
          tr$P5$avg_slope == P5$avg_slope, tr$P5b$avg_slope == P5b$avg_slope,
          abs(tr$P5$fe_plim - (P5$avg_slope + P5$kappa * Vs2 / Es2)) < 1e-12, tr$P5b$fe_plim == P5b$avg_slope)
for (b in c("P5", "P5b")) {
  d <- D[[b]]; u <- d[d$t == 1, ]
  stopifnot(nrow(d) == P[[b]]$N * P[[b]]$T, all(u$s_true %in% P[[b]]$s_types), min(d$y) > 0, min(d$x) > 0,
            max(abs(d$y - (P[[b]]$mu0 + d$b_true * d$x + d$alpha_true + d$eps_true))) < 1e-12)
}
u5 <- D$P5[D$P5$t == 1, ]; u5b <- D$P5b[D$P5b$t == 1, ]
stopifnot(cor(u5$b_true, u5$s_true^2) > 0.8, abs(cor(u5$b_true, u5$mu_true)) < 0.1,      # P5: slopes follow the within SD
          abs(cor(u5b$b_true, u5b$s_true)) < 0.1, cor(u5b$b_true, u5b$mu_true) > 0.8)    # P5b: slopes follow the level
stopifnot(abs(tr$P5$large$fe - tr$P5$fe_plim) < 0.01, abs(tr$P5$large$mg - tr$P5$avg_slope) < 0.01,
          abs(tr$P5b$large$fe - tr$P5$avg_slope) < 0.01)

## ---- P6: attrition on alpha_i versus attrition on the shock ------------------------------------
d <- D$P6; P6 <- P$P6
stopifnot(nrow(d) == P6$N * P6$T, tr$P6$beta == P6$beta, all(d$s_alpha[d$t == 1] == 1), all(d$s_eps[d$t == 1] == 1),
          max(abs(d$y - (P6$mu0 + P6$beta * d$x + d$alpha_true + d$eps_true))) < 1e-12)
absorbing <- function(s, id) all(tapply(s, id, function(v) all(diff(v) <= 0)))
stopifnot(absorbing(d$s_alpha, d$id), absorbing(d$s_eps, d$id))
stopifnot(all(d$y[d$s_eps == 1 & d$t > 1] >= P6$floor),                                  # stayers' months are above the floor
          all(tapply(seq_len(nrow(d)), d$id, function(k) {                               # exit month = first month below the floor
            v <- d$y[k] < P6$floor & d$t[k] > 1; s <- d$s_eps[k]
            if (!any(v)) all(s == 1) else all(s == as.integer(seq_along(k) < which(v)[1]))
          })))
ka <- mean(d$s_alpha); ke <- mean(d$s_eps)
stopifnot(ka > 0.6, ka < 0.85, ke > 0.6, ke < 0.85,
          mean(d$alpha_true[d$t == P6$T & d$s_alpha == 0]) < mean(d$alpha_true[d$t == P6$T & d$s_alpha == 1]) - 0.3,
          min(d$y) > 0, min(d$x) > 0)
stopifnot(abs(tr$P6$large$fe_alpha - P6$beta) < 0.01, tr$P6$large$fe_eps < P6$beta - 0.03)

## ---- P7: heteroskedastic, serially correlated errors; RE.1 holds, RE.3 fails ---------------------
d <- D$P7; P7 <- P$P7
stopifnot(nrow(d) == P7$N * P7$T, P7$R_mc >= 1000, P7$rho_e > 0, P7$alpha_slope > 0, P7$eps_slope > 0,
          max(abs(d$y - (P7$mu0 + P7$beta * d$x + d$alpha_true + d$eps_true))) < 1e-12)
xb7 <- tapply(d$x, d$id, mean); a7 <- d$alpha_true[d$t == 1]
stopifnot(abs(cor(a7, xb7)) < 0.1,                                                       # alpha_i mean-independent of X_i
          cor(abs(a7), abs(xb7 - P7$x_mean)) > 0.2)                                      # but its spread grows with |Xbar_i - 5|
