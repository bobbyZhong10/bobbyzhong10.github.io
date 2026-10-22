# ------------------------------------------------------------------
# tests/test_dgp.R -- checks on the simulated Vireo data and truth03.rds written by 01_dgp.R.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("mech_fns.R")
tr <- readRDS("truth03.rds"); P <- tr$params
V1 <- readRDS("vireo_V1.rds")
fin <- function(x, k = 1L) is.numeric(x) && length(x) == k && all(is.finite(x))

## ---- truth03 fields named on the chapter card --------------------------------------------
stopifnot(is.list(tr$V1$F), identical(tr$V1$F$family, "truncnorm"), fin(tr$V1$F$mean), fin(tr$V1$F$sd),
          identical(c(tr$V1$F$lower, tr$V1$F$upper), c(0, 3)),
          is.numeric(tr$V1$v), length(tr$V1$v) == nrow(V1), identical(tr$V1$v, V1$v_true),
          is.list(tr$V2$theta_dist), identical(names(tr$V2$theta_dist), c("regular", "irregular")),
          is.list(tr$V4$prefs), identical(names(tr$V4$prefs), c("prop", "recv")))
# seeds: integers, all different
sd <- unlist(P$seed)
stopifnot(length(sd) == 5L, !anyDuplicated(sd), all(sd == round(sd)))

## ---- the generators reproduce the saved data from the stored seeds ------------------------
# (compared as data.frames: a data.table read back from .rds carries a new internal self-reference pointer)
set.seed(P$seed$V1); stopifnot(identical(as.data.frame(gen_V1(P$V1)), as.data.frame(V1)))
set.seed(P$seed$V4); g4 <- gen_V4(P$V4)
stopifnot(identical(g4$prefs, tr$V4$prefs), identical(g4$common, tr$V4$common))

## ---- V0: hand block (no random numbers) ----------------------------------------------------
stopifnot(is.list(tr$V0), fin(tr$V0$vcg_rev), fin(tr$V0$pers_x), fin(tr$V0$ms_p_trade_sb), fin(tr$V0$mh_rent))

## ---- V1: ad-exchange auctions ----------------------------------------------------------------
P1 <- P$V1; D1 <- value_dist(tr$V1$F); rs <- tr$V1$r_star
stopifnot(identical(P1$n_set, 2:6), P1$n_auc_per_n * length(P1$n_set) >= 5000,
          identical(names(V1), c("auction", "n", "bidder", "v_true", "b_fp", "b_fp_res")),
          nrow(V1) == P1$n_auc_per_n * sum(P1$n_set), length(unique(V1$auction)) == P1$n_auc_per_n * length(P1$n_set))
na <- V1[bidder == 1L, .N, by = n]
stopifnot(identical(na$n, P1$n_set), all(na$N == P1$n_auc_per_n),
          all(V1[, .(ok = .N == n[1] && identical(bidder, seq_len(n[1]))), by = auction]$ok))
stopifnot(all(V1$v_true >= 0 & V1$v_true <= 3), all(V1$b_fp < V1$v_true), all(V1$b_fp > 0),
          identical(is.na(V1$b_fp_res), V1$v_true < rs), all(V1$b_fp_res[!is.na(V1$b_fp_res)] >= rs),
          all(V1$b_fp_res[!is.na(V1$b_fp_res)] <= V1$v_true[!is.na(V1$b_fp_res)]))
# within each n, the equilibrium bid is increasing in the value; spot check against fpa_bid
for (k in P1$n_set) { s <- V1[n == k][order(v_true)]; stopifnot(all(diff(s$b_fp) > 0)) }
chk <- V1[c(1, 7, 333, 5000, 12345, 22222, 39999), ]
stopifnot(max(abs(chk$b_fp - mapply(function(v, k) fpa_bid(v, k, D1$cdf), chk$v_true, chk$n))) < 1e-9)
# truth: r* solves Psi(r*) = 0; expected revenues by numerical integration
stopifnot(abs(virtual_value(rs, D1)) < 1e-9, abs(rs - optimal_reserve(D1)) < 1e-12,
          identical(tr$V1$rev_true$n, P1$n_set), all(tr$V1$rev_true$reserve > tr$V1$rev_true$none),
          all(diff(tr$V1$rev_true$none) > 0), all(diff(tr$V1$rev_true$reserve - tr$V1$rev_true$none) < 0))
# revenue equivalence in the truth: E[b(v_(1))] = E[v_(2)] for each n
for (k in c(2, 4, 6)) {
  bf <- fpa_bid_fun(k, D1)
  efp <- integrate(function(v) bf(v) * k * D1$cdf(v)^(k - 1) * D1$pdf(v), 0, 3, rel.tol = 1e-10)$value
  stopifnot(abs(efp - tr$V1$rev_true$none[tr$V1$rev_true$n == k]) < 1e-6)
}
stopifnot(identical(tr$V1$v_quantile$tau, P1$taus), max(abs(D1$cdf(tr$V1$v_quantile$v) - P1$taus)) < 1e-12)

## ---- V2: SaaS menus ------------------------------------------------------------------------------
P2 <- P$V2; th <- tr$V2$grid
stopifnot(length(th) == P2$n_grid, th[1] == P2$theta_lo, tail(th, 1) == P2$theta_hi)
dr <- type_density(tr$V2$theta_dist$regular, th); di <- type_density(tr$V2$theta_dist$irregular, th)
stopifnot(all(dr > 0), all(di > 0))                                                   # I.1: positive density on the support
mr <- optimal_menu(th, dr, P2$cost2); mi <- optimal_menu(th, di, P2$cost2)
stopifnot(all(diff(mr$psi) > 0), any(diff(mi$psi) < 0))                               # regular vs needs ironing
stopifnot(identical(names(tr$V2$two_type), c("theta_L", "theta_H", "mu", "cost2")))

## ---- V3: LEN parameters ----------------------------------------------------------------------------
P3 <- P$V3
stopifnot(fin(P3$A), fin(P3$c2), fin(P3$sigma2), P3$A > 0, P3$c2 > 0, P3$sigma2 > 0,
          identical(names(P3$signals), c("report", "peer", "noise")),
          P3$signals$report$k != 0, P3$signals$peer$k == 0, P3$signals$peer$rho != 0,
          P3$signals$noise$k == 0, P3$signals$noise$rho == 0, P3$N_sim >= 1000, P3$R_sim >= 100)

## ---- V4: the 200 x 200 matching market ------------------------------------------------------------
P4 <- P$V4; pf <- tr$V4$prefs
stopifnot(P4$n_workers == 200L, P4$n_firms == 200L, length(pf$prop) == 200L, length(pf$recv) == 200L,
          all(vapply(pf$prop, function(l) identical(sort(l), 1:200), logical(1))),
          all(vapply(pf$recv, function(l) identical(sort(l), 1:200), logical(1))),
          length(tr$V4$common$firm) == 200L, length(tr$V4$common$worker) == 200L)
# the common component: firms with a higher common quality sit higher on workers' lists on average
avg_rank_firm <- rowMeans(sapply(pf$prop, function(l) match(1:200, l)))
stopifnot(cor(avg_rank_firm, tr$V4$common$firm) < -0.5)
stopifnot(P4$n_trunc_prop + P4$n_swap_prop >= 1000, P4$n_trunc_recv >= 1000)
