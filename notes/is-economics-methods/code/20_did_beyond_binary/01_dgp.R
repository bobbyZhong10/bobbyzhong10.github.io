# ------------------------------------------------------------------
# 01_dgp.R -- Tessera simulation world for chapter 20
#             (DiD II: fuzzy, continuous and switching treatments)
#
# Tessera is a (fictional) food-delivery platform. One world, three blocks:
#   A  fuzzy take-up   : a courier incentive program opens in some cities;
#                        couriers enroll voluntarily (two groups, two periods);
#                        variant A' (tessera_A_rising.rds) lets enrollment also
#                        rise in the control cities
#   B  continuous dose : cities get different per-order bonuses (with stayers);
#                        tessera_B_had.rds is the same world with no stayers
#   C  switching dose  : a weekly city-level bonus tier in {0, 1, 2} that is
#                        switched on, raised, lowered, paused and resumed
# Every number in chapter 20 comes from these data and from truth20.rds.
#
# Outputs (in this directory):
#   tessera_A.rds         data.table: unit city g t d y type_true tau_true
#   tessera_A_rising.rds  data.table: same columns (variant A')
#   tessera_B.rds         data.table: unit t post g dose d y tau_true
#   tessera_B_had.rds     data.table: same columns, no zero-dose unit
#   tessera_C.rds         data.table: unit t g dose y tau_true
#   truth20.rds           list(A, A_rising, B, C, params): true parameters
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md): block A's g, t, d are the chapter's S_i, T_i, D_i
# (2.2); block B's dose is D_i with values d (2.3); block C's g is G_i, the
# first period in which unit i's treatment changes (0 = never switches; 2.2
# rule 6), and the horizon l is the chapter's ell (2.4). The field name
# truth20$C$delta follows the Task 10 brief; in the chapter's mathematics the
# l-period effect is NOT written delta_l (2.4 rule 5: delta is reserved for
# Wald-TC's delta_d in chapter 20).
#
# Seeds: block A (and A') 2001, block-A population check 2004, block B 2002,
# block C 2003.
# ------------------------------------------------------------------
library(data.table)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(code_dir, "fuzzy_did_fns.R"))

## =================================================================
## Block A: fuzzy take-up (de Chaisemartin & D'Haultfoeuille 2018)
## =================================================================
## Units are couriers nested in cities (cluster = city). Two periods:
## t = 0 is the month before the program opens in the new cities, t = 1 the
## month after. Outcome y: weekly completed orders (four-week average,
## continuous, so the CIC support condition, Assumption 8, holds).
##
## Groups (g = S_i):
##   g = 1 (80 cities): the program opens city-wide at t = 1. Before that an
##     invite-only beta had enrolled about 8% of couriers (always-takers);
##     about 46% enroll when it opens (switchers).
##   g = 0 (80 cities): the program opened here two years earlier; about 30%
##     of couriers are enrolled and enrollment does not change between the two
##     months (D(0) = D(1) for every courier: a stable control group).
## In both groups enrollment is self-selected on the gain beta_i and on the
## level a_i through a latent index (Assumption 3: D = 1{V >= v_gt}):
##   V = rho_beta z_beta + kappa_a z_a + s xi + v_city,
## with population thresholds giving the shares pA1, pS1 (g = 1) and pA0 (g = 0).
##
## Potential outcomes:
##   Y_it(0) = mu0 + c_city + a_i + gamma1 * t + eta_{city,t} + eps_it
##   Y_it(1) = Y_it(0) + tau_it,   tau_it = beta_i + lambda * t
## lambda: at t = 1 the platform raises the per-order bonus for every enrolled
## courier in every city, so every enrolled courier's effect rises by lambda.
##
## Which dCDH (2018) assumptions hold (Review of Economic Studies numbering):
##   Assumptions 1, 3                                      hold
##   Assumption 2 (stable control-group treatment rate)    holds: pA0 at t = 0, 1
##   Assumption 4 (common trends in Y(0))                  holds (common gamma1)
##   Assumption 4' (common trends within D(0) = d)         holds: Y(0) trends are
##     common to all couriers; the Y(1) trend among D(0) = 1 units is
##     gamma1 + lambda in both groups
##   Assumption 5 (LATE of units treated at T = 0 stable)  FAILS: lambda > 0
##   Assumptions 7, 8 (CIC model, common support)          hold: shift model,
##     h_d(u, t) = u + gamma1 * t + lambda * t * d is common to both groups
## Hence (Theorems 2, 3) W_TC = W_CIC = Delta. With a stable control rate,
## Theorem 1 gives W_DID = Delta under Assumptions 1, 3, 4, 5; here Assumption 5
## fails and, from DID_Y = pS1 Delta + (pA1 - pA0) lambda and DID_D = pS1,
##   W_DID = Delta + (pA1 - pA0) lambda / pS1,
## negative here because the control group has more incumbent enrollees
## (pA0 > pA1) whose effect rises by lambda.
##
## Truth (target parameter of dCDH 2018, Section 2.1):
##   Delta = E(Y_11(1) - Y_11(0) | S),  S = {D(0) < D(1), G = 1},
## the LATE of treatment-group switchers at t = 1.
##   truth20$A$late      sample version: mean of tau_i1 over the g = 1
##                       switchers in tessera_A (finite-population target)
##   truth20$A$late_pop  population value (closed form below)
##
## Variant A' (tessera_A_rising.rds): the same couriers, outcome draws and
## treatment group; in the control cities the incumbent 30% stay enrolled and a
## further pS0_rising = 10% join at t = 1 after a recruitment drive aimed at
## busy couriers: the control never-takers with the highest level index
## W = z_a + sw_noise * xi_w switch, exactly pS0_rising of all control couriers
## (rank rule; self-selection on level). Assumption 4'
## still holds, Assumption 2 fails. dCDH (2018) Section 3.1, Theorem 4: Delta
## is then only partially identified (the TC bounds are infinite because the
## outcome is unbounded), and the Wald-TC estimand equals Delta plus
##   - [pA1 pS0 / (pA0 + pS0)] [E(Y(1)|S0,1) - E(Y(1)|A0,1)] / pS1
##   - [(1 - pA1) pS0 / (pN0 + pS0)] [E(Y(0)|N0,1) - E(Y(0)|S0,1)] / pS1
## (A0, S0, N0: control always-, switchers, never-takers). The control
## switchers' delta_0 change mixes the time trend with the departure of
## high-level couriers from the D = 0 cell, so Wald-TC (and Wald-CIC) are
## biased upwards. truth20$A_rising$late equals truth20$A$late (same
## treatment group); its large-sample W_TC, W_CIC, W_DID and the two bias
## terms are stored for the chapter.
PA <- list(
  n_city     = 80,     # cities per group
  n_min      = 90,     # couriers per city drawn uniformly from n_min:n_max
  n_max      = 150,
  mu0        = 70,     # mean weekly orders without the program
  sd_city    = 8,      # city level c_city ~ N(shift_g1 * g, sd_city^2)
  shift_g1   = 3,      # program cities are slightly busier (levels only)
  sd_a       = 12,     # courier level a_i
  gamma1     = 4,      # common seasonal rise from t = 0 to t = 1
  sd_eta     = 0.6,    # city-by-period shock eta_{city,t}
  sd_eps     = 4,      # transitory courier noise eps_it
  mu_beta    = 6,      # courier gain beta_i ~ N(mu_beta, sd_beta^2) at the t = 0 bonus
  sd_beta    = 5,
  lambda     = 8,      # rise in every enrolled courier's effect at t = 1
  pA1        = 0.08,   # g = 1: always-takers (beta testers)
  pS1        = 0.46,   # g = 1: switchers (enroll when the program opens)
  pA0        = 0.30,   # g = 0: incumbent enrollees, the same couriers at t = 0 and t = 1
  rho_beta   = 0.7,    # latent index V = rho_beta z_beta + kappa_a z_a + s xi + v_city
  kappa_a    = 0.3,    #   (z_beta, z_a standardized gain and level; s = sqrt(1 - rho^2 - kappa^2))
  sd_vcity   = 0.15,   #   city shift in take-up, v_city ~ N(0, sd_vcity^2)
  pS0_rising = 0.10,   # variant A': control-group switchers (share of all control couriers)
  sw_noise   = 0.5     # variant A': level index W = z_a + sw_noise * xi_w
)

## gen_A: one draw of couriers and outcomes; returns block A and variant A'
gen_A <- function(P, n_city) {
  cities <- data.table(city = seq_len(2 * n_city), g = rep(1:0, each = n_city))
  cities[, n := sample(P$n_min:P$n_max, .N, replace = TRUE)]
  cities[, c_lvl := rnorm(.N, 0, P$sd_city) + P$shift_g1 * g]
  cities[, v_c := rnorm(.N, 0, P$sd_vcity)]
  cities[, eta0 := rnorm(.N, 0, P$sd_eta)]
  cities[, eta1 := rnorm(.N, 0, P$sd_eta)]
  cr <- cities[rep(seq_len(.N), n)]
  cr[, unit := .I]
  N <- nrow(cr)
  cr[, z_a := rnorm(N)]
  cr[, a := P$sd_a * z_a]
  cr[, z_beta := rnorm(N)]
  cr[, beta := P$mu_beta + P$sd_beta * z_beta]
  s_xi <- sqrt(1 - P$rho_beta^2 - P$kappa_a^2)
  cr[, V := P$rho_beta * z_beta + P$kappa_a * z_a + s_xi * rnorm(N) + v_c]
  sdV <- sqrt(1 + P$sd_vcity^2)
  thr_A1 <- qnorm(1 - P$pA1) * sdV
  thr_S1 <- qnorm(1 - P$pA1 - P$pS1) * sdV
  thr_A0 <- qnorm(1 - P$pA0) * sdV
  cr[, type_true := ifelse(g == 1,
                           ifelse(V >= thr_A1, "always", ifelse(V >= thr_S1, "switcher", "never")),
                           ifelse(V >= thr_A0, "always", "never"))]
  # variant A': the control never-takers with the highest level index switch at t = 1
  cr[, W := z_a + P$sw_noise * rnorm(N)]
  n_sw0 <- round(P$pS0_rising * sum(cr$g == 0))
  cand <- cr[g == 0 & type_true == "never", unit[order(-W)]]
  cr[, type_rising := type_true]
  cr[unit %in% cand[seq_len(n_sw0)], type_rising := "switcher"]
  pan <- cr[rep(seq_len(N), 2)]
  pan[, t := rep(0:1, each = N)]
  pan[, y0 := P$mu0 + c_lvl + a + P$gamma1 * t + ifelse(t == 1, eta1, eta0) + rnorm(.N, 0, P$sd_eps)]
  pan[, tau_true := beta + P$lambda * t]
  out <- function(type_col) {
    x <- copy(pan)
    x[, type_true := get(type_col)]
    x[, d := as.integer(type_true == "always" | (type_true == "switcher" & t == 1))]
    x[, y := y0 + d * tau_true]
    setkey(x, unit, t)
    x[, .(unit, city, g, t, d, y, type_true, tau_true)]
  }
  list(A = out("type_true"), A_rising = out("type_rising"), y0 = pan[, .(unit, t, y0)])
}

## closed-form population values for block A (V standardized:
## Corr(z_beta, V) = rho_beta / sdV; truncated-normal means)
popA <- with(PA, {
  sdV <- sqrt(1 + sd_vcity^2); r <- rho_beta / sdV
  u1 <- qnorm(1 - pA1); l1 <- qnorm(1 - pA1 - pS1); u0 <- qnorm(1 - pA0)
  Ez_S1 <- r * (dnorm(l1) - dnorm(u1)) / (pnorm(u1) - pnorm(l1))  # E(z_beta | switcher, g = 1)
  Ez_A1 <- r * dnorm(u1) / (1 - pnorm(u1))                        # E(z_beta | always-taker, g = 1)
  Ez_A0 <- r * dnorm(u0) / (1 - pnorm(u0))                        # E(z_beta | always-taker, g = 0)
  late <- mu_beta + sd_beta * Ez_S1 + lambda                      # Delta
  list(late = late,
       w_did = late + (pA1 - pA0) * lambda / pS1,
       bias_did = (pA1 - pA0) * lambda / pS1,
       att_always = c(g1_t0 = mu_beta + sd_beta * Ez_A1, g1_t1 = mu_beta + sd_beta * Ez_A1 + lambda,
                      g0_t0 = mu_beta + sd_beta * Ez_A0, g0_t1 = mu_beta + sd_beta * Ez_A0 + lambda))
})

set.seed(2001)
gA <- gen_A(PA, PA$n_city)
A  <- gA$A
AR <- gA$A_rising

type_shares <- function(X, pop) {
  s <- X[t == 0, .N, by = .(g, type_true)][, share := N / sum(N), by = g]
  if (!is.null(pop)) s[, share_pop := pop[paste0(type_true, g)]]
  setorder(s, -g, type_true)
  s[]
}
shares_A  <- type_shares(A, with(PA, c(always1 = pA1, switcher1 = pS1, never1 = 1 - pA1 - pS1,
                                       always0 = pA0, never0 = 1 - pA0)))
shares_AR <- type_shares(AR, NULL)

## W_TC bias terms of variant A' from cell means of potential outcomes at t = 1
## (formula in the block header; S0, A0, N0 are the control-group types)
tc_bias_terms <- function(X, y0) {
  Z <- merge(X[t == 1], y0[t == 1], by = c("unit", "t"))
  sh <- function(gg, ty) Z[g == gg, mean(type_true == ty)]
  pA1 <- sh(1, "always"); pS1 <- sh(1, "switcher")
  pA0 <- sh(0, "always"); pS0 <- sh(0, "switcher"); pN0 <- sh(0, "never")
  m1 <- function(ty) Z[g == 0 & type_true == ty, mean(y0 + tau_true)]   # E(Y(1) | type, g = 0, t = 1)
  m0 <- function(ty) Z[g == 0 & type_true == ty, mean(y0)]              # E(Y(0) | type, g = 0, t = 1)
  c(switchers_vs_always_Y1 = -pA1 * pS0 / (pA0 + pS0) * (m1("switcher") - m1("always")) / pS1,
    never_vs_switchers_Y0  = -(1 - pA1) * pS0 / (pN0 + pS0) * (m0("never") - m0("switcher")) / pS1)
}

## population check on a large auxiliary draw (1000 cities per group, not saved):
## the hand-coded estimators against the closed form (block A) and the
## large-sample W_TC of variant A'
set.seed(2004)
gB <- gen_A(PA, 1000)
est_all <- function(X) c(w_did = wald_did(X, "y", "d", "g", "t")$estimate,
                         w_tc  = wald_tc(X, "y", "d", "g", "t")$estimate,
                         w_cic = wald_cic(X, "y", "d", "g", "t")$estimate)
late_of <- function(X) X[g == 1 & t == 1 & type_true == "switcher", mean(tau_true)]
pop_check_A  <- c(list(n_city_per_group = 1000L, n_couriers = uniqueN(gB$A$unit), late = late_of(gB$A)),
                  as.list(est_all(gB$A)))
pop_check_AR <- c(list(n_city_per_group = 1000L, n_couriers = uniqueN(gB$A_rising$unit), late = late_of(gB$A_rising),
                       late_control = gB$A_rising[g == 0 & t == 1 & type_true == "switcher", mean(tau_true)]),
                  as.list(est_all(gB$A_rising)))
tc_bias_AR <- tc_bias_terms(gB$A_rising, gB$y0)
rm(gB)

truth_A <- list(
  late           = late_of(A),
  late_pop       = popA$late,
  lambda         = PA$lambda,
  att_always     = A[type_true == "always", .(att = mean(tau_true)), keyby = .(g, t)],
  att_always_pop = popA$att_always,
  shares         = shares_A,
  w_did_pop      = popA$w_did,
  bias_did_pop   = popA$bias_did,
  pop_check      = pop_check_A,
  definition = paste(
    "late = mean of tau_true (= Y(1) - Y(0) at t = 1) over the g = 1 switchers (D(0) = 0, D(1) = 1);",
    "dCDH (2018) Delta = E(Y_11(1) - Y_11(0) | S), S = {D(0) < D(1), G = 1}; late_pop is its closed form.",
    "Stable control rate (Assumption 2) and Assumption 4' (7) give W_TC = W_CIC = Delta (Theorems 2, 3);",
    "Assumption 5 fails, so w_did_pop = Delta + (pA1 - pA0) lambda / pS1.")
)
truth_A_rising <- list(
  late          = late_of(AR),
  late_pop      = popA$late,
  late_control  = AR[g == 0 & t == 1 & type_true == "switcher", mean(tau_true)],
  shares        = shares_AR,
  pop_check     = pop_check_AR,
  tc_bias_parts = tc_bias_AR,
  definition = paste(
    "Variant A': same couriers, outcomes and treatment group as block A, so late = truth20$A$late;",
    "control-group take-up rises by pS0_rising (never-takers with the highest level index switch).",
    "Assumption 2 fails: dCDH (2018) Theorem 4 only bounds Delta (TC bounds infinite for an unbounded outcome).",
    "pop_check: large-sample W_DID, W_TC, W_CIC; tc_bias_parts: the two W_TC bias terms from large-sample cell means.")
)

## =================================================================
## Block B: continuous dose (Callaway, Goodman-Bacon & Sant'Anna;
##          dCDH et al. stayers / heterogeneous adoption designs)
## =================================================================
## Units are cities; two periods t = 1 (pre) and t = 2 (post); post = 1{t = 2}.
## Outcome y: city average of weekly orders per courier.
## At t = 2, 74% of cities get a per-order bonus dose ~ U(0.2, 3.0) (USD per
## order); 26% get none (stayers, dose = 0). Small cities are more likely to
## get none (selection on levels: harmless under parallel trends).
##
## Potential outcomes: Y_it(d) = alpha_i + gamma2 * post + post * theta_i(d) + eps_it,
##   theta_i(d) = b_i * m(d),  m(d) = a1 * d - a2 * d^2   (concave dose response)
##   b_i = (1 + kappa (D_i - dbar)) * exp(sd_nu nu_i - sd_nu^2 / 2) for treated
##   cities (dbar = 1.6 = E(D | D > 0)), b_i = exp(sd_nu nu_i - sd_nu^2 / 2) for
##   the rest. Selection on gains: E(b | D = d) = 1 + kappa (d - dbar), so the
##   platform gives larger bonuses where couriers respond more.
## Parallel trends holds (common gamma2), so ATT(d|d) is identified; strong
## parallel trends fails (E(b | D = d) varies with d), so the slope of
## ATT(d|d) = ACRT(d|d) + kappa * m(d) (selection bias), not ACRT(d|d).
##
## Truth (truth20$B):
##   att_curve  data.frame on dose = 0.20, 0.25, ..., 3.00:
##     att         ATT(d|d)  = E[Y_2(d) - Y_2(0) | D = d] = (1 + kappa (d - dbar)) m(d)
##     att_strong  ATT(d)    = E[Y_2(d) - Y_2(0) | D > 0] = m(d)
##     acrt        ACRT(d|d) = E[theta_i'(d) | D = d]     = (1 + kappa (d - dbar)) m'(d)
##     acrt_strong ACRT(d)   = E[theta_i'(d) | D > 0]     = m'(d)
##     slope_att   d ATT(d|d) / dd = acrt + sel_bias,  sel_bias = kappa m(d)
##     Eb          E(b | D = d)
##   att_overall       mean of theta_i(D_i) over treated cities (sample);
##                     population E[theta(D) | D > 0] in att_overall_pop
##   acrt_overall      mean of b_i m'(D_i) over treated cities (ACRT^o, sample)
##   acrt_overall_pt_pop  E[slope_att(D) | D > 0]: what a dose-slope estimator
##                     identifies under parallel trends only
##   as, was           Average (Weighted Average) of Switchers' Slopes targeted
##                     by did_multiplegt_stat on tessera_B:
##                     as = mean of theta_i(D_i) / D_i, was = sum theta_i(D_i) / sum D_i
##   had_was           sum theta_i(D_i) / sum D_i over all cities of tessera_B_had,
##                     the WAS targeted by did_had (no untreated city)
## tessera_B_had: same cities, same alpha_i, b_i and eps_it; the 26% stayers
## instead get small doses ~ U(0.005, 0.25) (quasi-untreated); the other
## cities keep their dose.
PB <- list(
  N           = 500,    # cities
  share_zero  = 0.26,   # stayers in tessera_B (exact share: the 130 lowest s_i)
  sel_level   = 0.5,    # s_i = sel_level z_alpha + sqrt(1 - sel_level^2) xi: stayers have low s_i
  d_lo        = 0.2,    # treated dose ~ U(d_lo, d_hi), USD per order
  d_hi        = 3.0,
  kappa       = 0.3,    # selection on gains: E(b | D = d) = 1 + kappa (d - dbar)
  sd_nu       = 0.3,    # idiosyncratic gain heterogeneity (mean-one lognormal factor)
  a1          = 4,      # m(d) = a1 d - a2 d^2
  a2          = 0.6,
  mu_alpha    = 70,     # city level alpha_i ~ N(mu_alpha, sd_alpha^2)
  sd_alpha    = 8,
  gamma2      = 3,      # common change from t = 1 to t = 2
  sd_eps      = 1.5,    # city-period noise
  had_lo      = 0.005,  # tessera_B_had: stayers' dose ~ U(had_lo, had_hi)
  had_hi      = 0.25
)
PB$dbar <- (PB$d_lo + PB$d_hi) / 2
m_B  <- function(d) PB$a1 * d - PB$a2 * d^2
mp_B <- function(d) PB$a1 - 2 * PB$a2 * d
Eb_B <- function(d) 1 + PB$kappa * (d - PB$dbar)

set.seed(2002)
uB <- data.table(unit = seq_len(PB$N), z_alpha = rnorm(PB$N), xi = rnorm(PB$N))
uB[, alpha := PB$mu_alpha + PB$sd_alpha * z_alpha]
uB[, s := PB$sel_level * z_alpha + sqrt(1 - PB$sel_level^2) * xi]
uB[, stayer := rank(s) <= round(PB$share_zero * PB$N)]
uB[, dose_B := ifelse(stayer, 0, runif(.N, PB$d_lo, PB$d_hi))]
uB[, b := ifelse(stayer, 1, Eb_B(dose_B)) * exp(PB$sd_nu * rnorm(.N) - PB$sd_nu^2 / 2)]
uB[, dose_H := ifelse(stayer, runif(.N, PB$had_lo, PB$had_hi), dose_B)]
uB[, eps1 := rnorm(.N, 0, PB$sd_eps)]
uB[, eps2 := rnorm(.N, 0, PB$sd_eps)]

make_B <- function(dose_col) {
  x <- uB[rep(seq_len(.N), 2)]
  x[, t := rep(1:2, each = PB$N)]
  x[, post := as.integer(t == 2)]
  x[, dose := get(dose_col)]                    # constant within unit (contdid convention)
  x[, g := ifelse(dose > 0, 2L, 0L)]            # timing group: first treated period, 0 = never
  x[, d := dose * post]                         # treatment actually received in period t
  x[, tau_true := post * b * m_B(dose)]
  x[, y := alpha + PB$gamma2 * post + ifelse(t == 1, eps1, eps2) + tau_true]
  setkey(x, unit, t)
  x[, .(unit, t, post, g, dose, d, y, tau_true)]
}
B <- make_B("dose_B")
H <- make_B("dose_H")

grid_B <- seq(PB$d_lo, PB$d_hi, by = 0.05)
att_curve <- data.frame(
  dose = grid_B,
  att = Eb_B(grid_B) * m_B(grid_B),
  att_strong = m_B(grid_B),
  acrt = Eb_B(grid_B) * mp_B(grid_B),
  acrt_strong = mp_B(grid_B),
  slope_att = Eb_B(grid_B) * mp_B(grid_B) + PB$kappa * m_B(grid_B),
  sel_bias = PB$kappa * m_B(grid_B),
  Eb = Eb_B(grid_B)
)
## population averages over D ~ U(d_lo, d_hi): exact midpoint-rule quadrature
dq <- PB$d_lo + (PB$d_hi - PB$d_lo) * (seq_len(200000) - 0.5) / 200000
avg_D <- function(f) mean(f(dq))
trB <- uB[stayer == FALSE]
trH <- uB
truth_B <- list(
  att_overall          = B[post == 1 & dose > 0, mean(tau_true)],
  att_overall_pop      = avg_D(function(d) Eb_B(d) * m_B(d)),
  att_curve            = att_curve,
  acrt_overall         = trB[, mean(b * mp_B(dose_B))],
  acrt_overall_pop     = avg_D(function(d) Eb_B(d) * mp_B(d)),
  acrt_overall_pt_pop  = avg_D(function(d) Eb_B(d) * mp_B(d) + PB$kappa * m_B(d)),
  as                   = trB[, mean(b * m_B(dose_B) / dose_B)],
  as_pop               = avg_D(function(d) Eb_B(d) * m_B(d) / d),
  was                  = trB[, sum(b * m_B(dose_B)) / sum(dose_B)],
  was_pop              = avg_D(function(d) Eb_B(d) * m_B(d)) / PB$dbar,
  had_was              = trH[, sum(b * m_B(dose_H)) / sum(dose_H)],
  share_zero           = mean(B[post == 1, dose] == 0),
  definition = paste(
    "ATT(d|d) = E[Y_2(d) - Y_2(0) | D = d]; ATT(d) = E[Y_2(d) - Y_2(0) | D > 0];",
    "ACRT(d|d) = E[dY_2(l)/dl at l = d | D = d]; att_overall = mean realized effect over treated cities;",
    "acrt_overall = mean of b_i m'(D_i) over treated cities; acrt_overall_pt_pop = E[d ATT(d|d)/dd at D | D > 0];",
    "as/was = AOSS/WAOSS of did_multiplegt_stat on tessera_B; had_was = WAS of did_had on tessera_B_had.")
)

## =================================================================
## Block C: switching, non-binary treatment (dCDH 2026, DID_l)
## =================================================================
## Units are cities, weeks t = 1..30. Treatment dose in {0, 1, 2}: a weekly
## peak-bonus tier (0 none, 1 = USD 1 per order, 2 = USD 2 per order).
## Every city starts at 0 (status quo D_{g,1} = 0). 25% of cities never get a
## bonus; the others switch on at week g ~ U{6, ..., 20} (the chapter's G_i) to
## tier 1 (prob 0.6) or 2 (0.4). Afterwards the tier is reviewed weekly
## following a Markov chain (rows: from 0, 1, 2; columns: to 0, 1, 2), so it
## can be paused, resumed, raised and lowered.
##
## Potential outcomes for a treatment path d = (d_1, ..., d_30):
##   Y_gt(d) = alpha_g + trend * t + season * sin(2 pi t / 13)
##             + b_g * sum_{k=0}^{3} omega_k f(d_{t-k}) + eps_gt,
##   f(0) = 0, f(1) = 1, f(2) = 1.6 (diminishing returns in the tier),
##   omega = (1, 0.6, 0.4, 0.2): the effect builds up over four weeks of
##   exposure and decays over three weeks after a pause (carryover);
##   b_g = 2 + 2 (20 - g) / 14 + N(0, 0.5^2): early switchers gain more;
##   never-switchers b_g = 3 + N(0, 0.5^2) (irrelevant: never treated).
## tau_true = Y_gt(actual path) - Y_gt(0, ..., 0), the effect of the actual
## path relative to the status-quo path (no anticipation: 0 before week g).
## Static TWFE (y on dose with unit and week FE) is biased by the dynamics
## (paused cells with carryover act as controls), by the heterogeneity across
## switching cohorts, and by the concave tier response.
##
## Truth (truth20$C). Because every status quo is 0, every switcher is a
## switcher-in (S_g = +1) and no cell has both a larger and a lower dose than
## the status quo, so did_multiplegt_dyn drops no cell and DID_l estimates the
## effect of having been exposed to a weakly higher dose for l periods.
## With never-switchers observed until week 30 and g <= 23, every switcher
## enters every horizon l = 1..8, and
##   delta[l]      = mean over switchers of tau_true at week g - 1 + l
##                   (the effect DID_l estimates; equals did_multiplegt_dyn on
##                   the noise-free column tau_true, checked in tests/test_dgp.R)
##   dose_cum[l]   = mean over switchers of sum_{k=0}^{l-1} D_{g, g+k}
##   delta_norm[l] = delta[l] / dose_cum[l]   (normalized effect DID^n_l)
##   avg_total_effect = sum_l delta[l] / sum_l mean(D_{g, g-1+l})
##                   (the Av_tot_eff row of did_multiplegt_dyn with effects = 8)
##   avg_effect    = sum of tau_true / sum of dose over all cells with dose > 0:
##                   the true effect per unit of dose over treated cells
##                   (dose-weighted, the benchmark for the TWFE dose coefficient);
##   avg_effect_unweighted = mean of tau_true / dose over the same cells.
PC <- list(
  N          = 200,     # cities
  nT         = 30,      # weeks
  share_never = 0.25,   # never-switchers
  F_lo       = 6,       # first switch week g ~ U{F_lo, ..., F_hi}
  F_hi       = 20,
  p_first2   = 0.4,     # probability that the first non-zero tier is 2
  trans      = rbind(c(0.70, 0.25, 0.05),    # from 0 (paused): stay, resume at 1, resume at 2
                     c(0.10, 0.80, 0.10),    # from 1: pause, stay, raise
                     c(0.08, 0.12, 0.80)),   # from 2: pause, lower, stay
  f          = c(0, 1, 1.6),   # f(0), f(1), f(2)
  omega      = c(1.0, 0.6, 0.4, 0.2),
  b_lo       = 2,       # b_g = b_lo + (b_hi - b_lo) (F_hi - g) / (F_hi - F_lo) + N(0, sd_b^2)
  b_hi       = 4,
  sd_b       = 0.5,
  mu_alpha   = 70,      # alpha_g ~ N(mu_alpha, sd_alpha^2)
  sd_alpha   = 8,
  trend      = 0.1,     # common weekly trend
  season     = 1,       # amplitude of the 13-week cycle
  sd_eps     = 2        # city-week noise
)

set.seed(2003)
uC <- data.table(unit = seq_len(PC$N))
uC[, never := unit %in% sample(PC$N, round(PC$share_never * PC$N))]
uC[, g := ifelse(never, 0L, sample(PC$F_lo:PC$F_hi, .N, replace = TRUE))]
uC[, b := ifelse(never, (PC$b_lo + PC$b_hi) / 2,
                 PC$b_lo + (PC$b_hi - PC$b_lo) * (PC$F_hi - g) / (PC$F_hi - PC$F_lo)) + rnorm(.N, 0, PC$sd_b)]
uC[, alpha := rnorm(.N, PC$mu_alpha, PC$sd_alpha)]
Dm <- matrix(0L, PC$N, PC$nT)
for (i in which(!uC$never)) {
  Fg <- uC$g[i]
  Dm[i, Fg] <- if (runif(1) < PC$p_first2) 2L else 1L
  if (Fg < PC$nT) for (tt in (Fg + 1):PC$nT)
    Dm[i, tt] <- sample(0:2, 1, prob = PC$trans[Dm[i, tt - 1] + 1, ])
}
fD <- matrix(PC$f[Dm + 1], PC$N, PC$nT)
tau_m <- matrix(0, PC$N, PC$nT)
for (k in seq_along(PC$omega)) {
  lag <- k - 1
  cols <- (lag + 1):PC$nT
  tau_m[, cols] <- tau_m[, cols] + PC$omega[k] * fD[, cols - lag, drop = FALSE]
}
tau_m <- tau_m * uC$b
C <- data.table(unit = rep(seq_len(PC$N), PC$nT), t = rep(seq_len(PC$nT), each = PC$N),
                dose = as.vector(Dm), tau_true = as.vector(tau_m))
C <- merge(C, uC[, .(unit, g, alpha)], by = "unit")
C[, y := alpha + PC$trend * t + PC$season * sin(2 * pi * t / 13) + tau_true + rnorm(.N, 0, PC$sd_eps)]
setkey(C, unit, t)
C <- C[, .(unit, t, g, dose, y, tau_true)]

L <- 8
sw_C <- C[g > 0]
stopifnot(all(uC[g > 0, g] - 1 + L <= PC$nT))                  # every switcher enters every horizon
delta_C    <- sapply(seq_len(L), function(l) sw_C[t == g - 1 + l, mean(tau_true)])
n_sw_C     <- sapply(seq_len(L), function(l) sw_C[t == g - 1 + l, .N])
dose_cum_C <- sapply(seq_len(L), function(l) sw_C[t >= g & t <= g - 1 + l, .(s = sum(dose)), by = unit][, mean(s)])
dose_inc_C <- sapply(seq_len(L), function(l) sw_C[t == g - 1 + l, mean(dose)])
truth_C <- list(
  delta                 = delta_C,
  delta_norm            = delta_C / dose_cum_C,
  dose_cum              = dose_cum_C,
  dose_inc              = dose_inc_C,
  n_switchers           = n_sw_C,
  avg_total_effect      = sum(delta_C) / sum(dose_inc_C),
  avg_effect            = C[dose > 0, sum(tau_true) / sum(dose)],
  avg_effect_unweighted = C[dose > 0, mean(tau_true / dose)],
  share_switch_2plus    = C[, .(n = sum(diff(dose) != 0)), by = unit][, mean(n >= 2)],
  definition = paste(
    "delta[l] = mean over switchers of Y_{g,G+l-1}(actual path) - Y_{g,G+l-1}(status quo 0), l = 1..8",
    "(all switchers are switchers-in and all enter every horizon); delta_norm = delta / dose_cum;",
    "avg_total_effect = sum(delta) / sum(dose_inc); avg_effect = sum(tau_true) / sum(dose) over cells with dose > 0.")
)

## =================================================================
## Save
## =================================================================
truth20 <- list(A = truth_A, A_rising = truth_A_rising, B = truth_B, C = truth_C,
                params = list(seed = c(A = 2001L, A_pop_check = 2004L, B = 2002L, C = 2003L),
                              A = PA, B = PB, C = PC))
saveRDS(A,  file.path(code_dir, "tessera_A.rds"))
saveRDS(AR, file.path(code_dir, "tessera_A_rising.rds"))
saveRDS(B,  file.path(code_dir, "tessera_B.rds"))
saveRDS(H,  file.path(code_dir, "tessera_B_had.rds"))
saveRDS(C,  file.path(code_dir, "tessera_C.rds"))
saveRDS(truth20, file.path(code_dir, "truth20.rds"))

## ---- Report -----------------------------------------------------
rate <- function(X, gg, tt) mean(X$d[X$g == gg & X$t == tt])
cat("== Tessera block A (fuzzy take-up, stable control group) ==\n")
cat(sprintf("Couriers %d in %d cities; take-up g=1: %.3f -> %.3f, g=0: %.3f -> %.3f\n",
            uniqueN(A$unit), uniqueN(A$city), rate(A, 1, 0), rate(A, 1, 1), rate(A, 0, 0), rate(A, 0, 1)))
cat(sprintf("True LATE of g=1 switchers: sample %.4f, population %.4f\n", truth_A$late, truth_A$late_pop))
e_A <- est_all(A)
cat(sprintf("Estimates on tessera_A: Wald-DID %.4f, Wald-TC %.4f, Wald-CIC %.4f\n", e_A[1], e_A[2], e_A[3]))
cat(sprintf("Population W_DID %.4f (bias (pA1 - pA0) lambda / pS1 = %.4f)\n", popA$w_did, popA$bias_did))
cat(sprintf("Large-sample check (%d couriers): late %.4f, W_DID %.4f, W_TC %.4f, W_CIC %.4f\n",
            pop_check_A$n_couriers, pop_check_A$late, pop_check_A$w_did, pop_check_A$w_tc, pop_check_A$w_cic))
cat("\n== Variant A' (control-group take-up rises, switchers selected on level) ==\n")
cat(sprintf("take-up g=0: %.3f -> %.3f; true LATE %.4f (control switchers' LATE %.4f)\n",
            rate(AR, 0, 0), rate(AR, 0, 1), truth_A_rising$late, truth_A_rising$late_control))
e_AR <- est_all(AR)
cat(sprintf("Estimates on tessera_A_rising: Wald-DID %.4f, Wald-TC %.4f, Wald-CIC %.4f\n", e_AR[1], e_AR[2], e_AR[3]))
cat(sprintf("Large-sample check: late %.4f, W_DID %.4f, W_TC %.4f, W_CIC %.4f; W_TC bias terms %.4f + %.4f\n",
            pop_check_AR$late, pop_check_AR$w_did, pop_check_AR$w_tc, pop_check_AR$w_cic,
            tc_bias_AR[1], tc_bias_AR[2]))
cat("\n== Tessera block B (dose) ==\n")
cat(sprintf("Cities %d, stayers %.3f; ATT^o sample %.4f (pop %.4f); ACRT^o sample %.4f (pop %.4f); PT-only slope %.4f\n",
            PB$N, truth_B$share_zero, truth_B$att_overall, truth_B$att_overall_pop,
            truth_B$acrt_overall, truth_B$acrt_overall_pop, truth_B$acrt_overall_pt_pop))
cat(sprintf("AS %.4f (pop %.4f), WAS %.4f (pop %.4f), HAD WAS %.4f\n",
            truth_B$as, truth_B$as_pop, truth_B$was, truth_B$was_pop, truth_B$had_was))
cat("\n== Tessera block C (switching) ==\n")
cat(sprintf("Cities %d, weeks %d; never-switchers %d; switching at least twice %.3f\n",
            PC$N, PC$nT, sum(uC$never), truth_C$share_switch_2plus))
cat("delta_l (l = 1..8):", sprintf("%.4f", truth_C$delta), "\n")
cat("normalized:        ", sprintf("%.4f", truth_C$delta_norm), "\n")
cat(sprintf("avg_effect %.4f (unweighted %.4f); avg_total_effect %.4f\n",
            truth_C$avg_effect, truth_C$avg_effect_unweighted, truth_C$avg_total_effect))
cat("Saved tessera_A.rds, tessera_A_rising.rds, tessera_B.rds, tessera_B_had.rds, tessera_C.rds, truth20.rds\n")
