# ------------------------------------------------------------------
# 01_dgp.R -- the Plover world for chapter 1 (consumer and producer theory, welfare, uncertainty)
#
# Plover is a (fictional) content-subscription platform; its users are subscribers who split a monthly
# budget between Plover premium content (good 1, priced per hour) and everything else (good 2, the
# numeraire, p2 = 1). Money is in US dollars. Four blocks:
#   W0   one representative subscriber, by hand (inherited from the legacy chapter; no random numbers):
#        Cobb-Douglas a = 2/5, w = 100, p1 from 2 to 4; Slutsky, CV / EV / dCS, Gorman, the producer's
#        cost function, the square-root lottery, SOSD; plus the markup identity (new) and the Willig and
#        Vartia checks on the same numbers.
#   W1   premium price experiment (new): 10,000 subscribers with Stone-Geary preferences, heterogeneous
#        marginal budget shares beta_i and subsistence quantities (g1_i, g2_i), lognormal monthly budgets;
#        the premium price is randomised over {2, 2.5, 3, 3.5, 4} USD and usage is observed. Variant W1g:
#        the same subscribers with a common beta (Gorman holds).
#   W2   credit-bundle revealed-preference panel (new): 300 users x 8 random budgets; 70% choose the CES
#        optimum exactly, 30% stay on the budget line with multiplicative noise; Bronars uniform random
#        choosers on the same budgets as the power benchmark.
#   W3   add-on protection and risk aversion (new, deterministic): parameters only; 04_risk.R computes.
# Every number in chapter 1 comes from these data, truth01.rds and the res_*.rds files.
#
# Outputs (in this directory):
#   plover_W1.rds          data.frame: subscriber p1 w x1 beta_true g1_true g2_true
#   plover_W1g.rds         the same columns, beta_true = beta_mean for everyone
#   plover_W2.rds          data.frame: user t type_true alpha_true rho_true p1 p2 m x1 x2 x1_opt_true x2_opt_true
#   plover_W2_bronars.rds  data.frame: chooser user t p1 p2 m x1 x2
#   truth01.rds            list(W0, W1, W1g, W2, params)
# Columns ending in _true are simulation truth, not observable data.
#
# Notation (notation_v3.md 5.2): x1, x2, p1, p2 = 1, income w; CV, EV, dCS loss-positive (positive for a price
# rise); Stone-Geary beta_i, gamma_i = (g1_i, g2_i); eta = the income elasticity of good 1 (E_1w in the text);
# the representative consumer's CV is CV^rep, the subscriber average is the mean CV.
#
# Seeds (every simulated block its own set.seed at its start): W1 1101, W2 1102, W2 Bronars choosers 1103.
# W1g calls set.seed(1101) again on purpose: it is the same draw as W1 with beta_i replaced by beta_mean,
# so W1 and W1g differ only in the heterogeneity that breaks Gorman. The block parameters were calibrated
# on exploratory seeds (1 to 4, outside this directory) before these seeds were fixed; the seeds were not
# tuned afterwards.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("welfare_fns.R")

seed <- list(W1 = 1101, W2 = 1102, W2_bronars = 1103)

## =================================================================
## W0: one subscriber by hand (legacy 6.1 to 6.6 and 4.4; no random numbers)
## =================================================================
## Cobb-Douglas u = x1^(2/5) x2^(3/5), budget w = 100 USD a month, p2 = 1, premium price from 2 to 4 USD per hour.
P0 <- list(a = 2 / 5, w = 100, p2 = 1, p1_0 = 2, p1_1 = 4)
W0 <- local({
  a <- P0$a; w <- P0$w; p2 <- P0$p2; p1_0 <- P0$p1_0; p1_1 <- P0$p1_1
  x1 <- function(p1, m = w) a * m / p1
  x2 <- function(p2, m = w) (1 - a) * m / p2
  v <- function(p1, m = w) x1(p1, m)^a * x2(p2, m)^(1 - a)                    # indirect utility (level form)
  e <- function(p1, u) u * (p1 / a)^a * (p2 / (1 - a))^(1 - a)                # expenditure function
  h1 <- function(p1, u) a * u * (1 / a)^a * (p2 / (1 - a))^(1 - a) * p1^(a - 1)   # Shephard: de/dp1
  u0 <- v(p1_0); u1 <- v(p1_1)
  out <- list(a = a, w = w, p2 = p2, p1_0 = p1_0, p1_1 = p1_1,
              # 6.1 UMP
              x1_0 = x1(p1_0), x2_0 = x2(p2), x1_1 = x1(p1_1), x2_1 = x2(p2), u0 = u0, u1 = u1,
              spend_0 = p1_0 * x1(p1_0), spend_1 = p1_1 * x1(p1_1),
              # 6.2 EMP and Shephard; the compensated bundle at the new price (Hicks decomposition)
              e_check_0 = e(p1_0, u0), e_check_1 = e(p1_1, u1), h1_0 = h1(p1_0, u0),
              e_p1_u0 = e(p1_1, u0), h1_p1_u0 = h1(p1_1, u0), h2_p1_u0 = (1 - a) * e(p1_1, u0) / p2,
              # 6.3 Slutsky decomposition of dx1/dp1 at the initial point
              slutsky_total = -a * w / p1_0^2, slutsky_income = -x1(p1_0) * (a / p1_0),
              # 6.4 the three welfare numbers (loss positive)
              ev = w - e(p1_0, u1), dcs = a * w * log(p1_1 / p1_0), cv = e(p1_1, u0) - w)
  out$slutsky_subst <- out$slutsky_total - out$slutsky_income
  # Willig (eta = 1 for Cobb-Douglas) and Vartia on the same numbers (Ruling P13 check, see 02_welfare.R)
  out$eta <- 1
  out$willig <- willig_cv(out$dcs, out$eta, w)
  out$willig_ev <- willig_cv(out$dcs, -out$eta, w)
  out$vartia <- vartia_cv(function(p1, m) a * m / p1, p1_0, p1_1, w)
  # 6.6 Slutsky matrix at (p0, w) = ((2, 1), 100)
  x1s <- x1(p1_0); x2s <- x2(p2)
  S <- matrix(c(-a * w / p1_0^2 + x1s * (a / p1_0), 0 + x2s * (a / p1_0),
                0 + x1s * ((1 - a) / p2), -(1 - a) * w / p2^2 + x2s * ((1 - a) / p2)), nrow = 2, byrow = TRUE)
  out$S <- S
  out$S_sym <- S[1, 2] - S[2, 1]
  out$S_p1 <- (S %*% c(p1_0, p2))[1]
  out$S_eig <- zapsmall(eigen(S, symmetric = TRUE)$values)
  # 6.6 Gorman aggregation: two LES consumers, gamma^A = (5, 10), gamma^B = (15, 0), W = 300
  x1_les <- function(g1, g2, wi, b = 2 / 5) g1 + b * (wi - p1_0 * g1 - p2 * g2) / p1_0
  out$gorman_W <- 300
  out$x1_A_100 <- x1_les(5, 10, 100); out$x1_B_200 <- x1_les(15, 0, 200)
  out$x1_A_150 <- x1_les(5, 10, 150); out$x1_B_150 <- x1_les(15, 0, 150)
  out$X1 <- out$x1_A_100 + out$x1_B_200
  out$X1_redis <- out$x1_A_150 + out$x1_B_150
  out$X1_b07 <- x1_les(5, 10, 100, 0.4) + x1_les(15, 0, 200, 0.7)
  out$X1_b07_redis <- x1_les(5, 10, 150, 0.4) + x1_les(15, 0, 150, 0.7)
  # 6.5 producer: y = z1^(1/2) z2^(1/2), r = (4, 1), y = 10: c(r, y) = 2 y sqrt(r1 r2), Shephard z = dc/dr
  r1 <- 4; r2 <- 1; y <- 10
  out$r1 <- r1; out$r2 <- r2; out$y <- y
  out$cost <- 2 * y * sqrt(r1 * r2)
  out$z1 <- y * sqrt(r2 / r1); out$z2 <- y * sqrt(r1 / r2)
  out$cost_check <- r1 * out$z1 + r2 * out$z2
  # markup identity (new): output sold at p_y = 6 USD; z2 (freelance creator hours) is the flexible input.
  # mc = dc/dy; M = p / mc; the output elasticity of z2 is 1/2; its revenue share s_v = r2 z2 / (p y).
  # Cost minimisation r2 = mc df/dz2 gives M = E_yv / s_v.
  out$p_y <- 6                                                         # output price (w0_p_y)
  out$mc <- 2 * sqrt(r1 * r2)
  out$markup <- out$p_y / out$mc
  out$lerner <- (out$p_y - out$mc) / out$p_y
  out$elas_v <- 1 / 2
  out$share_v <- r2 * out$z2 / (out$p_y * y)
  out$markup_identity <- out$elas_v / out$share_v
  out$foc_check <- out$mc * 0.5 * y / out$z2                           # mc times the marginal product of z2 = r2
  # 6.5 uncertainty: u(x) = sqrt(x), 50-50 lottery {64, 144}
  out$lot_lo <- 64; out$lot_hi <- 144
  out$Ew <- 0.5 * 64 + 0.5 * 144
  out$Eu <- 0.5 * sqrt(64) + 0.5 * sqrt(144)
  out$ce <- out$Eu^2
  out$rp <- out$Ew - out$ce
  out$A_Ew <- 1 / (2 * out$Ew)
  out$varw <- 0.5 * (64 - out$Ew)^2 + 0.5 * (144 - out$Ew)^2
  out$rp_approx <- 0.5 * out$A_Ew * out$varw
  # 6.6 SOSD: F = {100, 108} against the mean-preserving spread G = {64, 144}
  out$mean_F <- 0.5 * 100 + 0.5 * 108; out$mean_G <- 0.5 * 64 + 0.5 * 144
  out$eu_F <- 0.5 * sqrt(100) + 0.5 * sqrt(108); out$eu_G <- 0.5 * sqrt(64) + 0.5 * sqrt(144)
  # 4.4 CRRA: local risk premium 0.5 (rho / E[w]) Var for a risk of variance 100
  out$crra_rho <- c(0.5, 1, 2); out$crra_var <- 100
  out$crra_rp_local <- 0.5 * (out$crra_rho / out$Ew) * out$crra_var
  # the two observed choices (p1, x) = (2, (20, 60)) and (4, (10, 60)) satisfy GARP (Afriat)
  out$garp_pass <- garp_check(rbind(c(p1_0, p2), c(p1_1, p2)), rbind(c(out$x1_0, out$x2_0), c(out$x1_1, out$x2_1)))$pass
  out
})
stopifnot(W0$ev < W0$dcs, W0$dcs < W0$cv, isTRUE(all.equal(W0$S[1, 2], W0$S[2, 1])), all(W0$S_eig <= 0),
          W0$X1 == W0$X1_redis, W0$eu_F > W0$eu_G, abs(W0$markup - W0$markup_identity) < 1e-12, abs(W0$foc_check - W0$r2) < 1e-12)

## =================================================================
## W1: premium price experiment (Stone-Geary subscribers)
## =================================================================
## 10,000 subscribers. Monthly budget w ~ lognormal, median 100 USD, sdlog 0.3 (as W0's
## w = 100 for the median subscriber). Marginal budget share of premium beta_i ~ Beta(7, 13): mean 0.35 (W0's
## a = 2/5 is the same order), sd 0.104, independent of w; subsistence premium hours g1_i ~ U(0.5, 1.5) and
## committed other spending g2_i ~ U(0, 10) USD, independent of everything. The ranges keep supernumerary
## income w - p1 g1 - g2 positive up to p1 = 12 USD (Ruling P7's out-of-support price) for every subscriber.
## The premium price is randomised: 2,000 subscribers at each of 2, 2.5, 3, 3.5, 4 USD per hour.
## Welfare question (as in W0): the loss from raising p1 from p0 = 2 to p1 = 4 USD; Ruling P7's out-of-support
## comparison: from 2 to p_out = 12 USD. W1g: the same draw with beta_i = 0.35 for everyone.
## beta independent of w is deliberate: the cross-sectional income slope of demand then equals the average
## individual income effect, so the analyst's Willig and Vartia corrections are not contaminated by
## composition; heterogeneity in beta still breaks Gorman (redistribution moves aggregate demand) and makes
## the representative consumer's CV miss the mean CV (Jensen gap; Hausman and Newey 2016: the cross-section
## cannot reveal it).
P1 <- list(N = 10000, p_grid = c(2, 2.5, 3, 3.5, 4), p2 = 1, p0 = 2, p1 = 4, p_out = 12,
           w_median = 100, w_sdlog = 0.3, beta_mean = 0.35, beta_phi = 20,
           g1_range = c(0.5, 1.5), g2_range = c(0, 10))
set.seed(seed$W1)
W1 <- gen_W1(P1)
set.seed(seed$W1)                                   # W1g: the same draw, common beta
W1g <- gen_W1(P1, common_beta = TRUE)

# Truth for one sample: per-subscriber welfare (closed forms), sample means, the population's representative
# consumer, and the per-capita demand curve.
truth_W1 <- function(d, P) {
  s <- sg_welfare(P$p0, P$p1, d$w, d$beta_true, d$g1_true, d$g2_true, P$p2)
  so <- sg_welfare(P$p0, P$p_out, d$w, d$beta_true, d$g1_true, d$g2_true, P$p2)
  willig_i <- willig_cv(s$dcs, s$eta, d$w)
  # Representative consumer of per-capita demand: the LES whose demand equals mean demand at mean income and
  # whose income effect is the mean individual income effect (b = mean beta). Its parameters solve
  # g1 (1 - b) = mean(g1_i (1 - beta_i)) and b (wbar - p2 g2) = mean(beta_i (w_i - p2 g2_i)).
  b <- mean(d$beta_true); wbar <- mean(d$w)
  g1r <- mean(d$g1_true * (1 - d$beta_true)) / (1 - b)
  g2r <- (wbar - mean(d$beta_true * (d$w - P$p2 * d$g2_true)) / b) / P$p2
  rep <- sg_welfare(P$p0, P$p1, wbar, b, g1r, g2r, P$p2)
  rep_out <- sg_welfare(P$p0, P$p_out, wbar, b, g1r, g2r, P$p2)
  x0 <- sg_demand(P$p0, d$w, d$beta_true, d$g1_true, d$g2_true, P$p2)
  pcov <- function(u, v) mean(u * v) - mean(u) * mean(v)
  pgrid <- seq(1.5, P$p_out, by = 0.05)
  list(cv_i = s$cv, ev_i = s$ev, dcs_i = s$dcs, eta_i = s$eta, willig_i = willig_i,
       cv12_i = so$cv,
       cv_mean = mean(s$cv), ev_mean = mean(s$ev), dcs_mean = mean(s$dcs), eta_mean = mean(s$eta),
       willig_mean = mean(willig_i), cv12_mean = mean(so$cv), ev12_mean = mean(so$ev), dcs12_mean = mean(so$dcs),
       cv_rep_pop = rep$cv, ev_rep_pop = rep$ev, cv12_rep_pop = rep_out$cv,
       rep = list(b = b, g1 = g1r, g2 = g2r, wbar = wbar),
       # second-order theory: mean CV - CV^rep ~ (1/2) (p1 - p0)^2 Cov(x1_i(p0), dx1_i/dw) with dx1_i/dw = beta_i / p0
       cov_x0_beta = pcov(x0, d$beta_true), cov_spend0_beta = pcov(P$p0 * x0, d$beta_true), cov_w_beta = pcov(d$w, d$beta_true),
       gap_second_order = 0.5 * (P$p1 - P$p0)^2 * pcov(x0, d$beta_true / P$p0),
       demand_curve = data.frame(p1 = pgrid, x1_mean = vapply(pgrid, function(p)
         mean(sg_demand(p, d$w, d$beta_true, d$g1_true, d$g2_true, P$p2)), numeric(1))))
}
truth_W1_list <- truth_W1(W1, P1)
truth_W1g_list <- truth_W1(W1g, P1)

## =================================================================
## W2: credit-bundle revealed-preference panel
## =================================================================
## 300 users (210 rational, 90 noisy, assigned at random), 8 monthly budgets each: a 60 USD credit and per-hour
## prices of premium (p1) and standard (p2) content drawn from U(1, 3) USD, so the budget lines cross often.
## CES utility with alpha ~ U(0.3, 0.7), rho ~ U(-1.5, 0.6) (elasticity of substitution 0.4 to 2.5). Noisy
## users: logit of the budget share of premium = optimal value + 0.5 N(0, 1) (the odds are multiplied by a
## lognormal factor). Bronars benchmark: 10 uniform random choosers on each user's budgets (3,000 in all).
P2 <- list(n_users = 300, n_budgets = 8, share_rational = 0.7, m = 60, price_range = c(1, 3),
           alpha_range = c(0.3, 0.7), rho_range = c(-1.5, 0.6), noise_sd = 0.5, n_bronars_per_user = 10)
set.seed(seed$W2)
W2 <- gen_W2(P2)
set.seed(seed$W2_bronars)
W2b <- gen_W2_bronars(W2, P2$n_bronars_per_user)
first <- W2[W2$t == 1, ]
truth_W2 <- list(type = first$type_true, alpha = first$alpha_true, rho = first$rho_true)

## =================================================================
## W3: add-on protection and risk aversion (deterministic)
## =================================================================
## Plover sells its reading tablet with an optional one-year protection plan. Replacement cost 150 USD,
## probability of a covered loss in the year 0.08 (expected loss 12 USD). Plan prices on offer: 15, 20, 24,
## 30 USD. The worked subscriber buys at 20 USD and declines at 24 USD. Rabin-type calibration: the 50-50 bet
## "lose 10, win 11" and bets that lose 100 USD.
P3 <- list(loss = 150, prob = 0.08, prices = c(15, 20, 24, 30), buy_at = 20, decline_at = 24,
           rabin_loss = 10, rabin_gain = 11, big_loss = 100)

## =================================================================
## Save
## =================================================================
saveRDS(W1, "plover_W1.rds")
saveRDS(W1g, "plover_W1g.rds")
saveRDS(W2, "plover_W2.rds")
saveRDS(W2b, "plover_W2_bronars.rds")
saveRDS(list(W0 = W0, W1 = truth_W1_list, W1g = truth_W1g_list, W2 = truth_W2,
             params = list(seed = seed, W0 = P0, W1 = P1, W2 = P2, W3 = P3)), "truth01.rds")

cat(sprintf("W0: EV %.4f  dCS %.4f  CV %.4f  Willig %.4f\n", W0$ev, W0$dcs, W0$cv, W0$willig))
cat(sprintf("W1:  mean CV %.4f  EV %.4f  dCS %.4f  CV^rep %.4f  CV 2 -> 12 %.4f\n", truth_W1_list$cv_mean,
            truth_W1_list$ev_mean, truth_W1_list$dcs_mean, truth_W1_list$cv_rep_pop, truth_W1_list$cv12_mean))
cat(sprintf("W1g: mean CV %.4f  EV %.4f  dCS %.4f  CV^rep %.4f\n", truth_W1g_list$cv_mean,
            truth_W1g_list$ev_mean, truth_W1g_list$dcs_mean, truth_W1g_list$cv_rep_pop))
cat(sprintf("W2: %d users x %d budgets (%d rational), %d Bronars choosers\n", P2$n_users, P2$n_budgets,
            sum(truth_W2$type == "rational"), max(W2b$chooser)))
cat("Saved plover_W1.rds, plover_W1g.rds, plover_W2.rds, plover_W2_bronars.rds, truth01.rds\n")
