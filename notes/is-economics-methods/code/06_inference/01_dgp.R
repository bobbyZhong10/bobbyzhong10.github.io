# ============================================================================
# 01_dgp.R  --  Halcyon food-delivery platform: inference DGP (three pieces)
# ----------------------------------------------------------------------------
# "Halcyon" is a fictional food-delivery platform. Three data pieces, one seed,
# each built to expose a different inference problem.
#
#  (B) WTP / value-of-time.  Users rate a completed order; satisfaction falls
#      with delivery TIME and with the delivery FEE:
#         sat = b0 + b_time*time + b_fee*fee + noise,   b_time<0, b_fee<0.
#      The value of one saved minute in dollars is the RATIO
#         WTP = b_time / b_fee  (a nonlinear function of two estimates)
#      => Delta method vs bootstrap. A weak-denominator variant (fee barely
#      varies => b_fee imprecise) makes the ratio's sampling distribution skew,
#      so the symmetric Delta CI misleads and percentile-t is more honest.
#
#  (A) City-clustered rollout.  A city-level policy D_c is switched on in half
#      of G cities; the outcome is measured at the user-SESSION level. Sessions
#      in a city share a common city shock a_c, so naive i.i.d. SEs treat
#      correlated sessions as independent and badly OVERSTATE precision.
#      => cluster-robust "sandwich" SE, the Moulton factor, the few-cluster
#      failure, wild cluster bootstrap, randomization inference.
#
#  (C) Daily time series.  A platform metric with serial correlation; i.i.d.
#      SEs understate => Newey-West HAC.
#
# All TRUE parameters are known. Seed fixed. Scripts use relative paths; run
# with the working directory set to this folder.
# ============================================================================

set.seed(626)

## ---- (B) WTP / value-of-time -----------------------------------------------
b0 <- 5.0; b_time <- -0.08; b_fee <- -0.40      # WTP truth = b_time/b_fee = 0.20
WTP_true <- b_time / b_fee
gen_wtp <- function(n, fee_lo, fee_hi) {
  time <- runif(n, 10, 45)                      # delivery time (minutes)
  fee  <- runif(n, fee_lo, fee_hi)              # delivery fee (dollars)
  sat  <- b0 + b_time*time + b_fee*fee + rnorm(n, 0, 1.5)
  data.frame(sat, time, fee)
}
wtp_main <- gen_wtp(3000, 1.0, 7.0)             # fee varies a lot: strong denom
wtp_weak <- gen_wtp(3000, 3.05, 3.95)           # fee varies little: weak denom

## ---- (A) city-clustered rollout --------------------------------------------
make_clusters <- function(G, n_c, tau, sigma_a = 0.30, sigma_e = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  a  <- rnorm(G, 0, sigma_a)                    # city shock (intraclass corr)
  D  <- rep(0:1, length.out = G)[sample(G)]     # half cities treated
  city <- rep(seq_len(G), each = n_c)
  y  <- 2.0 + tau*D[city] + a[city] + rnorm(G*n_c, 0, sigma_e)
  data.frame(y, D = D[city], city = factor(city))
}
G_main <- 40; n_c <- 200; tau_true <- 0.25
clus_main <- make_clusters(G_main, n_c, tau_true, seed = 6261)
# intraclass correlation implied by the variance components
rho  <- 0.30^2 / (0.30^2 + 1.0^2)
moulton <- sqrt(1 + (n_c - 1) * rho)

## ---- (C) daily time series with AR(1) errors -------------------------------
Tn <- 400; rho_ar <- 0.6; beta_ts <- 0.30
xt <- as.numeric(arima.sim(list(ar = 0.5), n = Tn))
u  <- as.numeric(arima.sim(list(ar = rho_ar), n = Tn, sd = 1.0))
yt <- 1.0 + beta_ts * xt + u
ts_dat <- data.frame(yt, xt)

saveRDS(list(wtp_main = wtp_main, wtp_weak = wtp_weak,
             clus_main = clus_main, ts_dat = ts_dat,
             truth = list(WTP = WTP_true, b_time = b_time, b_fee = b_fee,
                          tau = tau_true, rho = rho, moulton = moulton,
                          beta_ts = beta_ts, G = G_main, n_c = n_c),
             mkclus = make_clusters), "halcyon.rds")

cat(sprintf("WTP truth (b_time/b_fee) = %.4f\n", WTP_true))
cat(sprintf("clustered rollout: G=%d cities, n_c=%d, tau=%.2f\n", G_main, n_c, tau_true))
cat(sprintf("intraclass corr rho = %.4f ; Moulton factor = %.3f\n", rho, moulton))
cat(sprintf("time series: T=%d, AR(1) error rho=%.2f, beta=%.2f\n", Tn, rho_ar, beta_ts))
