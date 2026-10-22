# ------------------------------------------------------------------
# 03_personalized.R -- block C2: personalised pricing of Cirrus Pro, the output test, the limits of price
# discrimination.
#   OLS of log WTP on the features (training half, conventional standard errors); segments by quantiles of the
#   predicted index; each segment priced at its profit-maximising price on the training half; profit, consumer
#   surplus, total surplus and the share served on the evaluation half (per merchant per month, USD). Regimes:
#   uniform (K = 1), coarse (K = 4), fine (K = 50, "personalised"), and the full ladder of K.
#   Fine against uniform, merchant by merchant: who pays less, who gains, who is newly served.
#   BBM (Bergemann, Brooks and Morris 2015): the consumer-optimal segmentation of the evaluation half's discrete
#   WTP distribution, by construction (bbm_segment), and the surplus triangle it spans.
#   Output test (Schmalensee 1981; Varian 1985): linear segment demands p = 12 - q and p = 8 - q with marginal cost
#   2 (C0's legacy 6.2 segments; both served at one price), and p = 12 - q with p = 4 - q (the weak segment is
#   not served at one price).
#   Monte Carlo: 200 fresh populations and splits (seed 4104 at the start of the block), regimes uniform, coarse,
#   fine; mean and 2.5% to 97.5% range across replications (Ruling P23).
# Reads cirrus_C2.rds, truth04.rds. Writes res_pd.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("pricing_fns.R")
C2 <- readRDS("cirrus_C2.rds")
tr <- readRDS("truth04.rds")
P2 <- tr$params$C2
regimes <- names(P2$K)                                                  # uniform, coarse, fine

## ---- the ladder of segmentations on the single draw ---------------------------------------------------------
run <- pd_run(C2, P2$K_ladder, P2$cost)
lad <- run$ladder
row_of <- function(K) which(lad$K == K)
get <- function(col) setNames(vapply(P2$K, function(K) lad[[col]][row_of(K)], numeric(1)), regimes)

wev <- C2$wtp[!C2$train]
p_u <- run$prices[[as.character(P2$K[["uniform"]])]]
p_f <- run$prices[[as.character(P2$K[["fine"]])]]
buy_u <- wev >= p_u; buy_f <- wev >= p_f
cs_u <- (wev - p_u) * buy_u; cs_f <- (wev - p_f) * buy_f
by_merchant <- list(pay_less = mean(p_f < p_u), pay_more = mean(p_f > p_u),
                    better_off = mean(cs_f > cs_u), worse_off = mean(cs_f < cs_u),
                    newly_served = mean(buy_f & !buy_u), no_longer_served = mean(buy_u & !buy_f))

cf <- summary(run$fit)$coefficients
ols <- data.frame(term = rownames(cf), estimate = cf[, "Estimate"], se = cf[, "Std. Error"],
                  truth = c(P2$beta[["const"]], P2$beta[["sales"]], P2$beta[["tenure"]], P2$beta[["electronics"]], P2$beta[["home"]]),
                  row.names = NULL)

seg_c <- run$segments[[as.character(P2$K[["coarse"]])]]
seg_f <- run$segments[[as.character(P2$K[["fine"]])]]
prices <- list(uniform = run$segments[["1"]]$price, coarse = seg_c$price, fine = seg_f$price,
               coarse_n_train = seg_c$n_train, fine_n_train = seg_f$n_train)

## ---- BBM on the evaluation half ---------------------------------------------------------------------------------
vals <- sort(unique(wev))
probs <- as.numeric(table(factor(wev, levels = vals))) / length(wev)
bb <- bbm_segment(vals, probs)
unif_rev <- vapply(vals, function(p) p * sum(probs[vals >= p]), numeric(1))
bbm <- list(uniform_profit = bb$uniform_profit, uniform_price = vals[which.max(unif_rev)], efficient = bb$efficient_total,
            profit = bb$profit, total = bb$total, cs_max = bb$efficient_total - bb$uniform_profit,
            n_segments = length(bb$segments$weight), n_values = length(vals),
            largest_segment_weight = max(bb$segments$weight), min_price = min(bb$segments$price),
            max_price = max(bb$segments$price))
# the one-price regime's price was set on the training half; on the evaluation half it earns this much less than the
# evaluation half's best single price (the triangle's floor), so its point lies just below the triangle
bbm$uniform_gap <- bbm$uniform_profit - lad$profit[row_of(P2$K[["uniform"]])]

## ---- the output test with linear segment demands ---------------------------------------------------------------
C0 <- tr$C0
linear <- linear_pd(a = C0$seg_a, b = C0$b, cost = C0$cost)
open_a <- c(C0$seg_a[1], 4)                                            # the weak segment: p = 4 - q
linear_open <- linear_pd(a = open_a, b = C0$b, cost = C0$cost)

## ---- Monte Carlo: fresh populations and splits ------------------------------------------------------------------
set.seed(tr$params$seed$C2_mc)
draws <- t(vapply(seq_len(P2$R_mc), function(r) {
  d <- gen_C2(P2)
  rr <- pd_run(d, unname(P2$K), P2$cost)
  L <- rr$ladder
  pr <- setNames(L$profit, regimes); cs <- setNames(L$cs, regimes); ts <- setNames(L$ts, regimes); sv <- setNames(L$served, regimes)
  c(profit_uniform = pr[["uniform"]], profit_gain_coarse = pr[["coarse"]] / pr[["uniform"]] - 1,
    profit_gain_fine = pr[["fine"]] / pr[["uniform"]] - 1, profit_fine_over_coarse = pr[["fine"]] / pr[["coarse"]] - 1,
    cs_change_coarse = cs[["coarse"]] / cs[["uniform"]] - 1, cs_change_fine = cs[["fine"]] / cs[["uniform"]] - 1,
    ts_change_coarse = ts[["coarse"]] / ts[["uniform"]] - 1, ts_change_fine = ts[["fine"]] / ts[["uniform"]] - 1,
    served_change_coarse = sv[["coarse"]] - sv[["uniform"]], served_change_fine = sv[["fine"]] - sv[["uniform"]],
    pay_less_fine = mean(rr$prices[[as.character(P2$K[["fine"]])]] < rr$prices[["1"]]),
    r2 = summary(rr$fit)$r.squared)
}, numeric(12)))
mc_summary <- data.frame(stat = colnames(draws), mean = colMeans(draws),
                         lo = apply(draws, 2, stats::quantile, 0.025, names = FALSE),
                         hi = apply(draws, 2, stats::quantile, 0.975, names = FALSE),
                         share_positive = colMeans(draws > 0), row.names = NULL)

res <- list(K = P2$K, profit = get("profit"), cs = get("cs"), ts = get("ts"), served = get("served"),
            mean_price_paid = get("mean_price_paid"), ladder = lad, ols = ols, r2 = summary(run$fit)$r.squared,
            n_train = sum(C2$train), n_eval = sum(!C2$train), prices = prices, by_merchant = by_merchant,
            efficient = mean(wev), bbm = bbm, linear = linear, linear_open = linear_open,
            linear_a = C0$seg_a, linear_open_a = open_a,
            mc = list(R = P2$R_mc, draws = draws, summary = mc_summary))
saveRDS(res, "res_pd.rds")
print(lad, digits = 4, row.names = FALSE)
print(mc_summary, digits = 3, row.names = FALSE)
cat(sprintf("BBM: uniform profit %.3f, efficient %.3f, consumer-optimal CS %.3f, %d segments\n",
            bbm$uniform_profit, bbm$efficient, bbm$cs_max, bbm$n_segments))
cat("Saved res_pd.rds\n")
