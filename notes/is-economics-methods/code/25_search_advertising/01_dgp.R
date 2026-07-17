# ------------------------------------------------------------------
# 01_dgp.R -- "Lumen" search-platform log.
#
# N consumers each issue a query and see the SAME J products, but in
# a RANDOMIZED position order (the platform randomizes ranking). A
# consumer searches sequentially from the top, paying a search cost c
# to inspect each listing. Inspecting product j reveals its total
# utility u_ij = gamma - alpha * p_j + eps_ij, where eps ~ Gumbel(0,1),
# so a conditional logit over the TRUE choice set is correctly
# specified. Both price and the idiosyncratic match value are learned
# only on inspection. The consumer keeps the best option seen so far
# and stops optimally (Weitzman reservation rule): keep searching
# while the best inspected utility is below the reservation value z,
# where z solves  c = E[ max(u - z, 0) ].  At the stop the consumer
# buys the best inspected product if its utility exceeds the outside
# option (0), else buys nothing.
#
# Because position is randomized, WHICH products a consumer inspects
# is exogenous. This is the lever that (a) identifies the search cost
# c from how deep people search and (b) measures the causal value of
# prominence (top ranking), i.e. the "advertising" effect, and (c)
# lets a conditional logit on the INSPECTED set recover the true price
# coefficient, while the naive full-information logit (choice set =
# all J shown) is attenuated because consumers never inspect most of
# the cheap options lower down.
# ------------------------------------------------------------------
library(data.table)

set.seed(21)

## ---- dimensions --------------------------------------------------
N_cons <- 8000
J_prod <- 8

## ---- structural parameters --------------------------------------
gamma_true <- 2.90     # mean utility intercept
alpha_true <- 0.030    # price sensitivity (per dollar), logit scale
c_true     <- 0.40     # search cost per inspection (utils)
o0_true    <- 0.60     # outside-option location; u0_i = o0 + Gumbel(0,1)

## ---- product assortment (fixed prices) --------------------------
prices <- round(rnorm(J_prod, 100, 18), 1)
mean_u <- gamma_true - alpha_true * prices    # pre-noise mean utility

## ---- Gumbel(0,1) draw --------------------------------------------
rgum <- function(n) -log(-log(runif(n)))

## ---- reservation value z : solve c = E[max(u - z, 0)] -----------
## Expectation over a randomly drawn product's price and match shock.
Emax <- function(z, nsim = 300000) {
  pp <- sample(prices, nsim, replace = TRUE)
  uu <- gamma_true - alpha_true * pp + rgum(nsim)
  mean(pmax(uu - z, 0))
}
z_star <- uniroot(function(z) Emax(z) - c_true, c(-5, 20))$root

## ---- simulate search + purchase ---------------------------------
pos_mat  <- matrix(NA_integer_, N_cons, J_prod)   # product id at each position
insp_mat <- matrix(0L, N_cons, J_prod)            # inspected? by position
buy_prod <- integer(N_cons)                        # purchased product id (0 = none)
n_insp   <- integer(N_cons)

u0_all <- o0_true + rgum(N_cons)          # per-consumer outside option
for (i in 1:N_cons) {
  order_i <- sample(J_prod)              # randomized ranking: products by position
  pos_mat[i, ] <- order_i
  u0 <- u0_all[i]
  best_pu <- -Inf                        # best inspected product utility
  best_j  <- 0L
  k <- 0L
  for (r in 1:J_prod) {
    j <- order_i[r]
    u <- mean_u[j] + rgum(1)
    insp_mat[i, r] <- 1L
    k <- k + 1L
    if (u > best_pu) { best_pu <- u; best_j <- j }
    if (max(u0, best_pu) >= z_star) break   # stop once best open option clears z
  }
  n_insp[i]   <- k
  buy_prod[i] <- if (best_pu > u0) best_j else 0L   # buy only if best product beats outside
}

## ---- assemble a long consumer x position log --------------------
rows <- vector("list", N_cons)
for (i in 1:N_cons) {
  rows[[i]] <- data.table(
    cons = i,
    position  = 1:J_prod,
    product   = pos_mat[i, ],
    price     = prices[pos_mat[i, ]],
    inspected = insp_mat[i, ],
    bought    = as.integer(pos_mat[i, ] == buy_prod[i] & buy_prod[i] > 0)
  )
}
log_dt <- rbindlist(rows)

## ---- summary moments --------------------------------------------
buy_rate   <- mean(buy_prod > 0)
mean_depth <- mean(n_insp)
top_buy <- log_dt[position == 1, mean(bought)]
bot_buy <- log_dt[position == J_prod, mean(bought)]

cat("== Lumen search-platform DGP ==\n")
cat(sprintf("Consumers: %d, products: %d\n", N_cons, J_prod))
cat(sprintf("True gamma=%.2f, alpha=%.3f (logit scale), c=%.3f\n",
            gamma_true, alpha_true, c_true))
cat(sprintf("Reservation value z* = %.3f\n", z_star))
cat(sprintf("Prices: %s\n", paste(sort(prices), collapse = ", ")))
cat(sprintf("Buy rate: %.3f  |  mean products inspected: %.3f of %d\n",
            buy_rate, mean_depth, J_prod))
cat(sprintf("Prominence effect: P(buy|pos 1)=%.3f vs P(buy|pos %d)=%.3f (lift %.3f)\n",
            top_buy, J_prod, bot_buy, top_buy - bot_buy))

## ---- save --------------------------------------------------------
out_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
saveRDS(list(log = log_dt, n_insp = n_insp, buy_prod = buy_prod,
             pos_mat = pos_mat, prices = prices),
        file.path(out_dir, "lumen_log.rds"))
saveRDS(list(gamma = gamma_true, alpha = alpha_true, c = c_true, o0 = o0_true,
             z_star = z_star, prices = prices, buy_rate = buy_rate,
             mean_depth = mean_depth, top_buy = top_buy, bot_buy = bot_buy),
        file.path(out_dir, "truth25.rds"))
fwrite(log_dt, file.path(out_dir, "lumen_log.csv"))
cat("Saved lumen_log.rds / truth25.rds / lumen_log.csv\n")
