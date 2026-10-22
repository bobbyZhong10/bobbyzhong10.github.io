# ------------------------------------------------------------------
# 02b_contagion_r0.R -- block S1m: which coefficient does the seller equation's misspecification contaminate,
# and does that depend on the city plan's coupon-fee link?
#
# 02_system.R shows that in S1m (the seller equation omits coupon x buyers) 3SLS is inconsistent for the buyer
# coefficient n_sellers while 2SLS of the buyer equation is not. That contagion runs through the same
# corr(coupon, fee) = -0.7 as the 3SLS efficiency gain of S1: the coupon is the seller equation's excluded
# shifter whose moment condition is violated, and 3SLS passes the violation on to the buyer equation only for the
# coefficients that the coupon helps to identify. With corr(coupon, fee) = 0 the coupon says nothing about the
# first stage of n_sellers (3SLS = 2SLS for n_sellers, see pop_avar$se_ratio_independent in 01_dgp.R), so the
# n_sellers limit is clean and the contamination lands on the buyer equation's own coupon coefficient.
#
# This step draws two large S1m samples (one million markets each), one with corr(coupon, fee) = 0 and, for
# comparison, one with the city plan's -0.7 (a different draw from the plim_S1m of 01_dgp.R), and stores the
# 2SLS and 3SLS limits of the buyer coefficients n_sellers and coupon (truth 1.2 and 0.4), the 2SLS limit of the
# seller coefficient on buyers (3SLS = 2SLS there) and the over-identification statistic per market.
# Seeds (new, used nowhere else): 9113 for corr(coupon, fee) = 0 and 9114 for corr(coupon, fee) = -0.7.
#
# Output: res_contagion_r0.rds
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("iv_fns.R")

tr <- readRDS("truth09.rds")
P1m <- tr$params$S1m
N_large <- tr$params$N_large
seed_contagion <- list(r0 = 9113, r_minus0p7 = 9114)
stopifnot(!any(unlist(seed_contagion) %in% unlist(tr$params$seed)))          # new seeds, none reused

limits <- function(par, seed) {
  set.seed(seed)
  d <- gen_S1m(par, n = N_large)
  sy <- s1_system(d)
  tB <- tsls(sy$y$buyer, sy$X$buyer, sy$Z, vcov = "homo"); tS <- tsls(sy$y$seller, sy$X$seller, sy$Z, vcov = "homo")
  th <- three_sls(sy$y, sy$X, sy$Z)
  J <- hansen_j(sy$y$seller, sy$X$seller, sy$Z)
  list(r = par$r, seed = seed, N = nrow(d), corr_coupon_fee = cor(d$coupon, d$fee),
       n_sellers_tsls = tB$b[["n_sellers"]], n_sellers_three_sls = th$b[["buyer:n_sellers"]],
       coupon_tsls = tB$b[["coupon"]], coupon_three_sls = th$b[["buyer:coupon"]],
       seller_buyers_tsls = tS$b[["buyers"]], seller_buyers_three_sls = th$b[["seller:buyers"]],
       J_seller_per_n = J$J / nrow(d))
}
t0 <- proc.time()[["elapsed"]]
res <- list(seed = seed_contagion, N = N_large, truth = c(n_sellers = tr$S1$buyer[["n_sellers"]], coupon = tr$S1$buyer[["coupon"]]),
            r0 = limits(modifyList(P1m, list(r = 0)), seed_contagion$r0),
            r_minus0p7 = limits(P1m, seed_contagion$r_minus0p7))
saveRDS(res, "res_contagion_r0.rds")

for (k in c("r0", "r_minus0p7")) with(res[[k]], cat(sprintf(
  "S1m, corr(coupon, fee) = %.1f (sample %.3f, N = %d): n_sellers 2SLS %.4f, 3SLS %.4f | buyer coupon 2SLS %.4f, 3SLS %.4f | seller buyers 2SLS %.4f | J/n %.4f\n",
  r, corr_coupon_fee, N, n_sellers_tsls, n_sellers_three_sls, coupon_tsls, coupon_three_sls, seller_buyers_tsls, J_seller_per_n)))
cat(sprintf("Time %.1f s\nSaved res_contagion_r0.rds\n", proc.time()[["elapsed"]] - t0))
