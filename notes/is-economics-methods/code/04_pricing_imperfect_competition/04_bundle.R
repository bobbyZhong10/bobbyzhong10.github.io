# ------------------------------------------------------------------
# 04_bundle.R -- block C3: pure bundling of J information goods with independent U[0, 1] values (Bakos and
# Brynjolfsson 1999). For every J in {1, 2, 5, 20, 100}: the optimal bundle price (bundle_opt, zero marginal
# cost), revenue and deadweight loss per consumer, and per good: revenue per good (against the mean value 1/2)
# and the deadweight loss as a share of the efficient surplus (the mean bundle value). Selling the goods
# separately earns 1/4 per good and loses 1/4 of the surplus whatever J (J = 1 is that benchmark).
# Reads cirrus_C3.rds, truth04.rds. Writes res_bundle.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("pricing_fns.R")
C3 <- readRDS("cirrus_C3.rds")
tr <- readRDS("truth04.rds")
J <- tr$params$C3$J

rows <- lapply(seq_along(J), function(k) {
  v <- C3[, k]
  b <- bundle_opt(v)
  buy <- v >= b$price
  eff <- mean(v)                                                       # efficient surplus per consumer (cost 0)
  data.frame(J = J[k], price = b$price, rev = b$rev, dwl = b$dwl, share_buy = mean(buy),
             cs = sum(v[buy] - b$price) / length(v), efficient = eff,
             price_per_good = b$price / J[k], rev_per_good = b$rev / J[k], dwl_share = b$dwl / eff,
             cs_per_good = sum(v[buy] - b$price) / length(v) / J[k], value_per_good = eff / J[k],
             sd_bundle_value = sd(v))
})
tab <- do.call(rbind, rows)
res <- c(as.list(tab), list(table = tab, separate = list(rev_per_good = tab$rev_per_good[1], dwl_share = tab$dwl_share[1])))
saveRDS(res, "res_bundle.rds")
print(tab, digits = 4, row.names = FALSE)
cat("Saved res_bundle.rds\n")
