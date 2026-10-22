# ------------------------------------------------------------------
# 03_revealed_pref.R -- block W2: GARP and the Afriat efficiency index (CCEI) on the credit-bundle panel,
# with Bronars uniform random choosers on the same budgets as the power benchmark.
# For every user (and every Bronars chooser): garp_check on the 8 observations (prices p^t, chosen bundle
# x^t) and ccei. Pass rates and CCEI summaries by type. The power of the GARP test against uniform random
# choice is 1 - pass_bronars (Bronars 1987).
# Reads truth01.rds, plover_W2.rds, plover_W2_bronars.rds; writes res_rp.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("welfare_fns.R")
tr <- readRDS("truth01.rds")
W2 <- readRDS("plover_W2.rds"); W2b <- readRDS("plover_W2_bronars.rds")

test_panel <- function(d, id) {
  out <- lapply(split(d, d[[id]]), function(g) {
    P <- cbind(g$p1, g$p2); X <- cbind(g$x1, g$x2)
    gc <- garp_check(P, X)
    data.frame(id = g[[id]][1], user = g$user[1], pass = gc$pass, n_viol = gc$n_viol, ccei = ccei(P, X))
  })
  do.call(rbind, out)
}
users <- test_panel(W2, "user")
users$type <- tr$W2$type[users$user]
bron <- test_panel(W2b, "chooser")

rat <- users$type == "rational"; noi <- users$type == "noisy"
q <- function(x, p) unname(stats::quantile(x, p, type = 7))
rp <- list(n_users = nrow(users), n_rational = sum(rat), n_noisy = sum(noi), n_bronars = nrow(bron),
           n_budgets = max(W2$t),
           pass_rational = mean(users$pass[rat]), ccei_rational_min = min(users$ccei[rat]),
           pass_noisy = mean(users$pass[noi]), ccei_noisy_mean = mean(users$ccei[noi]),
           ccei_noisy_median = stats::median(users$ccei[noi]), ccei_noisy_min = min(users$ccei[noi]),
           ccei_noisy_p10 = q(users$ccei[noi], 0.1), share_noisy_ccei_ge_095 = mean(users$ccei[noi] >= 0.95),
           pass_bronars = mean(bron$pass), power = 1 - mean(bron$pass),
           ccei_bronars_mean = mean(bron$ccei), ccei_bronars_median = stats::median(bron$ccei),
           ccei_bronars_min = min(bron$ccei), ccei_bronars_p10 = q(bron$ccei, 0.1),
           share_bronars_ccei_ge_095 = mean(bron$ccei >= 0.95),
           pass_all_users = mean(users$pass), n_viol_noisy_mean_failers = mean(users$n_viol[noi & !users$pass]),
           users = users, bronars = bron)
saveRDS(rp, "res_rp.rds")

cat(sprintf("W2: %d users (%d rational, %d noisy) x %d budgets; %d Bronars choosers\n",
            rp$n_users, rp$n_rational, rp$n_noisy, rp$n_budgets, rp$n_bronars))
cat(sprintf("GARP pass rate: rational %.4f (min CCEI %.10f), noisy %.4f (median CCEI %.4f, min %.4f), Bronars %.4f (median CCEI %.4f)\n",
            rp$pass_rational, rp$ccei_rational_min, rp$pass_noisy, rp$ccei_noisy_median, rp$ccei_noisy_min,
            rp$pass_bronars, rp$ccei_bronars_median))
cat("Saved res_rp.rds\n")
