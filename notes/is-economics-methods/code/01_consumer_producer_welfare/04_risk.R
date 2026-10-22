# ------------------------------------------------------------------
# 04_risk.R -- block W3 (deterministic): add-on protection, CARA risk aversion and a Rabin-type calibration.
#
# Protection: a covered loss of L = 150 USD happens with probability q = 0.08. Under CARA u(x) = -exp(-A x)
# (wealth plays no role) a subscriber buys the plan at price pi iff pi <= (1/A) log(1 - q + q exp(A L)), the
# certainty-equivalent cost of the uninsured risk. Declining the plan is taking the bet "lose L - pi with
# probability q, keep pi otherwise", so the indifference coefficient is
#   A*(pi) = cara_reject_threshold(L - pi, pi, prob = q):
# buying at pi reveals A >= A*(pi), declining reveals A <= A*(pi). The worked subscriber buys at 20 and
# declines at 24 USD: A in [A*(20), A*(24)].
# Rabin-type calibration: A_rabin = cara_reject_threshold(10, 11) rejects the 50-50 bet "lose 10, win 11".
# A CARA agent rejects the 50-50 bet "lose B, win G" for every G as soon as exp(A B) / 2 >= 1, that is
# B >= log(2) / A. So any A above log(2) / 100 (and A_rabin is above it) rejects every 50-50 bet that loses
# 100 USD, however much it wins; and the protection buyer's lower bound A*(20) is above it too.
# Reads truth01.rds; writes res_risk.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("welfare_fns.R")
P3 <- readRDS("truth01.rds")$params$W3

ce_cost <- function(A, L, q) log(1 - q + q * exp(A * L)) / A          # certainty-equivalent cost of the risk
A_star <- vapply(P3$prices, function(pr) cara_reject_threshold(P3$loss - pr, pr, prob = P3$prob), numeric(1))
names(A_star) <- P3$prices
A_interval <- c(A_star[[as.character(P3$buy_at)]], A_star[[as.character(P3$decline_at)]])
A_rabin <- cara_reject_threshold(P3$rabin_loss, P3$rabin_gain)
A_big <- log(2) / P3$big_loss                                          # rejects every 50-50 bet losing 100 above this
expected_loss <- P3$prob * P3$loss

res <- list(expected_loss = expected_loss, prices = P3$prices, A_star = A_star,
            markup_over_expected_loss = P3$prices / expected_loss - 1,
            risk_premium_at_A_star = vapply(seq_along(A_star), function(j) ce_cost(A_star[j], P3$loss, P3$prob) - expected_loss, numeric(1)),
            A_interval = A_interval,
            A_rabin = A_rabin, A_big = A_big,
            loss_never_accepted_at_A_rabin = log(2) / A_rabin,          # 50-50 bets losing more than this are rejected whatever they win
            loss_never_accepted_at_A_lo = log(2) / A_interval[1],
            loss_never_accepted_at_A_hi = log(2) / A_interval[2],
            # the worked subscriber also rejects "lose 10, win 11" only if A >= A_rabin: true at the upper end of
            # the interval, not at the lower end
            rabin_bet_rejected_at_A_lo = as.numeric(A_interval[1] > A_rabin),
            rabin_bet_rejected_at_A_hi = as.numeric(A_interval[2] > A_rabin))
stopifnot(abs(res$risk_premium_at_A_star - (P3$prices - expected_loss)) < 1e-8)   # at A*(pi) the CE cost is exactly pi
saveRDS(res, "res_risk.rds")

cat(sprintf("W3: expected loss %.2f USD; A*(price): %s\n", expected_loss,
            paste(sprintf("%g -> %.6f", P3$prices, A_star), collapse = ", ")))
cat(sprintf("Buys at %g, declines at %g: A in [%.6f, %.6f]; A_rabin (lose 10, win 11) %.6f; log(2)/100 = %.6f\n",
            P3$buy_at, P3$decline_at, A_interval[1], A_interval[2], A_rabin, A_big))
cat(sprintf("At A_rabin every 50-50 bet losing more than %.4f USD is rejected whatever it wins\n", res$loss_never_accepted_at_A_rabin))
cat("Saved res_risk.rds\n")
