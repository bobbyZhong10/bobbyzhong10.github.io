# ------------------------------------------------------------------
# 05_matching.R -- block V4: deferred acceptance on Vireo's 200 x 200 talent market.
#
# 1. DA with workers proposing (worker-optimal stable matching) and with firms proposing (firm-optimal);
#    stability of both under the true lists; partner ranks on each side.
# 2. Manipulation checks (seed 3105), the other agents always truthful, outcomes judged by true lists:
#    workers (proposers): 1,000 random truncations and 1,000 random swaps of two entries; firms (receivers):
#    1,000 random truncations, and for every firm the standard misreport of Roth (1982) and Gale and Sotomayor
#    (1985): truncate the list just below the firm's partner in the firm-optimal stable matching.
#    A misreport is profitable when it leaves the agent with a partner higher on its true list.
# Output: res_match.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("mech_fns.R")
tr <- readRDS("truth03.rds"); P4 <- tr$params$V4; seed <- tr$params$seed
pp <- tr$V4$prefs$prop; pr <- tr$V4$prefs$recv
nw <- P4$n_workers; nf <- P4$n_firms

## ---- 1. both proposing directions ----------------------------------------------------------------
mP <- deferred_acceptance(pp, pr)                                 # worker -> firm, workers propose
mR <- invert_match(deferred_acceptance(pr, pp), nf)               # worker -> firm, firms propose
fP <- invert_match(mP, nf); fR <- invert_match(mR, nf)            # firm -> worker
wrP <- partner_rank(mP, pp); wrR <- partner_rank(mR, pp)
frP <- partner_rank(fP, pr); frR <- partner_rank(fR, pr)

## ---- 2. manipulation checks (seed 3105) --------------------------------------------------------------
set.seed(seed$V4_manip)
prop_dev <- function(i, li) {                                     # worker i reports li
  ppi <- pp; ppi[[i]] <- li
  m2 <- deferred_acceptance(ppi, pr)
  c(gain = wrP[i] - partner_rank(m2, pp)[i], changed = !identical(m2, mP))
}
recv_dev <- function(j, lj) {                                     # firm j reports lj
  prj <- pr; prj[[j]] <- lj
  f2 <- invert_match(deferred_acceptance(pp, prj), nf)
  c(gain = frP[j] - partner_rank(f2, pr)[j], to_receiver_optimal = !is.na(f2[j]) && f2[j] == fR[j])
}
prop_trunc <- t(vapply(seq_len(P4$n_trunc_prop), function(r) {
  i <- sample.int(nw, 1L); k <- sample.int(nf - 1L, 1L)
  c(i = i, k = k, prop_dev(i, truncate_pref(pp[[i]], k)))
}, numeric(4)))
prop_swap <- t(vapply(seq_len(P4$n_swap_prop), function(r) {
  i <- sample.int(nw, 1L); pos <- sample.int(nf, 2L)
  c(i = i, k = pos[1] * 1000 + pos[2], prop_dev(i, swap_pref(pp[[i]], pos)))
}, numeric(4)))
recv_rand <- t(vapply(seq_len(P4$n_trunc_recv), function(r) {
  j <- sample.int(nf, 1L); k <- sample.int(nw - 1L, 1L)
  c(j = j, k = k, recv_dev(j, truncate_pref(pr[[j]], k)))
}, numeric(4)))
recv_target <- t(vapply(seq_len(nf), function(j) {
  k <- match(fR[j], pr[[j]])
  c(j = j, k = k, recv_dev(j, truncate_pref(pr[[j]], k)))
}, numeric(4)))
prop_all <- rbind(prop_trunc, prop_swap)

## ---- save ----------------------------------------------------------------------------------------
m4 <- list(n_workers = nw, n_firms = nf, match_P = mP, match_R = mR,
           stable_P = is_stable(mP, pp, pr), stable_R = is_stable(mR, pp, pr),
           n_diff_partner = sum(mP != mR),
           worker_rank_P = wrP, worker_rank_R = wrR, firm_rank_P = frP, firm_rank_R = frR,
           share_worker_top1_P = mean(wrP == 1), share_worker_top1_R = mean(wrR == 1),
           share_firm_top1_P = mean(frP == 1), share_firm_top1_R = mean(frR == 1),
           n_prop_trunc = nrow(prop_trunc), n_prop_swap = nrow(prop_swap), n_prop_misreports = nrow(prop_all),
           proposer_gain_any = sum(prop_all[, "gain"] > 0), n_prop_changed = sum(prop_all[, "changed"] == 1),
           n_prop_worse = sum(prop_all[, "gain"] < 0), prop_max_gain = max(prop_all[, "gain"]),
           n_recv_random = nrow(recv_rand), receiver_gain_random = sum(recv_rand[, "gain"] > 0),
           receiver_gain_targeted = sum(recv_target[, "gain"] > 0),
           targeted_get_receiver_optimal = sum(recv_target[, "to_receiver_optimal"] == 1),
           targeted_mean_gain = mean(recv_target[recv_target[, "gain"] > 0, "gain"]),
           targeted_max_gain = max(recv_target[, "gain"]),
           random_recv_worse = sum(recv_rand[, "gain"] < 0),
           prop_misreports = prop_all, recv_random = recv_rand, recv_targeted = recv_target)
m4$stable_both <- m4$stable_P && m4$stable_R
# receiver_gain_any (the card's field) counts profitable firm reports over both checks, random plus targeted
# truncations, so a firm can count more than once; it is not a number of firms. receiver_gain_n_firms counts
# the distinct firms with at least one profitable report, receiver_gain_random_n_firms those found by the
# random truncations alone.
m4$receiver_gain_any <- m4$receiver_gain_random + m4$receiver_gain_targeted
firms_random <- unique(recv_rand[recv_rand[, "gain"] > 0, "j"])
m4$receiver_gain_n_firms <- length(union(firms_random, recv_target[recv_target[, "gain"] > 0, "j"]))
m4$receiver_gain_random_n_firms <- length(firms_random)
saveRDS(m4, "res_match.rds")

cat(sprintf("stable: workers propose %s, firms propose %s; %d of %d workers have different partners in the two matchings\n",
            m4$stable_P, m4$stable_R, m4$n_diff_partner, nw))
cat(sprintf("mean partner rank: workers %.3f (workers propose) vs %.3f (firms propose); firms %.3f vs %.3f\n",
            mean(wrP), mean(wrR), mean(frP), mean(frR)))
cat(sprintf("workers: %d misreports, %d change the matching, %d profitable, %d hurt the worker\n",
            m4$n_prop_misreports, m4$n_prop_changed, m4$proposer_gain_any, m4$n_prop_worse))
cat(sprintf("firms: %d random truncations, %d profitable; targeted truncation profitable for %d firms (mean gain %.2f places), %d of %d land on their firm-optimal partner\n",
            m4$n_recv_random, m4$receiver_gain_random, m4$receiver_gain_targeted, m4$targeted_mean_gain,
            m4$targeted_get_receiver_optimal, nf))
cat(sprintf("firms: %d profitable reports in all (receiver_gain_any) from %d distinct firms; the random ones from %d firms\n",
            m4$receiver_gain_any, m4$receiver_gain_n_firms, m4$receiver_gain_random_n_firms))
cat("Saved res_match.rds\n")
