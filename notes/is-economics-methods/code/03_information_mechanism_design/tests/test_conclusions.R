# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 3.
# The first block is the chapter card, verbatim (comments translated, Ruling P6). The blocks after it
# lock the conclusions the text will draw from the same results: revenue equivalence and the reserve
# gain for every bidder count, and their Monte Carlo ranges (Ruling P23); the GPV band; the shape of
# the regular and ironed menus; the three extra signals of the LEN contract; and the matching checks
# behind Dubins-Freedman and Roth (1982), including the brief's targeted receiver truncation.
# Run from the chapter code directory after 01 to 05 (tests/run_tests.R does this).
# ------------------------------------------------------------------

## ---- chapter card (verbatim) --------------------------------------------------------------
source("../../tools/R/chapter_kit.R")
a1 <- readRDS("res_auction.rds"); s2 <- readRDS("res_menu.rds"); k3 <- readRDS("res_len.rds"); m4 <- readRDS("res_match.rds")
stopifnot(abs(a1$rev_fp - a1$rev_sp) < 3 * a1$se_diff, a1$rev_reserve > a1$rev_sp)   # revenue equivalence; the reserve collects more
stopifnot(a1$gpv_med_relerr < 0.02)                                                  # interior quantiles 0.05 to 0.95
stopifnot(abs(s2$regular$q_top - s2$regular$q_eff_top) < 1e-8, all(diff(s2$regular$q) >= 0), s2$irregular$bunch_width > 0)
stopifnot(k3$b_with_signal > k3$b_without, k3$profit_with_signal > k3$profit_without)  # informativeness principle
stopifnot(m4$stable_both, m4$proposer_gain_any == 0, m4$receiver_gain_any > 0)

tr <- readRDS("truth03.rds")

## ---- V1: revenue equivalence and the reserve, bidder count by bidder count ---------------------
bn <- a1$by_n
stopifnot(identical(bn$n, 2:6), all(abs(bn$rev_fp - bn$rev_sp) < 3 * bn$se_diff),         # RE holds for every n
          all(abs(bn$rev_fp_res - bn$rev_sp_res) < 3 * bn$se_diff_res),                    # and with the reserve
          all(abs(bn$gain_reserve - (bn$rev_true_res - bn$rev_true)) < 3 * bn$se_gain),   # the gain matches the truth for every n
          bn$gain_reserve[1] > 10 * abs(bn$gain_reserve[5]),                                # the reserve matters in thin auctions
          all(diff(bn$rev_sp) > 0), all(diff(bn$shading) < 0))                             # more bidders: more revenue, less shading
# the sample means sit inside the Monte Carlo ranges around the truth, and the truth inside the ranges
mc <- a1$mc
stopifnot(mc$R >= 1000, covers(mc$diff_range, 0), covers(mc$diff_res_range, 0), mc$gain_range[[1]] > 0,
          abs(mc$diff_mean) < 3 * mc$diff_sd / sqrt(mc$R), mc$diff_sd / a1$se_diff > 0.8, mc$diff_sd / a1$se_diff < 1.25,
          covers(mc$rev_sp_range, tr$V1$rev_true_pooled), covers(mc$rev_fp_range, tr$V1$rev_true_pooled),
          covers(mc$rev_sp_res_range, tr$V1$rev_reserve_true_pooled), mc$share_within_3se > 0.98)
stopifnot(abs(a1$rev_sp - tr$V1$rev_true_pooled) < 3 * a1$se_rev_sp)
## ---- V1: GPV recovers the value quantiles; the Monte Carlo band covers the truth -----------------------------
gb <- mc$gpv_band
stopifnot(mean(gb$lo <= gb$v_true & gb$v_true <= gb$hi) > 0.9, max(abs(gb$mean - gb$v_true) / gb$v_true) < 0.02,
          # 0.03: just above the realised n = 2 median (0.0247, fewest bids); the Monte Carlo keeps no per-n error to bound it
          mc$gpv_relerr_range[[2]] < 0.02, all(a1$gpv_med_relerr_by_n < 0.03),
          all(a1$gpv$b < a1$gpv$v_true))                                                   # bids sit below values at every quantile

## ---- V2: the regular menu distorts below the top only; the irregular one bunches -------------------------
rg <- s2$regular; ir <- s2$irregular
stopifnot(all(rg$menu$q[rg$menu$theta < max(rg$menu$theta)] < rg$menu$q_eff[rg$menu$theta < max(rg$menu$theta)]),
          abs(rg$theta_exclusion - rg$theta_psi0) < 2 * diff(rg$menu$theta[1:2]), rg$bunch_width == 0,
          rg$ic_gain < 1e-9, ir$ic_gain < 1e-9, ir$ic_gain_naive > 1,                     # ironing restores incentive compatibility
          ir$profit_naive_if_truthful > ir$profit, ir$bunch_q > 0, ir$bunch_mass > 0.1,
          abs(ir$q_top - ir$q_eff_top) < 1e-8, all(diff(ir$q) >= 0))
# the printed bunch is the ironed one: read the served bunches off Psi_bar (maximal runs of one value with
# q > 0); on the widest, Psi_bar is the density-weighted mean of Psi and lies between Psi just outside and just
# inside each end; the expected virtual surplus with Psi_bar equals the profit; and the profit beats the
# monotone, IC menu that flattens Psi by its running maximum
mm <- ir$menu; Km <- nrow(mm); c2 <- tr$V2$cost2
bu <- Filter(function(ii) length(ii) > 1L && mm$q[ii[1]] > 0, split(seq_len(Km), cumsum(c(TRUE, diff(mm$psi_bar) != 0))))
ib <- bu[[which.max(vapply(bu, function(ii) diff(range(mm$theta[ii])), numeric(1)))]]
lo <- min(ib); hi <- max(ib); pb <- mm$psi_bar[lo]
cm <- pmax(cummax(mm$psi), 0) / c2; Ucm <- c(0, cumsum((cm[-1] + cm[-Km]) / 2 * diff(mm$theta)))
stopifnot(length(bu) == ir$n_bunches, ir$bunch_lo == mm$theta[lo], ir$bunch_hi == mm$theta[hi], ir$bunch_q == mm$q[lo],
          abs(pb - sum(mm$w[ib] * mm$psi[ib]) / sum(mm$w[ib])) < 1e-10,
          mm$psi[lo - 1L] <= pb, pb <= mm$psi[lo], mm$psi[hi] <= pb, pb <= mm$psi[hi + 1L],
          abs(sum(mm$w * (mm$psi_bar * mm$q - c2 * mm$q^2 / 2)) - ir$profit) < 1e-5,
          ir$profit > sum(mm$w * (mm$theta * cm - Ucm - c2 * cm^2 / 2)) + 1e-8)
tt <- s2$two_type
stopifnot(tt$profit < tt$profit_fb, tt$q_L < tt$q_L_fb, tt$q_H == tt$q_H_fb, tt$rent_H > 0)

## ---- V3: what an extra signal is worth depends on what it says about effort -------------------------------
cs <- k3$cases
stopifnot(identical(cs$case, c("without", "report", "peer", "noise")),
          cs$beta[cs$case == "report"] > cs$beta[cs$case == "without"],
          cs$beta[cs$case == "peer"] > cs$beta[cs$case == "without"],                         # relative performance evaluation
          abs(cs$beta[cs$case == "noise"] - cs$beta[cs$case == "without"]) < 1e-12,           # pure noise: zero weight
          abs(cs$b_signal[cs$case == "noise"]) < 1e-12, cs$b_signal[cs$case == "report"] < 0,
          cs$b_signal[cs$case == "peer"] < 0, all(abs(cs$profit - cs$beta / 2) < 1e-12))     # profit = beta / (2 c'') at the optimum
sm <- k3$sim
stopifnot(covers(sm$gain_range, k3$profit_with_signal - k3$profit_without),
          abs(sm$gain_mean - (k3$profit_with_signal - k3$profit_without)) < 3 * sm$gain_sd / sqrt(sm$R),
          abs(sm$ce_without_mean) < 0.01, abs(sm$ce_with_mean) < 0.01, sm$share_gain_positive > 0.9)

## ---- V4: DA is strategy-proof for proposers but not for receivers --------------------------------------
stopifnot(m4$n_prop_misreports >= 1000, m4$n_prop_changed > 0, m4$n_recv_random >= 1000,
          m4$receiver_gain_targeted == m4$n_diff_partner, m4$n_diff_partner > 0,
          m4$targeted_get_receiver_optimal == m4$n_firms,                                   # each truncating firm gets its receiver-optimal partner
          m4$receiver_gain_any == m4$receiver_gain_random + m4$receiver_gain_targeted,
          all(m4$worker_rank_P <= m4$worker_rank_R), all(m4$firm_rank_R <= m4$firm_rank_P),   # proposer-optimal, receiver-pessimal
          mean(m4$worker_rank_P) < mean(m4$worker_rank_R), mean(m4$firm_rank_R) < mean(m4$firm_rank_P))
# receiver_gain_any counts profitable reports, not firms. The firms that gain are exactly the firms whose
# partners differ in the two matchings: a truncated list leaves a matching that is stable under the true lists
# (or leaves the firm unmatched), so no truncation gives a firm more than its firm-optimal stable partner
rr <- m4$recv_random; rt <- m4$recv_targeted; gr <- rr[rr[, "gain"] > 0, , drop = FALSE]
diff_f <- which(m4$firm_rank_P != m4$firm_rank_R)                                     # firms with different partners
stopifnot(length(diff_f) == m4$n_diff_partner, m4$receiver_gain_n_firms == length(diff_f),
          m4$receiver_gain_random_n_firms == length(unique(gr[, "j"])), m4$receiver_gain_random_n_firms <= m4$receiver_gain_n_firms,
          setequal(union(gr[, "j"], rt[rt[, "gain"] > 0, "j"]), diff_f),
          all(gr[, "gain"] <= rt[match(gr[, "j"], rt[, "j"]), "gain"]))                 # random truncations never beat the targeted one
