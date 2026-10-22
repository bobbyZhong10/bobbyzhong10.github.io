# Chapter 3 running example: Vireo's ad slot auction and other mechanisms
# Every number can be computed by hand; this script transcribes the hand calculations and checks them.
# Run from this directory.

## 6.1 Efficient allocation and VCG payments -------------------------------
v   <- c(10, 6, 4)           # value per click (sorted)
ctr <- c(100, 40, 0)         # slot click rates (3rd entry = no slot)
# Efficient allocation: advertiser i -> slot i
# VCG payment = others' max value in its absence - others' realized value under the chosen allocation
pay1 <- (v[2]*ctr[1] + v[3]*ctr[2]) - (v[2]*ctr[2] + v[3]*ctr[3])
pay2 <- (v[1]*ctr[1] + v[3]*ctr[2]) - (v[1]*ctr[1] + v[3]*ctr[3])
VCG_rev <- pay1 + pay2
util1 <- v[1]*ctr[1] - pay1
util2 <- v[2]*ctr[2] - pay2
cat(sprintf("VCG: pay1=%d (%.1f/click) pay2=%d (%.1f/click) revenue=%d\n",
            pay1, pay1/ctr[1], pay2, pay2/ctr[2], VCG_rev))
cat(sprintf("     util1=%d util2=%d (both nonnegative -> IR holds; truthful dominance follows from the VCG theorem of Section 3.6)\n", util1, util2))
stopifnot(VCG_rev == 680, util1 >= 0, util2 >= 0)

## 6.2 GSP comparison ------------------------------------------------------
# Truthful bids: the top slot pays the second bid, the second slot pays the third bid
gsp_truth <- v[2]*ctr[1] + v[3]*ctr[2]
# Local envy-freeness check: advertiser 2 swapping with the one above (top slot at price 6) gets 0 <= its current 80;
# advertiser 3 moving up to the second slot (at price 4) gets 0 <= 0.
stopifnot((v[2] - v[2])*ctr[1] <= (v[2] - v[3])*ctr[2], (v[3] - v[3])*ctr[2] <= 0)
# Downward deviation check: advertiser 1 dropping to the second slot (at price 4) gets 240 <= its current 400; advertiser 2 exiting gets 0 <= 80.
# An upward deviation costs at least the swap price, so the swap checks above cover it. Hence truthful bidding is a locally envy-free equilibrium here.
stopifnot((v[1] - v[3])*ctr[2] <= (v[1] - v[2])*ctr[1], 0 <= (v[2] - v[3])*ctr[2])
# The lowest revenue locally envy-free (LEF) equilibrium reproduces VCG:
# second slot price = v3, top slot price p1_low makes advertiser 2 exactly indifferent between the two slots
p1_low  <- (v[2]*(ctr[1] - ctr[2]) + v[3]*ctr[2])/ctr[1]
gsp_low <- p1_low*ctr[1] + v[3]*ctr[2]
stopifnot(abs(gsp_low - VCG_rev) < 1e-9)
cat(sprintf("GSP truthful revenue=%d (also a locally envy-free equilibrium here); lowest LEF revenue=%.0f, reproduces VCG=%d\n",
            gsp_truth, gsp_low, VCG_rev))

## 6.3 Revenue equivalence and reserve price -------------------------------
# Two advertisers with independent uniform[0,1] values
n <- 2
spa <- (n - 1)/(n + 1)                 # second price: expected second highest value
fpa <- (n - 1)/n * n/(n + 1)           # first price: shading b(v)=(n-1)/n v, expected highest bid
cat(sprintf("revenue equivalence: second price=%.4f first price=%.4f\n", spa, fpa))
stopifnot(abs(spa - fpa) < 1e-12)
# Optimal reserve r* = 1/2 (solves r = (1-F)/f = 1-r); expected revenue with n=2
r <- 0.5
# one bidder above r pays r; both above pay the second highest value E[min | both >= r]
rev_reserve <- 2*r*(1 - r)*r + (1 - r)^2*(r + (1 - r)/3)
stopifnot(abs(rev_reserve - 5/12) < 1e-12)
cat(sprintf("optimal reserve r*=0.5, expected revenue=%.4f (vs %.4f without a reserve)\n", rev_reserve, spa))

## 3.1 Akerlof adverse selection check -------------------------------------
# Quality theta~U[0,1], seller values it at theta, buyer at 1.5 theta
# At price p the sellers are those with theta<=p, mean quality p/2, buyer WTP = 1.5*(p/2) = 0.75p < p
p_grid <- seq(0.1, 1, by = 0.1)
wtp <- 1.5 * (p_grid/2)
cat("Akerlof: WTP < price at every price ->", all(wtp < p_grid), " market unravels to theta=0\n")

## 3.4 Moral hazard check --------------------------------------------------
# High effort costs 5; success probability 0.8 (0.2 under low effort); limited liability w>=0
# IC: 0.6 (wH - wL) >= 5 -> wH-wL >= 25/3; wL=0 -> wH=25/3; rent = 0.8 wH - 5
wH <- 25/3; wL <- 0
rent <- 0.8*wH - 5
cat(sprintf("moral hazard: wH-wL>=%.3f, wH=%.3f, agent rent=%.3f\n", 25/3, wH, rent))
stopifnot(abs(0.6*(wH - wL) - 5) < 1e-9)

## 4.6 Comparative statics checks ------------------------------------------
# First-price equilibrium b(v)=(n-1)v/n: shading falls as bidder count rises.
n_grid <- 2:8; v <- 1
shading <- v - ((n_grid - 1)/n_grid)*v
cat("first-price shading by bidder count:", paste(round(shading, 3), collapse=", "), "\n")
stopifnot(all(diff(shading) < 0))
