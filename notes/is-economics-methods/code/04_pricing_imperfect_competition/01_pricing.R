# Chapter 4 running example: pricing and competition in the Cirrus market
# Every number can be computed by hand; this script transcribes the hand calculations and checks them.
# Run from this directory.

a <- 12; b <- 1; cost <- 2       # demand p = 12 - q, marginal cost 2

## 6.1 Single-price monopoly and the Lerner index --------------------------
qm <- (a - cost)/(2*b); pm <- a - b*qm; profit_m <- (pm - cost)*qm
qc <- (a - cost)/b                          # competitive quantity (p = c)
CS_m <- 0.5*(a - pm)*qm; DWL <- 0.5*(qc - qm)*(pm - cost)
lerner <- (pm - cost)/pm; inv_elas <- qm/pm  # |eps| = 1*p/q, inverse = q/p
cat(sprintf("monopoly: qm=%g pm=%g profit=%g CS=%g DWL=%g; Lerner=%.4f 1/|e|=%.4f\n",
            qm, pm, profit_m, CS_m, DWL, lerner, inv_elas))
stopifnot(abs(lerner - inv_elas) < 1e-9)

## 6.2 Third-degree price discrimination (inverse elasticity rule) ---------
q1 <- (12 - cost)/2; p1 <- 12 - q1
q2 <- (8  - cost)/2; p2 <- 8  - q2
cat(sprintf("third-degree PD: segment 1 p=%g |e|=%.4f markup=%.4f; segment 2 p=%g |e|=%.4f markup=%.4f\n",
            p1, p1/q1, (p1-cost)/p1, p2, p2/q2, (p2-cost)/p2))

## 6.3 Second-degree price discrimination (tiered menu = screening) --------
thH <- 10; thL <- 6; mu <- 0.5
# Single-package benchmark (within the two-type model): a single offer earns at most 25
single_both <- optimize(function(q) thL*q - q^2/2, c(0, 12), maximum = TRUE)$objective  # serve both types
single_high <- mu * (thH*thH - thH^2/2)      # serve only the high type: q=10, T=50, prob 1/2
uniform_bench <- max(single_both, single_high)
cat(sprintf("second-degree PD single-package benchmark: serve both=%.0f serve high only=%.0f -> ceiling=%.0f\n",
            single_both, single_high, uniform_bench))
# Uniform unit price (within the two-type model): each type buys q = theta - p
unit_rev <- function(p) p*(mu*pmax(thH - p, 0) + (1 - mu)*pmax(thL - p, 0))
unit_opt <- optimize(unit_rev, c(0, thH), maximum = TRUE)
cat(sprintf("second-degree PD uniform unit price: p=%.4f profit=%.4f\n", unit_opt$maximum, unit_opt$objective))
stopifnot(abs(unit_opt$maximum - 4) < 1e-3, abs(unit_opt$objective - 16) < 1e-6,
          unit_opt$objective < uniform_bench)
qH <- thH                                    # no distortion at the top
qL <- thL - (mu/(1 - mu))*(thH - thL)        # distortion at the bottom
rentH <- (thH - thL)*qL
TL <- thL*qL - qL^2/2                          # low type IR binds
TH <- (thH*qH - qH^2/2) - rentH                # high type IC binds
profit2 <- mu*TH + (1 - mu)*TL
profit_fb <- mu*thH^2/2 + (1 - mu)*thL^2/2     # full information (first degree)
uH_own <- thH*qH - qH^2/2 - TH
uH_mim <- thH*qL - qL^2/2 - TL
uL_own <- thL*qL - qL^2/2 - TL
cat(sprintf("second-degree PD: qL=%g qH=%g rentH=%g TL=%g TH=%g profit=%g (full information=%g)\n",
            qL, qH, rentH, TL, TH, profit2, profit_fb))
stopifnot(abs(uH_own - uH_mim) < 1e-9, abs(uL_own) < 1e-9, profit2 > uniform_bench)

## 6.4 Oligopoly competition -----------------------------------------------
q_c <- (a - cost)/(3*b); Q <- 2*q_c
p_cournot <- a - b*Q; profit_each <- (p_cournot - cost)*q_c
p_bertrand <- cost                             # undifferentiated price war
t <- 1; p_hotelling <- cost + t; profit_hot <- (p_hotelling - cost)*0.5
cat(sprintf("Cournot p=%.4f profit each=%.4f total=%.4f; Bertrand p=%g; Hotelling p=%g profit=%g\n",
            p_cournot, profit_each, 2*profit_each, p_bertrand, p_hotelling, profit_hot))
stopifnot(pm > p_cournot, p_cournot > p_bertrand)

## 6.5 Nonlinear pricing: two-part tariff and bundling ---------------------
c_ <- 2
q_star <- 12 - c_                 # efficient quantity at p = c
CS <- 0.5*(12 - c_)*q_star        # entire consumer surplus at this price
F  <- CS                          # fixed fee sweeps up all the surplus
profit_two_part <- F + (c_ - c_)*q_star
vA <- c(10, 2); vB <- c(2, 10)                          # two merchant types' valuations of A and B
revA <- max(2*sum(vA >= 2), 10*sum(vA >= 10))           # price 2 selling two vs price 10 selling one
revB <- max(2*sum(vB >= 2), 10*sum(vB >= 10))
rev_bundle <- 12*sum((vA + vB) >= 12)
cat(sprintf("two-part tariff: q=%g F=%g profit=%g (single price 25); bundling: separate=%g pure bundle=%g\n",
            q_star, F, profit_two_part, revA + revB, rev_bundle))
stopifnot(F == 50, 12*sum((vA + vB) >= 12) > revA + revB)

## 6.6 Rent path under Hotelling's rule ------------------------------------
r <- 0.05; rent0 <- 4; cext <- 2; t <- 0:2
rent <- rent0*(1 + r)^t; price <- cext + rent; pv <- rent/(1 + r)^t
cat("Hotelling rent:", paste(round(rent, 4), collapse = ", "),
    "price:", paste(round(price, 4), collapse = ", "),
    "present value:", paste(round(pv, 4), collapse = ", "), "\n")
stopifnot(all(abs(pv - 4) < 1e-9))

## 6.7 Vertical relations: double marginalization --------------------------
qm <- (12 - 2)/2; pm <- 12 - qm               # integrated monopoly
w_opt <- optimize(function(w) (w - 2)*(12 - w)/2, c(2, 12), maximum = TRUE)$maximum  # upstream chooses w
w <- 7; qd <- (12 - w)/2; pd <- (12 + w)/2    # separated: downstream MR=w -> q=(12-w)/2
joint_profit <- (w - 2)*qd + (pd - w)*qd
cat(sprintf("vertical: integrated p=%g profit=%g; separated w=%g (upstream optimum %.4f) retail price=%g joint profit=%g\n",
            pm, (pm - 2)*qm, w, w_opt, pd, joint_profit))
stopifnot(abs(w_opt - 7) < 1e-3, abs((w - 2)*qd + (pd - w)*qd - 18.75) < 1e-9)

## 4.8 Comparative statics checks ------------------------------------------
# Linear monopoly pass-through is 1/2; symmetric Cournot is n/(n+1).
n_grid <- 1:8
pass_through <- n_grid/(n_grid + 1)
cat("linear-demand pass-through, monopoly to Cournot n=8:",
    paste(round(pass_through, 3), collapse=", "), "\n")
stopifnot(abs(pass_through[1] - 0.5) < 1e-12, all(diff(pass_through) > 0))
