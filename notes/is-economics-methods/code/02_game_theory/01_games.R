# Chapter 2 running example: games between Aster and Boreal
# Every number can be computed by hand; this script transcribes the hand calculations and checks them.
# Run from this directory.

## 6.1 Pricing game: dominance and Nash equilibrium ------------------------
# Payoff matrix A[Aster action, Boreal action], entries are Aster's payoffs; symmetric game, B = t(A)
A <- matrix(c(10, 2, 14, 6), nrow = 2, byrow = TRUE,
            dimnames = list(c("High", "Low"), c("High", "Low")))
High_dominated <- all(A["Low", ] > A["High", ])   # Low strictly dominates High
cat("High strictly dominated by Low (Aster):", High_dominated, " -> NE=(Low,Low)=(6,6)\n")

## 6.2 Compatibility game: pure and mixed strategy equilibria --------------
# (S1,S1)=(3,1) (S2,S2)=(1,3) mismatch=(0,0); A prefers S1, B prefers S2
# Pure strategy NE: check each cell for a profitable unilateral deviation
# Mixed: A chooses S1 with prob p, B chooses S1 with prob q
p <- 3/4        # from B's indifference p*1 = (1-p)*3
q <- 1/4        # from A's indifference 3q = 1-q
coord_prob <- p*q + (1-p)*(1-q)
EU_A <- 3*q     # A's expected payoff (both pure strategies pay the same in equilibrium)
cat(sprintf("mixed NE: p=%.2f q=%.2f coordination prob=%.3f A expected=%.2f\n",
            p, q, coord_prob, EU_A))
stopifnot(abs(coord_prob - 3/8) < 1e-12)

## 6.3 Sequential game: backward induction and SPE -------------------------
# A moves first, B observes and matches (matching pays 1 or 3, both > 0); A foresees this -> picks S1, gets 3
aster_choice <- if (3 > 1) "S1" else "S2"
cat("SPE: Aster ->", aster_choice, "; Boreal matches -> payoffs (3, 1)\n")

## 6.4 Incomplete information: Bayesian Nash equilibrium -------------------
# Boreal: Normal type (prob beta) has dominant action Low; Committed type (1-beta) has dominant action High
# Aster's payoffs: against High, High=10/Low=7; against Low, High=3/Low=6
beta <- 2/5
EU_high <- beta*3 + (1 - beta)*10     # = 10 - 7 beta
EU_low  <- beta*6 + (1 - beta)*7      # = 7 - beta
aster_bne <- ifelse(EU_high >= EU_low, "High", "Low")
cat(sprintf("BNE beta=%.1f: EU_high=%.1f EU_low=%.1f -> Aster=%s (threshold beta*=0.5)\n",
            beta, EU_high, EU_low, aster_bne))

## 6.5 Repeated game: cooperation threshold --------------------------------
R <- 10; T <- 14; P <- 6
delta_star <- (T - R) / (T - P)
coop  <- function(d) R / (1 - d)
cheat <- function(d) T + d * P / (1 - d)
cat(sprintf("grim threshold delta*=%.2f; at the threshold coop=%.1f cheat=%.1f\n",
            delta_star, coop(delta_star), cheat(delta_star)))
stopifnot(abs(coop(delta_star) - cheat(delta_star)) < 1e-9)

## 4.1 Comparative statics checks ------------------------------------------
# Higher continuation value makes cooperation strictly easier.
deltas <- c(0.4, 0.5, 0.7, 0.9)
slack <- coop(deltas) - cheat(deltas)
cat("cooperation slack by delta:", paste(round(slack, 2), collapse=", "), "\n")
stopifnot(all(diff(slack) > 0), abs(slack[deltas == delta_star]) < 1e-9)

## 6.6 Bargaining: Nash solution and Rubinstein solution --------------------
# Surplus V = 100; Aster's share d1 + alpha * (V - d1 - d2)
V <- 100
nb <- function(d1, d2, alpha = 0.5) d1 + alpha*(V - d1 - d2)
cat(sprintf("Nash bargaining (Aster): sym=%.0f asym d=(20,0)=%.0f weighted alpha=0.6=%.0f\n",
            nb(0, 0), nb(20, 0), nb(20, 0, 0.6)))
stopifnot(nb(0, 0) == 50, nb(20, 0) == 60, abs(nb(20, 0, 0.6) - 68) < 1e-9)
# Rubinstein alternating offers: proposer 1/(1+delta), responder delta/(1+delta)
rub <- function(delta) c(mover = 1/(1+delta), responder = delta/(1+delta))
r80 <- round(100 * rub(0.80), 4); r99 <- round(100 * rub(0.99), 4)
cat(sprintf("Rubinstein x100: delta=0.80 -> %.4f / %.4f; delta=0.99 -> %.4f / %.4f\n",
            r80[1], r80[2], r99[1], r99[2]))
stopifnot(all(abs(r80 - c(55.5556, 44.4444)) < 1e-9),
          all(abs(r99 - c(50.2513, 49.7487)) < 1e-9))

## 6.7 Correlated equilibrium: a public coin --------------------------------
CE_A <- 0.5*3 + 0.5*1; CE_B <- 0.5*1 + 0.5*3; mixed_NE <- 3/4
cat(sprintf("public-coin CE: CE_A=%.2f CE_B=%.2f vs mixed NE=%.2f\n",
            CE_A, CE_B, mixed_NE))
stopifnot(CE_A == 2, CE_B == 2, mixed_NE == 0.75)

## 6.8 Global game: the unique tipping threshold ----------------------------
# Indifference of the marginal player: theta* + E[l] - 1 = 0 with E[l] = 1/2
theta_star <- 1 - 1/2
cat(sprintf("global-game threshold theta*=%.1f\n", theta_star))
stopifnot(theta_star == 0.5)
