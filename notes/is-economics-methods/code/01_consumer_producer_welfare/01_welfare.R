# Chapter 1 running example: Lumen's Cobb-Douglas consumer and welfare measures
# Every number can be computed by hand; this script transcribes the hand calculations and checks them.
# Run from this directory.

a  <- 2/5      # expenditure share of good 1
w  <- 100      # budget
p2 <- 1        # price of the outside good
p1_0 <- 2      # initial price of premium content
p1_1 <- 4      # price after the increase

## 6.1 UMP: Walrasian demand and indirect utility --------------------------
x1 <- function(p1) a * w / p1
x2 <- function(p2) (1 - a) * w / p2
u  <- function(p1) x1(p1)^a * x2(p2)^(1 - a)   # indirect utility v(p,w)

u0 <- u(p1_0); u1 <- u(p1_1)
cat("x1_0 =", x1(p1_0), " x2 =", x2(p2), " u0 =", round(u0, 4), "\n")
cat("x1_1 =", x1(p1_1), " u1 =", round(u1, 4), "\n")
cat("spending p1*x1 before increase =", p1_0 * x1(p1_0), " after increase =", p1_1 * x1(p1_1), "\n")

## 6.2 EMP: expenditure function and Hicksian demand (Shephard check) ------
e <- function(p1, uu) uu * (p1 / a)^a * (p2 / (1 - a))^(1 - a)
cat("e(p0,u0) (should be 100) =", round(e(p1_0, u0), 6), "\n")
cat("e(p1,u1) (should be 100) =", round(e(p1_1, u1), 6), "\n")
# Shephard: h1 = de/dp1, which at (p0,u0) should equal the Walrasian x1=20
h1_shephard <- function(p1, uu) a * uu * (1 / a)^a * (p2 / (1 - a))^(1 - a) * p1^(a - 1)
cat("h1(p0,u0) =", round(h1_shephard(p1_0, u0), 4), " (should be 20)\n")

## 6.3 Slutsky decomposition (own price, initial point) --------------------
total  <- -a * w / p1_0^2               # dx1/dp1
income <- -x1(p1_0) * (a / p1_0)        # -x1 * dx1/dw
subst  <- total - income                # substitution effect
cat("Slutsky: total =", total, " subst =", subst, " income =", income, "\n")

## 6.4 Three welfare numbers -----------------------------------------------
dCS <- 40 * log(2)                      # integral of Walrasian demand
EV  <- w - e(p1_0, u1)                  # valued at old prices
CV  <- e(p1_1, u0) - w                  # valued at new prices
cat("EV =", round(EV, 4), " dCS =", round(dCS, 4), " CV =", round(CV, 4), "\n")
stopifnot(EV < dCS, dCS < CV)           # ordering EV < dCS < CV

## 6.5 Producer duality example --------------------------------------------
# y = z1^.5 z2^.5, r=(4,1), y=10 -> c=2y sqrt(r1 r2); Shephard z_i
r1 <- 4; r2 <- 1; y <- 10
cost <- 2 * y * sqrt(r1 * r2)
z1 <- y * sqrt(r2 / r1); z2 <- y * sqrt(r1 / r2)
cat("cost =", cost, " z1 =", z1, " z2 =", z2, " r1*z1 + r2*z2 =", r1 * z1 + r2 * z2, "\n")

## 6.5 Uncertainty example -------------------------------------------------
# u(w)=sqrt(w), 50-50 lottery {64,144}
uw <- sqrt
Ew  <- 0.5 * 64 + 0.5 * 144
Eu  <- 0.5 * uw(64) + 0.5 * uw(144)
CE  <- Eu^2
rp  <- Ew - CE
A_Ew <- 1 / (2 * Ew)
varw <- 0.5 * (64 - Ew)^2 + 0.5 * (144 - Ew)^2
rp_approx <- 0.5 * A_Ew * varw
cat("E[w] =", Ew, " CE =", CE, " risk premium =", rp,
    " approx =", round(rp_approx, 4), "\n")

## 6.6 Reverse self-checks: integrability, Gorman aggregation, SOSD ---------
# Slutsky matrix at the initial point (p0, w) = ((2, 1), 100)
q1  <- p1_0
x1s <- a * w / q1; x2s <- (1 - a) * w / p2
S <- matrix(c(-a * w / q1^2 + x1s * (a / q1),  0 + x2s * (a / q1),
              0 + x1s * ((1 - a) / p2),        -(1 - a) * w / p2^2 + x2s * ((1 - a) / p2)),
            nrow = 2, byrow = TRUE)
print(S)
cat("symmetry S12 - S21 =", S[1, 2] - S[2, 1], " (S p)_1 =", (S %*% c(q1, p2))[1], "\n")
cat("eigenvalues of S:", zapsmall(eigen(S)$values), "\n")
stopifnot(isTRUE(all.equal(S[1, 2], S[2, 1])), all(zapsmall(eigen(S)$values) <= 0))

# Gorman aggregation: two LES consumers with common marginal budget share b1 = 2/5, W = 300
b1 <- 2/5
x1_les <- function(g1, g2, wi, b = b1) g1 + b * (wi - q1 * g1 - p2 * g2) / q1
X1       <- x1_les(5, 10, 100) + x1_les(15, 0, 200)
X1_redis <- x1_les(5, 10, 150) + x1_les(15, 0, 150)
cat("individual x1 at (100, 200):", x1_les(5, 10, 100), x1_les(15, 0, 200),
    " at (150, 150):", x1_les(5, 10, 150), x1_les(15, 0, 150), "\n")
cat("X1 =", X1, " X1_redis =", X1_redis, "\n")
X1_b07       <- x1_les(5, 10, 100, 0.4) + x1_les(15, 0, 200, 0.7)
X1_b07_redis <- x1_les(5, 10, 150, 0.4) + x1_les(15, 0, 150, 0.7)
cat("B's slope 0.7: orig =", X1_b07, " redis =", X1_b07_redis, "\n")
stopifnot(isTRUE(all.equal(X1, X1_redis)))

# SOSD: tight lottery F = {100, 108} vs mean-preserving spread G = {64, 144}, u = sqrt
mean_F <- 0.5 * 100 + 0.5 * 108; mean_G <- 0.5 * 64 + 0.5 * 144
EU_F <- 0.5 * sqrt(100) + 0.5 * sqrt(108); EU_G <- 0.5 * sqrt(64) + 0.5 * sqrt(144)
cat("mean_F =", mean_F, " mean_G =", mean_G, " EU_F =", round(EU_F, 4), " EU_G =", EU_G, "\n")
stopifnot(mean_F == mean_G, EU_F > EU_G)

## 4.4 Comparative statics checks ------------------------------------------
# CRRA u(w)=w^(1-rho)/(1-rho): local risk premium rises with rho and variance.
rho_grid <- c(0.5, 1, 2)
risk_var <- 100
rp_local <- 0.5 * (rho_grid / Ew) * risk_var
cat("local risk premia by CRRA rho:", paste(round(rp_local, 3), collapse=", "), "\n")
stopifnot(all(diff(rp_local) > 0))
