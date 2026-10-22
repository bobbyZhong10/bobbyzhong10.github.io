# ------------------------------------------------------------------
# 04_twosided.R -- Armstrong (2006) monopoly two-sided pricing,
# solved numerically for the "Cadence" platform to give the chapter
# a real worked number for the price-below-cost / subsidy-side point.
#
# Two sides: merchants (side 1) and cardholders (side 2). Linear
# membership demands with a cross-group externality:
#   n1 = a1 + b1 * (alpha1 * n2 - p1)
#   n2 = a2 + b2 * (alpha2 * n1 - p2)
# alpha1 = benefit a merchant gets per cardholder (large: merchants
# love reach); alpha2 = benefit a cardholder gets per merchant
# (small). Platform per-member costs f1, f2. Monopoly maximizes
#   Pi = n1*(p1 - f1) + n2*(p2 - f2).
# The Armstrong formula predicts p_i = f_i - alpha_j*n_j + markup, so
# the side that confers the larger externality on the other side gets
# the larger discount -- here cardholders end up subsidized. The markup
# n_i/b_i also changes with the platform's scale, so the gap between
# the platform price and the isolated price has two parts.
# ------------------------------------------------------------------
code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))

## ---- parameters --------------------------------------------------
a1 <- 0.6; b1 <- 0.5; alpha1 <- 3.0; f1 <- 1.0   # merchants
a2 <- 0.5; b2 <- 0.5; alpha2 <- 0.6; f2 <- 1.0   # cardholders

## Given prices, solve the membership fixed point (linear system):
##   n1 - b1*alpha1*n2 = a1 - b1*p1
##   -b2*alpha2*n1 + n2 = a2 - b2*p2
memberships <- function(p1, p2) {
  A <- matrix(c(1, -b2 * alpha2, -b1 * alpha1, 1), 2, 2)
  rhs <- c(a1 - b1 * p1, a2 - b2 * p2)
  solve(A, rhs)          # c(n1, n2)
}

profit <- function(p) {
  n <- memberships(p[1], p[2])
  n[1] * (p[1] - f1) + n[2] * (p[2] - f2)
}

opt <- optim(c(1, 1), function(p) -profit(p), method = "BFGS")
p_star <- opt$par
n_star <- memberships(p_star[1], p_star[2])

## Benchmark: naive one-sided monopoly that ignores the externality
## (sets alpha1 = alpha2 = 0 when pricing, i.e. treats each side in
## isolation), for contrast.
profit_naive_side <- function(p, side) {
  # each side priced as an isolated monopolist against its own demand
  if (side == 1) { n <- a1 + b1 * (-p); return(n * (p - f1)) }
  else           { n <- a2 + b2 * (-p); return(n * (p - f2)) }
}
p1_naive <- optimize(function(p) -profit_naive_side(p, 1), c(-10, 10))$minimum
p2_naive <- optimize(function(p) -profit_naive_side(p, 2), c(-10, 10))$minimum

## Decompose each price into Armstrong's three terms. With linear
## demand dn_i/du_i = b_i, so p_i = f_i - alpha_j * n_j + n_i / b_i.
## The isolated benchmark has no externality term, only the markup.
ext_star     <- c(alpha2 * n_star[2], alpha1 * n_star[1])  # alpha_j * n_j, side i = 1, 2
markup_star  <- n_star / c(b1, b2)                          # n_i / b_i
n_naive      <- c(a1 - b1 * p1_naive, a2 - b2 * p2_naive)  # isolated memberships
markup_naive <- n_naive / c(b1, b2)
z3 <- function(x) round(x, 3) + 0                           # avoid printing -0.000

cat("== Armstrong monopoly two-sided pricing (Cadence) ==\n")
cat(sprintf("Params: alpha1(merchant per cardholder)=%.1f, alpha2(cardholder per merchant)=%.1f\n",
            alpha1, alpha2))
cat(sprintf("Optimal prices:   merchants p1* = %.3f,  cardholders p2* = %.3f\n",
            p_star[1], p_star[2]))
cat(sprintf("Memberships:      n1* = %.3f,  n2* = %.3f\n", n_star[1], n_star[2]))
cat(sprintf("Externality discount on cardholder price alpha1*n1* = %.3f\n",
            alpha1 * n_star[1]))
cat(sprintf("One-side-blind benchmark prices: p1 = %.3f, p2 = %.3f (isolated monopoly prices, no subsidy)\n",
            p1_naive, p2_naive))
cat(sprintf("Isolated benchmark memberships: n1 = %.3f, n2 = %.3f (cardholder price = cost, zero demand)\n",
            z3(n_naive[1]), z3(n_naive[2])))
cat(sprintf("Externality term alpha_j*n_j*: merchants (alpha2*n2*) = %.3f, cardholders (alpha1*n1*) = %.3f\n",
            ext_star[1], ext_star[2]))
cat(sprintf("Markup n_i/b_i: platform merchants = %.3f, cardholders = %.3f; isolated merchants = %.3f, cardholders = %.3f\n",
            markup_star[1], markup_star[2], z3(markup_naive[1]), z3(markup_naive[2])))
cat(sprintf("Check f_i - alpha_j*n_j + n_i/b_i: merchants = %.3f, cardholders = %.3f\n",
            f1 - ext_star[1] + markup_star[1], f2 - ext_star[2] + markup_star[2]))
cat(sprintf("Total platform profit: %.3f\n", profit(p_star)))

saveRDS(list(p_star = p_star, n_star = n_star,
             p_naive = c(p1_naive, p2_naive),
             n_naive = n_naive, ext_star = ext_star,
             markup_star = markup_star, markup_naive = markup_naive,
             params = list(a1=a1,b1=b1,alpha1=alpha1,f1=f1,
                           a2=a2,b2=b2,alpha2=alpha2,f2=f2)),
        file.path(code_dir, "res24_twosided.rds"))
cat("Saved res24_twosided.rds\n")
