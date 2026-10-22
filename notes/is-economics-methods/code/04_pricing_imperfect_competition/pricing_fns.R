# ------------------------------------------------------------------
# pricing_fns.R -- functions for chapter 4 (pricing and imperfect competition).
#
# Chapter card interface:
#   cournot_n(a, b, c, n)              -> list(q, p, profit): symmetric n-firm Cournot with p = a - bQ and
#                                         constant marginal cost c (q and profit per firm)
#   passthrough_cournot(n)             -> n / (n + 1): dp/dc of a common cost shock, symmetric Cournot, linear demand
#   logit_share(p, delta, alpha)       -> shares s_j = exp(delta_j - alpha p_j) / (1 + sum_k exp(delta_k - alpha p_k))
#   bertrand_logit(delta, alpha, c)    -> Bertrand-Nash prices, one single-product firm per product (nleqslv)
#   joint_max_logit(delta, alpha, c)   -> prices that maximise the joint profit of all products (nleqslv)
#   price_leader(delta, alpha, c)      -> c(p1, p2): sequential pricing in a logit duopoly, product 1 leads
#   bundle_opt(vals)                   -> list(price, rev, dwl): optimal pure-bundle price at zero marginal cost
#   bbm_segment(values, probs)         -> list(segments, profit, uniform_profit, total, efficient_total): the
#                                         consumer-optimal segmentation of Bergemann, Brooks and Morris (2015)
# Helpers:
#   passthrough_monopoly(Q, dp, d2p)   1 / (2 + Q p''(Q) / p'(Q)), the monopoly pass-through
#   logit_foc, solve_logit             logit first-order conditions under conduct kappa, and their solution
#   passthrough_logit                  the equilibrium pass-through matrix dp_k / dc_j (implicit function theorem)
#   logit_profit                       (p - c) s per unit of market size
#   opt_price                          the profit-maximising single price on an empirical distribution of values
#   linear_pd                          uniform price against third-degree discrimination, linear segment demands
#   gen_C1, gen_C2, gen_C3             the data generators of blocks C1 to C3 (01_dgp.R)
#   c1_equilibria                      Bertrand and joint prices and shares in every C1 market
#   pd_fit, pd_segment_prices, pd_outcomes, pd_run   block C2: OLS prediction, segment prices, welfare
#
# Models, in the notation the text uses (notation_v3.md 5.5). Code argument names are not notation:
#   delta in the code is the quality xi_j of the text, c is marginal cost c_j (mc), kappa the conduct parameter.
#   Linear demand: inverse demand p = a - bQ, constant marginal cost c. Monopoly pass-through
#     PT = 1 / (2 + Q p''(Q) / p'(Q)): 1/2 for linear demand, e / (e - 1) for constant elasticity e,
#     1 for exponential demand p = lambda log(A / Q). Symmetric Cournot with linear demand: n / (n + 1).
#   Logit: u_ij = xi_j - alpha p_j + eps_ij, u_i0 = eps_i0. Firm j maximises pi_j + kappa sum_{k != j} pi_k
#     (kappa = 0: Bertrand-Nash; kappa = 1: joint profit maximisation). With m_j = p_j - c_j the first-order
#     condition divided by s_j is
#       1 - alpha (1 - s_j) m_j + kappa alpha sum_{k != j} s_k m_k = 0,
#     so Bertrand gives m_j = 1 / (alpha (1 - s_j)) and joint maximisation m_j = 1 / (alpha s_0) for every j.
#     Single-product monopoly pass-through is 1 - s. Pass-through matrix: PT[k, j] = dp_k / dc_j (rows are
#     prices, columns are costs); the text's PT_{jk} = dp_k / dc_j is PT[k, j] here.
#   Price leadership (duopoly): the follower (product 2) plays its Bertrand best response b(p1), whose slope
#     on the best response is b'(p1) = s1 s2 / (1 - s2); the leader (product 1) maximises
#     (p1 - c1) s1(p1, b(p1)), with first-order condition alpha m1 [(1 - s1) - s1 s2^2 / (1 - s2)] = 1.
#   Bundling and BBM: unit demand, zero marginal cost (information goods); a consumer buys when her value is at
#     least the price. bundle_opt's rev and dwl are per consumer.
#   BBM greedy construction: on the support S of the remaining mass, the extremal market x^S makes the seller
#     indifferent among all prices in S: Pr(v >= v_k) = min(S) / v_k. Put the largest weight on it that leaves no
#     negative mass, drop the values whose mass is used up, repeat. Each segment is priced at its lowest value
#     (the seller is indifferent, every buyer buys), so total surplus is efficient; the uniform monopoly price
#     stays in every segment's support, so the producer still earns exactly the uniform profit.
# Packages: stats, nleqslv (solve_logit, price_leader).
# ------------------------------------------------------------------

## =================================================================
## Linear demand: monopoly, Cournot, pass-through
## =================================================================

#' Symmetric n-firm Cournot with inverse demand p = a - bQ and constant marginal cost c.
#' Each firm's first-order condition p - c - b q = 0 gives q = (a - c) / (b (n + 1)), p = (a + n c) / (n + 1).
cournot_n <- function(a, b, c, n) {
  q <- (a - c) / (b * (n + 1))
  p <- (a + n * c) / (n + 1)
  list(q = q, p = p, profit = (p - c) * q)
}

#' dp/dc of a common cost shock in symmetric n-firm Cournot with linear demand.
passthrough_cournot <- function(n) n / (n + 1)

#' Monopoly pass-through dp/dc = 1 / (2 + Q p''(Q) / p'(Q)), evaluated at output Q with p'(Q) = dp, p''(Q) = d2p.
passthrough_monopoly <- function(Q, dp, d2p) 1 / (2 + Q * d2p / dp)

## =================================================================
## Logit demand: shares, first-order conditions, equilibria, pass-through
## =================================================================

#' Logit shares with an outside good of utility 0 (computed stably).
logit_share <- function(p, delta, alpha) {
  u <- delta - alpha * p
  m <- max(0, u)
  e <- exp(u - m)
  e / (exp(-m) + sum(e))
}

#' Profit per unit of market size, (p - c) s.
logit_profit <- function(p, delta, alpha, c) (p - c) * logit_share(p, delta, alpha)

#' First-order conditions (divided by s_j) when firm j maximises pi_j + kappa sum_{k != j} pi_k.
logit_foc <- function(p, delta, alpha, c, kappa = 0) {
  s <- logit_share(p, delta, alpha)
  m <- p - c
  sm <- s * m
  1 - alpha * (1 - s) * m + kappa * alpha * (sum(sm) - sm)
}

# Joint maximisation in closed form up to one scalar: every markup is m with alpha m = 1 + sum_j exp(delta_j - alpha c_j - alpha m).
.joint_markup <- function(delta, alpha, c) {
  f <- function(m) alpha * m - 1 - sum(exp(delta - alpha * c - alpha * m))
  hi <- 1 / alpha
  while (f(hi) < 0) hi <- 2 * hi
  stats::uniroot(f, c(0, hi), tol = 1e-14)$root
}

#' Equilibrium prices under conduct kappa (0 = Bertrand-Nash, 1 = joint maximisation), solved with nleqslv.
#' Stops unless the first-order conditions hold to 1e-11.
solve_logit <- function(delta, alpha, c, kappa = 0, start = NULL) {
  J <- length(delta)
  if (length(c) != J) stop("solve_logit: delta and c must have the same length", call. = FALSE)
  if (is.null(start)) {
    start <- c + 1 / alpha
    if (kappa > 0) start <- start + kappa * (c + .joint_markup(delta, alpha, c) - start)
  }
  sol <- nleqslv::nleqslv(start, function(p) logit_foc(p, delta, alpha, c, kappa), method = "Newton",
                          control = list(xtol = 1e-15, ftol = 1e-14, maxit = 500))
  if (!(sol$termcd %in% c(1, 2)) || max(abs(sol$fvec)) > 1e-11)
    stop("solve_logit: no convergence (termcd ", sol$termcd, ", max |foc| ", format(max(abs(sol$fvec))), ")", call. = FALSE)
  sol$x
}

#' Bertrand-Nash prices with one single-product firm per product: p_j - c_j = 1 / (alpha (1 - s_j)).
bertrand_logit <- function(delta, alpha, c) solve_logit(delta, alpha, c, kappa = 0)

#' Prices that maximise the joint profit of all products: p_j - c_j = 1 / (alpha s_0) for every j.
joint_max_logit <- function(delta, alpha, c) solve_logit(delta, alpha, c, kappa = 1)

#' Equilibrium pass-through matrix under conduct kappa: PT[k, j] = dp_k / dc_j, from the implicit function
#' theorem on logit_foc (analytic Jacobians), at the equilibrium prices solve_logit returns.
passthrough_logit <- function(delta, alpha, c, kappa = 0, p = NULL) {
  if (is.null(p)) p <- solve_logit(delta, alpha, c, kappa)
  J <- length(p)
  s <- logit_share(p, delta, alpha)
  m <- p - c
  Ds <- alpha * (outer(s, s) - diag(s, J))                     # Ds[k, l] = ds_k / dp_l
  one_s <- matrix(s, J, J, byrow = TRUE) - diag(s, J)          # [j, l] = s_l for l != j, 0 on the diagonal
  Jp <- alpha * diag(m, J) %*% Ds - alpha * diag(1 - s, J) +
    kappa * alpha * (matrix(drop(crossprod(m, Ds)), J, J, byrow = TRUE) - diag(m, J) %*% Ds) +
    kappa * alpha * one_s
  Jc <- alpha * diag(1 - s, J) - kappa * alpha * one_s
  -solve(Jp, Jc)
}

#' Sequential pricing in a logit duopoly: product 1 commits to its price first, product 2 best-responds.
#' Returns c(p1, p2). Solves the leader's and the follower's first-order conditions jointly (nleqslv),
#' starting from the Bertrand-Nash prices.
price_leader <- function(delta, alpha, c) {
  if (length(delta) != 2L || length(c) != 2L) stop("price_leader: a duopoly (two products) is required", call. = FALSE)
  f <- function(p) {
    s <- logit_share(p, delta, alpha)
    m <- p - c
    c(alpha * m[1] * ((1 - s[1]) - s[1] * s[2]^2 / (1 - s[2])) - 1,
      alpha * m[2] * (1 - s[2]) - 1)
  }
  sol <- nleqslv::nleqslv(bertrand_logit(delta, alpha, c), f, method = "Newton",
                          control = list(xtol = 1e-15, ftol = 1e-14, maxit = 500))
  if (!(sol$termcd %in% c(1, 2)) || max(abs(sol$fvec)) > 1e-11)
    stop("price_leader: no convergence (termcd ", sol$termcd, ")", call. = FALSE)
  sol$x
}

## =================================================================
## Unit demand: optimal single price, bundling, the limits of price discrimination
## =================================================================

#' The profit-maximising single price on an empirical distribution of values v (one unit each, marginal
#' cost `cost`): a consumer buys when v >= price, and the optimum is one of the values above cost.
#' Ties in profit go to the higher price.
opt_price <- function(v, cost = 0) {
  v <- v[v > cost]
  if (!length(v)) stop("opt_price: no value above cost", call. = FALSE)
  u <- sort(unique(v), decreasing = TRUE)
  n_ge <- cumsum(tabulate(match(v, u), length(u)))
  u[which.max((u - cost) * n_ge)]
}

#' Optimal pure bundle at zero marginal cost. vals: an N x J matrix of per-good values (the bundle value is
#' the row sum) or a vector of bundle values. rev = expected revenue per consumer, dwl = the value of the
#' bundle to the consumers it excludes, per consumer.
bundle_opt <- function(vals) {
  v <- if (is.matrix(vals)) rowSums(vals) else as.numeric(vals)
  price <- opt_price(v, 0)
  buy <- v >= price
  list(price = price, rev = price * mean(buy), dwl = sum(v[!buy]) / length(v))
}

#' Consumer-optimal segmentation of a discrete value distribution (Bergemann, Brooks and Morris 2015), built
#' by the greedy construction (see the header). values: strictly increasing, positive; probs: masses (>= 0).
#' segments = list(mass = segment x value matrix of masses (column sums = probs), weight = segment masses,
#' price = each segment's price, its lowest value). profit and total are the producer and total surplus with
#' every segment priced at its lowest value; uniform_profit is the best single price; efficient_total = sum(values probs).
bbm_segment <- function(values, probs) {
  K <- length(values)
  if (length(probs) != K || K == 0L) stop("bbm_segment: values and probs must have the same positive length", call. = FALSE)
  if (any(diff(values) <= 0) || values[1] <= 0) stop("bbm_segment: values must be positive and strictly increasing", call. = FALSE)
  if (any(probs < 0)) stop("bbm_segment: probs must be non-negative", call. = FALSE)
  tol <- 1e-14 * sum(probs)
  r <- probs
  mass <- matrix(0, 0, K)
  price <- numeric(0)
  while (any(r > tol)) {
    S <- which(r > tol)
    vS <- values[S]
    tail_S <- vS[1] / vS                                  # Pr(v >= v_k) in the extremal market on S
    x <- tail_S - c(tail_S[-1], 0)
    ratio <- r[S] / x
    w <- min(ratio)
    seg <- numeric(K)
    seg[S] <- w * x
    r[S] <- r[S] - seg[S]
    r[S[ratio <= w * (1 + 1e-12)]] <- 0                   # the values used up by this segment
    mass <- rbind(mass, seg, deparse.level = 0)
    price <- c(price, vS[1])
  }
  weight <- rowSums(mass)
  unif <- vapply(values, function(p) p * sum(probs[values >= p]), numeric(1))
  list(segments = list(mass = mass, weight = weight, price = price),
       profit = sum(price * weight),
       uniform_profit = max(unif),
       total = sum(mass %*% values),
       efficient_total = sum(values * probs))
}

#' Third-degree price discrimination with linear segment demands p = a_k - b_k q_k and marginal cost `cost`,
#' against the best single price. The uniform optimum is interior to one regime "the m segments with the
#' highest choke prices are served" (the profit function has convex kinks at the choke prices).
linear_pd <- function(a, b, cost) {
  K <- length(a)
  b <- rep_len(b, K)
  # discrimination: a monopoly price in every segment that is worth serving
  served_pd <- a > cost
  p_pd <- ifelse(served_pd, (a + cost) / 2, NA_real_)
  q_pd <- ifelse(served_pd, (a - p_pd) / b, 0)
  # one price: try every regime, keep the valid interior optimum with the highest profit
  ord <- order(a, decreasing = TRUE)
  best <- NULL
  for (m in seq_len(K)) {
    S <- ord[seq_len(m)]
    if (any(a[S] <= cost)) break
    p <- (sum(a[S] / b[S]) + cost * sum(1 / b[S])) / (2 * sum(1 / b[S]))
    lo <- if (m < K) max(a[ord[(m + 1):K]]) else -Inf
    if (p < min(a[S]) && p >= lo) {
      prof <- (p - cost) * sum((a[S] - p) / b[S])
      if (is.null(best) || prof > best$profit) best <- list(p = p, profit = prof, m = m)
    }
  }
  if (is.null(best)) stop("linear_pd: no segment is worth serving", call. = FALSE)
  p_u <- best$p
  q_u <- pmax(0, (a - p_u) / b)
  cs_u <- sum(b * q_u^2 / 2)
  cs_pd <- sum(b * q_pd^2 / 2)
  prof_pd <- sum(ifelse(served_pd, (p_pd - cost) * q_pd, 0))
  list(p_uniform = p_u, q_uniform = q_u, Q_uniform = sum(q_u), profit_uniform = best$profit, CS_uniform = cs_u,
       W_uniform = best$profit + cs_u, n_served_uniform = sum(q_u > 0),
       p_pd = p_pd, q_pd = q_pd, Q_pd = sum(q_pd), profit_pd = prof_pd, CS_pd = cs_pd, W_pd = prof_pd + cs_pd,
       n_served_pd = sum(q_pd > 0))
}

## =================================================================
## Data generators (01_dgp.R) and the block C1 equilibria
## =================================================================

#' Block C1: M local markets of the Cirrus-Garnet logit duopoly. Draws, in this order: cost shifters of
#' Cirrus and Garnet (uniform, +- cost_halfwidth around c_base), demand shocks xi1, xi2 (normal, sd xi_sd) added
#' to the qualities delta, and the intercept shock of the homogeneous-product variant (normal, sd hom_eta_sd).
#' Then solves every market's equilibria (c1_equilibria). Columns ending in _true are simulation truth.
gen_C1 <- function(P) {
  M <- P$M
  w1 <- stats::runif(M, -P$cost_halfwidth, P$cost_halfwidth)
  w2 <- stats::runif(M, -P$cost_halfwidth, P$cost_halfwidth)
  xi1 <- stats::rnorm(M, 0, P$xi_sd)
  xi2 <- stats::rnorm(M, 0, P$xi_sd)
  eta <- stats::rnorm(M, 0, P$hom_eta_sd)
  d <- data.frame(market = seq_len(M), c1 = P$c_base[1] + w1, c2 = P$c_base[2] + w2,
                  xi1_true = xi1, xi2_true = xi2, a_hom_true = P$hom_a + eta)
  cbind(d, c1_equilibria(d, P))
}

#' Bertrand-Nash (kappa = 0) and joint-maximising (kappa = 1) prices and shares in every C1 market, and the
#' homogeneous-product variant with inverse demand p = a_m - Q (slope hom_b = 1): Cournot
#' p = (a + c1 + c2) / 3, q_i = (a - 2 c_i + c_j) / 3; Bertrand with unequal costs: the low-cost firm sells
#' at the rival's cost, p = max(c1, c2).
c1_equilibria <- function(d, P) {
  if (P$hom_b != 1) stop("c1_equilibria: the homogeneous variant is written for slope 1", call. = FALSE)
  M <- nrow(d)
  out <- matrix(NA_real_, M, 8, dimnames = list(NULL, c("p1_bertrand", "p2_bertrand", "s1_bertrand", "s2_bertrand",
                                                        "p1_joint", "p2_joint", "s1_joint", "s2_joint")))
  for (m in seq_len(M)) {
    dm <- P$delta + c(d$xi1_true[m], d$xi2_true[m])
    cm <- c(d$c1[m], d$c2[m])
    pb <- bertrand_logit(dm, P$alpha, cm)
    pj <- joint_max_logit(dm, P$alpha, cm)
    out[m, ] <- c(pb, logit_share(pb, dm, P$alpha), pj, logit_share(pj, dm, P$alpha))
  }
  a <- d$a_hom_true
  data.frame(out,
             p_cournot_hom = (a + d$c1 + d$c2) / 3,
             q1_cournot_hom = (a - 2 * d$c1 + d$c2) / 3,
             q2_cournot_hom = (a - 2 * d$c2 + d$c1) / 3,
             p_bertrand_hom = pmax(d$c1, d$c2))
}

#' Block C2: N merchants deciding whether to subscribe to Cirrus Pro (unit demand, one month). Draws, in this
#' order: standardised log monthly sales, standardised tenure, category, the log-WTP noise, and the random
#' half split (exactly N / 2 in the training half). WTP in whole USD per month, at least 1:
#'   wtp = max(1, round(exp(index + eps))),  index = const + sales sales_z + tenure tenure_z + electronics 1{electronics}
#'   + home 1{home},  eps ~ N(0, sigma^2).
gen_C2 <- function(P) {
  N <- P$N
  sales_z <- stats::rnorm(N)
  tenure_z <- stats::rnorm(N)
  category <- factor(sample(P$cat_levels, N, replace = TRUE, prob = P$cat_probs), levels = P$cat_levels)
  b <- P$beta
  index <- b[["const"]] + b[["sales"]] * sales_z + b[["tenure"]] * tenure_z +
    b[["electronics"]] * (category == "electronics") + b[["home"]] * (category == "home")
  eps <- stats::rnorm(N, 0, P$sigma)
  wtp <- pmax(1, round(exp(index + eps)))
  train <- sample(rep(c(TRUE, FALSE), N / 2))
  data.frame(customer = seq_len(N), sales_z = sales_z, tenure_z = tenure_z, category = category, wtp = wtp,
             train = train, index_true = index, eps_true = eps)
}

#' Block C3: N consumers' values for J_max = max(J) goods, independent U[0, 1]; returns the N x length(J)
#' matrix of bundle values (the sum of the first J goods), columns J1, J2, ...
gen_C3 <- function(P) {
  V <- matrix(stats::runif(P$N * max(P$J)), P$N, max(P$J))
  B <- vapply(P$J, function(J) rowSums(V[, seq_len(J), drop = FALSE]), numeric(P$N))
  colnames(B) <- paste0("J", P$J)
  B
}

## =================================================================
## Block C2: prediction by OLS, segment prices, welfare on the evaluation half
## =================================================================

#' OLS of log WTP on the features in the training half; returns the fit and the predicted index for everyone.
pd_fit <- function(d) {
  fit <- stats::lm(log(wtp) ~ sales_z + tenure_z + category, data = d[d$train, ])
  list(fit = fit, index = unname(stats::predict(fit, newdata = d)))
}

#' K segments by quantiles of the training half's predicted index (type 7 quantiles; K = 1 is one price for
#' everyone); each segment's price is the profit-maximising price on its training customers' WTP.
pd_segment_prices <- function(index_train, wtp_train, K, cost = 0) {
  br <- if (K == 1) c(-Inf, Inf) else c(-Inf, stats::quantile(index_train, seq_len(K - 1) / K, names = FALSE), Inf)
  g <- findInterval(index_train, br, left.open = TRUE)
  list(breaks = br, price = vapply(seq_len(K), function(k) opt_price(wtp_train[g == k], cost), numeric(1)),
       n_train = tabulate(g, K))
}

#' Per-customer outcomes when each customer faces `price` (buy when wtp >= price).
pd_outcomes <- function(wtp, price, cost = 0) {
  buy <- wtp >= price
  profit <- mean((price - cost) * buy)
  cs <- mean((wtp - price) * buy)
  c(profit = profit, cs = cs, ts = profit + cs, served = mean(buy), mean_price_paid = sum(price * buy) / sum(buy))
}

#' The whole C2 pipeline on one data set: OLS on the training half, segment prices for every K in Ks set on
#' the training half, outcomes on the evaluation half. Returns the fit, the ladder (one row per K) and the
#' evaluation customers' prices under each K.
pd_run <- function(d, Ks, cost = 0) {
  f <- pd_fit(d)
  itr <- f$index[d$train]; iev <- f$index[!d$train]
  wtr <- d$wtp[d$train]; wev <- d$wtp[!d$train]
  prices <- list(); seg <- list()
  ladder <- t(vapply(Ks, function(K) {
    sp <- pd_segment_prices(itr, wtr, K, cost)
    p <- sp$price[findInterval(iev, sp$breaks, left.open = TRUE)]
    prices[[as.character(K)]] <<- p
    seg[[as.character(K)]] <<- sp
    pd_outcomes(wev, p, cost)
  }, numeric(5)))
  list(fit = f$fit, index = f$index, ladder = data.frame(K = Ks, ladder), prices = prices, segments = seg)
}
