# ------------------------------------------------------------------
# 02_estimate.R -- recovering price sensitivity and the search cost
# in the Lumen log, three ways.
#
#   (1) naive full-information conditional logit: choice set = all J
#       shown products. Attenuated, because consumers never inspect
#       most of the cheap options further down the list.
#   (2) inspected-set conditional logit: choice set = the products the
#       consumer actually inspected. OVER-corrects: the inspected set
#       is itself selected by the optimal stopping rule (a consumer who
#       meets a cheap high-match option early stops early), exaggerating
#       the price-purchase link.
#   (3) structural sequential-search estimator by simulated method of
#       moments (SMM): simulate the full search-and-purchase model at
#       candidate (gamma, alpha, c) and match data moments (mean search
#       depth, buy rate, mean purchased price, share stopping after one
#       inspection). Recovers the true price coefficient AND the search
#       cost. The outside-option location o0 is fixed at its normalized
#       value; common random numbers keep the objective smooth.
# ------------------------------------------------------------------
library(data.table)
library(survival)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
d      <- readRDS(file.path(code_dir, "lumen_log.rds"))
truth  <- readRDS(file.path(code_dir, "truth25.rds"))
log_dt <- d$log; setDT(log_dt)
prices <- d$prices; J <- length(prices); o0 <- truth$o0

buyers <- log_dt[, .(any = max(bought)), by = cons][any == 1, cons]

## ---- (1) naive full-information conditional logit ---------------
naive <- copy(log_dt[cons %in% buyers]); naive[, cs := cons]
m_naive <- clogit(bought ~ price + strata(cs), data = naive)
a_naive <- -coef(m_naive)["price"]; se_naive <- sqrt(vcov(m_naive)["price", "price"])

## ---- (2) inspected-set conditional logit ------------------------
insp <- copy(log_dt[cons %in% buyers & inspected == 1]); insp[, cs := cons]
m_insp <- clogit(bought ~ price + strata(cs), data = insp)
a_insp <- -coef(m_insp)["price"]; se_insp <- sqrt(vcov(m_insp)["price", "price"])

## ---- data moments -----------------------------------------------
n_insp   <- d$n_insp; buy <- as.integer(d$buy_prod > 0)
purch_pr <- prices[d$buy_prod[d$buy_prod > 0]]
m_data <- c(depth = mean(n_insp), buy = mean(buy),
            price = mean(purch_pr), p1 = mean(n_insp == 1))

## ---- structural SMM (common random numbers) ---------------------
set.seed(2100)
N <- 6000
ordL <- lapply(1:N, function(i) sample(J))
G    <- matrix(-log(-log(matrix(runif(N * J), N, J))), N, J)
U0   <- o0 - log(-log(runif(N)))
uzP  <- sample(prices, 150000, TRUE); uzG <- -log(-log(runif(150000)))

sim_moments <- function(g, a, cc) {
  mu <- g - a * prices
  z <- tryCatch(uniroot(function(z) mean(pmax(g - a * uzP + uzG - z, 0)) - cc,
                        c(-8, 25))$root, error = function(e) NA_real_)
  if (is.na(z)) return(rep(NA_real_, 4))
  dep <- integer(N); bj <- integer(N)
  for (i in 1:N) {
    ord <- ordL[[i]]; bp <- -Inf; jj <- 0L; k <- 0L
    for (r in 1:J) {
      j <- ord[r]; u <- mu[j] + G[i, r]; k <- k + 1L
      if (u > bp) { bp <- u; jj <- j }
      if (max(U0[i], bp) >= z) break
    }
    dep[i] <- k; bj[i] <- if (bp > U0[i]) jj else 0L
  }
  c(mean(dep), mean(bj > 0), mean(prices[bj[bj > 0]]), mean(dep == 1))
}

W <- diag(1 / m_data^2)
obj <- function(p) {
  ms <- sim_moments(p[1], p[2], exp(p[3]))
  if (any(is.na(ms))) return(1e6)
  dv <- ms - m_data; as.numeric(dv %*% W %*% dv)
}

starts <- list(c(2.5, 0.02, log(0.30)), c(3.2, 0.04, log(0.50)),
               c(2.9, 0.03, log(0.40)))
best <- NULL; bestv <- Inf
for (st in starts) {
  f <- optim(st, obj, method = "Nelder-Mead",
             control = list(maxit = 500, reltol = 1e-9))
  if (f$value < bestv) { bestv <- f$value; best <- f }
}
g_hat <- best$par[1]; a_hat <- best$par[2]; c_hat <- exp(best$par[3])
sc_dollars <- c_hat / a_hat          # search cost per inspection in dollars

cat("== Recovering price sensitivity (true alpha = 0.030) ==\n")
cat(sprintf("(1) naive full-info clogit alpha = %.4f (SE %.4f)  [attenuated]\n", a_naive, se_naive))
cat(sprintf("(2) inspected-set clogit  alpha = %.4f (SE %.4f)  [over-corrects]\n", a_insp, se_insp))
cat(sprintf("(3) structural SMM        alpha = %.4f\n", a_hat))
cat(sprintf("    structural SMM gamma = %.3f (true 2.90), c = %.3f (true %.3f)\n",
            g_hat, c_hat, truth$c))
cat(sprintf("    search cost per inspection = c/alpha = $%.2f\n", sc_dollars))
cat(sprintf("Data moments : depth %.3f, buy %.3f, price %.2f, p1 %.3f\n",
            m_data[1], m_data[2], m_data[3], m_data[4]))
cat(sprintf("SMM  moments : %s\n", paste(round(sim_moments(g_hat, a_hat, c_hat), 3), collapse = ", ")))
cat(sprintf("Prominence (advertising) lift from randomized ranking = %.3f\n",
            truth$top_buy - truth$bot_buy))

saveRDS(list(a_naive = a_naive, se_naive = se_naive,
             a_insp = a_insp, se_insp = se_insp,
             g_hat = g_hat, a_hat = a_hat, c_hat = c_hat, sc_dollars = sc_dollars,
             m_data = m_data, prom_lift = truth$top_buy - truth$bot_buy),
        file.path(code_dir, "res25_estimate.rds"))
cat("Saved res25_estimate.rds\n")
