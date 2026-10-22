# ------------------------------------------------------------------
# 05_policy.R -- Policy learning, Module 13.
# The tool has a per-merchant cost c (in log-GMV-growth units). The
# platform should deploy it only where tau(X) > c. We learn a shallow,
# deployable targeting rule (depth-2 policy tree) from doubly robust
# scores and compare its value to treat-all / treat-none / oracle.
# Writes res_policy.rds
# ------------------------------------------------------------------
set.seed(1313)
suppressMessages(library(policytree))

sc    <- readRDS("dr_scores.rds")
truth <- readRDS("truth.rds")
X <- sc$X; gamma <- sc$gamma; tau_true <- sc$tau_true
n <- nrow(X)

cost <- 0.10                       # per-merchant deployment cost
# net doubly robust reward matrix: subtract cost from the "treat" column
Gamma <- gamma
Gamma[, 2] <- Gamma[, 2] - cost

## ---- learn a depth-2 policy tree ---------------------------------
# restrict to the two economically meaningful covariates for a readable
# rule (size X1, tenure X2); grf CATE already told us these drive tau.
Xp <- X[, c(1, 2)]
colnames(Xp) <- c("size_X1", "tenure_X2")
ptree <- policy_tree(Xp, Gamma, depth = 2)
pi_hat <- predict(ptree, Xp)       # actions in {1,2} (1=control,2=treat)

## ---- policy values via doubly robust scores ----------------------
val <- function(action) mean(Gamma[cbind(1:n, action)])
v_none  <- val(rep(1L, n))
v_all   <- val(rep(2L, n))
v_tree  <- val(pi_hat)
# oracle: treat exactly where true tau > cost
pi_oracle <- ifelse(tau_true > cost, 2L, 1L)
v_oracle <- val(pi_oracle)

share_treated <- function(a) mean(a == 2L)

## ---- clean evaluation against the KNOWN tau (noise-free) ----------
# value gain over treat-none = E[(tau - c) * 1{treated}]; oracle exact.
tgain <- function(action) mean((tau_true - cost) * (action == 2L))
tg_none   <- tgain(rep(1L, n))
tg_all    <- tgain(rep(2L, n))
tg_tree   <- tgain(pi_hat)
tg_oracle <- tgain(pi_oracle)

res_policy <- list(
  cost = cost,
  value = c(none = v_none, all = v_all, tree = v_tree, oracle = v_oracle),
  share = c(none = 0, all = 1,
            tree = share_treated(pi_hat),
            oracle = share_treated(pi_oracle)),
  gain_tree_over_all  = v_tree - v_all,
  gain_tree_over_none = v_tree - v_none,
  regret_vs_oracle    = v_oracle - v_tree,
  agree_oracle = mean(pi_hat == pi_oracle),
  true_gain = c(none = tg_none, all = tg_all, tree = tg_tree,
                oracle = tg_oracle),
  tree = ptree
)
saveRDS(res_policy, "res_policy.rds")

cat(sprintf("Per-merchant cost c = %.3f\n\n", cost))
cat("Policy values (doubly robust, net of cost):\n")
cat(sprintf("  treat none   : %+.4f  (share treated %.2f)\n", v_none, 0))
cat(sprintf("  treat all    : %+.4f  (share treated %.2f)\n", v_all, 1))
cat(sprintf("  policy tree  : %+.4f  (share treated %.2f)\n",
            v_tree, share_treated(pi_hat)))
cat(sprintf("  oracle rule  : %+.4f  (share treated %.2f)\n",
            v_oracle, share_treated(pi_oracle)))
cat(sprintf("\n  tree gain over treat-all  : %+.4f\n", v_tree - v_all))
cat(sprintf("  tree gain over treat-none : %+.4f\n", v_tree - v_none))
cat(sprintf("  regret vs oracle          : %.4f\n", v_oracle - v_tree))
cat(sprintf("  agreement with oracle rule: %.3f\n", mean(pi_hat == pi_oracle)))
cat("\nValue GAIN over treat-none, evaluated at the KNOWN tau:\n")
cat(sprintf("  treat all   : %+.4f\n", tg_all))
cat(sprintf("  policy tree : %+.4f  (%.1f%% of oracle gain)\n",
            tg_tree, 100 * tg_tree / tg_oracle))
cat(sprintf("  oracle rule : %+.4f\n", tg_oracle))
cat("\nLearned policy tree:\n"); print(ptree)
