# =====================================================================
# 22_synthetic_control_factor_models / 06_design.R
# Block S4: designing a pilot with synthetic controls (Abadie and Zhao,
# "Synthetic controls for experimental design", arXiv 2108.02196v5,
# Sections 2 and 3) versus random assignment.
#
# Setting: J = 40 eligible metros, 24 history months (the paper's T_0),
# a 6-month commission pilot in months 25-30 in at most m = 5 metros.
# aster_S4.rds holds the untreated outcomes Y0 only; the analyst sees
# Y0 + tau * D_j in the pilot months (tau = 0.04 for every pilot metro).
# The target is the average effect over all 40 metros (the paper's tau_t
# with f_j = 1/J), here 0.04.
#
# Synthetic design (what is implemented): the paper's Constrained design,
# problem (6), with f_j = 1/J, predictors X_j = the 24 history outcomes
# (E = {1, ..., 24}, no covariates, no rescaling: one outcome in one
# unit), and 1 <= ||w||_0 <= 5:
#     min_{w, v} ||Xbar - sum_j w_j X_j||^2 + ||Xbar - sum_j v_j X_j||^2
#     s.t. w, v on the simplex, w_j v_j = 0, 1 <= ||w||_0 <= 5,
# Xbar = the 40-metro average history. For a given treated set S the
# problem splits into two SC problems (w on S, v on the other metros;
# the paper's Online Appendix OA.5, "enumeration"), each solved with
# sc_weights() (V = I). The paper enumerates every S; with 40 metros and
# up to 5 treated (about 760,000 sets) this script instead runs a local
# search over S (moves: swap one treated metro for a control, add one,
# drop one; best improvement until no move helps) from 100 random starts
# (seed in est_seeds.R) plus a greedy forward start, and keeps the best set found.
# The paper does not require exact optimality, only a feasible (w, v)
# with a close fit (its Section 2). Units with w_j > 0 get the pilot;
# the estimator is (7): tau_hat_t = sum_j w_j Y_jt - sum_j v_j Y_jt,
# averaged over the 6 pilot months.
# Not implemented: the paper's other formulations (8)-(11), the
# bias-corrected estimator (OA.1) and its permutation inference (19).
#
# Random assignment (the paper's RND, Section 5.2.4): 5 of the 40 metros
# drawn at random, a new draw in every replicate; estimators: the
# pilot-window difference in means (the paper's RND), and, as a stronger
# baseline, the difference in differences that also subtracts the
# history-period difference.
#
# Evaluation: the 500 replicate pilot windows Y0_pilot_reps (the same
# history, new demand shocks and noise in months 25-30). For each
# replicate the error is the 6-month average estimate minus 0.04;
# RMSE = sqrt(mean over replicates of error^2). The synthetic design is
# chosen once from the history (it cannot depend on the future).
# Writes res_design.rds. The wall time of the design search is printed below but not
# stored (a stored time would change the file's bytes on every rerun). Run from this
# directory:  Rscript 06_design.R
# =====================================================================

source("sc_fns.R")
source("est_seeds.R")
S4 <- readRDS("aster_S4.rds"); tr <- readRDS("truth22.rds")
m4 <- S4$meta
J <- m4$N; T_pre <- m4$T_pre; m_max <- m4$m_pilot; R <- m4$R; tau <- tr$S4$tau
stopifnot(tau == m4$tau)
X  <- t(S4$Y0[, 1:T_pre])                  # 24 x 40, metros in columns
xb <- rowMeans(X)                           # Xbar with f_j = 1/J
cores <- max(1L, min(8L, parallel::detectCores() - 1L))

# ---- the synthetic design: problem (6) by local search over treated sets ----
fit_set <- function(S) {
  w <- sc_weights(xb, X[, S, drop = FALSE], diag(T_pre))
  C <- setdiff(seq_len(J), S)
  v <- sc_weights(xb, X[, C, drop = FALSE], diag(T_pre))
  ww <- numeric(J); ww[S] <- w
  vv <- numeric(J); vv[C] <- v
  list(S = sort(S), w = ww, v = vv,
       obj = sum((xb - drop(X %*% ww))^2) + sum((xb - drop(X %*% vv))^2))
}
local_search <- function(S) {
  cur <- fit_set(S)
  repeat {
    S <- cur$S; out <- setdiff(seq_len(J), S)
    moves <- c(
      unlist(lapply(S, function(i) lapply(out, function(o) sort(c(setdiff(S, i), o)))), recursive = FALSE),
      if (length(S) < m_max) lapply(out, function(o) sort(c(S, o))),
      if (length(S) > 1) lapply(S, function(i) setdiff(S, i)))
    fits <- lapply(moves, fit_set)
    objs <- vapply(fits, `[[`, numeric(1), "obj")
    if (min(objs) < cur$obj - 1e-14) cur <- fits[[which.min(objs)]] else break
  }
  cur
}
# greedy forward start: add the metro that lowers the objective most, up to m_max
greedy <- integer(0)
for (k in seq_len(m_max)) {
  cand <- setdiff(seq_len(J), greedy)
  objs <- vapply(cand, function(o) fit_set(c(greedy, o))$obj, numeric(1))
  greedy <- c(greedy, cand[which.min(objs)])
}
set.seed(est_seeds22[["s4_search_starts"]])
n_starts <- 100L
starts <- c(list(sort(greedy)), lapply(seq_len(n_starts), function(k) sort(sample.int(J, m_max))))
t0 <- proc.time()[["elapsed"]]
found <- parallel::mclapply(starts, local_search, mc.cores = cores)
search_sec <- proc.time()[["elapsed"]] - t0
objs <- vapply(found, `[[`, numeric(1), "obj")
best <- found[[which.min(objs)]]
treated <- which(best$w > 1e-8)
if (length(treated) < length(best$S)) best <- fit_set(treated)   # drop zero-weight members (never raises obj)
treated <- which(best$w > 1e-8)
controls <- which(best$v > 1e-8)
n_starts_at_best <- sum(abs(objs - min(objs)) < 1e-12)

# ---- evaluation over the 500 replicate pilot windows ----------------------
reps <- S4$Y0_pilot_reps                   # R x J x T_pilot (untreated outcomes)
contrast <- best$w - best$v                # sum_j (w_j - v_j) Y_jt
err_syn <- vapply(seq_len(R), function(r) mean(drop(contrast %*% reps[r, , ])), numeric(1))
# (tau enters as sum_j w_j tau = tau, so estimate - tau = the untreated contrast)
est_syn <- tau + err_syn

set.seed(est_seeds22[["s4_random_draws"]])
draws <- lapply(seq_len(R), function(r) sort(sample.int(J, m_max)))
hist_mean <- rowMeans(S4$Y0[, 1:T_pre])
err_dim <- err_did <- numeric(R)
for (r in seq_len(R)) {
  D <- draws[[r]]
  pil <- rowMeans(reps[r, , ])                           # pilot-window mean of Y0, by metro
  err_dim[r] <- mean(pil[D]) - mean(pil[-D])
  err_did[r] <- err_dim[r] - (mean(hist_mean[D]) - mean(hist_mean[-D]))
}
rmse <- function(e) sqrt(mean(e^2))

# fit diagnostics of the chosen design (history months)
fit_tr <- drop(X %*% best$w); fit_co <- drop(X %*% best$v)
res <- list(
  rmse_synthetic_design = rmse(err_syn),
  rmse_random_design    = rmse(err_dim),            # the paper's RND: random assignment + difference in means
  rmse_random_did       = rmse(err_did),            # random assignment + difference in differences
  bias_synthetic = mean(err_syn), sd_synthetic = sd(err_syn),
  bias_random = mean(err_dim), sd_random = sd(err_dim), bias_random_did = mean(err_did), sd_random_did = sd(err_did),
  ratio_random_to_synthetic = rmse(err_dim) / rmse(err_syn), ratio_random_did_to_synthetic = rmse(err_did) / rmse(err_syn),
  est_rep1 = est_syn[1], est_syn = est_syn, err_syn = err_syn, err_random = err_dim, err_random_did = err_did,
  synthetic = list(treated = S4$unit_names[treated], w = setNames(best$w[treated], S4$unit_names[treated]),
                   controls = S4$unit_names[controls], v = setNames(best$v[controls], S4$unit_names[controls]),
                   obj = best$obj, n_treated = length(treated), n_controls = length(controls),
                   rmse_fit_treated = sqrt(mean((xb - fit_tr)^2)), rmse_fit_control = sqrt(mean((xb - fit_co)^2)),
                   rmse_fit_gap = sqrt(mean((fit_tr - fit_co)^2)),
                   path_treated = fit_tr, path_control = fit_co, xbar = xb,
                   n_starts = length(starts), n_starts_at_best = n_starts_at_best,
                   greedy = S4$unit_names[sort(greedy)], greedy_obj = fit_set(greedy)$obj),
  random_draws = draws,
  tau = tau, reps = R, m_pilot = m_max, n_metros = J, n_history = T_pre, n_pilot_months = m4$T_pilot
)
saveRDS(res, "res_design.rds")

f4 <- function(v) sprintf("%.4f", v)
cat("================ S4: pilot design (tau =", tau, ",", R, "replicate pilot windows) ================\n")
cat("synthetic design (Abadie-Zhao (6), m <= 5): treated", paste(res$synthetic$treated, collapse = ", "), "\n")
cat("  treated weights", paste(f4(res$synthetic$w), collapse = " "), "|", res$synthetic$n_controls, "control metros\n")
cat("  history fit RMSE: treated", f4(res$synthetic$rmse_fit_treated), "control", f4(res$synthetic$rmse_fit_control),
    "| best of", length(starts), "starts (", n_starts_at_best, "reached it ) in", sprintf("%.1f s", search_sec), "\n")
cat("RMSE synthetic design        ", f4(res$rmse_synthetic_design), " (bias", f4(res$bias_synthetic), ")\n")
cat("RMSE random + diff in means  ", f4(res$rmse_random_design), " (bias", f4(res$bias_random), ")\n")
cat("RMSE random + DiD            ", f4(res$rmse_random_did), " (bias", f4(res$bias_random_did), ")\n")
cat("replicate 1 (the stored pilot window): synthetic-design estimate", f4(res$est_rep1), "\n")
