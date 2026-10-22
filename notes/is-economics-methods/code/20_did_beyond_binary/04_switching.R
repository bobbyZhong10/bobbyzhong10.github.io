# ------------------------------------------------------------------
# 04_switching.R -- Block C (switching, non-binary treatment): static TWFE
# against the dynamic DID_l estimators of de Chaisemartin & D'Haultfoeuille.
#
#   (1) fixest::feols          static TWFE: y on the bonus tier, city and week FE,
#                              clustered by city
#   (2) TwoWayFEWeights        weights the TWFE coefficient puts on the treated cells
#   (3) DIDmultiplegtDYN       did_multiplegt_dyn(): DID_l (l = 1..8), placebos (1..4),
#                              the average total effect, and the normalised effects
#
# Data: tessera_C.rds (200 cities x 30 weeks; unit = city, `dose` = tier 0/1/2).
# Truth: truth20.rds$C (delta, delta_norm, avg_total_effect, avg_effect).
#
# ---- Return-object mapping (TwoWayFEWeights 2.1.0, DIDmultiplegtDYN 2.4.0) ------
# twowayfeweights() ($dat_result data.frame(T, G, weight); $nr_plus, $nr_minus,
#   $nr_weights, $sum_plus, $sum_minus, $beta, $sensibility, $sensibility2, $tot_cells)
#   -> res$weights: counts, sums, share_negative = nr_minus / nr_weights (by count),
#      the TWFE coefficient rebuilt from the weights, tables by weeks since the first
#      switch and by tier. For a non-binary D the weights apply to the per-unit-of-dose
#      effects tau_gt / D_gt of the cells with D > 0 (dCDH 2020).
# did_multiplegt_dyn() returns a list of class "did_multiplegt_dyn":
#   $results$Effects    8 x 8 matrix, rows Effect_1..Effect_8 (= DID_l; with
#                       normalized = TRUE the normalised effect DID^n_l), columns
#                       Estimate, SE, LB CI, UB CI, N, Switchers, N.w, Switchers.w
#   $results$Placebos   4 x 8 matrix, rows Placebo_1..Placebo_4
#   $results$ATE        1 x 8 matrix, row Av_tot_eff: average total effect per unit of
#                       treatment (identical in the plain and the normalised call)
#   $results$delta_D_avg_total   average number of periods over which the dose accumulates
#   $results$p_jointeffects, p_jointplacebo   joint nullity tests (H0: all effects, resp. all
#                       placebos, are zero), computed by the package as 1 - pchisq(b' V^- b, k).
#                       The subtraction returns exactly 0 when the true p-value is below about
#                       1e-16 (floating point), so p_jointeffects == 0 means "p < 1e-16", not
#                       "p = 0"; we store it as reported and flag it (p_jointeffects_underflow).
#                       We cannot rebuild the package's statistic from $coef$vcov: the placebo
#                       block of that matrix gives a different chi-square (2.09 against the
#                       package's 5.03), so the package's own p-values are used throughout.
#   $coef$b, $coef$vcov stacked Effect_1..Effect_8, Placebo_1..Placebo_4 and covariance
#   -> res$dyn$l1..l8 = list(est, se, lo, hi, n, switchers); res$dyn$effects (data.frame);
#      res$dyn$placebos; res$dyn$avg_total; res$dyn$norm (the same for the normalised call);
#      res$dyn$p_jointeffects, res$dyn$p_jointplacebo (as the package reports them).
#
# Horizon convention: DID_l compares the outcome l-1 weeks after a city's first switch
# (week g) with the last week before it (g - 1), against not-yet-switched cities with the
# same status quo (0); Placebo_l compares g - 1 - l with g - 1. The chapter's l is
# `ell`; the stored truth truth20$C$delta[l] is the true effect DID_l estimates.
# did_multiplegt_dyn needs polars attached first.
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(fixest); library(TwoWayFEWeights)
  library(polars); library(DIDmultiplegtDYN)
})

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
C  <- readRDS(file.path(code_dir, "tessera_C.rds"))
tr <- readRDS(file.path(code_dir, "truth20.rds"))
df <- as.data.frame(C[, .(unit, t, y, dose)])

warn_log <- character()
with_log <- function(expr) withCallingHandlers(expr, warning = function(w) {
  warn_log <<- c(warn_log, conditionMessage(w)); invokeRestart("muffleWarning") })
timer <- function(label, expr) { t0 <- Sys.time(); v <- force(expr)
  cat(sprintf("  [%s: %.1f s]\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs")))); v }

## ================================================================
## (1) static TWFE
## ================================================================
m_twfe <- feols(y ~ dose | unit + t, data = C, cluster = ~unit)
m_orac <- feols(tau_true ~ dose | unit + t, data = C, cluster = ~unit)   # TWFE on the noise-free effect column
ci <- confint(m_twfe)["dose", ]
twfe <- list(est = coef(m_twfe)[["dose"]], se = se(m_twfe)[["dose"]], lo = ci[[1]], hi = ci[[2]],
             p = pvalue(m_twfe)[["dose"]], df = degrees_freedom(m_twfe, "t"),
             oracle_est = coef(m_orac)[["dose"]], n_obs = nobs(m_twfe))

## ================================================================
## (2) weights
## ================================================================
tw <- with_log(twowayfeweights(df, Y = "y", G = "unit", T = "t", D = "dose", type = "feTR", summary_measures = TRUE))
cells <- merge(as.data.table(tw$dat_result)[, .(unit = G, t = T, weight)], C[, .(unit, t, g, dose, tau_true)], by = c("unit", "t"))
wt <- cells[dose > 0]
wt[, event_time := t - g]                         # weeks since the city's first switch (g > 0 for every treated cell)
stopifnot(nrow(wt) == tw$tot_cells, abs(sum(wt$weight) - 1) < 1e-8, sum(wt$weight < 0) == tw$nr_minus,
          all(wt$g > 0), all(cells$weight[cells$dose == 0] == 0))
implied <- wt[, sum(weight * tau_true / dose)]
## The package weight is eps_gt D_gt / sum(eps D), eps = residual of D on unit and week FE. The weights cover
## the cells with D > 0 only; cells that are paused (D = 0) but still carry an effect (carryover, tau_true > 0)
## enter the coefficient as "controls". Exactly, for the noise-free column,
##   oracle = sum_{D>0} w_gt tau_gt / D_gt  +  sum_{D=0} eps_gt tau_gt / sum(eps D).
eps_all <- resid(feols(dose ~ 1 | unit + t, data = C))
C[, eps := eps_all]
Dtot <- C[, sum(eps * dose)]
carry <- C[dose == 0, sum(eps * tau_true)] / Dtot
stopifnot(max(abs(wt$weight - merge(wt[, .(unit, t)], C[, .(unit, t, eps, dose)], by = c("unit", "t"))[, eps * dose / Dtot])) < 1e-10,
          abs(implied + carry - twfe$oracle_est) < 1e-8, abs(tw$beta - twfe$est) < 1e-8)
by_event <- wt[, .(n_cells = .N, n_negative = sum(weight < 0), sum_weight = sum(weight),
                   sum_negative_weight = sum(weight[weight < 0])), keyby = event_time]
by_dose  <- wt[, .(n_cells = .N, n_negative = sum(weight < 0), sum_weight = sum(weight),
                   sum_negative_weight = sum(weight[weight < 0]), mean_effect_per_dose = mean(tau_true / dose)), keyby = dose]
weights <- list(n_cells = tw$tot_cells, n_positive = tw$nr_plus, n_negative = tw$nr_minus,
                share_negative = tw$nr_minus / tw$nr_weights,
                sum_positive = tw$sum_plus, sum_negative = tw$sum_minus,
                min_sd_beta_delta0 = tw$sensibility, min_sd_opposite_sign = tw$sensibility2,
                implied_beta = implied,                       # sum of weight x per-dose effect over cells with D > 0
                carryover_term = carry,                       # contribution of paused cells that still carry an effect
                contrib_positive = wt[weight > 0, sum(weight * tau_true / dose)],
                contrib_negative = wt[weight < 0, sum(weight * tau_true / dose)],
                mean_effect_per_dose = wt[, mean(tau_true / dose)],                     # unweighted, over treated cells
                dose_weighted_effect_per_dose = wt[, sum(tau_true) / sum(dose)],       # = truth20$C$avg_effect
                by_event = as.data.frame(by_event), by_dose = as.data.frame(by_dose))
stopifnot(abs(weights$dose_weighted_effect_per_dose - tr$C$avg_effect) < 1e-10,
          abs(weights$mean_effect_per_dose - tr$C$avg_effect_unweighted) < 1e-10)

## ================================================================
## (3) did_multiplegt_dyn
## ================================================================
cat("== did_multiplegt_dyn ==\n")
L <- 8; P <- 4
run_dyn <- function(normalized) with_log(did_multiplegt_dyn(df, outcome = "y", group = "unit", time = "t", treatment = "dose",
                                                           effects = L, placebo = P, cluster = "unit",
                                                           normalized = normalized, graph_off = TRUE))
dyn_r  <- timer("DID_l", run_dyn(FALSE))
dynn_r <- timer("normalised DID^n_l", run_dyn(TRUE))

tidy <- function(mat, l) data.frame(l = l, est = mat[, "Estimate"], se = mat[, "SE"], lo = mat[, "LB CI"], hi = mat[, "UB CI"],
                                    n = mat[, "N"], switchers = mat[, "Switchers"], row.names = NULL)
package_dyn <- function(r) {
  E <- tidy(r$results$Effects, 1:L); Pl <- tidy(r$results$Placebos, 1:P); A <- r$results$ATE
  stopifnot(nrow(E) == L, nrow(Pl) == P, nrow(A) == 1)
  out <- list(effects = E, placebos = Pl,
              avg_total = list(est = A[1, "Estimate"], se = A[1, "SE"], lo = A[1, "LB CI"], hi = A[1, "UB CI"],
                               n = A[1, "N"], switchers = A[1, "Switchers"]),
              delta_D_avg_total = r$results$delta_D_avg_total,
              p_jointeffects = r$results$p_jointeffects, p_jointplacebo = r$results$p_jointplacebo,
              p_jointeffects_underflow = isTRUE(r$results$p_jointeffects == 0),   # 1 - pchisq() hit 0: p < 1e-16
              df_jointeffects = L, df_jointplacebo = P)
  stopifnot(is.finite(out$p_jointeffects), is.finite(out$p_jointplacebo))
  for (l in 1:L) out[[paste0("l", l)]] <- as.list(E[l, c("est", "se", "lo", "hi", "n", "switchers")])
  for (l in 1:P) out[[paste0("pl", l)]] <- as.list(Pl[l, c("est", "se", "lo", "hi", "n", "switchers")])
  out
}
dyn <- package_dyn(dyn_r)
dyn$norm <- package_dyn(dynn_r)
dyn$vcov <- dyn_r$coef$vcov
stopifnot(abs(dyn$norm$avg_total$est - dyn$avg_total$est) < 1e-10)   # Av_tot_eff does not depend on `normalized`

res_switch <- list(twfe = twfe, weights = weights, dyn = dyn,
                   versions = c(TwoWayFEWeights = as.character(packageVersion("TwoWayFEWeights")),
                                DIDmultiplegtDYN = as.character(packageVersion("DIDmultiplegtDYN")),
                                fixest = as.character(packageVersion("fixest"))),
                   warnings = warn_log)
saveRDS(res_switch, file.path(code_dir, "res_switch.rds"))

## ---- report -------------------------------------------------------------------------------
cat(sprintf("\n== Block C truth: avg_effect %.4f (unweighted %.4f); avg total effect %.4f ==\n",
            tr$C$avg_effect, tr$C$avg_effect_unweighted, tr$C$avg_total_effect))
cat(sprintf("static TWFE %.4f (se %.4f) [%.4f, %.4f]; noise-free TWFE %.4f; gap to avg_effect %.2f se\n",
            twfe$est, twfe$se, twfe$lo, twfe$hi, twfe$oracle_est, (twfe$est - tr$C$avg_effect) / twfe$se))
cat(sprintf("weights: %d positive, %d negative of %d treated cells (share %.4f); sums %.4f / %.4f\n",
            weights$n_positive, weights$n_negative, weights$n_cells, weights$share_negative,
            weights$sum_positive, weights$sum_negative))
cat(sprintf("  noise-free TWFE %.4f = weighted per-dose effects over treated cells %.4f + carryover term %.4f\n",
            twfe$oracle_est, weights$implied_beta, weights$carryover_term))
cat(sprintf("  contributions to beta: positive-weight cells %.4f, negative-weight cells %.4f\n",
            weights$contrib_positive, weights$contrib_negative))
show <- function(E, truth) { E$truth <- truth; E$covers <- E$lo <= truth & truth <= E$hi; print(E[, c("l", "est", "se", "lo", "hi", "truth", "covers")], digits = 4, row.names = FALSE) }
cat("\nDID_l vs delta_l:\n"); show(dyn$effects, tr$C$delta)
cat("\nnormalised DID^n_l vs delta_norm:\n"); show(dyn$norm$effects, tr$C$delta_norm)
cat("\nplacebos:\n"); print(dyn$placebos, digits = 4, row.names = FALSE)
cat(sprintf("\njoint placebo test: p %.4f (df %d); joint effects test: p %g%s (df %d)\n",
            dyn$p_jointplacebo, P, dyn$p_jointeffects, if (dyn$p_jointeffects_underflow) " (underflow: p < 1e-16)" else "", L))
cat(sprintf("average total effect %.4f (se %.4f) [%.4f, %.4f] vs truth %.4f; delta_D %.3f\n",
            dyn$avg_total$est, dyn$avg_total$se, dyn$avg_total$lo, dyn$avg_total$hi, tr$C$avg_total_effect, dyn$delta_D_avg_total))
cat("\nPackage warnings:\n"); if (length(warn_log)) print(table(warn_log)) else cat("none\n")
cat("Saved res_switch.rds\n")
