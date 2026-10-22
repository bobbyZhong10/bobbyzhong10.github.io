# ------------------------------------------------------------------
# 07_dcdh.R -- de Chaisemartin & D'Haultfoeuille diagnostics and estimator
# on the Kestrel Marketplace panel.
#
#   (1) TwoWayFEWeights::twowayfeweights(type = "feTR"):
#       the weights that the static TWFE coefficient puts on the
#       treated (g,t) cells; negative weights are the diagnostic.
#   (2) DIDmultiplegtDYN::did_multiplegt_dyn(): the dynamic
#       DID_l estimators (effects = 12) and placebos (placebo = 6),
#       clustered by city.
#
# Return-object mapping (packages inspected 2026-10-20: TwoWayFEWeights
# 2.1.0, DIDmultiplegtDYN 2.4.0).
#
# twowayfeweights() returns a list:
#   $dat_result  data.table with columns G, T, weight for EVERY (G,T) cell
#                (G = group = city_id, T = month); untreated cells have
#                weight 0, so we keep the D == 1 cells (tot_cells of them).
#   $nr_plus, $nr_minus, $nr_weights     number of positive / negative /
#                all weights (cells with a nonzero weight)
#   $sum_plus, $sum_minus                sums of positive / negative weights
#   $beta                                the static TWFE coefficient
#   $sensibility, $sensibility2          Corollary 1 summary measures
#                (dCDH 2020a): min sd of the effects compatible with
#                beta_fe and Delta_TR = 0, and with an effect of the sign
#                opposite to beta_fe in every cell.
#   -> res_dcdh$twfe_weights : data.frame, one row per treated cell
#                (city_id, month, cohort, event_time, weight)
#   -> res_dcdh$twfe_summary : the scalars above plus shares and
#                breakdowns by cohort and by event time
#
# did_multiplegt_dyn() returns a list of class "did_multiplegt_dyn":
#   $results$Effects   12 x 8 matrix, rows Effect_1..Effect_12 (= DID_l),
#                      columns Estimate, SE, LB CI, UB CI, N, Switchers,
#                      N.w, Switchers.w
#   $results$Placebos  6 x 8 matrix, rows Placebo_1..Placebo_6, same columns
#   $results$ATE       1 x 8 matrix, "Av_tot_eff": average cumulative
#                      (total) effect per treatment unit
#   $results$delta_D_avg_total   average number of periods over which the
#                      effect accumulates (6.5 = mean of 1..12)
#   $results$p_jointeffects / p_jointplacebo   joint nullity p-values
#   $coef$b, $coef$vcov   stacked estimates and their covariance matrix
#   $N_switchers_effect_l number of switchers behind DID_l
#   -> res_dcdh$dyn : list(effects, placebos, avg_total, p_jointeffects,
#                      p_jointplacebo, delta_D_avg_total, vcov)
#
# Horizon convention: DID_l compares the outcome l-1 periods after a
# city's first treated month F_g with the last pre-treatment month F_g - 1
# (so DID_1 is the impact effect, event time e = 0, and DID_l <-> e = l-1).
# Placebo_l compares F_g - 1 - l with F_g - 1, i.e. event time e = -1 - l
# in an event study whose reference period is e = -1.
#
# did_multiplegt_dyn needs the polars package attached beforehand, else
# it stops with "object 'pl' not found".
# ------------------------------------------------------------------
library(data.table)
library(TwoWayFEWeights)
library(polars)            # must be attached before DIDmultiplegtDYN
library(DIDmultiplegtDYN)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dt <- readRDS(file.path(code_dir, "kestrel_panel.rds"))
df <- as.data.frame(dt[, .(city_id, month, lgmv, D)])

## collect warnings so they can be reported, but still show them
warn_log <- character()
with_log <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    warn_log <<- c(warn_log, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

## ---- (1) TWFE weights --------------------------------------------
tw <- with_log(twowayfeweights(df, Y = "lgmv", G = "city_id", T = "month", D = "D",
                               type = "feTR", summary_measures = TRUE))
cat("== TWFE weights (type feTR) ==\n")
print(tw)

cells <- merge(as.data.frame(tw$dat_result),
               as.data.frame(dt[, .(G = city_id, T = month, D, cohort = g)]),
               by = c("G", "T"))
twfe_weights <- cells[cells$D == 1, ]
twfe_weights <- data.frame(city_id = twfe_weights$G, month = twfe_weights$T,
                           cohort = twfe_weights$cohort,
                           event_time = twfe_weights$T - twfe_weights$cohort,
                           weight = twfe_weights$weight)
twfe_weights <- twfe_weights[order(twfe_weights$city_id, twfe_weights$month), ]
rownames(twfe_weights) <- NULL
stopifnot(nrow(twfe_weights) == tw$tot_cells,
          abs(sum(twfe_weights$weight) - 1) < 1e-8,
          sum(twfe_weights$weight < 0) == tw$nr_minus,
          # untreated cells carry zero weight
          all(cells$weight[cells$D == 0] == 0))

by_cohort <- as.data.table(twfe_weights)[, .(n_cells = .N,
                                             n_negative = sum(weight < 0),
                                             sum_weight = sum(weight),
                                             sum_negative_weight = sum(weight[weight < 0])),
                                         keyby = cohort]
by_event <- as.data.table(twfe_weights)[, .(n_cells = .N,
                                            n_negative = sum(weight < 0),
                                            sum_weight = sum(weight),
                                            sum_negative_weight = sum(weight[weight < 0])),
                                        keyby = event_time]
cat("\nTWFE weights by cohort:\n"); print(by_cohort)
cat("\nTWFE weights by event time (selected):\n")
print(by_event[event_time %in% c(0, 3, 6, 9, 12, 15, 18, 21, 23)])

twfe_summary <- list(
  type            = tw$type,
  beta            = tw$beta,
  n_cells         = tw$tot_cells,
  n_positive      = tw$nr_plus,
  n_negative      = tw$nr_minus,
  share_negative  = tw$nr_minus / tw$nr_weights,   # by count
  sum_positive    = tw$sum_plus,
  sum_negative    = tw$sum_minus,
  min_sd_beta_delta0   = tw$sensibility,           # sensibility
  min_sd_opposite_sign = tw$sensibility2,          # sensibility2
  by_cohort       = as.data.frame(by_cohort),
  by_event        = as.data.frame(by_event)
)
cat(sprintf("\nNegative weights: %d of %d cells (%.1f%%), sum %.4f; positive sum %.4f\n",
            twfe_summary$n_negative, twfe_summary$n_cells, 100 * twfe_summary$share_negative,
            twfe_summary$sum_negative, twfe_summary$sum_positive))

## ---- (2) dCDH dynamic estimators ---------------------------------
dyn <- with_log(did_multiplegt_dyn(df, outcome = "lgmv", group = "city_id", time = "month",
                                   treatment = "D", effects = 12, placebo = 6,
                                   cluster = "city_id", graph_off = TRUE))
print(dyn)

tidy_mat <- function(mat, l, e) {
  out <- data.frame(l = l, e = e,
                    estimate = mat[, "Estimate"], se = mat[, "SE"],
                    lb = mat[, "LB CI"], ub = mat[, "UB CI"],
                    n = mat[, "N"], switchers = mat[, "Switchers"],
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
E <- dyn$results$Effects; P <- dyn$results$Placebos; A <- dyn$results$ATE
stopifnot(nrow(E) == 12, nrow(P) == 6, nrow(A) == 1)
effects  <- tidy_mat(E, l = 1:12, e = 0:11)
placebos <- tidy_mat(P, l = 1:6,  e = -(2:7))
avg_total <- list(estimate = A[1, "Estimate"], se = A[1, "SE"],
                  lb = A[1, "LB CI"], ub = A[1, "UB CI"],
                  n = A[1, "N"], switchers = A[1, "Switchers"])

dyn_out <- list(effects = effects, placebos = placebos, avg_total = avg_total,
                p_jointeffects = dyn$results$p_jointeffects,
                p_jointplacebo = dyn$results$p_jointplacebo,
                delta_D_avg_total = dyn$results$delta_D_avg_total,
                vcov = dyn$coef$vcov)

cat("\nDID_l effects (l -> e = l-1):\n"); print(effects, digits = 4)
cat("\nPlacebos (l -> e = -1-l):\n");     print(placebos, digits = 4)
cat(sprintf("\nAverage total effect: %.4f (SE %.4f); joint p effects %.4f, placebos %.4f\n",
            avg_total$estimate, avg_total$se, dyn_out$p_jointeffects, dyn_out$p_jointplacebo))

if (length(warn_log)) {
  cat("\nPackage warnings:\n"); print(table(warn_log))
} else cat("\nNo package warnings.\n")

saveRDS(list(twfe_weights = twfe_weights, twfe_summary = twfe_summary, dyn = dyn_out,
             versions = c(TwoWayFEWeights = as.character(packageVersion("TwoWayFEWeights")),
                          DIDmultiplegtDYN = as.character(packageVersion("DIDmultiplegtDYN"))),
             warnings = warn_log),
        file.path(code_dir, "res_dcdh.rds"))
cat("Saved res_dcdh.rds\n")
