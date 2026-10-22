# ------------------------------------------------------------------
# 03_dose.R -- Block B (continuous dose): what each estimator recovers.
#
#   (1) contdid::cont_did        ATT(d|d) curve, overall ATT, "ACRT" curve and
#                                overall "ACRT" (Callaway, Goodman-Bacon & Sant'Anna)
#   (2) fixest::feols            binarised DiD, and TWFE of y on the dose received
#                                (unit and time FE, clustered by city)
#   (3) TwoWayFEWeights          weights behind the TWFE dose coefficient
#   (4) DIDHAD::did_had          heterogeneous adoption design on tessera_B_had
#                                (no untreated city; quasi-untreated cities as controls)
#   (5) DIDmultiplegtSTAT        AOSS and WAOSS on tessera_B (stayers = dose 0)
#
# Data: tessera_B.rds (500 cities, t = 1, 2; `dose` constant within city, `d` =
# dose * post the dose received) and tessera_B_had.rds. Truth: truth20.rds$B.
#
# ---- Return-object mappings (packages inspected 2026-10-20: contdid 0.1.1,
# ---- DIDHAD 2.0.0, DIDmultiplegtSTAT 1.0.0, TwoWayFEWeights 2.1.0, fixest 0.14.2)
#
# cont_did() returns an object of class "dose_obj" (a list):
#   $dose            the grid dvals at which the curves are evaluated
#   $att.d, $att.d_se, $att.d_crit.val     ATT(d|d) at dvals, its se and the
#                    critical value (qnorm(.975) = 1.96 when cband = FALSE)
#   $acrt.d, $acrt.d_se, $acrt.d_crit.val  derivative of the fitted curve at dvals
#   $overall_att, $overall_att_se          average of ATT(d|d) over the treated
#   $overall_acrt, $overall_acrt_se        average derivative over the treated
#   (the *_inffunc entries are influence functions; pte_params echoes the call)
# -> res$cont$curve       data.frame(dose, est, se, lo, hi, truth, truth_strong)
#    res$cont$acrt_curve  data.frame(dose, est, se, lo, hi, truth_slope, truth_acrt, ...)
#    res$cont$overall     list(est, se, lo, hi)  [lo, hi = est -/+ qnorm(.975) se]
#    res$cont$acrt_overall list(est, se, lo, hi)
# The `level` and `slope` calls return the same point estimates; each is read
# for its own target (ATT from `level`, ACRT from `slope`).
#
# Three facts about contdid 0.1.1 that matter for the chapter (all checked here):
#  (a) The fitted curve is a global B-spline (degree 3, no knots = a cubic
#      polynomial in the dose) of the treated cities' outcome change on the
#      dose, minus the untreated cities' mean change. ATT(d|d) is identified
#      under parallel trends alone; ACRT(d|d) (the causal response) additionally
#      needs strong parallel trends.
#  (b) The grid basis is built with boundary knots = range(dvals), whereas the
#      fit uses range(treated doses). If dvals does not span the treated doses
#      exactly, att.d and acrt.d are evaluated on a different basis and are
#      WRONG (the default dvals, the 10th to 99th percentiles, gives a curve up
#      to 1.16 below the cubic fit). The overall ATT and ACRT are unaffected.
#      We therefore pass dvals that start at the smallest and end at the
#      largest treated dose, and check the curve against an lm() replication.
#  (c) The standard errors come from a multiplier bootstrap, so they depend on
#      the seed (+/- 8% at 999 draws); we use 9999 draws and a fixed seed. The
#      pointwise se of the curves agree with a city bootstrap of the
#      replication; the overall ACRT se does not: the package's is smaller by
#      about n / n_treated = 1.35 (see $cont$check), so $cont$acrt_overall_boot
#      carries the city-bootstrap interval, which the chapter should quote.
# The "ACRT" curve of contdid is the derivative of the ATT(d|d) fit, i.e.
# d ATT(d|d)/dd = ACRT(d|d) + selection bias; the causal ACRT(d|d) is NOT
# identified by parallel trends alone (truth20$B$acrt_overall vs acrt_overall_pt_pop).
#
# did_had() returns list(args, results, plot):
#   $results$resmat  matrix, row Effect_1, columns Estimate, SE, LB.CI, UB.CI,
#                    N, BW, N.BW, T, p.val, ID: Estimate = WAS of the treated per unit
#                    dose at period 2, BW the bandwidth below which cities count as
#                    quasi-untreated, N.BW the number of such cities, T and p.val the QUG
#                    test (H0: the lowest dose is zero, i.e. quasi-untreated cities exist)
#   $results$yatchew_test  (only with yatchew = TRUE) sigma2_lin, sigma2_diff, T_hr,
#                    p-value, N: heteroskedasticity-robust Yatchew test of linearity of
#                    E(outcome change | dose); under linearity TWFE is unbiased for the WAS
# -> res$had   list(est, se, lo, hi, n, bw, n_bw, qug_stat, qug_p, yatchew = list(...), twfe)
#
# did_multiplegt_stat() returns list(args, results); results$table is a 6 x 6 matrix with
# rows AOSS, aoss_2, WAOSS, waoss_2, IVWAOSS, ivwaoss_2 (the *_2 rows repeat the
# estimate for the single period pair) and columns Estimate, SE, LB CI, UB CI, Switchers,
# Stayers; results$aoss_vs_waoss (aoss_vs_waoss = TRUE) is a 1 x 6 matrix, the test
# H0: AOSS = WAOSS (Estimate, SE, LB CI, UB CI, t stat., pval.); its SE column is
# sqrt(N) times the se implied by CI and t (see the note at the extraction).
# -> res$stat  list(aoss, waoss, diff)
#
# twowayfeweights() (type "feTR"): $dat_result data.frame(T, G, weight) for every city-period
# cell, $nr_plus / $nr_minus / $nr_weights counts, $sum_plus / $sum_minus, $beta (= the TWFE
# coefficient), $sensibility / $sensibility2 (dCDH 2020 Corollary 1). For a non-binary D the
# weights apply to the per-unit-of-dose effects tau_gt / D_gt of the treated cells.
# -> res$twfe$weights
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(fixest); library(splines2)
  library(contdid); library(DIDHAD); library(DIDmultiplegtSTAT); library(TwoWayFEWeights)
})

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
B  <- readRDS(file.path(code_dir, "tessera_B.rds"))
H  <- readRDS(file.path(code_dir, "tessera_B_had.rds"))
tr <- readRDS(file.path(code_dir, "truth20.rds"))
PB <- tr$params$B

BITERS_CONT <- 9999L   # multiplier-bootstrap draws inside cont_did
B_CITY      <- 999L    # city bootstrap of the lm() replication (check on the package's se)
SEED_CONT   <- 20301L
SEED_BOOT   <- 20302L

warn_log <- character(); msg_log <- character()
with_log <- function(expr) withCallingHandlers(expr,
  warning = function(w) { warn_log <<- c(warn_log, conditionMessage(w)); invokeRestart("muffleWarning") },
  message = function(m) { msg_log <<- c(msg_log, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })
timer <- function(label, expr) { t0 <- Sys.time(); v <- force(expr)
  cat(sprintf("  [%s: %.1f s]\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs")))); v }
ci_norm <- function(est, se) list(est = est, se = se, lo = est - qnorm(0.975) * se, hi = est + qnorm(0.975) * se)

## true curves at any dose (closed forms of 01_dgp.R; checked below against truth20$B$att_curve)
Eb <- function(d) 1 + PB$kappa * (d - PB$dbar)
m  <- function(d) PB$a1 * d - PB$a2 * d^2
mp <- function(d) PB$a1 - 2 * PB$a2 * d
truth_at <- function(d) data.frame(
  truth = Eb(d) * m(d),                      # ATT(d|d)
  truth_strong = m(d),                       # ATT(d)
  truth_slope = Eb(d) * mp(d) + PB$kappa * m(d),   # d ATT(d|d)/dd (what parallel trends identifies)
  truth_acrt = Eb(d) * mp(d),                # ACRT(d|d)
  truth_acrt_strong = mp(d))                 # ACRT(d)
cur <- tr$B$att_curve
chk <- truth_at(cur$dose)
stopifnot(max(abs(chk$truth - cur$att)) < 1e-10, max(abs(chk$truth_strong - cur$att_strong)) < 1e-10,
          max(abs(chk$truth_slope - cur$slope_att)) < 1e-10, max(abs(chk$truth_acrt - cur$acrt)) < 1e-10)

## ================================================================
## (1) contdid
## ================================================================
cat("== (1) contdid ==\n")
d_pos <- B[post == 1 & dose > 0, dose]
n_all <- uniqueN(B$unit); n_treated <- length(d_pos)
dv <- sort(unique(c(min(d_pos), seq(0.25, 2.95, by = 0.05), max(d_pos))))   # spans the treated doses exactly
stopifnot(min(dv) == min(d_pos), max(dv) == max(d_pos))
df_B <- as.data.frame(B)
run_cont <- function(target, degree = 3, num_knots = 0, biters = BITERS_CONT) {
  set.seed(SEED_CONT)
  with_log(cont_did(yname = "y", dname = "dose", gname = "g", tname = "t", idname = "unit", data = df_B,
                    target_parameter = target, aggregation = "dose", treatment_type = "continuous",
                    dvals = dv, degree = degree, num_knots = num_knots,
                    control_group = "nevertreated", biters = biters, cband = FALSE))
}
rc_lev <- timer("cont_did level", run_cont("level"))
rc_slp <- timer("cont_did slope", run_cont("slope"))
stopifnot(abs(rc_lev$overall_att - rc_slp$overall_att) < 1e-10, abs(rc_lev$overall_acrt - rc_slp$overall_acrt) < 1e-10,
          max(abs(rc_lev$att.d - rc_slp$att.d)) < 1e-10, max(abs(rc_lev$acrt.d - rc_slp$acrt.d)) < 1e-10,
          max(abs(rc_lev$dose - dv)) < 1e-12)

tru <- truth_at(dv)
curve <- data.frame(dose = dv, est = unname(rc_lev$att.d), se = rc_lev$att.d_se)
curve$lo <- curve$est - rc_lev$att.d_crit.val * curve$se
curve$hi <- curve$est + rc_lev$att.d_crit.val * curve$se
curve$truth <- tru$truth; curve$truth_strong <- tru$truth_strong
acrt_curve <- data.frame(dose = dv, est = unname(rc_slp$acrt.d), se = rc_slp$acrt.d_se)
acrt_curve$lo <- acrt_curve$est - rc_slp$acrt.d_crit.val * acrt_curve$se
acrt_curve$hi <- acrt_curve$est + rc_slp$acrt.d_crit.val * acrt_curve$se
acrt_curve$truth_slope <- tru$truth_slope; acrt_curve$truth_acrt <- tru$truth_acrt
acrt_curve$truth_acrt_strong <- tru$truth_acrt_strong

fit_stats <- function(est, truth, lo, hi) list(rmse = sqrt(mean((est - truth)^2)), mean_abs_err = mean(abs(est - truth)),
  max_abs_err = max(abs(est - truth)), mean_signed_err = mean(est - truth), share_covered = mean(lo <= truth & truth <= hi))
fit_att   <- fit_stats(curve$est, curve$truth, curve$lo, curve$hi)
fit_acrt_slope <- fit_stats(acrt_curve$est, acrt_curve$truth_slope, acrt_curve$lo, acrt_curve$hi)
fit_acrt_true  <- fit_stats(acrt_curve$est, acrt_curve$truth_acrt,  acrt_curve$lo, acrt_curve$hi)

## replication by lm() (fact (b)) and city bootstrap of the replication (fact (c))
wide <- dcast(B, unit + dose ~ t, value.var = "y"); setnames(wide, c("1", "2"), c("y1", "y2")); wide[, dy := y2 - y1]
bnd <- range(d_pos)
rep_est <- function(w, grid = dv) {
  tt <- w[dose > 0]; u <- w[dose == 0]
  bs <- as.data.frame(bSpline(tt$dose, degree = 3, knots = numeric(0), Boundary.knots = bnd)); bs$dy <- tt$dy
  fit <- lm(dy ~ ., data = bs); b <- coef(fit)
  G  <- bSpline(grid, degree = 3, knots = numeric(0), Boundary.knots = bnd)
  dG <- dbs(grid, degree = 3, knots = numeric(0), Boundary.knots = bnd)
  dd <- dbs(tt$dose, degree = 3, knots = numeric(0), Boundary.knots = bnd)
  c(att_o = mean(fitted(fit)) - mean(u$dy), acrt_o = mean(dd %*% b[-1]),
    drop(cbind(1, G) %*% b) - mean(u$dy), drop(dG %*% b[-1]))
}
e0 <- rep_est(wide); K <- length(dv)
rep_err <- c(att_o = abs(e0[[1]] - rc_lev$overall_att), acrt_o = abs(e0[[2]] - rc_lev$overall_acrt),
             att_d = max(abs(e0[3:(2 + K)] - rc_lev$att.d)), acrt_d = max(abs(e0[(3 + K):(2 + 2 * K)] - rc_lev$acrt.d)))
stopifnot(all(rep_err < 1e-9))
set.seed(SEED_BOOT)
bt <- timer("city bootstrap of the lm replication", t(replicate(B_CITY, rep_est(wide[sample(.N, replace = TRUE)]))))
se_b <- apply(bt, 2, sd)
q_att  <- quantile(bt[, 1], c(0.025, 0.975), names = FALSE); q_acrt <- quantile(bt[, 2], c(0.025, 0.975), names = FALSE)
se_ratio_curve_att  <- median(se_b[3:(2 + K)] / curve$se)
se_ratio_curve_acrt <- median(se_b[(3 + K):(2 + 2 * K)] / acrt_curve$se)

overall <- ci_norm(rc_lev$overall_att, rc_lev$overall_att_se)
acrt_overall <- ci_norm(rc_slp$overall_acrt, rc_slp$overall_acrt_se)
acrt_overall_boot <- list(est = e0[[2]], se = se_b[[2]], lo = q_acrt[1], hi = q_acrt[2])
overall_boot <- list(est = e0[[1]], se = se_b[[1]], lo = q_att[1], hi = q_att[2])

## robustness of the overall ATT / ACRT and of the curve to the spline specification
specs <- data.frame(degree = c(1, 2, 3, 3), num_knots = c(0, 0, 0, 1))
robust <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  r <- run_cont("level", specs$degree[i], specs$num_knots[i], biters = 999L)
  data.frame(degree = specs$degree[i], num_knots = specs$num_knots[i], att = r$overall_att, att_se = r$overall_att_se,
             acrt = r$overall_acrt, acrt_se = r$overall_acrt_se,
             curve_rmse = sqrt(mean((r$att.d - tru$truth)^2)), curve_max_abs_err = max(abs(r$att.d - tru$truth)))
}))

cont <- list(curve = curve, acrt_curve = acrt_curve, overall = overall, acrt_overall = acrt_overall,
             acrt_overall_boot = acrt_overall_boot, overall_boot = overall_boot,
             fit = list(att = fit_att, acrt_vs_slope = fit_acrt_slope, acrt_vs_true_acrt = fit_acrt_true),
             robust = robust,
             check = list(replication_error = rep_err, se_ratio_curve_att = se_ratio_curve_att,
                          se_ratio_curve_acrt = se_ratio_curve_acrt,
                          se_ratio_overall_acrt = se_b[[2]] / rc_slp$overall_acrt_se,
                          se_ratio_overall_att = se_b[[1]] / rc_lev$overall_att_se,
                          n_over_n_treated = n_all / n_treated, B_city = B_CITY),
             spec = list(degree = 3, num_knots = 0, control_group = "nevertreated", biters = BITERS_CONT,
                         seed = SEED_CONT, n_dvals = K, n = n_all, n_treated = n_treated, n_untreated = n_all - n_treated))

## ================================================================
## (2) binarised DiD and TWFE dose regression; (3) TWFE weights
## ================================================================
cat("== (2) fixest, (3) TwoWayFEWeights ==\n")
B[, treat_post := as.integer(post == 1 & dose > 0)]
m_bin  <- feols(y ~ treat_post | unit + t, data = B, cluster = ~unit)
m_twfe <- feols(y ~ d | unit + t, data = B, cluster = ~unit)
m_orac <- feols(tau_true ~ d | unit + t, data = B, cluster = ~unit)      # TWFE on the noise-free effect column
stopifnot(abs(coef(m_bin)[["treat_post"]] - rc_lev$overall_att) < 1e-9)  # collapsing the dose to 0/1 reproduces the overall ATT
ci_t <- function(m, v) { ci <- confint(m)[v, ]; list(est = coef(m)[[v]], se = se(m)[[v]], lo = ci[[1]], hi = ci[[2]],
                                                   p = pvalue(m)[[v]], df = degrees_freedom(m, "t")) }
binary <- ci_t(m_bin, "treat_post")
tw_est <- ci_t(m_twfe, "d")

tw <- with_log(twowayfeweights(df_B[, c("unit", "t", "y", "d")], Y = "y", G = "unit", T = "t", D = "d",
                               type = "feTR", summary_measures = TRUE))
wcell <- merge(as.data.table(tw$dat_result)[, .(unit = G, t = T, weight)], B[, .(unit, t, d, dose, tau_true)], by = c("unit", "t"))
wt <- wcell[d > 0]
stopifnot(nrow(wt) == tw$tot_cells, abs(sum(wt$weight) - 1) < 1e-8, sum(wt$weight < 0) == tw$nr_minus,
          all(wcell$weight[wcell$d == 0] == 0))
implied <- wt[, sum(weight * tau_true / d)]                              # sum of weight x per-unit-dose effect
stopifnot(abs(implied - coef(m_orac)[["d"]]) < 1e-8, abs(tw$beta - tw_est$est) < 1e-8)   # the weights reproduce the oracle TWFE
weights <- list(
  n_cells = tw$tot_cells, n_positive = tw$nr_plus, n_negative = tw$nr_minus,
  share_negative = tw$nr_minus / tw$nr_weights,                          # by count
  sum_positive = tw$sum_plus, sum_negative = tw$sum_minus,
  min_sd_beta_delta0 = tw$sensibility, min_sd_opposite_sign = tw$sensibility2,
  implied_beta = implied,
  contrib_positive = wt[weight > 0, sum(weight * tau_true / d)], contrib_negative = wt[weight < 0, sum(weight * tau_true / d)],
  max_negative_weight_dose = wt[weight < 0, max(d)], min_positive_weight_dose = wt[weight > 0, min(d)])
twfe <- c(tw_est, list(oracle_est = coef(m_orac)[["d"]], weights = weights))

## ================================================================
## (4) DIDHAD on tessera_B_had
## ================================================================
cat("== (4) DIDHAD ==\n")
dfH <- as.data.frame(H)
hd <- timer("did_had", with_log(did_had(dfH, outcome = "y", group = "unit", time = "t", treatment = "d",
                                        effects = 1, placebo = 0, yatchew = TRUE, graph_off = TRUE)))
rm_ <- hd$results$resmat["Effect_1", ]; yt <- hd$results$yatchew_test["Effect_1", ]
stopifnot(length(yt) == 5, grepl("lin", names(yt)[1]), grepl("diff", names(yt)[2]), names(yt)[3] == "T_hr", names(yt)[5] == "N")
mH <- feols(y ~ d | unit + t, data = H, cluster = ~unit)
had <- list(est = rm_[["Estimate"]], se = rm_[["SE"]], lo = rm_[["LB.CI"]], hi = rm_[["UB.CI"]], n = rm_[["N"]],
            bw = rm_[["BW"]], n_bw = rm_[["N.BW"]], qug_stat = rm_[["T"]], qug_p = rm_[["p.val"]],
            yatchew = list(sigma2_lin = unname(yt[1]), sigma2_diff = unname(yt[2]), T_hr = unname(yt[3]),    # columns by position:
                           p = unname(yt[4]), n = unname(yt[5])),                                            # sigma2_lin, sigma2_diff, T_hr, p-value, N
            twfe = ci_t(mH, "d"), min_dose = H[post == 1, min(dose)], share_below_025 = H[post == 1, mean(dose < 0.25)])

## ================================================================
## (5) DIDmultiplegtSTAT on tessera_B
## ================================================================
cat("== (5) DIDmultiplegtSTAT ==\n")
st <- timer("did_multiplegt_stat", with_log(did_multiplegt_stat(df_B[, c("unit", "t", "y", "d")], Y = "y", ID = "unit", Time = "t", D = "d",
                                                              estimator = c("aoss", "waoss"), aoss_vs_waoss = TRUE)))
tab <- st$results$table
pick <- function(row) { r <- tab[row, ]; list(est = r[["Estimate"]], se = r[["SE"]], lo = r[["LB CI"]], hi = r[["UB CI"]],
                                             switchers = r[["Switchers"]], stayers = r[["Stayers"]]) }
## The SE column of the "Diff." row is sqrt(N) times the standard error implied by its own CI and
## t statistic (2.858 = 0.1278 * sqrt(500)), unlike the AOSS and WAOSS rows. We store the implied se,
## est / t, and check it against the CI half-width.
dv_ <- st$results$aoss_vs_waoss["Diff.", ]
se_diff <- abs(dv_[["Estimate"]] / dv_[["t stat."]])
stopifnot(abs(se_diff - (dv_[["UB CI"]] - dv_[["LB CI"]]) / (2 * qnorm(0.975))) < 1e-4,
          abs(dv_[["SE"]] / sqrt(nrow(wide)) - se_diff) < 1e-4)
stat <- list(aoss = pick("AOSS"), waoss = pick("WAOSS"),
             diff = list(est = dv_[["Estimate"]], se = se_diff, lo = dv_[["LB CI"]], hi = dv_[["UB CI"]],
                         t = dv_[["t stat."]], p = dv_[["pval."]], se_column_as_printed = dv_[["SE"]]))

res_dose <- list(cont = cont, binary = binary, twfe = twfe, had = had, stat = stat,
                 versions = c(contdid = as.character(packageVersion("contdid")), DIDHAD = as.character(packageVersion("DIDHAD")),
                              DIDmultiplegtSTAT = as.character(packageVersion("DIDmultiplegtSTAT")),
                              TwoWayFEWeights = as.character(packageVersion("TwoWayFEWeights")), fixest = as.character(packageVersion("fixest"))),
                 warnings = warn_log, messages = msg_log)
saveRDS(res_dose, file.path(code_dir, "res_dose.rds"))

## ---- report -------------------------------------------------------------------------------
f3 <- function(x) sprintf("%.3f", x)
ci_s <- function(r) sprintf("%s [%s, %s] (se %s)", f3(r$est), f3(r$lo), f3(r$hi), f3(r$se))
cat(sprintf("\n== Block B truth: ATT overall %.3f (pop %.3f); ACRT^o %.3f; PT-only slope %.3f; AS %.3f; WAS %.3f; HAD WAS %.3f ==\n",
            tr$B$att_overall, tr$B$att_overall_pop, tr$B$acrt_overall, tr$B$acrt_overall_pt_pop, tr$B$as, tr$B$was, tr$B$had_was))
cat("contdid overall ATT   ", ci_s(overall), "  covers:", overall$lo <= tr$B$att_overall && tr$B$att_overall <= overall$hi, "\n")
cat("contdid overall ATT (city bootstrap of replication) ", ci_s(overall_boot), "\n")
cat("contdid overall 'ACRT'", ci_s(acrt_overall), "\n")
cat("contdid overall 'ACRT' (city bootstrap)             ", ci_s(acrt_overall_boot), "\n")
cat(sprintf("  ATT(d|d) curve vs truth: rmse %.3f, mean signed err %.3f, max abs err %.3f, pointwise coverage %.2f\n",
            fit_att$rmse, fit_att$mean_signed_err, fit_att$max_abs_err, fit_att$share_covered))
cat(sprintf("  'ACRT' curve vs slope of ATT(d|d): rmse %.3f (coverage %.2f); vs causal ACRT(d|d): rmse %.3f (coverage %.2f)\n",
            fit_acrt_slope$rmse, fit_acrt_slope$share_covered, fit_acrt_true$rmse, fit_acrt_true$share_covered))
cat(sprintf("  se check: city bootstrap / package: overall ATT %.2f, overall ACRT %.2f, curve ATT %.2f, curve ACRT %.2f (n/n1 = %.2f)\n",
            cont$check$se_ratio_overall_att, cont$check$se_ratio_overall_acrt, se_ratio_curve_att, se_ratio_curve_acrt, n_all / n_treated))
cat("robustness to the spline specification:\n"); print(robust, digits = 4, row.names = FALSE)
cat("binarised DiD (dose > 0)", ci_s(binary), "\n")
cat(sprintf("TWFE on dose           %s ; noise-free TWFE %.3f ; p %.3g\n", ci_s(tw_est), twfe$oracle_est, tw_est$p))
cat(sprintf("  weights: %d positive, %d negative of %d cells (share %.3f, sum %.4f / %.4f); implied beta %.4f\n",
            weights$n_positive, weights$n_negative, weights$n_cells, weights$share_negative,
            weights$sum_positive, weights$sum_negative, weights$implied_beta))
cat("HAD (tessera_B_had)    ", ci_s(had), "  covers:", had$lo <= tr$B$had_was && tr$B$had_was <= had$hi, "\n")
cat(sprintf("  QUG T = %.3f (p %.3f); Yatchew T_hr %.3f (p %.3f); TWFE on tessera_B_had %s\n",
            had$qug_stat, had$qug_p, had$yatchew$T_hr, had$yatchew$p, ci_s(had$twfe)))
cat("AOSS                   ", ci_s(stat$aoss), "  covers:", stat$aoss$lo <= tr$B$as && tr$B$as <= stat$aoss$hi, "\n")
cat("WAOSS                  ", ci_s(stat$waoss), "  covers:", stat$waoss$lo <= tr$B$was && tr$B$was <= stat$waoss$hi, "\n")
cat(sprintf("AOSS - WAOSS test: %.3f (se %.3f), p %.4f\n", stat$diff$est, stat$diff$se, stat$diff$p))
cat("\nPackage warnings:\n"); print(table(warn_log)); cat("Package messages:\n"); print(table(msg_log))
cat("Saved res_dose.rds\n")
