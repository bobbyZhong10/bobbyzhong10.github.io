# =====================================================================
# 22_synthetic_control_factor_models / 04_hull.R
# Block S2: the fee reform in Aster-Capital, the metro with the highest
# intrinsic level (above every donor's), which puts it outside the donor
# convex hull. Same layout as S1: 20 donors, T_0 = 24 pre-months, post = months 25..36.
#
#   SC        sc_weights() on the 24 pre-period outcomes, V = I (the
#             recipe 03_s1_modern.R uses for S1's pre-RMSPE)
#   ASC       augsynth ridge, scm = TRUE (Ben-Michael, Feller, Rothstein
#             2021): SC weights plus a ridge correction that may go
#             negative, which is how ASC extrapolates outside the hull
#   SDID      synthdid (unit intercept absorbs the level gap), placebo se
#   DiD       synthdid::did_estimate (= the 2x2 DiD), for reference
#   penalized penalized_sc() (Abadie and L'Hour 2021, objective
#   SC        ||X1 - X0 w||^2 + pen * sum_j w_j ||X1 - X0_j||^2 on the
#             simplex, V = I). The penalty is chosen by Abadie and
#             L'Hour's Section 4.1 rule: leave-one-out cross-validation on
#             the post-period outcomes of the untreated units. Each donor
#             in turn plays the treated unit, its penalized SC is built
#             from the other 19 donors on the 24 pre-months, and the loss
#             is the sum over donors and the 12 post months of squared
#             prediction errors; pen minimises it over the grid
#             10^seq(-4, 1, by = 0.25) (21 values; pen = 0 is excluded
#             because the unpenalised problem need not have a unique
#             solution). Their Section 4.2 rule (holdout on the treated
#             unit's own pre-period, train months 1-16, validate 17-24) is
#             computed as a sensitivity check and stored, not used.
#             Uniqueness: penalized_sc has no starting point (the QP solver
#             is deterministic), so the weights are recomputed under 20
#             random orderings of the donors (est_seeds22) and the largest
#             absolute weight difference is stored.
# Writes res_hull.rds. Run from this directory:  Rscript 04_hull.R
# =====================================================================

suppressWarnings(suppressMessages({
  library(synthdid)
  library(augsynth)
}))
source("sc_fns.R")
source("est_seeds.R")

S2 <- readRDS("aster_S2.rds"); tr <- readRDS("truth22.rds")
Y  <- S2$Y
n_pre <- S2$meta$n_pre                       # 24
Tt <- ncol(Y); post <- (n_pre + 1):Tt
y1 <- Y[1, ]; Y0 <- t(Y[-1, ])
donors <- S2$unit_names[-1]
z <- qnorm(0.975)
gap_stats <- function(w) {
  g <- y1 - drop(Y0 %*% w)
  list(est = mean(g[post]), pre_rmspe = sqrt(mean(g[1:n_pre]^2)), gap = g)
}

# ---- SC on outcomes, V = I ------------------------------------------
w_sc <- sc_weights(y1[1:n_pre], Y0[1:n_pre, ], diag(n_pre))
gs <- gap_stats(w_sc)
sc <- list(est = gs$est, pre_rmspe = gs$pre_rmspe, w = w_sc, n_nonzero = sum(w_sc > 1e-6), gap = gs$gap,
           synth = drop(Y0 %*% w_sc))

# ---- ASC (augsynth ridge) -------------------------------------------
pan <- subset(S2$panel, select = -unit)      # augsynth breaks on an integer column named unit
asc_fit <- suppressMessages(augsynth(log_gmv ~ treat_it, metro, month, pan, progfunc = "ridge", scm = TRUE))
set.seed(est_seeds22[["s2_asc_summary"]])
asc_sum <- suppressMessages(summary(asc_fit))
asc_w <- setNames(drop(asc_fit$weights), rownames(asc_fit$weights))[donors]
att_path <- asc_sum$att
asc_gap <- att_path$Estimate[match(1:Tt, att_path$Time)]
asc <- list(est = asc_sum$average_att$Estimate, p_conformal = asc_sum$average_att$p_val, w = asc_w,
            n_neg = sum(asc_w < 0), min_w = min(asc_w), max_w = max(asc_w), sum_w = sum(asc_w),
            ridge_lambda = asc_fit$lambda, pre_rmspe = sqrt(mean(asc_gap[1:n_pre]^2)), gap = asc_gap,
            synth = y1 - asc_gap)

# ---- SDID and DiD (synthdid) ----------------------------------------
Ysd <- rbind(Y[-1, ], Y[1, , drop = FALSE]); N0 <- nrow(Ysd) - 1
sdid_fit <- synthdid_estimate(Ysd, N0, n_pre)
set.seed(est_seeds22[["s2_sdid_placebo"]])
se_sdid <- sqrt(as.numeric(vcov(sdid_fit, method = "placebo", replications = 500)))
w_sdid <- attr(sdid_fit, "weights")
# SDID's implied path: treated minus (omega-weighted donors + intercept), the intercept being
# the lambda-weighted pre-period level difference (synthdid_curves logic, written out)
om <- w_sdid$omega; la <- w_sdid$lambda
synth_sdid_raw <- drop(t(Y[-1, ]) %*% om)
shift <- sum(la * y1[1:n_pre]) - sum(la * synth_sdid_raw[1:n_pre])
sdid <- list(est = as.numeric(sdid_fit), se = se_sdid, ci = as.numeric(sdid_fit) + c(-1, 1) * z * se_sdid,
             omega = setNames(om, donors), lambda = la, n_omega_nonzero = sum(om > 1e-6),
             n_lambda_nonzero = sum(la > 1e-6), level_shift = shift, synth = synth_sdid_raw + shift)
did <- list(est = as.numeric(synthdid::did_estimate(Ysd, N0, n_pre)))

# ---- penalized SC: penalty by leave-one-out CV (Abadie and L'Hour 2021, Section 4.1) ----
pen_grid <- 10^seq(-4, 1, by = 0.25)
X0_all <- Y[-1, ]                                       # 20 x 36, donors in rows
loo_loss <- vapply(pen_grid, function(p) {
  sum(vapply(seq_len(nrow(X0_all)), function(i) {
    others <- X0_all[-i, , drop = FALSE]
    wi <- penalized_sc(X0_all[i, 1:n_pre], t(others[, 1:n_pre]), p)$w
    sum((X0_all[i, post] - drop(t(others[, post]) %*% wi))^2)
  }, numeric(1)))
}, numeric(1))
pen_loo <- pen_grid[which.min(loo_loss)]
# sensitivity: Section 4.2, holdout on the treated unit's pre-period (train 1-16, validate 17-24)
train <- 1:16; valid <- 17:n_pre
ho_loss <- vapply(pen_grid, function(p) {
  wi <- penalized_sc(y1[train], Y0[train, ], p)$w
  sum((y1[valid] - drop(Y0[valid, ] %*% wi))^2)
}, numeric(1))
pen_holdout <- pen_grid[which.min(ho_loss)]

fit_pen <- function(p) {
  ps <- penalized_sc(y1[1:n_pre], Y0[1:n_pre, ], p)
  g <- gap_stats(ps$w)
  list(pen = p, w = ps$w, n_nonzero = ps$n_nonzero, est = g$est, pre_rmspe = g$pre_rmspe, gap = g$gap,
       synth = drop(Y0 %*% ps$w))
}
pen <- fit_pen(pen_loo)
# uniqueness across donor orderings
set.seed(est_seeds22[["s2_pen_orderings"]])
n_orderings <- 20L
diffs <- vapply(seq_len(n_orderings), function(k) {
  o <- sample.int(ncol(Y0))
  wk <- penalized_sc(y1[1:n_pre], Y0[1:n_pre, o], pen$pen)$w
  max(abs(wk[donors] - pen$w[donors]))
}, numeric(1))
pen$max_diff_across_starts <- max(diffs)
pen$n_orderings <- n_orderings
pen$rule <- paste("Abadie and L'Hour (2021) Section 4.1: leave-one-out CV on the 12 post-period outcomes of the",
                  "20 untreated donors (sum of squared individual prediction errors), grid 10^seq(-4, 1, 0.25)")
pen$grid <- pen_grid; pen$loo_loss <- loo_loss
pen_holdout_fit <- fit_pen(pen_holdout)
pen$holdout <- list(pen = pen_holdout, est = pen_holdout_fit$est, n_nonzero = pen_holdout_fit$n_nonzero,
                    pre_rmspe = pen_holdout_fit$pre_rmspe, w = pen_holdout_fit$w, loss = ho_loss)
# the penalty path the chapter can describe (donors used at each pen on the grid)
pen$path <- data.frame(pen = pen_grid, n_nonzero = vapply(pen_grid, function(p) fit_pen(p)$n_nonzero, numeric(1)),
                       est = vapply(pen_grid, function(p) fit_pen(p)$est, numeric(1)))

res <- list(sc = sc, asc = asc, sdid = sdid, did = did, pen = pen, truth = tr$S2,
            level_gap = S2$meta$level_gap, hull_l1_gap = S2$meta$hull_l1_gap,
            treated_outside_hull = S2$meta$treated_outside_hull, y1 = y1, n_pre = n_pre)
saveRDS(res, "res_hull.rds")

f4 <- function(v) sprintf("%.4f", v)
err <- function(e) abs(e - tr$S2$att_post)
cat("================ S2: treated metro outside the donor hull (truth", f4(tr$S2$att_post), ") ================\n")
cat("SC (V = I):      ", f4(sc$est), " error", f4(err(sc$est)), " pre-RMSPE", f4(sc$pre_rmspe), " donors", sc$n_nonzero, "\n")
cat("ASC (ridge):     ", f4(asc$est), " error", f4(err(asc$est)), " negative weights", asc$n_neg, " min", f4(asc$min_w), "\n")
cat("SDID:            ", f4(sdid$est), " error", f4(err(sdid$est)), " 95% CI [", f4(sdid$ci[1]), ",", f4(sdid$ci[2]), "]\n")
cat("DiD (uniform):   ", f4(did$est), "\n")
cat("penalized SC:     pen =", format(pen$pen, digits = 4), "(LOO CV)", f4(pen$est), " donors", pen$n_nonzero,
    " pre-RMSPE", f4(pen$pre_rmspe), " max diff over", n_orderings, "orderings", format(pen$max_diff_across_starts, digits = 3), "\n")
cat("  sensitivity: holdout rule pen =", format(pen_holdout, digits = 4), "->", f4(pen$holdout$est), "with",
    pen$holdout$n_nonzero, "donors\n")
