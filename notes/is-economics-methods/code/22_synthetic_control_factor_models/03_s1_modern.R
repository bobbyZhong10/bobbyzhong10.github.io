# =====================================================================
# 22_synthetic_control_factor_models / 03_s1_modern.R
# Block S1 (single-metro fee reform, treated metro inside the donor hull)
# with the modern estimators and inference. T_0 = 24 pre-months (the
# chapter's T_0), post = months 25..36, 20 donors.
#
#   SC on outcomes   sc_weights() (sc_fns.R), V = I on the 24 pre-period
#                    outcomes; the pre-RMSPE recipe 04_hull.R uses on S2
#   SDID             synthdid::synthdid_estimate (AAHIW 2021); interval
#                    est -/+ 1.96 x placebo se (synthdid's jackknife and
#                    bootstrap return NA with one treated unit)
#   DiD, two ways    synthdid::did_estimate (SDID with uniform unit and
#                    time weights) vs the legacy 2x2 DiD, plus TWFE: one
#                    number three ways (the ledger 19->22 identity)
#   SC (synthdid)    synthdid::sc_estimate (no unit intercept, no time
#                    weights), with placebo se
#   ASC              augsynth ridge (Ben-Michael, Feller, Rothstein 2021)
#   conformal        conformal_p() (Chernozhukov, Wuthrich, Zhu 2021,
#                    arXiv 1712.09089v10): test of H0: no effect in months
#                    25-36; SC fitted on all 36 months under H0 (sc_fns.R
#                    header), statistic S_1. Reported scheme: iid
#                    permutations (Ruling S8). Rationale, from the paper and
#                    the DGP rather than the test: CWZ Section 2 says the
#                    choice of permutation set does not matter for exact
#                    validity when the residuals are exchangeable, and that
#                    the set of all permutations (Pi_all, "iid permutations")
#                    gives more precise p-values than the T moving blocks;
#                    their Section 3.1 uses Pi_all under Assumption 2.1 (iid
#                    shocks) and moving blocks only under the weaker
#                    stationary, mixing Assumption 2.2. S1 meets 2.1 exactly:
#                    the treated metro's factor structure is the w_true
#                    combination of the anchor donors (01_dgp.R), so its
#                    untreated outcome minus that combination is
#                    u_t = e_1t - sum_j w_true_j e_jt with e iid N(0, 0.02^2),
#                    an iid sequence. CWZ report both schemes in their
#                    application and sample 10,000 iid permutations; here
#                    9999 random ones plus the identity. augsynth's own
#                    conformal inference (the ASC p below) also defaults to
#                    type = "iid". The moving-block p-value is computed and
#                    stored next to it, and the two differ: with T = 36 the
#                    moving-block grid is 1/36, and the fit under H0 (which
#                    must also absorb the post-period effect) leaves serially
#                    correlated residuals, so neighbouring shifts give
#                    statistics close to the observed one. At the 10% level
#                    the iid test rejects and the moving-block test does not;
#                    tests/test_conclusions.R locks both results, and the
#                    chapter prints both with the assumption each needs.
#   scpi             prediction intervals for each post month (Cattaneo,
#                    Feng, Titiunik 2021), simplex weights, sims = 200,
#                    Gaussian out-of-sample bounds (CI.all.gaussian);
#                    the effect interval is y1_t minus the interval for
#                    the synthetic outcome
#   placebo ratio    mspe_ratio_perm() with sc_fit (hand-written check of
#                    the tidysynth permutation test of 02_core.R)
# Writes res_s1_ext.rds. Run from this directory:  Rscript 03_s1_modern.R
# =====================================================================

suppressWarnings(suppressMessages({
  library(synthdid)
  library(augsynth)
  library(scpi)
  library(fixest)
}))
source("sc_fns.R")
source("est_seeds.R")

S1 <- readRDS("aster_S1.rds"); tr <- readRDS("truth22.rds")
Y  <- S1$Y
n_pre <- S1$meta$T0 - 1                    # 24, the chapter's T_0
Tt <- ncol(Y); post <- (n_pre + 1):Tt
y1 <- Y[1, ]; Y0 <- t(Y[-1, ])            # Tt x 20, donors in columns
colnames(Y0) <- S1$unit_names[-1]          # aster_S1's Y has no dimnames: name the donors
z  <- qnorm(0.975)

# ---- SC on outcomes, V = I (the comparable pre-RMSPE recipe) ---------
w_sc <- sc_weights(y1[1:n_pre], Y0[1:n_pre, ], diag(n_pre))
gap_sc <- y1 - drop(Y0 %*% w_sc)
sc_hand <- list(est = mean(gap_sc[post]), pre_rmspe = sqrt(mean(gap_sc[1:n_pre]^2)),
                gap_last = gap_sc[Tt], w = w_sc, n_nonzero = sum(w_sc > 1e-6))

# ---- legacy naive 2x2 DiD, TWFE and synthdid's DiD -------------------
pre1 <- 1:n_pre
naive_did <- (mean(y1[post]) - mean(y1[pre1])) - (mean(colMeans(Y[-1, post])) - mean(colMeans(Y[-1, pre1])))
stopifnot(abs(naive_did - S1$meta$did_est) < 1e-12)      # two codings of one formula
twfe <- unname(coef(feols(log_gmv ~ treat_it | unit + month, S1$panel)))

Ysd <- rbind(Y[-1, ], Y[1, , drop = FALSE])   # synthdid wants controls first, treated last
N0 <- nrow(Ysd) - 1
did_sd  <- synthdid::did_estimate(Ysd, N0, n_pre)
did_uniform <- as.numeric(did_sd)

# ---- SDID and synthdid's SC, placebo variance ------------------------
sdid_fit <- synthdid_estimate(Ysd, N0, n_pre)
sc_sd    <- sc_estimate(Ysd, N0, n_pre)
set.seed(est_seeds22[["s1_sdid_placebo"]])
se_sdid  <- sqrt(as.numeric(vcov(sdid_fit, method = "placebo", replications = 500)))
set.seed(est_seeds22[["s1_sc_placebo"]])
se_scsd  <- sqrt(as.numeric(vcov(sc_sd, method = "placebo", replications = 500)))
set.seed(est_seeds22[["s1_did_placebo"]])
se_did   <- sqrt(as.numeric(vcov(did_sd, method = "placebo", replications = 500)))
w_sdid   <- attr(sdid_fit, "weights")
sdid <- list(est = as.numeric(sdid_fit), se = se_sdid,
             ci = as.numeric(sdid_fit) + c(-1, 1) * z * se_sdid,
             omega = setNames(w_sdid$omega, rownames(Y)[-1]), lambda = w_sdid$lambda,
             n_omega_nonzero = sum(w_sdid$omega > 1e-6), n_lambda_nonzero = sum(w_sdid$lambda > 1e-6),
             placebo_reps = 500L)
sc_synthdid <- list(est = as.numeric(sc_sd), se = se_scsd, ci = as.numeric(sc_sd) + c(-1, 1) * z * se_scsd)
did_placebo <- list(est = did_uniform, se = se_did, ci = did_uniform + c(-1, 1) * z * se_did)

# ---- ASC (augsynth, ridge) -----------------------------------------
# augsynth's internal `unit` clashes with the panel's integer column `unit`: drop it
pan <- subset(S1$panel, select = -unit)
asc_fit <- suppressMessages(augsynth(log_gmv ~ treat_it, metro, month, pan, progfunc = "ridge", scm = TRUE))
set.seed(est_seeds22[["s1_asc_summary"]])
asc_sum <- suppressMessages(summary(asc_fit))
asc_w <- setNames(drop(asc_fit$weights), rownames(asc_fit$weights))
asc <- list(est = asc_sum$average_att$Estimate, p_conformal = asc_sum$average_att$p_val,
            w = asc_w, ridge_lambda = asc_fit$lambda, n_neg = sum(asc_w < 0), min_w = min(asc_w))

# ---- conformal inference (CWZ 2021) ----------------------------------
w_h0 <- sc_weights(y1, Y0, diag(Tt))             # fitted on all 36 months under H0
y0hat_h0 <- drop(Y0 %*% w_h0)
set.seed(est_seeds22[["s1_conformal_iid"]])
p_iid  <- conformal_p(y1, y0hat_h0, T0 = n_pre, scheme = "iid", q = 1, n_perm = 9999)
p_mb   <- conformal_p(y1, y0hat_h0, T0 = n_pre, scheme = "moving_block", q = 1)
p_mb2  <- conformal_p(y1, y0hat_h0, T0 = n_pre, scheme = "moving_block", q = 2)
conformal <- list(p = p_iid, scheme = "iid", n_perm = 9999L, q = 1,
                  p_moving_block = p_mb, p_moving_block_q2 = p_mb2, n_shifts = Tt)

# ---- scpi prediction intervals ---------------------------------------
df <- S1$panel; df$metro <- as.character(df$metro)
sdat <- scdata(df = df, id.var = "metro", time.var = "month", outcome.var = "log_gmv",
               period.pre = 1:n_pre, period.post = post, unit.tr = "Aster-Metro",
               unit.co = S1$unit_names[-1], constant = FALSE, cointegrated.data = FALSE)
set.seed(est_seeds22[["s1_scpi"]])
sres <- scpi(data = sdat, w.constr = list(name = "simplex"), sims = 200, e.method = "gaussian",
             cores = 1, verbose = FALSE)
ci_y0 <- sres$inference.results$CI.all.gaussian    # interval for the synthetic (untreated) outcome
ci_in <- sres$inference.results$CI.in.sample
scpi_res <- list(lo = unname(y1[post] - ci_y0[, 2]), hi = unname(y1[post] - ci_y0[, 1]),
                 est = unname(y1[post] - drop(sres$est.results$Y.post.fit)),
                 lo_in = unname(y1[post] - ci_in[, 2]), hi_in = unname(y1[post] - ci_in[, 1]),
                 w = setNames(drop(sres$est.results$w), S1$unit_names[-1]),
                 sims = 200L, u_alpha = sres$inference.results$u.alpha, e_alpha = sres$inference.results$e.alpha)
scpi_res$covered <- scpi_res$lo <= tr$S1$tau_path & tr$S1$tau_path <= scpi_res$hi
scpi_res$covered_in <- scpi_res$lo_in <= tr$S1$tau_path & tr$S1$tau_path <= scpi_res$hi_in

# ---- hand-written placebo permutation (ADH convention) ---------------
mr <- mspe_ratio_perm(Y, treated = 1, T0 = n_pre, fit_fun = sc_fit)

res <- list(
  naive_did = naive_did, did_uniform = did_uniform, twfe = twfe,
  did_identity_gap = did_uniform - naive_did, twfe_identity_gap = twfe - naive_did,
  did_placebo = did_placebo,
  sdid = sdid, sc_synthdid = sc_synthdid, asc = asc,
  sc_hand = sc_hand, sc_pre_rmspe = sc_hand$pre_rmspe,
  conformal = conformal, scpi = scpi_res,
  perm_hand = list(ratio = mr$ratio, p = mr$p, rank = mr$rank, n_units = mr$n_units),
  truth = tr$S1
)
saveRDS(res, "res_s1_ext.rds")

f3 <- function(v) sprintf("%.3f", v); f4 <- function(v) sprintf("%.4f", v)
cat("================ S1: modern estimators (truth", f3(tr$S1$att_post), ") ================\n")
cat("DiD three ways: legacy 2x2", sprintf("%.10f", naive_did), "| synthdid did_estimate", sprintf("%.10f", did_uniform),
    "| TWFE", sprintf("%.10f", twfe), "\n")
cat("  gaps: did_estimate - 2x2 =", format(res$did_identity_gap, digits = 3), "; TWFE - 2x2 =",
    format(res$twfe_identity_gap, digits = 3), "\n")
cat("SC (sc_weights, V = I):", f4(sc_hand$est), " pre-RMSPE", f4(sc_hand$pre_rmspe), " donors", sc_hand$n_nonzero, "\n")
cat("SC (synthdid):", f4(sc_synthdid$est), " placebo se", f4(se_scsd), "\n")
cat("SDID:", f4(sdid$est), " placebo se", f4(se_sdid), " 95% CI [", f4(sdid$ci[1]), ",", f4(sdid$ci[2]), "]\n")
cat("ASC (ridge):", f4(asc$est), " conformal p", f3(asc$p_conformal), " negative weights", asc$n_neg, "\n")
cat("conformal p (S_1, SC under H0 on 36 months): iid", f4(p_iid), "(9999 perms) | moving block", f4(p_mb),
    "| moving block q = 2", f4(p_mb2), "\n")
cat("scpi (Gaussian) intervals cover the true effect in", sum(scpi_res$covered), "of 12 post months\n")
cat("hand placebo (sc_fit, ADH pools): ratio", f3(mr$ratio), " rank", mr$rank, "of", mr$n_units, " p", f3(mr$p), "\n")
