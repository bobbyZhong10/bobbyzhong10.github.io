# =====================================================================
# 22_synthetic_control_factor_models / 05c_stag_scpi.R
# Block S3 (aster_S3.rds): scpi prediction intervals for the staggered rollout
# (Phase 2 carry-over, Ruling S12). One run on the panel 05_staggered.R uses:
# the 15 adopters in three waves of 5 (months 29, 33, 37) and, as donors, the
# 65 never-treated metros (scdataMulti's pool is the never-treated units, the
# 20 contract metros among them; SDID by wave uses the same pool).
#   scdataMulti(effect = "time")  the average effect on the treated by event
#                                 month over the horizons every adopter reaches:
#                                 event months 0..11 (scpi's horizons 1..12),
#                                 each adopter weighted equally, i.e. the same
#                                 180 treated metro-months that multisynth's
#                                 Average row uses (s3_multisynth_n_cells)
#   scpi()                        simplex weights fitted for each adopter on its
#                                 own pre-period (months 1..g-1: 28, 32 or 36
#                                 months), V = "separate", no constant, sims =
#                                 200, Gaussian out-of-sample bounds
#                                 (CI.all.gaussian), u.alpha = e.alpha = 0.05
#                                 (coverage at least 90%), cores = 1, seed
#                                 est_seeds22[["s3_scpi"]]
# The effect interval at event month e is the adopters' mean outcome minus the
# interval for the synthetic (untreated) mean, as in 03_s1_modern.R. scpi gives
# one interval per event month and none for the 12-month average, so the summary
# is the mean of the 12 effects (est) and how many of the 12 intervals cover the
# truth 0.06 (n_covered).
# It runs in its own R process with its own seed and writes only its own file,
# so no random stream of the other scripts moves. About 2 minutes on this machine.
# Writes res_stag_scpi.rds. Run from this directory:  Rscript 05c_stag_scpi.R
# =====================================================================

suppressWarnings(suppressMessages(library(scpi)))
source("est_seeds.R")

S3 <- readRDS("aster_S3.rds"); tr <- readRDS("truth22.rds")
p <- S3$panel; m <- S3$meta; truth <- tr$S3$att
df <- data.frame(metro = as.character(p$metro), month = p$month, log_gmv = p$log_gmv, treat_it = p$treat_it)

t0 <- proc.time()[["elapsed"]]
sdat <- scdataMulti(df, id.var = "metro", time.var = "month", outcome.var = "log_gmv", treatment.var = "treat_it",
                    effect = "time", constant = FALSE, cointegrated.data = FALSE, verbose = FALSE)
sp <- sdat$specs
adopters <- S3$unit_names[is.finite(m$cohort)]; never <- S3$unit_names[is.infinite(m$cohort)]
# the sample must be 05_staggered.R's: every adopter treated, every never-treated metro a donor
stopifnot(setequal(sp$treated.units, adopters), length(adopters) == 15L,
          all(vapply(sp$donors.list, function(d) setequal(d, never), logical(1))))
n_h <- nrow(sdat$P)                                   # horizons every adopter reaches

set.seed(est_seeds22[["s3_scpi"]])
sres <- scpi(data = sdat, w.constr = list(name = "simplex"), sims = 200, e.method = "gaussian",
             cores = 1, verbose = FALSE)
secs <- proc.time()[["elapsed"]] - t0

# the adopters' mean outcome at event month e = 0..n_h-1 (month cohort + e), equal weights
coh <- m$cohort; ad <- which(is.finite(coh))
y1 <- vapply(seq_len(n_h) - 1L, function(e) mean(S3$Y[cbind(ad, coh[ad] + e)]), numeric(1))
y0 <- drop(sres$est.results$Y.post.fit)
ci <- sres$inference.results$CI.all.gaussian; cin <- sres$inference.results$CI.in.sample
stopifnot(length(y0) == n_h, nrow(ci) == n_h, nrow(cin) == n_h)
by_e <- data.frame(e = seq_len(n_h) - 1L, est = y1 - y0, lo = unname(y1 - ci[, 2]), hi = unname(y1 - ci[, 1]),
                   lo_in = unname(y1 - cin[, 2]), hi_in = unname(y1 - cin[, 1]))
by_e$covered <- by_e$lo <= truth & truth <= by_e$hi
by_e$covered_in <- by_e$lo_in <= truth & truth <= by_e$hi_in

# each adopter's pre-period fit (RMSPE of its own SC over months 1..g-1)
A <- drop(sdat$A); Afit <- drop(sres$est.results$Y.pre.fit)
unit_of <- sub("\\..*$", "", rownames(sdat$A))
pre_rmspe <- tapply(A - Afit, unit_of, function(r) sqrt(mean(r^2)))
# Theorem 7's q^pool: the adopters' mean pre-period gap over each one's last 12 pre months
# (event months -12..-1); a gap the separate fits leave on the average path carries into the post period
t_of <- suppressWarnings(as.numeric(sub("^.*\\.", "", rownames(sdat$A))))
stopifnot(!anyNA(t_of))
gap_u <- lapply(split(seq_along(A), unit_of), function(ix) { ix <- ix[order(t_of[ix])]; tail(A[ix] - Afit[ix], 12) })
stopifnot(all(lengths(gap_u) == 12L))
pool_gap <- rowMeans(do.call(cbind, gap_u))
w <- sres$est.results$w; w_unit <- sub("\\..*$", "", names(w))
n_w <- tapply(w, w_unit, function(x) sum(x > 0.01))

res <- list(est = mean(by_e$est), by_event = by_e, n_months = n_h, n_cells = length(ad) * n_h,
            n_covered = sum(by_e$covered), n_covered_in = sum(by_e$covered_in),
            n_below_truth = sum(by_e$hi < truth), n_above_truth = sum(by_e$lo > truth),
            lo_min = min(by_e$lo), hi_max = max(by_e$hi), mean_length = mean(by_e$hi - by_e$lo),
            n_treated = length(ad), n_donors = length(never), pre_rmspe = pre_rmspe, n_w_above_0p01 = n_w,
            pool_gap = pool_gap, pre_gap_mean = mean(pool_gap), q_pool12 = sqrt(mean(pool_gap^2)),
            w = w, truth = truth,
            settings = list(effect = "time", w_constr = "simplex", V = "separate", constant = FALSE,
                            cointegrated_data = FALSE, e_method = "gaussian", sims = 200L, cores = 1L,
                            seed = est_seeds22[["s3_scpi"]], u_alpha = sres$inference.results$u.alpha,
                            e_alpha = sres$inference.results$e.alpha,
                            failed_sims = sum(sres$inference.results$failed.sims)),
            versions = c(scpi = as.character(utils::packageVersion("scpi"))))
saveRDS(res, "res_stag_scpi.rds")

f4 <- function(v) sprintf("%.4f", v)
cat(sprintf("================ S3: scpi (scdataMulti, effect = \"time\"), truth %s, %.0f s ================\n", f4(truth), secs))
cat(sprintf("  %d adopters, %d never-treated donors, %d event months (%d treated metro-months)\n",
            res$n_treated, res$n_donors, n_h, res$n_cells))
for (i in seq_len(n_h))
  cat(sprintf("  e = %2d  %s  [%s, %s]  covers truth: %s\n", by_e$e[i], f4(by_e$est[i]), f4(by_e$lo[i]), f4(by_e$hi[i]),
              by_e$covered[i]))
cat("  mean of the", n_h, "effects", f4(res$est), "| intervals covering the truth:", res$n_covered, "of", n_h,
    "(in-sample only:", res$n_covered_in, ") | bounds", f4(res$lo_min), "to", f4(res$hi_max),
    "| mean length", f4(res$mean_length), "\n")
cat("  adopters' mean pre-period gap, event months -12..-1: mean", f4(res$pre_gap_mean), "rms", f4(res$q_pool12), "\n")
cat("  adopters' pre-period RMSPE: mean", f4(mean(pre_rmspe)), "max", f4(max(pre_rmspe)),
    "| failed simulations", res$settings$failed_sims, "\n")
