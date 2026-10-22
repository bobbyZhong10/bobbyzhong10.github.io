# ------------------------------------------------------------------
# stag_fns.R -- the block S3 estimator recipes, shared by 05_staggered.R
# (the stored draw, seed seeds22[["S3"]]) and 05b_stag_screen.R (the 24-seed screen).
#
# Each fit_*() takes an S3 panel (aster_S3.rds$panel or gen_S3(seed)$data$panel:
# columns unit, month, log_gmv, cohort (Inf = never treated), treat_it) and returns
# list(est, ci = c(lo, hi), ...) for the average effect on the treated cells (the
# effect is constant in S3, so every aggregation targets the same 0.06). Every
# bootstrap is seeded inside the function, so a call is deterministic and does not
# depend on the session's random stream or on the order of calls.
#
#   fit_twfe(p)                       fixest::feols(log_gmv ~ treat_it | unit + month),
#                                     clustered by metro; 95% confint()
#   fit_cs(p, control_group, biters, seed)
#                                     did::att_gt (est_method "reg", cohort Inf -> 0,
#                                     multiplier bootstrap, no uniform band) then
#                                     aggte(type = "simple"); est -/+ qnorm(.975) se
#   fit_fect(p, method, nboots, seed, ...)
#                                     fect::fect(log_gmv ~ treat_it, index = c("unit",
#                                     "month"), force = "two-way", se = TRUE, parallel =
#                                     FALSE); interval = est.avg CI.lower, CI.upper
#                                     (fect's bootstrap). "ife": CV over r in 0..4.
#                                     "mc": CV over lambda (fect defaults).
#   fit_mc_refit(p, lambda, nboots, seed)
#                                     fect MC at a given lambda with CV = FALSE, tol =
#                                     1e-7, max.iteration = 50000 (fect runs its CV, and
#                                     its final MC fit after CV, at tol >= 1e-3; Task 9
#                                     report 7.1, Ruling S6a)
#   fit_multisynth(p, seed)           augsynth::multisynth(log_gmv ~ treat_it, unit, month)
#                                     with defaults; summary() (bootstrap, seeded); row
#                                     Level "Average", Time NA
#   fit_sdid_cohort(p, B, seed, cores)
#                                     SDID per adoption cohort (AAHIW 2021, Section 8):
#                                     cohort g's units plus the never-treated units, all
#                                     48 months, pre = months 1..g-1; the cohort estimates
#                                     are averaged with weights equal to their shares of
#                                     treated metro-months. Interval: est -/+ 1.96 se, se
#                                     from a unit bootstrap of the whole procedure (AAHIW
#                                     Algorithm 2 applied to the aggregate: resample the 80
#                                     metros with replacement, refit every cohort present,
#                                     reweight; redraw when no treated or no never-treated
#                                     metro is drawn). The B index draws are made up front
#                                     from `seed`, so the result does not depend on `cores`.
# ------------------------------------------------------------------

.st_quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(expr))
  out
}

fit_twfe <- function(p) {
  m <- fixest::feols(log_gmv ~ treat_it | unit + month, data = p, cluster = ~unit)
  ci <- stats::confint(m, level = 0.95)
  list(est = unname(stats::coef(m)[["treat_it"]]), se = unname(fixest::se(m)[["treat_it"]]),
       ci = c(ci[1, 1], ci[1, 2]))
}

fit_cs <- function(p, control_group = "nevertreated", biters = 999, seed = 1) {
  p$gname <- ifelse(is.finite(p$cohort), p$cohort, 0)
  set.seed(seed)
  # att_gt warns that it cannot return the pre-test Wald statistic (singular covariance of
  # the pre-period ATT(g,t)); that pre-test is not used here
  gt <- withCallingHandlers(
    did::att_gt("log_gmv", "month", "unit", gname = "gname", data = p, control_group = control_group,
                est_method = "reg", bstrap = TRUE, cband = FALSE, biters = biters),
    warning = function(w) if (grepl("pre-test Wald statistic", conditionMessage(w))) invokeRestart("muffleWarning"))
  ag <- did::aggte(gt, type = "simple", bstrap = TRUE, biters = biters, cband = FALSE)
  est <- ag$overall.att; se <- ag$overall.se
  list(est = est, se = se, ci = est + c(-1, 1) * stats::qnorm(0.975) * se)
}

.fect_out <- function(f) {
  ea <- f$est.avg
  list(est = unname(ea[1, "ATT.avg"]), se = unname(ea[1, "S.E."]),
       ci = c(unname(ea[1, "CI.lower"]), unname(ea[1, "CI.upper"])))
}

fit_fect <- function(p, method = c("ife", "mc"), nboots = 100, seed = 1) {
  method <- match.arg(method)
  args <- list(log_gmv ~ treat_it, data = p, index = c("unit", "month"), method = method, force = "two-way",
               CV = TRUE, se = TRUE, nboots = nboots, parallel = FALSE, seed = seed)
  if (method == "ife") args$r <- c(0, 4)
  f <- .st_quiet(do.call(fect::fect, args))
  out <- .fect_out(f)
  if (method == "ife") out$r_cv <- unname(f$r.cv)
  if (method == "mc") { out$lambda_cv <- unname(f$lambda.cv); out$lambda_norm <- unname(f$lambda.cv / max(f$eigen.all)) }
  out$nboots <- nboots
  out
}

fit_mc_refit <- function(p, lambda, nboots = 100, seed = 1, tol = 1e-7, max_iter = 50000) {
  f <- .st_quiet(fect::fect(log_gmv ~ treat_it, data = p, index = c("unit", "month"), method = "mc",
                            force = "two-way", CV = FALSE, lambda = lambda, se = TRUE, nboots = nboots,
                            parallel = FALSE, seed = seed, tol = tol, max.iteration = max_iter))
  out <- .fect_out(f)
  out$lambda <- lambda; out$tol <- tol; out$max_iter <- max_iter; out$nboots <- nboots
  out$niter <- unname(f$niter)
  out
}

fit_multisynth <- function(p, seed = 1) {
  ms <- .st_quiet(augsynth::multisynth(log_gmv ~ treat_it, unit, month, p))
  set.seed(seed)
  sm <- .st_quiet(summary(ms))
  a <- sm$att
  row <- a[a$Level == "Average" & is.na(a$Time), ]
  list(est = row$Estimate, se = row$Std.Error, ci = c(row$lower_bound, row$upper_bound),
       nu = ms$nu, by_cohort = a[a$Level != "Average" & is.na(a$Time), c("Level", "Estimate")])
}

# SDID for each adoption cohort, aggregated by treated metro-months
.sdid_agg <- function(Y, cohort) {
  never <- which(is.infinite(cohort))
  gs <- sort(unique(cohort[is.finite(cohort)]))
  Tt <- ncol(Y)
  est <- numeric(0); cells <- numeric(0)
  for (g in gs) {
    trt <- which(cohort == g)
    e <- synthdid::synthdid_estimate(Y[c(never, trt), , drop = FALSE], length(never), g - 1)
    est[as.character(g)] <- as.numeric(e)
    cells[as.character(g)] <- length(trt) * (Tt - g + 1)
  }
  list(est = sum(cells * est) / sum(cells), by_cohort = est, cells = cells)
}

fit_sdid_cohort <- function(p, B = 500, seed = 1, cores = 1) {
  N <- max(p$unit); Tt <- max(p$month)
  Y <- matrix(NA_real_, N, Tt); Y[cbind(p$unit, p$month)] <- p$log_gmv
  cohort <- as.vector(tapply(p$cohort, p$unit, function(x) x[1]))
  point <- .sdid_agg(Y, cohort)
  set.seed(seed)
  draws <- vector("list", B); b <- 0L
  while (b < B) {
    ind <- sort(sample.int(N, N, replace = TRUE))
    if (any(is.finite(cohort[ind])) && any(is.infinite(cohort[ind]))) { b <- b + 1L; draws[[b]] <- ind }
  }
  boot <- unlist(parallel::mclapply(draws, function(ind) .sdid_agg(Y[ind, , drop = FALSE], cohort[ind])$est,
                                    mc.cores = cores))
  se <- stats::sd(boot)
  list(est = point$est, se = se, ci = point$est + c(-1, 1) * stats::qnorm(0.975) * se,
       by_cohort = point$by_cohort, cells = point$cells, B = B)
}

# Staleness guard for the cached S3 screen (Task 10 review, Important 3). The screen is
# recomputed only on demand (05b_stag_screen.R, RUN_SLOW=1), so an edit to a recipe above,
# to gen_S3() or a package upgrade could leave it describing an old recipe while res_stag.rds
# moves on. Returns a character vector of problems, empty when the cache is consistent:
#   * the screen's row for the stored seed must reproduce the point estimates of res_stag.rds
#     (twfe, cs, ife, mc_default, mc, multisynth) within tol; point estimates do not depend on
#     the number of bootstrap draws, so nboots 100 (screen) vs 500 (res_stag) does not matter;
#   * the package versions recorded with the screen must be the installed ones.
screen_staleness <- function(screen, st, stored_seed, tol = 1e-8) {
  if (is.null(screen)) return("there is no screen")
  row <- screen$per_seed[screen$per_seed$seed == stored_seed, , drop = FALSE]
  if (nrow(row) != 1L) return(paste("the screen has no row for the stored seed", stored_seed))
  probs <- character(0)
  for (k in c("twfe", "cs", "ife", "mc_default", "mc", "multisynth")) {
    a <- row[[paste0(k, "_est")]]; b <- st[[k]]$est
    # a missing column or a non-numeric entry must be reported, never pass vacuously
    if (!(is.numeric(a) && length(a) == 1L && is.finite(a) && is.numeric(b) && length(b) == 1L && is.finite(b))) {
      probs <- c(probs, sprintf("%s estimate for seed %d: missing or non-finite in the screen row or in res_stag.rds", k, stored_seed))
      next
    }
    if (!isTRUE(abs(a - b) < tol))
      probs <- c(probs, sprintf("%s estimate for seed %d: screen %.10g, res_stag.rds %.10g", k, stored_seed, a, b))
  }
  vers <- screen$versions
  if (length(vers) == 0L || is.null(names(vers))) {
    probs <- c(probs, "the screen records no package versions, so it cannot be checked against the installed ones")
  } else {
    inst <- vapply(names(vers), function(k) as.character(utils::packageVersion(k)), character(1))
    bad <- names(inst)[inst != vers[names(inst)]]
    if (length(bad))
      probs <- c(probs, paste0(bad, " ", inst[bad], " is installed but the screen was computed with ", vers[bad]))
  }
  probs
}
