# ------------------------------------------------------------------
# 08_honestdid.R -- Rambachan-Roth (2023) sensitivity analysis on the
# Callaway-Sant'Anna dynamic event study of the Kestrel panel.
#
# Event study:  did::att_gt(control_group = "notyettreated") followed by
#   did::aggte(type = "dynamic", min_e = -6, max_e = 12, bstrap = FALSE)
# so that the influence function is available. HonestDiD wants the
# pre-period coefficients to be LEVEL differences relative to the reference
# period e = -1, so att_gt uses base_period = "universal" (04_modern.R uses
# the default "varying" base period). Post-treatment ATT(g,t) are identical
# under both base periods (the base is g-1 either way), so the post-period
# event-study coefficients below coincide with cs_es in res_modern.rds; only
# the pre-period coefficients and the standard errors (analytic here,
# multiplier bootstrap in 04) differ.
#
# Covariance: V = IF' IF / n^2 with IF the n x K influence-function matrix
# of the K event-study coefficients (n = 60 cities), as in the HonestDiD
# example for the did package. The reference period e = -1 (coefficient 0,
# influence function 0) is dropped, leaving 5 pre-period coefficients
# (e = -6..-2) and 13 post-period coefficients (e = 0..12).
#
# Targets: theta(e = 6) and theta(e = 12), each a single post-period
# coefficient (l_vec = basis vector).
#   SD: DeltaSD(M), smoothness: second differences of the differential
#       trend are at most M (outcome units: log points). FLCI.
#   RM: DeltaRM(Mbar), relative magnitudes: post-period changes of the
#       differential trend are at most Mbar times the largest pre-period
#       change. Conditional least-favorable (C-LF) test inversion.
#
# ---- Display grids (the tables in the chapter) ----------------------
#   Mvec    = 0, 0.0001, 0.0002, 0.0005, 0.001, 0.0015, 0.002
#   Mbarvec = 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1
# (a grid is extended by doubling its last value if it never reaches a set
# that contains 0).
#
# ---- RM search grid (fix 1) -----------------------------------------
# The C-LF confidence set is the range of accepted points on a theta
# search grid. HonestDiD's default grid is +/- 20 se(theta-hat) around 0,
# which cuts off the e = 12 upper bounds at Mbar = 0.75 and 1 (they sat on
# the grid limit, 0.2827, and are really 0.352 and 0.426). We therefore pass
# an explicit grid: half-width L = 0.6 and step 0.0002 (gridPoints =
# 2 * ceiling(L / 0.0002) + 1 = 6001, odd so that theta = 0 is a grid
# point). L is widened by a factor 1.5 until no stored RM row has a bound on
# the grid limit; the `truncated` column of $e6 and $e12 records the check
# (all FALSE). RM bounds are inner approximations on the grid, accurate to
# the step 0.0002: quote at most three decimals.
#
# ---- Breakdown values (fix 2) ----------------------------------------
# Breakdown = the smallest M (Mbar) on a fine scan grid whose robust set
# contains 0 (lb <= 0 <= ub); the scan runs from 0 to the first display-grid
# value that contains 0. Scan steps (= the resolution of the stored value):
# SD 1e-5, RM 0.005. The fixed-length-CI (FLCI) behind DeltaSD is NOT
# monotone in M at fine scale for e = 6: the lower bound jumps by up to about
# 0.001 between neighbouring M (lb at M = 0.000155: +0.00006; 0.00016:
# -0.00007; 0.000165: +0.00071; 0.00017: -0.00035; 0.000175: +0.00044;
# 0.00017656: +0.00056; 0.00017666: -0.00054; 0.000185: +0.00039), so the
# sets are not nested at that scale and a bisection landed on a noise-driven
# crossing (0.000177). The smallest scanned M is stable: 0.00016 on the 1e-5
# scan, 0.000157 on a 1e-6 re-scan of the last 1e-5 (`breakdown_fine`). The
# e = 12 SD sets and all RM sets are smooth. `breakdown_reversals` counts
# scanned M above the breakdown whose set again excludes 0, and `sd_jitter`
# records the SD bound noise: the largest absolute deviation of lb (ub) from
# the straight line through the five evaluations at M * (1, 1 -/+ 0.5%,
# 1 -/+ 1%), over the displayed M > 0 (a smooth function would give about 0).
# Quote SD bounds to two decimals only; quote breakdowns to the stated
# resolution.
#
# The scan evaluates RM with a 7-point theta grid (same step, centred on 0,
# so it shares the points near 0 with the display grid; the decision "0
# accepted" does not depend on the other grid points) in parallel with
# forked workers; results are checked against the display tables at every
# M that lies on both. Below the breakdown the 7-point grid often accepts
# no point at all, which HonestDiD reports as lb = Inf, ub = -Inf with
# "no non-missing arguments to min/max" warnings; they are expected and
# tallied, as are "CI is open at one of the endpoints" warnings for sets
# that touch the 7-point grid edge.
#
# Output res_honest.rds:
#   $e6, $e12   data.frame(method, M, lb, ub, point, ci_method, grid_limit,
#               truncated); `point` is the original CS estimate theta-hat(e)
#   $breakdown  list(e6 = list(SD, RM), e12 = ...): smallest scanned M
#   $breakdown_resolution  same shape: the scan step
#   $breakdown_fine        same shape: SD 1e-6 re-scan result (RM: NA)
#   $breakdown_grid        same shape: smallest display-grid M containing 0
#   $breakdown_reversals   same shape: scanned M > breakdown excluding 0
#   $breakdown_set  data.frame(target, method, M, lb, ub, point, prev_M):
#               sets at the breakdown computed with the display settings
#   $orig       conventional CIs from HonestDiD::constructOriginalCS
#   $es, $beta, $V   coefficients, SEs, estimate vector and covariance used
#   $meta       grids, steps, RM grid, jitter, versions, warning tally, ...
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(did)
  library(HonestDiD)
  library(parallel)
})

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dt <- readRDS(file.path(code_dir, "kestrel_panel.rds"))
m_modern <- readRDS(file.path(code_dir, "res_modern.rds"))

t_start  <- Sys.time()
alpha    <- 0.05
targets  <- c(e6 = 6, e12 = 12)
Mvec     <- c(0, 0.0001, 0.0002, 0.0005, 0.001, 0.0015, 0.002)
Mbarvec  <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1)
step_scan     <- c(SD = 1e-5, RM = 0.005)   # breakdown scan steps = resolution
step_sd_fine  <- 1e-6
rm_step_theta <- 0.0002                      # RM theta-grid step
rm_L0         <- 0.6                         # initial theta-grid half-width
jitter_rel    <- c(-0.01, -0.005, 0.005, 0.01)
tol0          <- 1e-9                        # theta = 0 is a grid point up to rounding
n_cores       <- max(1L, min(12L, detectCores() - 2L))

## ---- Warnings are muffled but tallied (also from forked workers) -----
## HonestDiD (0.2.8) triggers a benign R >= 4.5 deprecation warning
## ("Recycling array of length 1 in vector-array arithmetic"), and for the
## RM sets "CI is open at one of the endpoints" when a bound is on the grid.
warn_tally <- integer()
add_warns <- function(w) for (msg in w)
  warn_tally[msg] <<- (if (msg %in% names(warn_tally)) warn_tally[[msg]] else 0L) + 1L
capture <- function(expr) {
  w <- character()
  v <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = v, warns = w)
}
par_map <- function(X, f) {                  # forked lapply; stops on any error
  out <- mclapply(X, function(x) capture(f(x)), mc.cores = n_cores)
  for (o in out) {
    if (!is.list(o) || inherits(o$value, "error"))
      stop("worker failed: ", if (is.list(o)) conditionMessage(o$value) else as.character(o))
    add_warns(o$warns)
  }
  lapply(out, `[[`, "value")
}

## ---- Event study with the influence function ---------------------
cs_u <- att_gt(yname = "lgmv", tname = "month", idname = "city_id", gname = "g",
               data = dt, control_group = "notyettreated", base_period = "universal",
               clustervars = "city_id", bstrap = FALSE, cband = FALSE)
es <- aggte(cs_u, type = "dynamic", min_e = -6, max_e = 12, bstrap = FALSE, cband = FALSE)

IF   <- es$inf.function$dynamic.inf.func.e
keep <- es$egt != -1                      # drop the reference period
stopifnot(sum(abs(IF[, !keep])) == 0, es$att.egt[!keep] == 0)
IFk  <- IF[, keep]
n    <- nrow(IFk)
V    <- crossprod(IFk) / n / n
beta <- es$att.egt[keep]
ev   <- es$egt[keep]
stopifnot(max(abs(sqrt(diag(V)) - es$se.egt[keep])) < 1e-12,   # V reproduces aggte SEs
          all(diff(ev) > 0))
npre  <- sum(ev < 0)
npost <- sum(ev >= 0)
stopifnot(npre == 5, npost == 13, identical(ev, as.numeric(c(-6:-2, 0:12))))

## post-period coefficients equal the CS dynamic estimates saved by 04
cmp <- merge(data.table(e = ev, est = beta), m_modern$cs_es, by = "e")
stopifnot(max(abs(cmp[e >= 0, est - estimate])) < 1e-10)

es_used <- data.frame(e = ev, estimate = beta, se = sqrt(diag(V)))
cat("== CS dynamic event study used by HonestDiD (universal base, e = -1 omitted) ==\n")
print(es_used, digits = 4)

basis <- function(e0) matrix(as.numeric(seq_len(npost) == e0 + 1), ncol = 1)
if (exists(".basisVector", envir = asNamespace("HonestDiD")))
  stopifnot(identical(unname(basis(6)),
                      unname(HonestDiD:::.basisVector(index = 7, size = npost))))

## ---- Evaluators: one robust set for one M, one target ------------------
sd_eval <- function(M, tg)
  createSensitivityResults(betahat = beta, sigma = V, numPrePeriods = npre,
                           numPostPeriods = npost, Mvec = M, l_vec = basis(tg), alpha = alpha)
grid_points <- function(L) 2 * ceiling(L / rm_step_theta - 1e-9) + 1   # odd: theta = 0 is a grid point
rm_eval <- function(Mbar, tg, L)             # theta grid: step rm_step_theta, centred on 0
  createSensitivityResults_relativeMagnitudes(
    betahat = beta, sigma = V, numPrePeriods = npre, numPostPeriods = npost,
    bound = "deviation from parallel trends", method = "C-LF", Mbarvec = Mbar,
    l_vec = basis(tg), alpha = alpha,
    grid.lb = -L, grid.ub = L, gridPoints = grid_points(L))
rm_L_scan <- 3 * rm_step_theta               # 7-point grid with the display grid's points near 0
has0 <- function(lb, ub) lb <= tol0 & ub >= -tol0

rows_of <- function(lst, method) {
  data.table(M = vapply(lst, function(r) if (method == "SD") r$M else r$Mbar, 0),
             lb = vapply(lst, function(r) r$lb, 0), ub = vapply(lst, function(r) r$ub, 0),
             ci_method = vapply(lst, function(r) r$method, ""))
}

## ---- Display tables ---------------------------------------------------
display_sd <- function(tg) {
  grid <- Mvec
  repeat {
    tab <- rows_of(par_map(grid, function(M) sd_eval(M, tg)), "SD")
    if (any(has0(tab$lb, tab$ub))) break
    grid <- c(grid, 2 * max(grid))
    stopifnot(length(grid) < 20)
  }
  tab
}
display_rm <- function(tg) {
  grid <- Mbarvec; L <- rm_L0
  repeat {
    tab <- rows_of(par_map(grid, function(Mb) rm_eval(Mb, tg, L)), "RM")
    trn <- tab$lb <= -L + tol0 | tab$ub >= L - tol0
    if (any(has0(tab$lb, tab$ub)) && !any(trn)) break
    if (any(trn)) L <- L * 1.5
    if (!any(has0(tab$lb, tab$ub))) grid <- c(grid, 2 * max(grid))
    stopifnot(L < 6, length(grid) < 16)
  }
  tab[, truncated := trn]
  attr(tab, "L") <- L
  tab
}

## ---- Scan for the breakdown ------------------------------------------
scan_breakdown <- function(tg, method, hi) {
  step <- step_scan[[method]]
  Ms   <- round(seq(0, hi, by = step), 10)
  stopifnot(abs(tail(Ms, 1) - hi) < 1e-9)
  f <- if (method == "SD") function(M) sd_eval(M, tg) else function(M) rm_eval(M, tg, rm_L_scan)
  tab <- rows_of(par_map(Ms, f), method)
  hit <- has0(tab$lb, tab$ub)
  stopifnot(!hit[1], any(hit))               # M = 0 excludes 0 and the scan reaches 0
  first <- which(hit)[1]
  list(first = Ms[first], step = step, tab = tab, hit = hit,
       reversals = sum(!hit[first:length(hit)]))
}

res_tab <- list(); disp <- list(); scans <- list(); bd <- list(); bd_res <- list()
bd_fine <- list(); bd_grid <- list(); bd_rev <- list(); orig <- list(); jit <- list()
bd_set <- list(); rm_L <- list()

for (nm in names(targets)) {
  tg    <- targets[[nm]]
  point <- beta[ev == tg]
  oc    <- constructOriginalCS(betahat = beta, sigma = V, numPrePeriods = npre,
                               numPostPeriods = npost, l_vec = basis(tg), alpha = alpha)
  orig[[nm]] <- list(lb = as.numeric(oc$lb), ub = as.numeric(oc$ub))
  cat(sprintf("\n== Target %s: theta-hat = %.4f, conventional CI [%.4f, %.4f] ==\n",
              nm, point, orig[[nm]]$lb, orig[[nm]]$ub))

  t0 <- Sys.time()
  tab_sd <- display_sd(tg); tab_rm <- display_rm(tg)
  rm_L[[nm]] <- attr(tab_rm, "L")
  cat(sprintf("display sets: %.0f s; RM theta grid half-width %.2f (%d points)\n",
              as.numeric(Sys.time() - t0, units = "secs"), rm_L[[nm]], grid_points(rm_L[[nm]])))

  ## smallest display-grid value containing 0, and the fine scans
  bd_grid[[nm]] <- list(SD = tab_sd[has0(lb, ub)][1, M], RM = tab_rm[has0(lb, ub)][1, M])
  t0 <- Sys.time()
  sc <- list(SD = scan_breakdown(tg, "SD", bd_grid[[nm]]$SD),
             RM = scan_breakdown(tg, "RM", bd_grid[[nm]]$RM))
  scans[[nm]] <- sc
  for (me in c("SD", "RM")) {
    bd[[nm]][[me]]     <- sc[[me]]$first
    bd_res[[nm]][[me]] <- sc[[me]]$step
    bd_rev[[nm]][[me]] <- sc[[me]]$reversals
  }
  ## consistency with the display tables at every M on both grids
  for (me in c("SD", "RM")) {
    d  <- if (me == "SD") tab_sd else tab_rm
    st <- sc[[me]]$tab
    for (i in seq_len(nrow(d))) {
      j <- which(abs(st$M - d$M[i]) < 1e-9)
      if (length(j) == 1) {
        stopifnot(has0(d$lb[i], d$ub[i]) == has0(st$lb[j], st$ub[j]))
        if (me == "SD") stopifnot(abs(d$lb[i] - st$lb[j]) < 1e-9, abs(d$ub[i] - st$ub[j]) < 1e-9)
      }
    }
  }
  ## SD: 1e-6 re-scan of the last 1e-5 before the breakdown
  Mf <- round(seq(max(0, bd[[nm]]$SD - step_scan[["SD"]]), bd[[nm]]$SD, by = step_sd_fine), 10)
  ft <- rows_of(par_map(Mf, function(M) sd_eval(M, tg)), "SD")
  bd_fine[[nm]] <- list(SD = Mf[which(has0(ft$lb, ft$ub))[1]], RM = NA_real_)
  cat(sprintf("scans: %.0f s; SD breakdown %.5f (step %g; 1e-6 re-scan %.6f; reversals %d); RM breakdown %.3f (step %g; reversals %d)\n",
              as.numeric(Sys.time() - t0, units = "secs"),
              bd[[nm]]$SD, step_scan[["SD"]], bd_fine[[nm]]$SD, bd_rev[[nm]]$SD,
              bd[[nm]]$RM, step_scan[["RM"]], bd_rev[[nm]]$RM))

  ## SD bound noise: range of lb and ub over M(1 +/- 0.5%, 1%) at each displayed M > 0
  jin <- CJ(i = which(tab_sd$M > 0), r = c(0, jitter_rel))
  jr  <- rows_of(par_map(seq_len(nrow(jin)), function(k)
                 sd_eval(tab_sd$M[jin$i[k]] * (1 + jin$r[k]), tg)), "SD")
  jin[, `:=`(lb = jr$lb, ub = jr$ub)]
  jj <- jin[, {
    Mi <- tab_sd$M[i[1]] * (1 + r)
    .(M = tab_sd$M[i[1]],
      lb_resid = max(abs(resid(lm(lb ~ Mi)))), ub_resid = max(abs(resid(lm(ub ~ Mi)))))
  }, by = i]
  jit[[nm]] <- list(by_M = as.data.frame(jj[, .(M, lb_resid, ub_resid)]),
                    max_lb_resid = max(jj$lb_resid), max_ub_resid = max(jj$ub_resid))

  ## the sets at the breakdown, with the display settings (RM: full theta grid)
  st_sd <- sd_eval(bd[[nm]]$SD, tg)
  st_rm <- rm_eval(bd[[nm]]$RM, tg, rm_L[[nm]])
  bd_set[[nm]] <- data.frame(target = nm, method = c("SD", "RM"),
                             M = c(bd[[nm]]$SD, bd[[nm]]$RM),
                             lb = c(st_sd$lb, st_rm$lb), ub = c(st_sd$ub, st_rm$ub),
                             point = point,
                             prev_M = c(bd[[nm]]$SD - step_scan[["SD"]], bd[[nm]]$RM - step_scan[["RM"]]))
  stopifnot(all(has0(bd_set[[nm]]$lb, bd_set[[nm]]$ub)))

  mk <- function(tab, method, L) data.frame(
    method = method, M = tab$M, lb = tab$lb, ub = tab$ub, point = point, ci_method = tab$ci_method,
    grid_limit = if (method == "RM") L else NA_real_,
    truncated = if (method == "RM") tab$truncated else FALSE, stringsAsFactors = FALSE)
  res_tab[[nm]] <- rbind(mk(tab_sd, "SD", NA), mk(tab_rm, "RM", rm_L[[nm]]))
  rownames(res_tab[[nm]]) <- NULL
  cat("Robust confidence sets:\n"); print(res_tab[[nm]], digits = 4)
  cat(sprintf("SD bound noise (max |deviation| from a line over M +/- 0.5..1%%): lb %.5f, ub %.5f\n",
              jit[[nm]]$max_lb_resid, jit[[nm]]$max_ub_resid))
}
bd_set <- do.call(rbind, bd_set); rownames(bd_set) <- NULL

## the point estimate lies in every robust set, no RM bound sits on the grid limit
for (nm in names(res_tab)) with(res_tab[[nm]], {
  stopifnot(all(lb <= point & point <= ub), !any(truncated))
})
cat("\nBreakdown sets (display settings):\n"); print(bd_set, digits = 4)

if (length(warn_tally)) {
  cat("\nPackage warnings (tally):\n"); print(warn_tally)
}
runtime <- as.numeric(Sys.time() - t_start, units = "secs")
cat(sprintf("\nTotal runtime %.0f s (%d forked workers)\n", runtime, n_cores))

saveRDS(list(e6 = res_tab$e6, e12 = res_tab$e12,
             breakdown = bd, breakdown_resolution = bd_res, breakdown_fine = bd_fine,
             breakdown_grid = bd_grid, breakdown_reversals = bd_rev, breakdown_set = bd_set,
             orig = orig, es = es_used, beta = beta, V = V,
             meta = list(alpha = alpha, targets = targets, Mvec = Mvec, Mbarvec = Mbarvec,
                         step_scan = step_scan, step_sd_fine = step_sd_fine,
                         rm_grid = list(half_width = unlist(rm_L), step = rm_step_theta,
                                        gridPoints = vapply(rm_L, grid_points, 1)),
                         sd_jitter = jit, npre = npre, npost = npost, n_clusters = n,
                         base_period = "universal", control_group = "notyettreated",
                         versions = c(HonestDiD = as.character(packageVersion("HonestDiD")),
                                      did = as.character(packageVersion("did"))),
                         warnings = warn_tally, runtime_sec = runtime, n_cores = n_cores)),
        file.path(code_dir, "res_honest.rds"))
cat("Saved res_honest.rds\n")
