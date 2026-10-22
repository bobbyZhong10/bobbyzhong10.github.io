# tests/test_new.R
# New chapter-19 results: dCDH, HonestDiD, did2s / etwfe, numbers.json.
# Run from the code directory:  Rscript tests/test_new.R   (about 1 minute: it
# recomputes the robust sets at the stored HonestDiD breakdown values)
suppressPackageStartupMessages({ library(HonestDiD); library(parallel) })
m <- readRDS("res_modern.rds")
d <- readRDS("res_dcdh.rds"); f <- readRDS("res_family.rds"); h <- readRDS("res_honest.rds")

## ---- dCDH ----------------------------------------------------------
stopifnot(sum(d$twfe_weights$weight < 0) >= 1, abs(sum(d$twfe_weights$weight) - 1) < 1e-8)
stopifnot(nrow(d$dyn$effects) == 12, all(is.finite(d$dyn$effects$estimate)))
stopifnot(nrow(d$dyn$placebos) == 6, all(is.finite(d$dyn$placebos$estimate)))
# DID_l is the effect l-1 periods after adoption: it equals the CS dynamic
# (not-yet-treated) estimate at e = l - 1 (the chapter's teaching point)
ef <- d$dyn$effects
stopifnot(all(ef$e == ef$l - 1))
for (l in 1:12)
  stopifnot(abs(ef$estimate[ef$l == l] - m$cs_es$estimate[m$cs_es$e == l - 1]) < 1e-10)

## ---- did2s / etwfe vs BJS ------------------------------------------
# Gardner two-step and BJS imputation give the same overall ATT; so does
# etwfe with not-yet-treated controls (marginaleffects numerics: ~2e-8)
stopifnot(abs(f$did2s$est - m$imp$est) < 1e-6)
stopifnot(is.finite(f$etwfe$est), is.finite(f$etwfe$se), f$etwfe$se > 0)
stopifnot(abs(f$etwfe$est - m$imp$est) < 1e-6)
stopifnot(abs(f$etwfe_never$est - m$sa_att$est) < 1e-6)

## ---- HonestDiD: stored tables ---------------------------------------
has0 <- function(lb, ub) lb <= 1e-9 & ub >= -1e-9   # theta = 0 is a grid point up to rounding
for (nm in c("e6", "e12")) {
  x <- h[[nm]]
  stopifnot(is.data.frame(x), all(c("method", "M", "lb", "ub", "point", "truncated") %in% names(x)),
            all(x$method %in% c("SD", "RM")), all(c("SD", "RM") %in% x$method),
            all(x$lb <= x$point & x$point <= x$ub))
  # fix 1: no stored bound sits on the RM search-grid limit
  stopifnot(!any(x$truncated), all(abs(x$lb[x$method == "RM"]) < x$grid_limit[x$method == "RM"] - 1e-6),
            all(abs(x$ub[x$method == "RM"]) < x$grid_limit[x$method == "RM"] - 1e-6))
  # M = 0 (Mbar = 0): the robust set excludes 0
  for (me in c("SD", "RM")) stopifnot(!has0(x$lb[x$method == me & x$M == 0], x$ub[x$method == me & x$M == 0]))
}
# external anchors (independent rerun of e = 12 with grid.lb = -0.6, grid.ub = 0.6):
# true upper bounds 0.3515 (Mbar 0.75) and 0.4259 (Mbar 1), not the old 0.2827 cap
e12rm <- h$e12[h$e12$method == "RM", ]
stopifnot(abs(e12rm$ub[e12rm$M == 0.75] - 0.3515) < 0.002, abs(e12rm$ub[e12rm$M == 1] - 0.4259) < 0.002)

## ---- HonestDiD: breakdown values -----------------------------------
# breakdown = smallest scanned M whose set contains 0; the set at the breakdown
# contains 0, the set one scan step below it does not, and so does none of a
# sample of smaller M (the sets are NOT monotone in M for e = 6 SD at fine scale)
stopifnot(!is.null(h$breakdown), !is.null(h$breakdown_resolution), !is.null(h$breakdown_set))
for (nm in c("e6", "e12")) for (me in c("SD", "RM")) {
  bd <- h$breakdown[[nm]][[me]]; rs <- h$breakdown_resolution[[nm]][[me]]
  stopifnot(is.finite(bd), bd > 0, is.finite(rs), rs > 0)
  s <- h$breakdown_set[h$breakdown_set$target == nm & h$breakdown_set$method == me, ]
  stopifnot(nrow(s) == 1, abs(s$M - bd) < 1e-12, has0(s$lb, s$ub), abs(s$prev_M - (bd - rs)) < 1e-12)
}
# external anchor: the first e = 6 SD M with 0 inside is about 0.000157 (not 0.000177)
stopifnot(abs(h$breakdown$e6$SD - 0.000157) < 2e-5)
stopifnot(abs(h$breakdown$e12$SD - 0.000853) < 2e-5)

# independent recomputation from the stored beta and V
npre <- h$meta$npre; npost <- h$meta$npost; rg <- h$meta$rm_grid
basis <- function(e0) matrix(as.numeric(seq_len(npost) == e0 + 1), ncol = 1)
eval_set <- function(nm, me, M) {
  tg <- h$meta$targets[[nm]]
  if (me == "SD")
    createSensitivityResults(betahat = h$beta, sigma = h$V, numPrePeriods = npre,
                             numPostPeriods = npost, Mvec = M, l_vec = basis(tg), alpha = h$meta$alpha)
  else {
    L <- rg$half_width[[nm]]
    createSensitivityResults_relativeMagnitudes(
      betahat = h$beta, sigma = h$V, numPrePeriods = npre, numPostPeriods = npost,
      bound = "deviation from parallel trends", method = "C-LF", Mbarvec = M,
      l_vec = basis(tg), alpha = h$meta$alpha, grid.lb = -L, grid.ub = L,
      gridPoints = 2 * ceiling(L / rg$step - 1e-9) + 1)
  }
}
tasks <- list()
for (nm in c("e6", "e12")) for (me in c("SD", "RM")) {
  bd <- h$breakdown[[nm]][[me]]; rs <- h$breakdown_resolution[[nm]][[me]]
  below <- if (me == "SD") {                       # up to 16 scan points below the breakdown, incl. bd - rs
    k <- round(bd / rs) - 1
    unique(c(round(seq(0, k, length.out = min(k + 1, 16)), 0), k)) * rs
  } else c(bd - rs, bd / 2)
  tasks <- c(tasks, list(list(nm = nm, me = me, M = bd, expect = TRUE)),
             lapply(below, function(M) list(nm = nm, me = me, M = M, expect = FALSE)))
}
out <- mclapply(tasks, function(tk) {
  r <- suppressWarnings(eval_set(tk$nm, tk$me, tk$M))
  list(tk = tk, lb = r$lb, ub = r$ub)
}, mc.cores = max(1L, min(12L, detectCores() - 2L)))
for (o in out) {
  if (inherits(o, "try-error")) stop("recomputation failed: ", as.character(o))
  got <- has0(o$lb, o$ub)
  if (got != o$tk$expect)
    stop(sprintf("%s %s: set at M = %.8g is [%.5f, %.5f]; contains 0 = %s, expected %s",
                 o$tk$nm, o$tk$me, o$tk$M, o$lb, o$ub, got, o$tk$expect))
  if (o$tk$expect) {                               # recomputation reproduces the stored set
    s <- h$breakdown_set[h$breakdown_set$target == o$tk$nm & h$breakdown_set$method == o$tk$me, ]
    stopifnot(abs(o$lb - s$lb) < 1e-8, abs(o$ub - s$ub) < 1e-8)
  }
}

## ---- numbers.json: flat {name: number}, prefixed by script ----------
nj <- jsonlite::fromJSON("numbers.json")
stopifnot(is.list(nj), length(nj) > 50, all(vapply(nj, function(z) is.numeric(z) && length(z) == 1, TRUE)))
need <- c("twfe_static_est", "dcdh_dyn_l1_est", "honest_e6_sd_breakdown", "honest_e6_sd_breakdown_resolution",
          "honest_e12_rm_breakdown_resolution", "family_did2s_est", "family_etwfe_est", "modern_imp_est",
          "dgp_true_att_simple")
stopifnot(all(need %in% names(nj)))
for (nm in c("e6", "e12")) for (me in c("sd", "rm")) {
  stopifnot(abs(nj[[sprintf("honest_%s_%s_breakdown", nm, me)]] - h$breakdown[[nm]][[toupper(me)]]) < 1e-12,
            abs(nj[[sprintf("honest_%s_%s_breakdown_resolution", nm, me)]] -
                h$breakdown_resolution[[nm]][[toupper(me)]]) < 1e-12)
  # every RM bound in numbers.json is a true bound, away from the search-grid limit
  lim <- nj[[sprintf("honest_%s_rm_grid_limit", nm)]]
  keys <- grep(sprintf("^honest_%s_rm_M.*_(lb|ub)$", nm), names(nj), value = TRUE)
  stopifnot(length(keys) == 16, all(abs(unlist(nj[keys])) < lim - 1e-6))
}
cat("test_new: OK\n")
