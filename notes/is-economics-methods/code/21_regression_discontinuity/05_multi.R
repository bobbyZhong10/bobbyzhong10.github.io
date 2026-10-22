# =====================================================================
# 21_regression_discontinuity / 05_multi.R
# Block R4 (aster_R4.rds): badge cutoffs 70 (vintage), 80 (home), 90 (electronics) with
# effects 0.40, 0.25, 0.10; parallel Y(0) curves across categories.
#   by_cutoff     rdmc (Cattaneo, Titiunik and Vazquez-Bare 2020): rdrobust at each cutoff on
#                 its own sample and MSE-optimal h; robust intervals
#   weights       ruling R3: pooled_cutoff_weights() at one common bandwidth (rdmc's pooled
#                 h), the sample analog of the Lemma 1 weights f(c | c) P(C = c) of Cattaneo,
#                 Keele, Titiunik and Vazquez-Bare (2016)
#   pooled        the weighted estimator sum(weights * by_cutoff$est) (the Lemma 1 identity
#                 holds exactly); interval from the bias-corrected estimates and robust
#                 variances of the independent cutoff samples, as rdmc's own "weighted" row
#   pooled_norm   rdmc's pooled row: rdrobust on the normalized score x - c_i, what applied
#                 work runs; it targets sum(weights * tau_c) (ruling R3 test)
#   extrap        Cattaneo, Keele, Titiunik and Vazquez-Bare (2021): the vintage effect at the
#                 interior score 75 from the parallel-Y(0) (constant bias) assumption,
#                   tau_70(75) = E[Y | 75, vintage, treated] - E[Y | 75, home, control]
#                                - (E[Y | 70, vintage, control] - E[Y | 70, home, control]),
#                 four triangular local-linear fits (h = 4), stacked sandwich standard error
# Block R4g (aster_R4g.rds): same-day zone = within 25 km of the hub; six boundary points.
#   geo           rdms (Keele and Titiunik 2015 style): rdrobust on the distance to each
#                 boundary point, signed by the zone; pooled row on the signed distance to
#                 the zone boundary (xnorm = 25 - r)
# Output: res_multi.rds. Run from this directory:  Rscript 05_multi.R   (a few seconds)
# =====================================================================

suppressWarnings(suppressMessages({
  library(rdrobust)
  library(rdmulti)
}))

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || here == "") here <- "."
source(file.path(here, "rd_fns.R"))
R4 <- readRDS(file.path(here, "aster_R4.rds"))
G  <- readRDS(file.path(here, "aster_R4g.rds"))
tr <- readRDS(file.path(here, "truth21.rds"))

quiet <- function(expr) { out <- NULL; invisible(utils::capture.output(out <- suppressWarnings(expr))); out }
z975 <- qnorm(0.975)
EXTRAP_H <- 4

## ---- R4: per-cutoff, weights, pooled -------------------------------------------------------------
cuts <- sort(unique(R4$cutoff))
stopifnot(identical(cuts, sort(R4$meta$cutoffs)), length(cuts) == 3L)
o <- quiet(rdmc(Y = R4$y, X = R4$x, C = R4$cutoff))
# rdmc columns: 1..3 the cutoffs (sorted), 4 its "weighted" row, 5 the pooled row (V_cl, CI_cl
# carry no column names, so columns are taken by position)
kc <- seq_along(cuts); k_w <- length(cuts) + 1L; k_p <- length(cuts) + 2L
stopifnot(identical(colnames(o$Coefs), c(as.character(kc), "weighted", "pooled")))
by_cutoff <- data.frame(cutoff = cuts, category = R4$meta$categories$category[match(cuts, R4$meta$categories$cutoff)],
                        est = as.numeric(o$Coefs[1, kc]), est_bc = as.numeric(o$B[1, kc]),
                        se = sqrt(as.numeric(o$V_cl[1, kc])), se_rob = sqrt(as.numeric(o$V[1, kc])),
                        h = as.numeric(o$H["left", kc]), n_left = as.integer(o$Nh["left", kc]),
                        n_right = as.integer(o$Nh["right", kc]), stringsAsFactors = FALSE)
by_cutoff$ci <- lapply(kc, function(k) as.numeric(o$CI[, k]))
by_cutoff$truth <- as.numeric(tr$R4$tau_c[as.character(cuts)])
by_cutoff$n_sample <- vapply(cuts, function(cc) sum(R4$cutoff == cc), integer(1))
# rdmc's per-cutoff fit is rdrobust on the cutoff's own sample: check one of them
chk <- rdrobust(R4$y[R4$cutoff == cuts[2]], R4$x[R4$cutoff == cuts[2]], c = cuts[2])
stopifnot(abs(chk$coef["Conventional", 1] - by_cutoff$est[2]) < 1e-10)

h_common <- as.numeric(o$H["left", k_p])
w <- pooled_cutoff_weights(R4$cutoff, R4$x, h = h_common)
stopifnot(identical(names(w), as.character(cuts)))
pooled_est <- sum(w * by_cutoff$est)
pooled_bc <- sum(w * by_cutoff$est_bc)
pooled_se_rob <- sqrt(sum(w^2 * by_cutoff$se_rob^2))
pooled <- list(est = pooled_est, se = sqrt(sum(w^2 * by_cutoff$se^2)), est_bc = pooled_bc, se_rob = pooled_se_rob,
               ci = c(pooled_bc - z975 * pooled_se_rob, pooled_bc + z975 * pooled_se_rob),
               target = sum(w * tr$R4$tau_c), h_common = h_common)
pooled_norm <- list(est = as.numeric(o$Coefs[1, k_p]), est_bc = as.numeric(o$B[1, k_p]),
                    se_rob = sqrt(as.numeric(o$V[1, k_p])), ci = as.numeric(o$CI[, k_p]),
                    h = h_common, n_left = as.integer(o$Nh["left", k_p]), n_right = as.integer(o$Nh["right", k_p]),
                    target = sum(w * tr$R4$tau_c))
# rdmc's own "weighted" row: weights are each cutoff's share of the effective sample N_h at
# its own bandwidth, not the Lemma 1 weights
rdmc_weighted <- list(est = as.numeric(o$Coefs[1, k_w]), ci = as.numeric(o$CI[, k_w]),
                      weights = stats::setNames(as.numeric(o$W[1, ]), as.character(cuts)))
# the same weighted estimator with every cutoff fitted at the common bandwidth
est_common <- vapply(cuts, function(cc) {
  s <- R4$cutoff == cc
  rdrobust(R4$y[s], R4$x[s], c = cc, h = h_common)$coef["Conventional", 1]
}, numeric(1))
by_cutoff_common <- list(est = stats::setNames(est_common, as.character(cuts)), pooled = sum(w * est_common), h = h_common)

## ---- R4: extrapolation to score 75 (CKTV 2021) ---------------------------------------------------------
# local-linear fit at x0 as a linear smoother: est = sum(l * y); residuals for the sandwich
ll_smoother <- function(idx, x0, h) {
  xi <- R4$x[idx]; yi <- R4$y[idx]
  kw <- pmax(1 - abs(xi - x0) / h, 0); keep <- kw > 0
  X <- cbind(1, xi[keep] - x0); W <- kw[keep]
  A <- solve(crossprod(X, W * X))
  l <- as.numeric((A %*% t(X * W))[1, ])                     # weights on y of the intercept
  beta <- A %*% crossprod(X, W * yi[keep])
  e <- as.numeric(yi[keep] - X %*% beta)
  n <- sum(keep)
  list(units = idx[keep], l = l, e = e, est = sum(l * yi[keep]), dof = n / (n - 2), n = n)
}
extrap_fit <- function(h) {
  ce <- R4$meta$extrap_cutoff; cr <- R4$meta$ref_cutoff; x0 <- R4$meta$x_extrap   # 70, 80, 75
  V <- R4$cutoff == ce; Hm <- R4$cutoff == cr
  parts <- list(a  = ll_smoother(which(V & R4$x >= ce), x0, h),     # vintage, treated, at 75
                b  = ll_smoother(which(Hm & R4$x < cr), x0, h),     # home, control, at 75
                cl = ll_smoother(which(V & R4$x < ce), ce, h),      # vintage, control, at 70 (boundary)
                ch = ll_smoother(which(Hm & R4$x < cr), ce, h))     # home, control, at 70 (interior)
  sgn <- c(a = 1, b = -1, cl = -1, ch = 1)
  est <- sum(sgn * vapply(parts, `[[`, numeric(1), "est"))
  # stacked sandwich: each unit's influence summed over the fits it enters (b and ch share
  # home controls between 71 and 74, so their covariance is kept)
  infl <- numeric(length(R4$y))
  for (nm in names(parts)) {
    p <- parts[[nm]]
    infl[p$units] <- infl[p$units] + sgn[[nm]] * sqrt(p$dof) * p$l * p$e
  }
  se <- sqrt(sum(infl^2))
  se_indep <- sqrt(sum(vapply(parts, function(p) p$dof * sum((p$l * p$e)^2), numeric(1))))
  list(est = est, se = se, ci = c(est - z975 * se, est + z975 * se), se_indep = se_indep, h = h,
       parts = vapply(parts, `[[`, numeric(1), "est"), n = vapply(parts, `[[`, numeric(1), "n"),
       naive_cross = parts$a$est - parts$b$est,                      # no constant-bias correction
       bias_at_70 = parts$cl$est - parts$ch$est)                     # Y(0) gap vintage - home at 70
}
extrap <- extrap_fit(EXTRAP_H)
extrap$truth <- tr$R4$tau_extrap
extrap$sens <- lapply(c(h3 = 3, h5 = 5), function(h) { e <- extrap_fit(h); e[c("est", "se", "ci", "h")] })

## ---- R4g: geographic RD --------------------------------------------------------------------------------
bp <- G$meta$boundary
r_hub <- sqrt(G$x1^2 + G$x2^2)
xnorm <- G$meta$radius - r_hub                   # signed distance to the zone boundary, > 0 inside
g <- quiet(rdms(Y = G$y, X = G$x1, C = bp$b1, X2 = G$x2, zvar = G$d, C2 = bp$b2, xnorm = xnorm))
# rdms columns: 1..6 the boundary points, 7 the pooled row on xnorm (taken by position)
kb <- seq_len(nrow(bp)); kb_p <- nrow(bp) + 1L
stopifnot(identical(colnames(g$Coefs), c(as.character(kb), "pooled")))
geo <- data.frame(point = seq_len(nrow(bp)), angle_deg = bp$angle_deg, b1 = bp$b1, b2 = bp$b2,
                  est = as.numeric(g$Coefs[1, kb]), est_bc = as.numeric(g$B[1, kb]),
                  se_rob = sqrt(as.numeric(g$V[1, kb])), h = as.numeric(g$H["left", kb]),
                  n_left = as.integer(g$Nh["left", kb]), n_right = as.integer(g$Nh["right", kb]))
geo$ci <- lapply(kb, function(k) as.numeric(g$CI[, k]))
geo$truth <- as.numeric(tr$R4g$tau_b)
geo$covers <- mapply(function(ci, v) ci[1] <= v && v <= ci[2], geo$ci, geo$truth)
# pooled along the whole boundary: directions are uniform, so the seller density along the
# circle is uniform and the target is the mean effect over the circle, tau_center + tau_east
# E[cos] = tau_center (0.20); the mean over the six equally spaced points is the same number
geo_pooled <- list(est = as.numeric(g$Coefs[1, kb_p]), ci = as.numeric(g$CI[, kb_p]),
                   h = as.numeric(g$H["left", kb_p]), n_left = as.integer(g$Nh["left", kb_p]),
                   n_right = as.integer(g$Nh["right", kb_p]), target = G$meta$tau_center)
stopifnot(abs(mean(tr$R4g$tau_b) - geo_pooled$target) < 1e-12)
geo_naive <- mean(G$y[G$d == 1]) - mean(G$y[G$d == 0])     # inside vs outside, whole metro

res <- list(by_cutoff = by_cutoff, weights = w, pooled = pooled, pooled_norm = pooled_norm,
            simple_mean = mean(by_cutoff$est), rdmc_weighted = rdmc_weighted, by_cutoff_common = by_cutoff_common,
            extrap = extrap, geo = geo, geo_pooled = geo_pooled, geo_naive = geo_naive,
            truth = list(tau_c = tr$R4$tau_c, omega = tr$R4$omega, tau_pooled = tr$R4$tau_pooled,
                         tau_extrap = tr$R4$tau_extrap, tau_b = tr$R4g$tau_b, simple_mean = mean(tr$R4$tau_c)),
            sample = list(N4 = length(R4$x), n_by_cutoff = stats::setNames(by_cutoff$n_sample, as.character(cuts)),
                          N4g = length(G$y), n_inside = sum(G$d), n_outside = sum(1 - G$d),
                          n_within4 = t(vapply(seq_len(nrow(bp)), function(k) {
                            dk <- sqrt((G$x1 - bp$b1[k])^2 + (G$x2 - bp$b2[k])^2) <= 4
                            c(inside = sum(dk & G$d == 1), outside = sum(dk & G$d == 0))
                          }, numeric(2)))))

saveRDS(res, file.path(here, "res_multi.rds"))

## ---- report ---------------------------------------------------------------------------------------------
ci_txt <- function(ci) sprintf("[%.3f, %.3f]", ci[1], ci[2])
cat("================ R4: multi-cutoff ================\n")
for (i in seq_len(nrow(by_cutoff)))
  cat(sprintf("cutoff %d (%s): %.4f robust %s  h %.2f  eff. n %d / %d  truth %.3f\n", by_cutoff$cutoff[i],
              by_cutoff$category[i], by_cutoff$est[i], ci_txt(by_cutoff$ci[[i]]), by_cutoff$h[i], by_cutoff$n_left[i],
              by_cutoff$n_right[i], by_cutoff$truth[i]))
cat(sprintf("Lemma 1 sample weights at common h %.2f: %s (population %s)\n", h_common,
            paste(sprintf("%.3f", w), collapse = ", "), paste(sprintf("%.3f", tr$R4$omega), collapse = ", ")))
cat(sprintf("pooled weighted estimator: %.4f robust %s  target %.4f  (simple mean of estimates %.4f, of truths %.3f)\n",
            pooled$est, ci_txt(pooled$ci), pooled$target, res$simple_mean, mean(tr$R4$tau_c)))
cat(sprintf("pooled normalized score  : %.4f robust %s  h %.2f\n", pooled_norm$est, ci_txt(pooled_norm$ci), h_common))
cat(sprintf("rdmc weighted row (N_h weights %s): %.4f\n", paste(sprintf("%.3f", rdmc_weighted$weights), collapse = ", "),
            rdmc_weighted$est))
cat(sprintf("weighted estimator at the common h: %.4f\n", by_cutoff_common$pooled))
cat(sprintf("extrapolation tau_70(75): %.4f  %s  se %.4f (independent-parts se %.4f), h %d; truth %.3f\n", extrap$est,
            ci_txt(extrap$ci), extrap$se, extrap$se_indep, EXTRAP_H, extrap$truth))
cat(sprintf("  naive cross-category gap at 75 %.4f; Y(0) gap at 70 %.4f; h 3: %.4f, h 5: %.4f\n", extrap$naive_cross,
            extrap$bias_at_70, extrap$sens$h3$est, extrap$sens$h5$est))
cat("================ R4g: geographic RD ================\n")
for (i in seq_len(nrow(geo)))
  cat(sprintf("point %d (%3d deg): %.3f robust %s  h %.1f  n %d / %d  truth %.2f  covers %s\n", i, geo$angle_deg[i],
              geo$est[i], ci_txt(geo$ci[[i]]), geo$h[i], geo$n_left[i], geo$n_right[i], geo$truth[i], geo$covers[i]))
cat(sprintf("share covered %.3f; pooled boundary %.3f %s (target %.2f); naive inside - outside %.3f\n",
            mean(geo$covers), geo_pooled$est, ci_txt(geo_pooled$ci), geo_pooled$target, geo_naive))
