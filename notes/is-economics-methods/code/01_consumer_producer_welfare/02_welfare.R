# ------------------------------------------------------------------
# 02_welfare.R -- block W1 (and W1g): welfare from the premium price experiment, and the Willig check.
#
# The analyst sees (p1, w, x1) for every subscriber, fits three demand forms by OLS (fit_demand_forms in
# welfare_fns.R: linear, log-linear, and the reciprocal-price form with an income term) and, for each form,
# computes for the price rise from 2 to 4 USD, averaging over the subscribers' own budgets w_i:
#   dcs     area under the fitted demand                                    (dcs_numeric)
#   willig  dCS_i + eta_i dCS_i^2 / (2 w_i), eta_i the form's income elasticity at (2, w_i)   (willig_cv)
#   cv, ev  exact CV and EV of the fitted demand (Vartia ODE; EV along the reversed path)   (vartia_cv)
#   cv_rep  the representative consumer's CV: the fitted demand at mean income treated as one consumer
#   cv_out  the CV of raising p1 from 2 to 12 USD (Ruling P7: far outside the 2 to 4 USD support)
# The headline numbers use the income-term form, the only one whose functional form contains the Stone-Geary
# conditional mean when beta is independent of income:
#   w1$dcs, w1$willig, w1$vartia, w1$cv_rep (W1), w1$cv_rep_g (W1g).
# spread_in_range = (max - min) / mean of the three forms' CV from 2 to 4 USD; spread_outside the same from
# 2 to 12 USD. willig_check (Ruling P13): the card's Willig formula against the exact CV of W0 and of every
# W1 / W1g subscriber (true dCS_i and eta_i), with the opposite sign and with coefficient 1 for contrast.
# Reads truth01.rds, plover_W1.rds, plover_W1g.rds; writes res_welfare.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("welfare_fns.R")
tr <- readRDS("truth01.rds")
P1 <- tr$params$W1
W1 <- readRDS("plover_W1.rds"); W1g <- readRDS("plover_W1g.rds")
pgrid <- tr$W1$demand_curve$p1

analyse_variant <- function(d, P) {
  fm <- fit_demand_forms(d)
  w <- d$w
  forms <- lapply(fm, function(f) {
    dcs_i <- dcs_numeric(f$x1_fun, P$p0, P$p1, w)
    eta_i <- f$eta_fun(P$p0, w)
    list(coef = f$coef, r2 = f$r2, smear = if (is.null(f$smear)) NA_real_ else f$smear,
         dcs = mean(dcs_i), eta = mean(eta_i),
         willig = mean(willig_cv(dcs_i, eta_i, w)), willig_ev = mean(willig_cv(dcs_i, -eta_i, w)),
         cv = mean(vartia_cv(f$x1_fun, P$p0, P$p1, w)), ev = mean(-vartia_cv(f$x1_fun, P$p1, P$p0, w)),
         cv_rep = vartia_cv(f$x1_fun, P$p0, P$p1, mean(w)),
         dcs_out = mean(dcs_numeric(f$x1_fun, P$p0, P$p_out, w)),
         cv_out = mean(vartia_cv(f$x1_fun, P$p0, P$p_out, w)),
         x1_at_p0_mean = mean(f$x1_fun(P$p0, w)), x1_at_pout_mean = mean(f$x1_fun(P$p_out, w)),
         demand_curve = vapply(pgrid, function(p) mean(f$x1_fun(p, w)), numeric(1)))
  })
  cl <- forms$linear$coef; wbar <- mean(w)
  forms$linear$choke_at_wbar <- -(cl[["a"]] + cl[["c_w"]] * wbar) / (cl[["b_p1"]] + cl[["d_p1w"]] * wbar)
  forms$loglinear$price_elasticity <- forms$loglinear$coef[["b_logp1"]]
  arms <- data.frame(p1 = sort(unique(d$p1)), x1_mean = as.numeric(tapply(d$x1, d$p1, mean)),
                     n = as.integer(table(d$p1)), w_mean = as.numeric(tapply(d$w, d$p1, mean)))
  cv <- vapply(forms, `[[`, numeric(1), "cv"); cv_out <- vapply(forms, `[[`, numeric(1), "cv_out")
  list(forms = forms, arms = arms, wbar = wbar, cv = cv, cv_out = cv_out,
       spread_in_range = diff(range(cv)) / mean(cv), spread_outside = diff(range(cv_out)) / mean(cv_out))
}

res_W1 <- analyse_variant(W1, P1)
res_W1g <- analyse_variant(W1g, P1)

## ---- Ruling P13: the Willig formula against exact CVs --------------------------------------------
W0 <- tr$W0
removed <- function(approx, exact, base) 1 - abs(approx - exact) / abs(base - exact)   # share of the base error removed
willig_W0 <- list(cv = W0$cv, ev = W0$ev, dcs = W0$dcs, eta = W0$eta,
                  willig = willig_cv(W0$dcs, W0$eta, W0$w),
                  willig_minus = willig_cv(W0$dcs, -W0$eta, W0$w),               # opposite sign: the first-order EV
                  willig_coef1 = W0$dcs + W0$eta * W0$dcs^2 / W0$w)               # coefficient 1 instead of 1/2
willig_W0$share_removed <- removed(willig_W0$willig, W0$cv, W0$dcs)
willig_W0$share_removed_ev <- removed(willig_W0$willig_minus, W0$ev, W0$dcs)
willig_truth <- function(t, d) {
  wc1 <- t$dcs_i + t$eta_i * t$dcs_i^2 / d$w
  wm <- willig_cv(t$dcs_i, -t$eta_i, d$w)
  list(willig_mean = mean(t$willig_i), willig_coef1_mean = mean(wc1), willig_minus_mean = mean(wm),
       share_removed = removed(mean(t$willig_i), t$cv_mean, t$dcs_mean),
       share_removed_ev = removed(mean(wm), t$ev_mean, t$dcs_mean),
       mean_abs_err_dcs = mean(abs(t$dcs_i - t$cv_i)), mean_abs_err_willig = mean(abs(t$willig_i - t$cv_i)),
       max_rel_err_willig = max(abs(t$willig_i - t$cv_i) / t$cv_i),
       min_share_removed_i = min(1 - abs(t$willig_i - t$cv_i) / abs(t$dcs_i - t$cv_i)),
       all_closer = all(abs(t$willig_i - t$cv_i) < abs(t$dcs_i - t$cv_i)))
}
willig_check <- list(W0 = willig_W0, W1 = willig_truth(tr$W1, W1), W1g = willig_truth(tr$W1g, W1g))

## ---- headline numbers (income-term form) ------------------------------------------------------------
inc <- res_W1$forms$income
res <- list(form_headline = "income",
            dcs = inc$dcs, willig = inc$willig, vartia = inc$cv, ev = inc$ev, willig_ev = inc$willig_ev,
            cv_rep = inc$cv_rep, cv_rep_g = res_W1g$forms$income$cv_rep,
            spread_in_range = res_W1$spread_in_range, spread_outside = res_W1$spread_outside,
            W1 = res_W1, W1g = res_W1g, willig_check = willig_check, pgrid = pgrid)
saveRDS(res, "res_welfare.rds")

cat(sprintf("W1 truth: mean CV %.4f, mean dCS %.4f, mean EV %.4f, CV 2 -> 12 %.4f\n",
            tr$W1$cv_mean, tr$W1$dcs_mean, tr$W1$ev_mean, tr$W1$cv12_mean))
for (v in c("W1", "W1g")) {
  r <- res[[v]]
  cat(sprintf("%s:\n", v))
  for (nm in names(r$forms)) {
    f <- r$forms[[nm]]
    cat(sprintf("  %-9s dCS %.4f  Willig %.4f  CV %.4f  EV %.4f  CV^rep %.4f  CV 2 -> 12 %.4f  eta %.4f\n",
                nm, f$dcs, f$willig, f$cv, f$ev, f$cv_rep, f$cv_out, f$eta))
  }
  cat(sprintf("  spread in range %.4f, outside %.4f\n", r$spread_in_range, r$spread_outside))
}
cat(sprintf("Headline (income-term form): dCS %.4f, Willig %.4f, Vartia %.4f vs true mean CV %.4f;\n",
            res$dcs, res$willig, res$vartia, tr$W1$cv_mean))
cat(sprintf("  |Willig - CV| / |dCS - CV| = %.4f; CV^rep W1 %.4f (true mean %.4f), W1g %.4f (true mean %.4f)\n",
            abs(res$willig - tr$W1$cv_mean) / abs(res$dcs - tr$W1$cv_mean), res$cv_rep, tr$W1$cv_mean,
            res$cv_rep_g, tr$W1g$cv_mean))
cat(sprintf("Willig check (P13): W0 removes %.4f of the CS error; W1 truth removes %.4f; W1g %.4f\n",
            willig_W0$share_removed, willig_check$W1$share_removed, willig_check$W1g$share_removed))
cat("Saved res_welfare.rds\n")
