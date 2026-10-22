# ------------------------------------------------------------------
# 06_figures.R -- chapter 6 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R
# and 20_did_beyond_binary/05_figures.R.
#   fig_06_censoring.svg   blocks C1 and C3: the fee coefficient beta_fee (PH scale) from the naive
#                          regressions on censored tenure, the censored-data Weibull MLE, the grouped
#                          (monthly cloglog) MLE and the infeasible regression on latent tenure, with
#                          95% intervals, against the truth
#   fig_06_qmle.svg        block C2: for a 60-dollar subscriber without referral, the true hazard and
#                          survival of the log-logistic mixture (and the two segments' hazards), the
#                          Weibull pseudo-true and QMLE curves, and the Kaplan-Meier estimate
#   fig_06_coverage.svg    blocks C2 and C4: Monte Carlo coverage of nominal 95% intervals by standard error
# Reads truth06.rds, sorrel_C2.rds, res_cens.rds, res_qmle.rds, res_grouped.rds, res_twostep.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("ml_fns.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(survival) })
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

tr <- readRDS("truth06.rds")
C2 <- readRDS("sorrel_C2.rds")
c1 <- readRDS("res_cens.rds")
c2 <- readRDS("res_qmle.rds")
c3 <- readRDS("res_grouped.rds")
c4 <- readRDS("res_twostep.rds")

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"; col_grey <- "grey55"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name, height = 4.8) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = height, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}

## ---- Figure (a): blocks C1 and C3, beta_fee under censoring ------------------------------------------
lab_a <- c(obs = "OLS, censored spells as exits",
           unc = "OLS, uncensored spells only",
           mle = "Weibull MLE, censored-data likelihood",
           grp = "Grouped-duration MLE (months 1 to 24)",
           lat = "OLS on latent tenure t* (infeasible)")
grp_a <- c(obs = "Ignores censoring", unc = "Ignores censoring", mle = "Uses the likelihood", grp = "Uses the likelihood",
           lat = "Infeasible benchmark")
fa <- data.frame(key = names(lab_a),
                 est = c(c1$ols_observed$est, c1$ols_uncensored$est, c1$mle$est, c3$fee$est, c1$ols_latent$est),
                 lo = c(c1$ols_observed$ci[1], c1$ols_uncensored$ci[1], c1$mle$ci[1], c3$fee$ci[1], c1$ols_latent$ci[1]),
                 hi = c(c1$ols_observed$ci[2], c1$ols_uncensored$ci[2], c1$mle$ci[2], c3$fee$ci[2], c1$ols_latent$ci[2]))
fa$label <- factor(lab_a[fa$key], levels = rev(unname(lab_a)))
fa$group <- factor(grp_a[fa$key], levels = unique(grp_a))
p1 <- ggplot(fa, aes(x = est, y = label, colour = group, shape = group)) +
  geom_vline(xintercept = tr$C1$beta, linewidth = 0.5, colour = col_true, linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(size = 2.4) +
  annotate("text", x = tr$C1$beta, y = 5.45, label = sprintf("True value %.2f ", tr$C1$beta), hjust = 1.05, size = 3.2, colour = col_true) +
  scale_colour_manual(values = setNames(c(col_orange, col_blue, col_grey), levels(fa$group))) +
  scale_shape_manual(values = setNames(c(17, 16, 15), levels(fa$group))) +
  scale_x_continuous(breaks = seq(0.006, 0.026, 0.004)) +
  coord_cartesian(ylim = c(0.8, 5.6)) +
  labs(x = "Fee coefficient beta_fee in the hazard (log hazard per USD), 95% CI", y = NULL,
       colour = NULL, shape = NULL) +
  base_theme
save_fig(p1, "fig_06_censoring.svg", height = 4.2)

## ---- Figure (b): block C2, the misspecified Weibull ------------------------------------------------------
P2 <- tr$params$C2
ts <- exp(seq(log(0.1), log(36), length.out = 400))
ty <- list(pi = plogis(P2$g0), a_d = P2$med_deal, a_r = P2$med_reg, p_d = P2$p_deal, p_r = P2$p_reg)   # fee 60, not referred
hz_ll <- function(t, a, p) (p / t) * (t / a)^p / (1 + (t / a)^p)
wb_h <- function(par, t) exp(par[1] + 60 * par[2]) * exp(par[4]) * t^(exp(par[4]) - 1)
wb_S <- function(par, t) exp(-exp(par[1] + 60 * par[2]) * t^exp(par[4]))
panel_b <- c(h = "Cancellation hazard (per month, log scale)", S = "Survival: share still subscribed")
lab_b <- c(true = "True mixture", deal = "Deal-seekers' segment", reg = "Regulars' segment",
           pseudo = "Weibull, pseudo-true value", qmle = "Weibull QMLE (sample)", km = "Kaplan-Meier (sample)")
sel <- C2$fee_usd == 60 & C2$referred == 0
km <- survfit(Surv(t_months, delta) ~ 1, data = C2[sel, ])
cb <- rbind(
  data.frame(panel = panel_b[["h"]], t = ts, value = c2_hazard(ts, 60, 0, P2), series = lab_b[["true"]]),
  data.frame(panel = panel_b[["h"]], t = ts, value = hz_ll(ts, ty$a_d, ty$p_d), series = lab_b[["deal"]]),
  data.frame(panel = panel_b[["h"]], t = ts, value = hz_ll(ts, ty$a_r, ty$p_r), series = lab_b[["reg"]]),
  data.frame(panel = panel_b[["h"]], t = ts, value = wb_h(tr$C2$pseudo_vec, ts), series = lab_b[["pseudo"]]),
  data.frame(panel = panel_b[["h"]], t = ts, value = wb_h(c2$qmle$par, ts), series = lab_b[["qmle"]]),
  data.frame(panel = panel_b[["S"]], t = ts, value = c2_surv(ts, 60, 0, P2), series = lab_b[["true"]]),
  data.frame(panel = panel_b[["S"]], t = ts, value = wb_S(tr$C2$pseudo_vec, ts), series = lab_b[["pseudo"]]),
  data.frame(panel = panel_b[["S"]], t = ts, value = wb_S(c2$qmle$par, ts), series = lab_b[["qmle"]]))
cb <- cb[cb$panel != panel_b[["h"]] | (cb$value > 0.004 & cb$value < 3), ]
kmd <- data.frame(panel = panel_b[["S"]], t = c(0, km$time), value = c(1, km$surv), series = lab_b[["km"]])
cb$series <- factor(cb$series, levels = unname(lab_b)); kmd$series <- factor(kmd$series, levels = unname(lab_b))
sc_col <- scale_colour_manual(values = setNames(c(col_true, col_orange, col_green, col_blue, col_blue, col_grey), unname(lab_b)),
                              drop = FALSE, limits = unname(lab_b))
sc_lty <- scale_linetype_manual(values = setNames(c("solid", "dotted", "dotdash", "solid", "dashed", "solid"), unname(lab_b)),
                                drop = FALSE, limits = unname(lab_b))
xlab_b <- "Tenure (months), 60-dollar plan, no referral"
hd <- rbind(cb[cb$panel == panel_b[["h"]], ],                          # an empty Kaplan-Meier row: its legend key only
            data.frame(panel = panel_b[["h"]], t = 1, value = NA_real_, series = factor(lab_b[["km"]], levels = unname(lab_b))))
p2h <- ggplot(hd, aes(x = t, y = value, colour = series, linetype = series)) +
  geom_line(linewidth = 0.75, na.rm = TRUE) +
  scale_y_log10(breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2), labels = c("0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1", "2")) +
  scale_x_continuous(breaks = seq(0, 36, 6)) + sc_col + sc_lty +
  labs(title = panel_b[["h"]], x = xlab_b, y = NULL, colour = NULL, linetype = NULL) +
  base_theme + theme(plot.title = element_text(size = 10))
p2s <- ggplot() +
  geom_step(data = kmd, aes(x = t, y = value, colour = series, linetype = series), linewidth = 0.4) +
  geom_line(data = cb[cb$panel == panel_b[["S"]], ], aes(x = t, y = value, colour = series, linetype = series), linewidth = 0.75) +
  scale_x_continuous(breaks = seq(0, 36, 6)) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + sc_col + sc_lty +
  labs(title = panel_b[["S"]], x = xlab_b, y = NULL, colour = NULL, linetype = NULL) +
  base_theme + theme(plot.title = element_text(size = 10)) + guides(colour = "none", linetype = "none")
# one legend, drawn by the hazard panel with every series (drop = FALSE), placed under both panels
p2h <- p2h + guides(colour = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE))
p2 <- (p2h | p2s) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(p2, "fig_06_qmle.svg")

## ---- Figure (c): blocks C2 and C4, Monte Carlo coverage --------------------------------------------------
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
panel_c <- c(C2 = sprintf("Block C2: Weibull QMLE of the shape kappa_* (%s samples of %s)", fmt(c2$mc$R), fmt(c2$mc$n)),
             C4 = sprintf("Block C4: two-step MLE of beta_u (%s samples of %s)", fmt(c4$mc$R), fmt(c4$mc$n)))
fc <- rbind(
  data.frame(panel = panel_c[["C2"]], se = "Model-based (Hessian)", cover = c2$mc$cover_model, kind = "Assumes the model or ignores step 1"),
  data.frame(panel = panel_c[["C2"]], se = "Outer product of scores (OPG)", cover = c2$mc$cover_opg, kind = "Assumes the model or ignores step 1"),
  data.frame(panel = panel_c[["C2"]], se = "Sandwich", cover = c2$mc$cover_sandwich, kind = "Valid"),
  data.frame(panel = panel_c[["C4"]], se = "Naive (Hessian, index as data)", cover = c4$mc$cover_naive, kind = "Assumes the model or ignores step 1"),
  data.frame(panel = panel_c[["C4"]], se = "Naive sandwich (index as data)", cover = c4$mc$cover_naive_sand, kind = "Assumes the model or ignores step 1"),
  data.frame(panel = panel_c[["C4"]], se = "Two-step corrected", cover = c4$mc$cover_corrected, kind = "Valid"),
  data.frame(panel = panel_c[["C4"]], se = "Oracle (true index)", cover = c4$mc$cover_oracle, kind = "Infeasible benchmark"))
fc$panel <- factor(fc$panel, levels = unname(panel_c))
fc$se <- factor(fc$se, levels = rev(unique(fc$se)))
fc$kind <- factor(fc$kind, levels = c("Assumes the model or ignores step 1", "Valid", "Infeasible benchmark"))
# Monte Carlo band around 0.95: +/- qnorm(0.975) sqrt(0.95 0.05 / R) for each block's R
bands <- data.frame(panel = factor(unname(panel_c), levels = unname(panel_c)),
                    half = qnorm(0.975) * sqrt(0.95 * 0.05 / c(c2$mc$R, c4$mc$R)))
p3 <- ggplot(fc, aes(x = cover, y = se, colour = kind, shape = kind)) +
  geom_rect(data = bands, aes(xmin = 0.95 - half, xmax = 0.95 + half, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = col_true, alpha = 0.08) +
  geom_vline(xintercept = 0.95, linewidth = 0.4, colour = col_true) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", cover)), hjust = 1.35, size = 3.1, show.legend = FALSE) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = setNames(c(col_orange, col_blue, col_grey), levels(fc$kind))) +
  scale_shape_manual(values = setNames(c(17, 16, 15), levels(fc$kind))) +
  scale_x_continuous(limits = c(0.45, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = "Coverage of the nominal 95% interval", y = NULL, colour = NULL, shape = NULL) +
  base_theme + theme(strip.text = element_text(face = "plain", hjust = 0))
save_fig(p3, "fig_06_coverage.svg", height = 5)
