# ------------------------------------------------------------------
# 05_figures.R -- Chapter 20 figures (SVG), same style as 19_did_event_studies/05_figures.R.
#   fig_20_fuzzy.svg      block A: Wald-DID, Wald-TC, Wald-CIC with 95% city-bootstrap CIs
#                         against the true LATE, with and without a stable control take-up rate
#   fig_20_dose.svg       block B: estimated ATT(d|d) and the slope of the fit (contdid) against truth
#   fig_20_switching.svg  block C: DID_l and placebos (did_multiplegt_dyn) against the true effects
#   fig_20_twfe.svg       block C: static TWFE against the dCDH average total effect, per unit of dose
# Event-time convention as in chapter 19: DID_l sits at e = l - 1 (weeks since the first
# switch), placebo l at e = -1 - l, and the reference week e = -1 is normalised to 0.
# ------------------------------------------------------------------
library(data.table)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
fig_dir  <- normalizePath(file.path(code_dir, "..", "..", "assets", "fig"))

tr <- readRDS(file.path(code_dir, "truth20.rds"))
f  <- readRDS(file.path(code_dir, "res_fuzzy.rds"))
dz <- readRDS(file.path(code_dir, "res_dose.rds"))
sw <- readRDS(file.path(code_dir, "res_switch.rds"))

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}

## ---- Figure (a): block A --------------------------------------------
lab_est <- c(did = "Wald-DID", tc = "Wald-TC", cic = "Wald-CIC")
lab_var <- c("Stable control take-up (block A)", "Rising control take-up (variant A')")
fz <- rbindlist(lapply(names(lab_est), function(k) rbind(
  data.table(estimator = lab_est[[k]], variant = lab_var[1], est = f[[k]]$est, lo = f[[k]]$lo, hi = f[[k]]$hi),
  data.table(estimator = lab_est[[k]], variant = lab_var[2], est = f$rising[[k]]$est, lo = f$rising[[k]]$lo, hi = f$rising[[k]]$hi))))
fz[, estimator := factor(estimator, levels = rev(unname(lab_est)))]
fz[, variant := factor(variant, levels = lab_var)]
pd <- position_dodge(width = 0.5)
p1 <- ggplot(fz, aes(x = est, y = estimator, colour = variant, shape = variant)) +
  geom_vline(xintercept = tr$A$late, linewidth = 0.5, colour = col_true, linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", position = pd, width = 0, linewidth = 0.45) +
  geom_point(position = pd, size = 2.2) +
  annotate("text", x = tr$A$late, y = 3.52, label = "True LATE", hjust = -0.1, vjust = 0, size = 3.2, colour = col_true) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), lab_var)) +
  scale_shape_manual(values = setNames(c(16, 17), lab_var)) +
  scale_x_continuous(breaks = seq(10, 24, 2)) +   # label the Wald-DID region (about 12) as well as 16 and 22
  labs(x = "Estimate of the LATE of switchers (weekly orders, 95% city-bootstrap CI)", y = NULL,
       colour = NULL, shape = NULL) +
  base_theme
save_fig(p1, "fig_20_fuzzy.svg")

## ---- Figure (b): block B ---------------------------------------------
cv <- as.data.table(dz$cont$curve); ac <- as.data.table(dz$cont$acrt_curve)
lab_panel <- c("ATT(d | d): effect on weekly orders", "Slope of the ATT(d | d) curve: orders per USD")
lab_s <- c("contdid estimate (95% CI)", "True: what parallel trends identifies", "True causal response ACRT(d | d)")
b_est <- rbind(
  data.table(panel = lab_panel[1], dose = cv$dose, est = cv$est, lo = cv$lo, hi = cv$hi),
  data.table(panel = lab_panel[2], dose = ac$dose, est = ac$est, lo = ac$lo, hi = ac$hi))
b_tru <- rbind(
  data.table(panel = lab_panel[1], dose = cv$dose, series = lab_s[2], value = cv$truth),
  data.table(panel = lab_panel[2], dose = ac$dose, series = lab_s[2], value = ac$truth_slope),
  data.table(panel = lab_panel[2], dose = ac$dose, series = lab_s[3], value = ac$truth_acrt))
b_est[, panel := factor(panel, levels = lab_panel)]; b_tru[, panel := factor(panel, levels = lab_panel)]
b_tru[, series := factor(series, levels = lab_s)]
b_est[, series := factor(lab_s[1], levels = lab_s)]
p2 <- ggplot() +
  geom_ribbon(data = b_est, aes(x = dose, ymin = lo, ymax = hi, fill = series), alpha = 0.18) +
  geom_line(data = b_est, aes(x = dose, y = est, colour = series, linetype = series), linewidth = 0.7) +
  geom_line(data = b_tru, aes(x = dose, y = value, colour = series, linetype = series), linewidth = 0.7) +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = setNames(c(col_blue, col_true, col_orange), lab_s), drop = FALSE) +
  scale_fill_manual(values = setNames(c(col_blue, col_true, col_orange), lab_s), drop = FALSE, guide = "none") +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed"), lab_s), drop = FALSE) +
  labs(x = "Dose d (USD per order)", y = NULL, colour = NULL, linetype = NULL) +
  base_theme + theme(strip.text = element_text(face = "plain"))
save_fig(p2, "fig_20_dose.svg")

## ---- Figure (c): block C event study ------------------------------------
E <- as.data.table(sw$dyn$effects)[, .(e = l - 1L, estimate = est, lo, hi)]
Pl <- as.data.table(sw$dyn$placebos)[, .(e = -1L - l, estimate = est, lo, hi)]
ref <- data.table(e = -1L, estimate = 0, lo = NA_real_, hi = NA_real_)
es <- rbind(Pl, ref, E)[order(e)]
es[, series := "did_multiplegt_dyn: DID_l and placebos"]
truth_es <- rbind(data.table(e = c(-5:-2, -1L), estimate = 0),
                  data.table(e = 0:7, estimate = tr$C$delta))[, `:=`(lo = NA_real_, hi = NA_real_, series = "True effect")]
cs <- rbind(es, truth_es)
cs[, series := factor(series, levels = c("True effect", "did_multiplegt_dyn: DID_l and placebos"))]
pd3 <- position_dodge(width = 0.4)
p3 <- ggplot(cs, aes(x = e, y = estimate, colour = series, shape = series)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = -0.5, linewidth = 0.3, colour = "grey40", linetype = "dashed") +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = pd3, width = 0, linewidth = 0.4, na.rm = TRUE) +
  geom_line(data = cs[series == "True effect"], aes(group = series), linewidth = 0.5, alpha = 0.7) +
  geom_point(position = pd3, size = 1.8, na.rm = TRUE) +
  scale_colour_manual(values = c("True effect" = col_true, "did_multiplegt_dyn: DID_l and placebos" = col_blue)) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(breaks = -5:7) +
  labs(x = "Event time (weeks since the first switch)", y = "Effect on weekly orders", colour = NULL, shape = NULL) +
  base_theme
save_fig(p3, "fig_20_switching.svg")

## ---- Figure (d): block C, per unit of dose ---------------------------------
tw_lab <- c("Static TWFE (feols)", "Static TWFE on noise-free effects", "Average total effect (did_multiplegt_dyn)")
st <- data.table(
  estimator = factor(tw_lab, levels = rev(tw_lab)),
  est = c(sw$twfe$est, sw$twfe$oracle_est, sw$dyn$avg_total$est),
  lo  = c(sw$twfe$lo, NA, sw$dyn$avg_total$lo),
  hi  = c(sw$twfe$hi, NA, sw$dyn$avg_total$hi))
refs <- data.table(x = c(tr$C$avg_total_effect, tr$C$avg_effect),
                   truth = factor(c("True average total effect", "True dose-weighted effect per unit of dose"),
                                  levels = c("True average total effect", "True dose-weighted effect per unit of dose")))
p4 <- ggplot() +
  geom_vline(data = refs, aes(xintercept = x, linetype = truth), linewidth = 0.5, colour = col_true) +
  geom_errorbar(data = st, aes(y = estimator, xmin = lo, xmax = hi, colour = estimator), orientation = "y", width = 0,
                linewidth = 0.45, na.rm = TRUE) +
  geom_point(data = st, aes(y = estimator, x = est, shape = estimator, colour = estimator), size = 2.4) +
  scale_colour_manual(values = setNames(c(col_orange, col_orange, col_blue), tw_lab), guide = "none") +
  scale_shape_manual(values = setNames(c(17, 2, 16), tw_lab), guide = "none") +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  labs(x = "Effect per unit of bonus tier (weekly orders; 95% CI where available)", y = NULL, linetype = NULL) +
  base_theme
save_fig(p4, "fig_20_twfe.svg")
