# ------------------------------------------------------------------
# 09_figures.R -- chapter 11 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R
# and 10_linear_panel_fe_re/09_figures.R (replaces the legacy 04_figures.R).
#   fig_11_estimators.svg     blocks D1 (rho = 0.70) and D1p (rho = 0.90), the inherited samples: pooled
#                             OLS, FE, the two FE bias corrections, Anderson-Hsiao, Arellano-Bond one- and
#                             two-step, Blundell-Bond with 95% intervals, the probability limits of pooled
#                             OLS and FE, and the truth (redrawn at 9 in wide; the legacy figure showed
#                             four bars for D1 only)
#   fig_11_nickell.svg        block D1T (the inherited Monte Carlo) and D5: FE against the months per merchant
#                             T_obs, with the exact Nickell limit, its first factor, the order-1/T term (the
#                             "nickell_approx" row of the inherited table; in the figure it is labelled as the
#                             order-1/T bias that Hahn-Kuersteiner remove: readers never see "legacy"),
#                             Arellano-Bond, and the two bias corrections (the inherited right panel,
#                             rho -> 1, is now the D1p panel of fig_11_estimators)
#   fig_11_proliferation.svg  block D2: one-step system GMM and its Hansen J p-value against the number of
#                             instruments r, every lag against collapsed, with OLS in levels, FE and the truth;
#                             the region r >= N of the J panel is shaded: there J = N exactly and the p-value
#                             depends on the degrees of freedom alone (the J test has no information)
# Reads truth11.rds, res_dynamic.rds, res_extra.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
source("dyn_fns.R")
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

tr <- readRDS("truth11.rds"); dy <- readRDS("res_dynamic.rds"); x <- readRDS("res_extra.rds")

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name, width = 8, height = 4.8) {
  ggsave(file.path(fig_dir, name), p, width = width, height = height, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}

## ---- Figure (a): blocks D1 and D1p, the estimators ----------------------------------------------
lab_e <- c(pool = "Pooled OLS", fe = "Fixed effects (within)", hk = "FE + Hahn-Kuersteiner correction",
           spj = "FE, half-panel jackknife (point only)", ah = "Anderson-Hsiao (instrument Y[t-2])",
           ab1 = "Arellano-Bond, one-step", ab2 = "Arellano-Bond, two-step (Windmeijer SE)",
           bb = "Blundell-Bond, two-step (Windmeijer SE)")
panel_a <- c(D1 = "Block D1: rho = 0.70 (500 merchants, 8 months)", D1p = "Block D1p: rho = 0.90 (same design)")
est <- do.call(rbind, lapply(names(panel_a), function(b) {
  z <- x[[b]]
  data.frame(panel = panel_a[[b]], key = names(lab_e),
             est = c(z$pool$est, z$fe$est, z$hk, z$spj, z$ah$est, z$ab_onestep$est, z$ab$est, z$bb$est),
             lo = c(z$pool$ci[1], z$fe$ci[1], z$hk_ci[1], NA, z$ah$ci[1], z$ab_onestep$ci[1], z$ab$ci[1], z$bb$ci[1]),
             hi = c(z$pool$ci[2], z$fe$ci[2], z$hk_ci[2], NA, z$ah$ci[2], z$ab_onestep$ci[2], z$ab$ci[2], z$bb$ci[2]))
}))
plim <- do.call(rbind, lapply(names(panel_a), function(b) {
  data.frame(panel = panel_a[[b]], key = c("pool", "fe"), value = c(tr[[b]]$pooled_plim, tr[[b]]$fe_plim))
}))
truth_a <- data.frame(panel = unname(panel_a), value = c(tr$rho, tr$rho_persist))
for (nm in c("est", "plim", "truth_a")) {
  z <- get(nm); z$panel <- factor(z$panel, levels = unname(panel_a))
  if (!is.null(z$key)) z$estimator <- factor(lab_e[z$key], levels = rev(unname(lab_e)))
  assign(nm, z)
}
lab_s <- c(e = "Estimate, 95% CI (clustered by merchant; GMM: robust or Windmeijer)",
           p = "Probability limit as the number of merchants grows (pooled OLS; FE: exact Nickell formula)")
est$series <- lab_s[["e"]]; plim$series <- lab_s[["p"]]
p1 <- ggplot() +
  geom_vline(data = truth_a, aes(xintercept = value), linewidth = 0.5, colour = col_true, linetype = "dashed") +
  geom_errorbar(data = est[!is.na(est$lo), ], aes(y = estimator, xmin = lo, xmax = hi, colour = series), orientation = "y",
                width = 0, linewidth = 0.5) +
  geom_point(data = est, aes(y = estimator, x = est, colour = series, shape = series), size = 2.2) +
  geom_point(data = plim, aes(y = estimator, x = value, colour = series, shape = series), size = 2.6, stroke = 0.9) +
  geom_text(data = est, aes(y = estimator, x = Inf, label = sprintf("%.3f", est)), colour = col_blue, hjust = 1.1, size = 2.9) +
  geom_text(data = truth_a, aes(x = value, y = 8.62, label = sprintf("True rho = %.2f", value)), hjust = -0.06, size = 3,
            colour = col_true) +
  facet_wrap(~panel, scales = "free_x") +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), unname(lab_s))) +
  scale_shape_manual(values = setNames(c(16, 5), unname(lab_s))) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.2))) +          # room on the right for the column of estimates
  scale_y_discrete(limits = rev(unname(lab_e)), expand = expansion(add = c(0.6, 1.1))) +
  labs(x = "Estimate of the persistence rho of monthly GMV (the number at the right of each row is the estimate)",
       y = NULL, colour = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  base_theme + theme(strip.text = element_text(face = "plain"))
save_fig(p1, "fig_11_estimators.svg", width = 10, height = 5.4)

## ---- Figure (b): blocks D1T and D5, Nickell bias and its corrections against T_obs -------------------
rho <- tr$rho
Tg <- seq(3, 33, by = 0.25)                                    # months per merchant (T = Tg - 1 regression periods)
lab_c <- c(ex = "FE probability limit: exact Nickell formula",
           ff = "First factor only: rho - (1 + rho)/(T - 1), T = T_obs - 1",
           o1 = "Order-1/T bias removed by Hahn-Kuersteiner: rho - (1 + rho)/T")
curves <- rbind(data.frame(T_obs = Tg, value = rho + nickell_bias(rho, Tg - 1), series = lab_c[["ex"]]),
                data.frame(T_obs = Tg[Tg >= 4], value = rho + nickell_first_factor(rho, Tg[Tg >= 4] - 1), series = lab_c[["ff"]]),
                data.frame(T_obs = Tg, value = rho - (1 + rho) / (Tg - 1), series = lab_c[["o1"]]))
curves$series <- factor(curves$series, levels = unname(lab_c))
lab_p <- c(fe = "FE, Monte Carlo mean (40 samples of 300 merchants)", ab = "Arellano-Bond, same Monte Carlo",
           hk = "FE + Hahn-Kuersteiner, Monte Carlo mean (1,000 samples of 500)", spj = "Half-panel jackknife, same Monte Carlo")
Tl <- as.numeric(colnames(dy$nick)); m5 <- x$D5
pts <- rbind(data.frame(T_obs = Tl, value = dy$nick["fe", ], series = lab_p[["fe"]]),
             data.frame(T_obs = Tl, value = dy$nick["ab", ], series = lab_p[["ab"]]),
             data.frame(T_obs = m5$T_obs, value = m5$hk, series = lab_p[["hk"]]),
             data.frame(T_obs = m5$T_obs, value = m5$spj, series = lab_p[["spj"]]))
pts$series <- factor(pts$series, levels = unname(lab_p))
p2 <- ggplot() +
  geom_hline(yintercept = rho, linewidth = 0.5, colour = col_true, linetype = "dashed") +
  annotate("text", x = 24, y = rho + 0.03, label = "True rho = 0.70", hjust = 0.5, size = 3.1, colour = col_true) +
  geom_line(data = curves, aes(x = T_obs, y = value, linetype = series), colour = col_blue, linewidth = 0.6) +
  geom_point(data = pts, aes(x = T_obs, y = value, colour = series, shape = series), size = 2.4) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "dotted"), unname(lab_c))) +
  scale_colour_manual(values = setNames(c(col_blue, col_true, col_orange, col_green), unname(lab_p))) +
  scale_shape_manual(values = setNames(c(16, 17, 15, 18), unname(lab_p))) +
  scale_x_continuous(breaks = c(4, 8, 16, 32), limits = c(3, 33)) +
  coord_cartesian(ylim = c(-0.05, 0.8)) +
  labs(x = "Months observed per merchant, T_obs (regression periods T = T_obs - 1)",
       y = "Estimate of rho (mean over samples) or limit", linetype = NULL, colour = NULL, shape = NULL) +
  guides(linetype = guide_legend(nrow = 3, order = 1), colour = guide_legend(nrow = 2, order = 2),
         shape = guide_legend(nrow = 2, order = 2)) +
  base_theme + theme(legend.box = "vertical")
save_fig(p2, "fig_11_nickell.svg", width = 8.5, height = 5.6)

## ---- Figure (c): block D2, instrument proliferation ----------------------------------------------------
D2 <- x$D2; pa <- D2$path
panel_c <- c(est = "One-step system GMM estimate of rho (95% CI)", J = "Hansen J p-value")
lab_k <- c(unc = "Every lag of each month (uncollapsed)", col = "Collapsed (one instrument per lag)")
pa$series <- factor(ifelse(pa$collapse, lab_k[["col"]], lab_k[["unc"]]), levels = unname(lab_k))
long <- rbind(data.frame(panel = panel_c[["est"]], r = pa$r, value = pa$est, lo = pa$lo, hi = pa$hi, series = pa$series),
              data.frame(panel = panel_c[["J"]], r = pa$r, value = pa$J_p, lo = NA, hi = NA, series = pa$series))
long$panel <- factor(long$panel, levels = unname(panel_c))
ref <- data.frame(panel = factor(panel_c[["est"]], levels = unname(panel_c)),
                  value = c(tr$D2$rho, D2$ols_transformed, D2$fe),
                  label = c(sprintf("True rho = %.2f", tr$D2$rho), sprintf("OLS in levels %.3f", D2$ols_transformed),
                            sprintf("FE %.3f", D2$fe)))
refJ <- data.frame(panel = factor(panel_c[["J"]], levels = unname(panel_c)), value = 0.05, label = "p = 0.05")
nline <- data.frame(panel = factor(unname(panel_c), levels = unname(panel_c)), r = D2$N)
jzone <- data.frame(panel = factor(panel_c[["J"]], levels = unname(panel_c)), xmin = D2$N, xmax = 140,
                    label = "J \u2261 N for r \u2265 N:\nthe p-value depends on\nthe degrees of freedom only")
p3 <- ggplot(long, aes(x = r, y = value, colour = series, shape = series)) +
  geom_vline(data = nline, aes(xintercept = r), linewidth = 0.4, colour = "grey55", linetype = "dotdash") +
  geom_hline(data = ref, aes(yintercept = value), linewidth = 0.4, colour = col_true, linetype = c("dashed", "dotted", "dotted")) +
  geom_text(data = ref, aes(x = 140, y = value, label = label), inherit.aes = FALSE, hjust = 1, vjust = -0.4, size = 2.9,
            colour = col_true) +
  geom_hline(data = refJ, aes(yintercept = value), linewidth = 0.4, colour = col_true, linetype = "dotted") +
  geom_text(data = refJ, aes(x = 140, y = value, label = label), inherit.aes = FALSE, hjust = 1, vjust = -0.4, size = 2.9,
            colour = col_true) +
  geom_text(data = nline[1, ], aes(x = r, y = 0.465, label = sprintf("r = N = %d", r)), inherit.aes = FALSE,
            hjust = -0.08, size = 2.9, colour = "grey45") +                # beside the line, below every interval
  geom_rect(data = jzone, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = "grey70", alpha = 0.25) +                              # r >= N: J = N, the p-value depends on df only
  geom_text(data = jzone, aes(x = (xmin + xmax) / 2, y = 0.27, label = label), inherit.aes = FALSE, size = 3,
            colour = "grey25", lineheight = 0.95) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.45, na.rm = TRUE) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.2) +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), unname(lab_k))) +
  scale_shape_manual(values = setNames(c(16, 17), unname(lab_k))) +
  scale_x_continuous(limits = c(0, 140), breaks = seq(0, 140, 20)) +
  labs(x = sprintf("Number of instruments r (block D2: %d key accounts, %d months; lag depth 2 to all)", D2$N, D2$T),
       y = NULL, colour = NULL, shape = NULL) +
  base_theme + theme(strip.text = element_text(face = "plain"))
save_fig(p3, "fig_11_proliferation.svg", width = 9.5, height = 4.8)
