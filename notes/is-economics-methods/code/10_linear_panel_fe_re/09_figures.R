# ------------------------------------------------------------------
# 09_figures.R -- chapter 10 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R
# and 20_did_beyond_binary/05_figures.R (replaces the legacy 04_figures.R).
#   fig_10_estimators.svg    block P1 (inherited static panel): pooled OLS, between, RE, FE, FD and
#                            Chamberlain MD with 95% intervals, their probability limits and the truth
#                            (redrawn at 8 in wide; the legacy figure showed pooled, RE and FE only)
#   fig_10_hausman_size.svg  block P7: Monte Carlo rejection rate of the classical and the robust
#                            (Mundlak) Hausman tests against the nominal level, RE.1 true, RE.3 false
#   fig_10_slopes.svg        blocks P5 and P5b: merchant slopes by within-merchant promotion SD, the
#                            FE weight of each type, the average slope, the FE probability limit and
#                            the FE estimate
# Reads truth10.rds, res_static.rds, res_extra.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

tr <- readRDS("truth10.rds"); st <- readRDS("res_static.rds"); x <- readRDS("res_extra.rds")

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name, height = 4.8) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = height, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}

## ---- Figure (a): block P1, static estimators -------------------------------------------------
P1 <- x$P1
lab_e <- c(pool = "Pooled OLS", be = "Between", re = "Random effects", fe = "Fixed effects (within)",
           fd = "First difference", ch = "Chamberlain (minimum distance)")
est <- data.frame(
  estimator = unname(lab_e),
  est = c(P1$pool$est, P1$be$est, P1$re$est, P1$fe$est, P1$fd$est, P1$chamberlain$beta),
  lo = c(P1$pool$ci[1], P1$be$ci[1], P1$re$ci[1], P1$fe$ci[1], P1$fd$ci[1], P1$chamberlain$ci[1]),
  hi = c(P1$pool$ci[2], P1$be$ci[2], P1$re$ci[2], P1$fe$ci[2], P1$fd$ci[2], P1$chamberlain$ci[2]))
plim <- data.frame(estimator = unname(lab_e),
                   value = c(tr$P1$plim_pols, tr$P1$plim_be, tr$P1$plim_re, rep(tr$P1$beta, 3)))
est$estimator <- factor(est$estimator, levels = rev(unname(lab_e)))
plim$estimator <- factor(plim$estimator, levels = rev(unname(lab_e)))
lab_s <- c(e = "Estimate, 95% CI (clustered by merchant; between: one row per merchant)",
           p = "Probability limit as the number of merchants grows")
est$series <- lab_s[["e"]]; plim$series <- lab_s[["p"]]
p1 <- ggplot() +
  geom_vline(xintercept = tr$P1$beta, linewidth = 0.5, colour = col_true, linetype = "dashed") +
  geom_errorbar(data = est, aes(y = estimator, xmin = lo, xmax = hi, colour = series), orientation = "y",
                width = 0, linewidth = 0.5) +
  geom_point(data = est, aes(y = estimator, x = est, colour = series, shape = series), size = 2.3) +
  geom_point(data = plim, aes(y = estimator, x = value, colour = series, shape = series), size = 2.6, stroke = 0.9) +
  geom_text(data = est, aes(y = estimator, x = est, label = sprintf("%.3f", est)), colour = col_blue,
            vjust = -0.85, size = 3) +
  annotate("text", x = tr$P1$beta, y = 6.45, label = "True beta = 0.5", hjust = -0.08, size = 3.2, colour = col_true) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), unname(lab_s))) +
  scale_shape_manual(values = setNames(c(16, 5), unname(lab_s))) +
  scale_x_continuous(breaks = seq(0.4, 1.6, 0.2), limits = c(0.4, 1.7)) +
  labs(x = "Estimate of beta (inherited static Nimbus panel: 600 merchants, 6 months)",
       y = NULL, colour = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  base_theme
save_fig(p1, "fig_10_estimators.svg")

## ---- Figure (b): block P7, size of the Hausman tests ----------------------------------------
cv <- x$P7$curve
lab_h <- c(c = "Classical Hausman (FE vs RE, non-robust variances)", r = "Robust Hausman (Mundlak, clustered by merchant)")
sz <- rbind(data.frame(nominal = cv$nominal, rate = cv$classic, series = lab_h[["c"]]),
            data.frame(nominal = cv$nominal, rate = cv$robust, series = lab_h[["r"]]))
sz$series <- factor(sz$series, levels = unname(lab_h))
band <- data.frame(nominal = cv$nominal, lo = pmax(0, cv$nominal - qnorm(0.975) * sqrt(cv$nominal * (1 - cv$nominal) / x$P7$R)),
                   hi = cv$nominal + qnorm(0.975) * sqrt(cv$nominal * (1 - cv$nominal) / x$P7$R))
y_top <- ceiling(20 * max(sz$rate, band$hi)) / 20
p2 <- ggplot() +
  geom_ribbon(data = band, aes(x = nominal, ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, colour = col_true) +
  geom_line(data = sz, aes(x = nominal, y = rate, colour = series), linewidth = 0.7) +
  geom_point(data = sz, aes(x = nominal, y = rate, colour = series, shape = series), size = 2.2) +
  annotate("text", x = 0.115, y = 0.03, label = "Black line: rejection rate = nominal level\nGrey band: Monte Carlo 95% band around it",
           size = 3.1, colour = col_true, hjust = 0) +
  scale_colour_manual(values = setNames(c(col_orange, col_blue), unname(lab_h))) +
  scale_shape_manual(values = setNames(c(17, 16), unname(lab_h))) +
  scale_x_continuous(breaks = cv$nominal, labels = sprintf("%g%%", 100 * cv$nominal)) +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, y_top, 0.05), labels = function(v) sprintf("%g%%", 100 * v)) +
  labs(x = "Nominal level of the test", y = sprintf("Rejection rate in %s samples (H0 true)", format(x$P7$R, big.mark = ",")),
       colour = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  base_theme
save_fig(p2, "fig_10_hausman_size.svg")

## ---- Figure (c): blocks P5 and P5b, random slopes ------------------------------------------
panel <- c(P5 = "Block P5: slopes rise with the within-merchant SD",
           P5b = "Block P5b: slopes follow the promotion level only")
pts <- do.call(rbind, lapply(names(panel), function(b) {
  u <- x[[b]]$unit; u$id <- seq_len(nrow(u)); u <- u[u$id %% 4 == 0, ]          # every 4th merchant
  u$rk <- ave(u$id, u$s, FUN = function(v) (rank(v) - 0.5) / length(v))
  data.frame(panel = panel[[b]], s = u$s + (u$rk - 0.5) * 0.28, b = u$b)
}))
typ <- do.call(rbind, lapply(names(panel), function(b) {
  bt <- x[[b]]$by_type
  data.frame(panel = panel[[b]], s = bt$s, mean_b = bt$mean_b, weight = bt$weight_fe)
}))
ref <- do.call(rbind, lapply(names(panel), function(b) {
  avg <- tr[[b]]$avg_slope; pl <- tr[[b]]$fe_plim
  data.frame(panel = panel[[b]], series = c("Average slope E[b_i]", "FE probability limit (within-variance weighted)"),
             value = c(avg, pl))
}))
fe <- do.call(rbind, lapply(names(panel), function(b) {
  f <- x[[b]]$fe; data.frame(panel = panel[[b]], est = f$est, lo = f$ci[1], hi = f$ci[2])
}))
lev_panel <- unname(panel)
pts$panel <- factor(pts$panel, levels = lev_panel); typ$panel <- factor(typ$panel, levels = lev_panel)
ref$panel <- factor(ref$panel, levels = lev_panel); fe$panel <- factor(fe$panel, levels = lev_panel)
ref$series <- factor(ref$series, levels = c("Average slope E[b_i]", "FE probability limit (within-variance weighted)"))
p3 <- ggplot() +
  geom_rect(data = fe, aes(ymin = lo, ymax = hi), xmin = 0.2, xmax = 2.3, fill = col_blue, alpha = 0.15) +
  geom_point(data = pts, aes(x = s, y = b), colour = "grey60", size = 0.7, alpha = 0.6) +
  geom_hline(data = ref, aes(yintercept = value, linetype = series, colour = series), linewidth = 0.6) +
  geom_point(data = typ, aes(x = s, y = mean_b, size = weight), colour = col_blue, alpha = 0.9) +
  geom_text(data = typ, aes(x = s, y = -0.25, label = sprintf("FE weight %.0f%%", 100 * weight)),
            size = 2.8, colour = col_blue) +
  facet_wrap(~panel) +
  scale_colour_manual(values = c(col_true, col_orange), name = NULL) +
  scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
  scale_size_area(max_size = 6, guide = "none") +
  scale_x_continuous(breaks = c(0.5, 1, 1.5, 2), labels = c("50", "100", "150", "200"), limits = c(0.2, 2.3)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(-0.3, NA)) +
  labs(x = "Within-merchant SD of monthly promotion spend (USD)",
       y = "Slope b_i (USD 1,000 of GMV per USD 100 of promotion)",
       caption = "Grey: true slopes b_i of every 4th merchant; blue dots: type averages (area = FE weight); blue band: FE 95% CI") +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  base_theme + theme(strip.text = element_text(face = "plain"), plot.caption = element_text(size = 8.5, colour = col_true))
save_fig(p3, "fig_10_slopes.svg", height = 5)
