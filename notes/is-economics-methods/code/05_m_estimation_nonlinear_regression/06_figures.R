# ------------------------------------------------------------------
# 06_figures.R -- chapter 5 figures (SVG), same style as 20_did_beyond_binary/05_figures.R.
#   fig_05_log1p_units.svg   block M1: OLS coefficient of log(1 + y) on Spotlight against the unit in
#                            which sales are recorded, with its population counterpart, against PPML
#                            (the same in every unit) and the true log(E[y(1)] / E[y(0)])
#   fig_05_misspec.svg       block M2m: the true mean, the NLS and the PPML pseudo-true exponential
#                            means and every second sample point (log scale: an exponential mean is a straight line)
#   fig_05_coverage.svg      blocks M4 and M3: Monte Carlo coverage of nominal 95% intervals; model-based
#                            or naive against sandwich or two-step corrected standard errors
# Reads truth05.rds, larkspur_M2m.rds, res_logs.rds, res_nls.rds, res_twostep.rds, res_mc.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

tr <- readRDS("truth05.rds")
M2m <- readRDS("larkspur_M2m.rds")
lg <- readRDS("res_logs.rds")
nl <- readRDS("res_nls.rds")
ts <- readRDS("res_twostep.rds")
mc <- readRDS("res_mc.rds")

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}

## ---- Figure (a): block M1, log(1 + y) against the unit ----------------------------------
sw <- lg$sweep
lab <- c(ols = "OLS of log(1 + y) on Spotlight (95% CI)", pop = "Population ATE on log(1 + y)",
         ppml = "PPML, any unit (95% CI)", truth = "True log(E[y(1)] / E[y(0)])")
brk <- 10^(-3:4)
brk_lab <- c("USD 1,000", "USD 100", "USD 10", "USD 1", "10 cents", "1 cent", "0.1 cent", "0.01 cent")
named <- sw[sapply(c(1e-3, 1, 100), function(a) which.min(abs(log10(sw$a) - log10(a)))), ]
lines_a <- rbind(data.frame(a = sw$a, value = sw$est, series = lab[["ols"]]),
                 data.frame(a = sw$a, value = sw$ate_pop, series = lab[["pop"]]),
                 data.frame(a = range(sw$a), value = lg$ppml$usd$est, series = lab[["ppml"]]),
                 data.frame(a = range(sw$a), value = tr$M1$log_ratio, series = lab[["truth"]]))
lines_a$series <- factor(lines_a$series, levels = unname(lab))
p1 <- ggplot() +
  annotate("rect", xmin = min(sw$a), xmax = max(sw$a), ymin = lg$ppml$usd$ci[1], ymax = lg$ppml$usd$ci[2],
           fill = col_orange, alpha = 0.12) +
  geom_ribbon(data = sw, aes(x = a, ymin = lo, ymax = hi), fill = col_blue, alpha = 0.18) +
  geom_line(data = lines_a, aes(x = a, y = value, colour = series, linetype = series), linewidth = 0.7) +
  geom_point(data = named, aes(x = a, y = est), colour = col_blue, size = 2.2) +
  geom_text(data = named, aes(x = a, y = est, label = sprintf("%.3f", est)), colour = col_blue,
            vjust = -0.9, hjust = 0.9, size = 3.2) +
  scale_x_log10(breaks = brk, labels = brk_lab) +
  scale_colour_manual(values = setNames(c(col_blue, col_blue, col_orange, col_true), unname(lab))) +
  scale_linetype_manual(values = setNames(c("solid", "dotted", "dashed", "solid"), unname(lab))) +
  labs(x = "Unit in which 30-day sales are recorded (log scale)", y = "Effect of Spotlight (coefficient on D)",
       colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  base_theme
save_fig(p1, "fig_05_log1p_units.svg")

## ---- Figure (b): block M2m, misspecified mean ---------------------------------------------
P2m <- tr$params$M2m
xs <- seq(-3 * P2m$sd_x, 3 * P2m$sd_x, length.out = 301)
lab2 <- c(true = "True mean E[y | x]", nls = "NLS pseudo-true exp(b0* + b1* x)",
          ppml = "PPML pseudo-true exp(b0* + b1* x)")
curves <- rbind(
  data.frame(x = xs, value = exp(P2m$a0 + P2m$a1 * xs + P2m$a2 * xs^2), series = lab2[["true"]]),
  data.frame(x = xs, value = exp(tr$M2m$pseudo_vec[1] + tr$M2m$pseudo_vec[2] * xs), series = lab2[["nls"]]),
  data.frame(x = xs, value = exp(tr$M2m$pseudo_ppml_vec[1] + tr$M2m$pseudo_ppml_vec[2] * xs), series = lab2[["ppml"]]))
curves$series <- factor(curves$series, levels = unname(lab2))
# every second merchant (the rows are i.i.d., so this is a random half, about 1,000 points): keeps the
# SVG vector and under inline_figures.py's 300 KB RASTER_THRESHOLD
pts <- M2m[seq(1, nrow(M2m), by = 2), ]
pts <- pts[abs(pts$x) <= 3 * P2m$sd_x, ]
p2 <- ggplot() +
  geom_point(data = pts, aes(x = x, y = y), colour = "grey60", alpha = 0.25, size = 0.8) +
  geom_line(data = curves, aes(x = x, y = value, colour = series, linetype = series), linewidth = 0.8) +
  scale_y_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8), labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8")) +
  scale_colour_manual(values = setNames(c(col_true, col_blue, col_orange), unname(lab2))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed"), unname(lab2))) +
  labs(x = "x = log(catalogue size / 40 listings)", y = "Monthly sales (USD thousands, log scale)",
       colour = NULL, linetype = NULL) +
  base_theme
save_fig(p2, "fig_05_misspec.svg")

## ---- Figure (c): blocks M4 and M3, coverage -----------------------------------------------
# strip labels are plotmath (facet labeller = label_parsed) so that the M3 axis can name beta[1] as the text does
panel <- c(M4 = 'atop("Block M4: Poisson PML on overdispersed orders", "(x: variance-to-mean ratio E[Var(y | x)] / E[y])")',
           M3 = 'atop("Block M3: two-step NLS on an estimated traffic index", "(x: true traffic elasticity " * beta[1] * ")")')
lab3 <- c(bad = "Model-based (M4) or naive (M3) standard error",
          good = "Sandwich (M4) or two-step corrected (M3) standard error")
g4 <- mc$grid; g3 <- ts$grid
cv <- rbind(
  data.frame(panel = panel[["M4"]], x = g4$avg_cond_var_to_mean, cover = g4$cover_model, series = lab3[["bad"]]),
  data.frame(panel = panel[["M4"]], x = g4$avg_cond_var_to_mean, cover = g4$cover_sandwich, series = lab3[["good"]]),
  data.frame(panel = panel[["M3"]], x = g3$theta1, cover = g3$cover_naive, series = lab3[["bad"]]),
  data.frame(panel = panel[["M3"]], x = g3$theta1, cover = g3$cover_corrected, series = lab3[["good"]]))
cv$panel <- factor(cv$panel, levels = unname(panel)); cv$series <- factor(cv$series, levels = unname(lab3))
band <- qnorm(0.975) * sqrt(0.95 * 0.05 / min(c(g4$R, g3$R)))       # Monte Carlo band around 0.95 at 1000 replications
p3 <- ggplot(cv, aes(x = x, y = cover, colour = series, shape = series)) +
  geom_hline(yintercept = 0.95, linewidth = 0.4, colour = col_true) +
  geom_hline(yintercept = 0.95 + c(-1, 1) * band, linewidth = 0.3, colour = col_true, linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~panel, scales = "free_x", labeller = label_parsed) +
  scale_colour_manual(values = setNames(c(col_orange, col_blue), unname(lab3))) +
  scale_shape_manual(values = setNames(c(17, 16), unname(lab3))) +
  scale_x_continuous(breaks = function(lim) if (max(lim) > 2) c(1, 2, 4, 6, 8) else seq(0, 1.2, 0.2)) +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
  labs(x = NULL, y = "Coverage of the nominal 95% interval", colour = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  base_theme + theme(strip.text = element_text(face = "plain"))
save_fig(p3, "fig_05_coverage.svg")
