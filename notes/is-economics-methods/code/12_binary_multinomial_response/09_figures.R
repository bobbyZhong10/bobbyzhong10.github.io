# ------------------------------------------------------------------
# 09_figures.R -- chapter 12 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R
# and 20_did_beyond_binary/05_figures.R.
#   fig_12_lpm_logit.svg      block B1: share adopting in 40 engagement bins, the LPM line (above 1 and
#                             below 0 at the ends), the logit MLE curve and the true probability
#   fig_12_cross_partial.svg  block B2: each firm's estimated cross-partial effect d^2 P / dx1 dx2 against
#                             its latent index (every fifth firm), the sample average (negative) and the
#                             "marginal effect of the interaction term" beta_jj' phi(index) (positive)
#   fig_12_fe_bias.svg        block B5: Monte Carlo distributions (1,000 panels per T) of the FE logit MLE
#                             with firm dummies, its analytical bias correction (alpaca::biasCorr: Fernandez-Val
#                             and Weidner 2016, L = 0; with firm effects only it coincides with Hahn and Newey
#                             2004) and the conditional MLE, T = 4 (a 4-month DGP, pooled-logit limit dgp_B5_T4_*)
#                             and T = 8
# Reads helix_data.rds, res_helix.rds, truth12.rds, res_extra.rds, res_panel.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

dat <- readRDS("helix_data.rds"); h <- readRDS("res_helix.rds"); tr <- readRDS("truth12.rds")
x <- readRDS("res_extra.rds"); pn <- readRDS("res_panel.rds")

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name, height = 4.8) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = height, device = svglite::svglite)
  cat(file.path(fig_dir, name), sprintf("(%.0f KB)", file.size(file.path(fig_dir, name)) / 1024), "\n")
}

## ---- Figure (a): block B1, LPM against logit ---------------------------------------------------------
br <- quantile(dat$x, seq(0, 1, length.out = 41))
bin <- cut(dat$x, br, include.lowest = TRUE)
bins <- data.frame(x = tapply(dat$x, bin, mean), share = tapply(dat$D, bin, mean))
xs <- seq(min(dat$x), max(dat$x), length.out = 301)
lab <- c(lpm = "LPM (OLS) fitted value", logit = "Logit MLE", truth = "True Pr(Y = 1 | X)")
curves <- rbind(data.frame(x = xs, value = h$lpm$intercept + h$lpm$slope * xs, series = lab[["lpm"]]),
                data.frame(x = xs, value = plogis(h$logit$beta0 + h$logit$beta1 * xs), series = lab[["logit"]]),
                data.frame(x = xs, value = plogis(tr$B1$beta[1] + tr$B1$beta[2] * xs), series = lab[["truth"]]))
curves$series <- factor(curves$series, levels = unname(lab))
mx <- data.frame(x = max(dat$x), y = h$lpm$max_phat)
mn <- data.frame(x = min(dat$x), y = h$lpm$min_phat)
p1 <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = Inf, fill = col_orange, alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = col_orange, alpha = 0.08) +
  geom_hline(yintercept = c(0, 1), linewidth = 0.3, colour = col_true) +
  geom_point(data = bins, aes(x = x, y = share), colour = "grey45", size = 1.8) +
  geom_line(data = curves, aes(x = x, y = value, colour = series, linetype = series), linewidth = 0.8) +
  geom_point(data = rbind(mx, mn), aes(x = x, y = y), colour = col_orange, size = 2.2) +
  geom_text(data = mx, aes(x = x, y = y, label = sprintf("%.3f", y)), colour = col_orange, hjust = 1.2, size = 3.3) +
  geom_text(data = mn, aes(x = x, y = y, label = sprintf("%.3f", y)), colour = col_orange, hjust = -0.25, size = 3.3) +
  scale_colour_manual(values = setNames(c(col_orange, col_blue, col_true), unname(lab))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed"), unname(lab))) +
  scale_y_continuous(breaks = seq(-0.5, 1.25, 0.25)) +
  labs(x = "Engagement X (z-score of past API activity)", y = "Probability of adopting Pro",
       colour = NULL, linetype = NULL, caption = "Points: share adopting in 40 equal-count engagement bins (n = 4,000 firms)") +
  base_theme + theme(plot.caption = element_text(colour = "grey40", size = 9))
save_fig(p1, "fig_12_lpm_logit.svg")

## ---- Figure (b): block B2, cross-partial effects --------------------------------------------------------
B2 <- x$B2
keep <- seq(1, length(B2$index), by = 5)
pts <- data.frame(index = B2$index[keep], cross = B2$cross_i[keep])
us <- seq(min(B2$index), max(B2$index), length.out = 301)
lab2 <- c(firm = "Cross-partial effect of one firm", avg = "Average cross-partial effect",
          naive = "Interaction coefficient times phi(index)")
lines2 <- rbind(data.frame(index = range(us), value = B2$avg_cross, series = lab2[["avg"]]),
                data.frame(index = us, value = B2$coef_interaction * dnorm(us), series = lab2[["naive"]]))
lines2$series <- factor(lines2$series, levels = unname(lab2)[2:3])
p2 <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = col_true) +
  geom_point(data = pts, aes(x = index, y = cross), colour = "grey55", alpha = 0.35, size = 0.8) +
  geom_line(data = lines2, aes(x = index, y = value, colour = series, linetype = series), linewidth = 0.9) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), unname(lab2)[2:3])) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), unname(lab2)[2:3])) +
  labs(x = "Estimated latent index X'beta (probit)", y = expression(partialdiff^2 * P / (partialdiff * x[1] * partialdiff * x[2])),
       colour = NULL, linetype = NULL,
       caption = "Grey points: estimated cross-partial effect of every fifth firm (800 of 4,000)") +
  base_theme + theme(plot.caption = element_text(colour = "grey40", size = 9))
save_fig(p2, "fig_12_cross_partial.svg")

## ---- Figure (c): block B5, incidental parameters ---------------------------------------------------------
lab3 <- c(fe = "FE logit MLE with firm dummies", bc = "Bias-corrected FE (Fernandez-Val and Weidner)", cl = "Conditional MLE (FE logit)")
dd <- do.call(rbind, lapply(names(pn$mc), function(k) {
  M <- pn$mc[[k]]$draws
  rbind(data.frame(T = k, est = M[, "fe"], series = lab3[["fe"]]),
        data.frame(T = k, est = M[, "bc"], series = lab3[["bc"]]),
        data.frame(T = k, est = M[, "cl"], series = lab3[["cl"]]))
}))
dd$panel <- factor(ifelse(dd$T == "T4", "T = 4 months", "T = 8 months"), levels = c("T = 4 months", "T = 8 months"))
dd$series <- factor(dd$series, levels = unname(lab3))
plim <- data.frame(panel = factor(c("T = 4 months", "T = 8 months"), levels = levels(dd$panel)),
                   value = c(tr$B5$T4$pooled_plim, tr$B5$T8$pooled_plim))
p3 <- ggplot(dd, aes(x = est, colour = series, linetype = series)) +
  geom_density(linewidth = 0.8, n = 256, bw = "nrd0", key_glyph = "path") +
  geom_vline(xintercept = tr$B5$beta, linewidth = 0.5, colour = col_true) +
  geom_vline(data = plim, aes(xintercept = value), linewidth = 0.4, colour = col_true, linetype = "dotted") +
  facet_wrap(~panel) +
  scale_colour_manual(values = setNames(c(col_orange, col_green, col_blue), unname(lab3))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "solid"), unname(lab3))) +
  labs(x = "Estimate of beta (true value 0.8: solid line; pooled logit limit: dotted line)", y = "Density",
       colour = NULL, linetype = NULL,
       caption = "1,000 simulated panels of N = 1,000 firms per panel length") +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  base_theme + theme(plot.caption = element_text(colour = "grey40", size = 9), strip.text = element_text(face = "plain"))
save_fig(p3, "fig_12_fe_bias.svg")
