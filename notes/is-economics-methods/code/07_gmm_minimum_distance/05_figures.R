# ------------------------------------------------------------------
# 05_figures.R -- chapter 7 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R
# and 20_did_beyond_binary/05_figures.R.
#   fig_07_md_markets.svg   blocks G1 and G1h: the 16 market price slopes with 95% intervals (point size: the
#                           optimal MD weight), the MD estimate and its interval, the truth (theta in G1, theta_m
#                           in G1h), and the J statistic of each panel
#   fig_07_gmm_coverage.svg block G2: Monte Carlo coverage of nominal 95% intervals for the elasticity against
#                           the number of regions, for 2SLS, identity-weighted, two-step (conventional and
#                           Windmeijer SE), iterated and continuously updated GMM
#   fig_07_s_set.svg        blocks G2 and G3: the S statistic S_n(theta) of the sample against theta, its
#                           chi^2 critical value, the S confidence set and the two-step Wald interval
# Reads truth07.rds, res_md.rds, res_gmm_mc.rds, res_weak.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

tr <- readRDS("truth07.rds")
md <- readRDS("res_md.rds")
g2 <- readRDS("res_gmm_mc.rds")
g3 <- readRDS("res_weak.rds")

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
col_purple <- "#CC79A7"; col_sky <- "#56B4E9"; col_amber <- "#E69F00"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name, height = 4.8) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = height, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}

## ---- Figure (a): blocks G1 and G1h, minimum distance across markets ------------------------------
panel <- c(G1 = sprintf("G1: one elasticity in every market\nJ = %.1f (df %d, p = %.2f)", md$G1$J, md$G1$df, md$G1$J_p),
           G1h = sprintf("G1h: elasticity differs across markets\nJ = %.1f (df %d, p = %.1e)", md$G1h$J, md$G1h$df, md$G1h$J_p))
ord <- md$G1$slopes$market[order(md$G1$slopes$weight)]            # same order in both panels: largest weight on top
sl <- rbind(cbind(md$G1$slopes[, c("market", "pi_hat", "lo", "hi", "weight")], panel = panel[["G1"]], truth = tr$G1$theta),
            cbind(md$G1h$slopes[, c("market", "pi_hat", "lo", "hi", "weight")], panel = panel[["G1h"]], truth = tr$G1h$theta_m))
sl$market <- factor(sl$market, levels = levels(ord)[as.integer(ord)])
sl$panel <- factor(sl$panel, levels = unname(panel))
mdl <- data.frame(panel = factor(unname(panel), levels = unname(panel)), est = c(md$G1$theta, md$G1h$theta),
                  lo = c(md$G1$ci[["lo"]], md$G1h$ci[["lo"]]), hi = c(md$G1$ci[["hi"]], md$G1h$ci[["hi"]]))
lab_t <- c(md = "MD estimate and 95% CI (inverse-variance weighted mean)",
           tr = "Truth: theta (G1), theta_m (G1h)")
p1 <- ggplot(sl, aes(y = market)) +
  geom_rect(data = mdl, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = col_blue, alpha = 0.15) +
  geom_vline(data = mdl, aes(xintercept = est, linetype = lab_t[["md"]]), colour = col_blue, linewidth = 0.6) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.4, colour = "grey35") +
  geom_point(aes(x = pi_hat, size = weight), colour = "grey15") +
  geom_point(aes(x = truth, shape = lab_t[["tr"]]), colour = col_orange, size = 2.4, stroke = 0.9) +
  facet_wrap(~panel) +
  scale_shape_manual(values = setNames(4, lab_t[["tr"]])) +
  scale_linetype_manual(values = setNames("solid", lab_t[["md"]])) +
  scale_size_area(max_size = 3.6, breaks = c(0.02, 0.08, 0.16), name = "MD weight") +
  coord_cartesian(xlim = c(-4.5, 1.5)) +
  labs(x = "Price slope of log weekly lessons booked (reduced form pi_m, 95% CI)", y = NULL, shape = NULL, linetype = NULL) +
  guides(linetype = guide_legend(order = 1), shape = guide_legend(order = 2), size = guide_legend(order = 3)) +
  base_theme + theme(strip.text = element_text(face = "plain"), legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
save_fig(p1, "fig_07_md_markets.svg", height = 6)

## ---- Figure (b): block G2, coverage against the number of regions ---------------------------------
lab_e <- c(one_step = "2SLS (one step), robust SE", identity = "Identity weight, robust SE",
           two_step_naive = "Two-step GMM, conventional SE", windmeijer = "Two-step GMM, Windmeijer SE",
           iterated = "Iterated GMM, conventional SE", cue = "CUE, conventional SE")
gr <- g2$grid
cv <- do.call(rbind, lapply(names(lab_e), function(k) data.frame(n = gr$n, cover = gr[[paste0("cover_", k)]], series = lab_e[[k]])))
cv$series <- factor(cv$series, levels = unname(lab_e))
band <- qnorm(0.975) * sqrt(0.95 * 0.05 / min(gr$R))                # Monte Carlo band around 0.95 at the smallest R
p2 <- ggplot(cv, aes(x = n, y = cover, colour = series, shape = series, linetype = series)) +
  geom_hline(yintercept = 0.95, linewidth = 0.4, colour = col_true) +
  geom_hline(yintercept = 0.95 + c(-1, 1) * band, linewidth = 0.3, colour = col_true, linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_x_log10(breaks = gr$n) +
  scale_y_continuous(limits = c(0.84, 0.98), breaks = seq(0.84, 0.98, 0.02)) +
  scale_colour_manual(values = setNames(c(col_green, col_amber, col_orange, col_blue, col_purple, col_sky), unname(lab_e))) +
  scale_shape_manual(values = setNames(c(15, 4, 17, 16, 18, 1), unname(lab_e))) +
  scale_linetype_manual(values = setNames(c("solid", "dotdash", "solid", "solid", "dashed", "dashed"), unname(lab_e))) +
  labs(x = "Number of regions n (log scale; 11 moments, 2 parameters)",
       y = "Coverage of the nominal 95% interval", colour = NULL, shape = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  base_theme
save_fig(p2, "fig_07_gmm_coverage.svg")

## ---- Figure (c): blocks G2 and G3, the S statistic and the S set -------------------------------------
pan <- c(G2 = sprintf("G2: 10 price experiments, first-stage F = %.0f", g2$sample$first_stage_F),
         G3 = sprintf("G3: 4 penny-sized perturbations, first-stage F = %.1f", g3$sample$first_stage_F))
s2 <- g2$sample$s_curve; s3 <- g3$sample$s_curve
curves <- rbind(data.frame(panel = pan[["G2"]], theta = s2$theta, S = s2$S, in_set = s2$in_set),
                data.frame(panel = pan[["G3"]], theta = s3$theta, S = s3$S, in_set = s3$in_set))
curves$panel <- factor(curves$panel, levels = unname(pan))
crit <- data.frame(panel = factor(unname(pan), levels = unname(pan)),
                   crit = c(qchisq(0.95, g2$sample$s_df), g3$sample$crit),
                   lab = c(sprintf("chi-square(%d)\n95%% critical value", g2$sample$s_df), sprintf("chi-square(%d)\n95%% critical value", g3$sample$df)))
wald <- data.frame(panel = factor(unname(pan), levels = unname(pan)),
                   lo = c(g2$sample$est$two_step_naive$ci[["lo"]], g3$sample$two_step$ci[["lo"]]),
                   hi = c(g2$sample$est$two_step_naive$ci[["hi"]], g3$sample$two_step$ci[["hi"]]),
                   est = c(g2$sample$est$two_step_naive$est, g3$sample$two_step$est))
sset <- data.frame(panel = factor(unname(pan), levels = unname(pan)),
                   lo = c(g2$sample$s_set[["lo"]], g3$sample$s_set[["lo"]]), hi = c(g2$sample$s_set[["hi"]], g3$sample$s_set[["hi"]]))
lab_c <- c(s = "S statistic S_n(theta)", set = "S confidence set (95%)", wald = "Two-step GMM Wald interval (95%)",
           tr = "True elasticity -1.6")
p3 <- ggplot() +
  geom_rect(data = sset, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = lab_c[["set"]]), alpha = 0.18) +
  geom_hline(data = crit, aes(yintercept = crit), colour = col_true, linewidth = 0.4, linetype = "dashed") +
  # label at the right end, below the line: there the S curve lies above the critical value (in both panels S > crit
  # for theta beyond the upper S-set bound), so the curve does not cross the text
  geom_text(data = crit, aes(x = Inf, y = crit, label = lab), hjust = 1.05, vjust = 1.25, size = 3, lineheight = 0.9, colour = col_true) +
  geom_line(data = curves, aes(x = theta, y = S, colour = lab_c[["s"]]), linewidth = 0.7) +
  geom_errorbar(data = wald, aes(xmin = lo, xmax = hi, y = -2.5, colour = lab_c[["wald"]]), orientation = "y",
                width = 2.5, linewidth = 0.7) +
  geom_point(data = wald, aes(x = est, y = -2.5, colour = lab_c[["wald"]]), size = 2) +
  geom_vline(data = data.frame(x = tr$G2$theta), aes(xintercept = x, colour = lab_c[["tr"]]), linewidth = 0.5, linetype = "dotted") +
  facet_wrap(~panel, scales = "free") +
  scale_colour_manual(values = setNames(c(col_blue, col_orange, col_true), unname(lab_c[c("s", "wald", "tr")])),
                      breaks = unname(lab_c[c("s", "wald", "tr")])) +
  scale_fill_manual(values = setNames(col_blue, lab_c[["set"]])) +
  labs(x = "Elasticity theta", y = "S statistic", colour = NULL, fill = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = c("solid", "solid", "dotted"))), fill = guide_legend(order = 2)) +
  base_theme + theme(strip.text = element_text(face = "plain"))
save_fig(p3, "fig_07_s_set.svg")
