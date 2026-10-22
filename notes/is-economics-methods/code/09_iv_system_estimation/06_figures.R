# ------------------------------------------------------------------
# 06_figures.R -- chapter 9 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R
# and 20_did_beyond_binary/05_figures.R.
#   fig_09_simultaneity.svg  block S1 (Working 1927): markets in the plane of sellers and buyers net of the
#                            buyer-side shifters and market size (Frisch-Waugh-Lovell), every third market;
#                            the OLS fit through the cloud, the true buyer equation, and the means of five
#                            fee groups, which line up along the 2SLS slope: the seller-side shifter traces
#                            out the buyer equation
#   fig_09_contagion.svg     blocks S1 and S1m: Monte Carlo sampling distributions of the buyer coefficient
#                            n_sellers under 2SLS and 3SLS, with both equations right and with the seller
#                            equation missing coupon x buyers
#   fig_09_weak_iv.svg       block S3: Monte Carlo coverage of the nominal 95% Wald and Anderson-Rubin sets and
#                            the median 2SLS estimate (with its interquartile range) against the mean
#                            first-stage F
# Reads truth09.rds, bramble_S1.rds, res_system.rds, res_weak.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("iv_fns.R")
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

tr <- readRDS("truth09.rds")
S1 <- readRDS("bramble_S1.rds")
s1 <- readRDS("res_system.rds")
s3 <- readRDS("res_weak.rds")
b <- tr$S1$buyer[["n_sellers"]]

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"; col_purple <- "#CC79A7"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), sprintf("(%.0f KB)", file.size(file.path(fig_dir, name)) / 1024), "\n")
}

## ---- Figure (a): block S1, the identification picture ---------------------------------------
W <- cbind(1, S1$ad_budget, S1$coupon, S1$log_pop)
xt <- drop(partial_out(S1$n_sellers, W)); yt <- drop(partial_out(S1$buyers, W)); ft <- drop(partial_out(S1$fee, W))
grp <- cut(ft, quantile(ft, seq(0, 1, 0.2)), include.lowest = TRUE, labels = FALSE)
cent <- data.frame(x = tapply(xt, grp, mean), y = tapply(yt, grp, mean))
pts <- data.frame(x = xt, y = yt)[seq(1, nrow(S1), by = 3), ]        # every third market (500 points): small vector SVG
lab1 <- c(ols = sprintf("OLS fit through the markets (slope %.2f)", s1$ols$buyer_est),
          tsls = sprintf("2SLS with the fee as instrument (slope %.2f)", s1$tsls$buyer_est),
          truth = sprintf("True buyer equation (slope %.1f)", b))
xr <- range(pts$x)
lines1 <- do.call(rbind, lapply(names(lab1), function(k) {
  sl <- switch(k, ols = s1$ols$buyer_est, tsls = s1$tsls$buyer_est, truth = b)
  data.frame(x = xr, y = sl * xr, series = lab1[[k]])
}))
lines1$series <- factor(lines1$series, levels = unname(lab1))
p1 <- ggplot() +
  geom_point(data = pts, aes(x = x, y = y), colour = "grey60", alpha = 0.35, size = 0.9) +
  geom_line(data = lines1, aes(x = x, y = y, colour = series, linetype = series), linewidth = 0.75) +
  geom_point(data = cent, aes(x = x, y = y), colour = col_orange, fill = "white", shape = 21, size = 3, stroke = 1.1) +
  annotate("text", x = cent$x[1], y = cent$y[1], label = "Means of five fee groups", hjust = -0.08, vjust = 1.8,
           size = 3.2, colour = col_orange) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange, col_true), unname(lab1))) +
  scale_linetype_manual(values = setNames(c("dashed", "solid", "solid"), unname(lab1))) +
  labs(x = "Active sellers (hundreds), net of ad budget, coupon and market size",
       y = "Active buyers (thousands), net of the same", colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  base_theme
save_fig(p1, "fig_09_simultaneity.svg")

## ---- Figure (b): blocks S1 and S1m, contagion ------------------------------------------------
panel <- c(S1 = "Block S1: both equations right", S1m = "Block S1m: seller equation omits coupon x buyers")
lab2 <- c(tsls = "2SLS, buyer equation alone", three = "3SLS, both equations jointly")
dens <- do.call(rbind, lapply(names(panel), function(k) {
  dr <- s1$mc$draws[[k]]
  rbind(data.frame(panel = panel[[k]], est = dr$tsls, series = lab2[["tsls"]]),
        data.frame(panel = panel[[k]], est = dr$three_sls, series = lab2[["three"]]))
}))
dens$panel <- factor(dens$panel, levels = unname(panel)); dens$series <- factor(dens$series, levels = unname(lab2))
p2 <- ggplot(dens, aes(x = est, colour = series, linetype = series)) +
  geom_vline(xintercept = b, linewidth = 0.4, colour = col_true, linetype = "dashed") +
  geom_density(linewidth = 0.75, n = 256, adjust = 1.1, key_glyph = "path") +
  facet_wrap(~panel) +
  scale_x_continuous(breaks = seq(0.25, 1.5, 0.25)) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), unname(lab2))) +
  scale_linetype_manual(values = setNames(c("solid", "solid"), unname(lab2))) +
  labs(x = sprintf("Estimate of the sellers' effect on buyers (true value %.1f, dashed); %d Monte Carlo samples each",
                   b, s1$mc$S1$R), y = "Density", colour = NULL, linetype = NULL) +
  base_theme + theme(strip.text = element_text(face = "plain"), plot.margin = margin(5.5, 14, 5.5, 5.5))
save_fig(p2, "fig_09_contagion.svg")

## ---- Figure (c): block S3, weak instrument -----------------------------------------------------
g <- s3$grid
pan3 <- c(cov = "Coverage of the nominal 95% set", med = "Median 2SLS estimate (interquartile range shaded)")
lab3 <- c(wald = "Wald interval (HC0)", ar = "Anderson-Rubin set", med = "2SLS median",
          truth = "True value", ols = "OLS probability limit")
cv <- rbind(data.frame(panel = pan3[["cov"]], F = g$mean_F, value = g$cover_wald, series = lab3[["wald"]]),
            data.frame(panel = pan3[["cov"]], F = g$mean_F, value = g$cover_ar, series = lab3[["ar"]]),
            data.frame(panel = pan3[["med"]], F = g$mean_F, value = g$median_est, series = lab3[["med"]]),
            data.frame(panel = pan3[["med"]], F = g$mean_F, value = b, series = lab3[["truth"]]),
            data.frame(panel = pan3[["med"]], F = g$mean_F, value = g$plim_ols, series = lab3[["ols"]]))
cv$panel <- factor(cv$panel, levels = unname(pan3)); cv$series <- factor(cv$series, levels = unname(lab3))
band <- qnorm(0.975) * sqrt(0.95 * 0.05 / min(g$R))                 # Monte Carlo band around 0.95 at 1000 replications
ref <- data.frame(panel = factor(pan3[["cov"]], levels = unname(pan3)), y = c(0.95, 0.95 - band, 0.95 + band),
                  type = c("solid", "dotted", "dotted"))
rib <- data.frame(panel = factor(pan3[["med"]], levels = unname(pan3)), F = g$mean_F, lo = g$q25_est, hi = g$q75_est)
p3 <- ggplot() +
  geom_hline(data = ref, aes(yintercept = y, linetype = I(type)), linewidth = 0.35, colour = col_true) +
  geom_vline(xintercept = 10, linewidth = 0.35, colour = "grey55", linetype = "dotdash") +
  geom_ribbon(data = rib, aes(x = F, ymin = lo, ymax = hi), fill = col_orange, alpha = 0.15) +
  geom_line(data = cv, aes(x = F, y = value, colour = series, linetype = series), linewidth = 0.7) +
  geom_point(data = cv[cv$series %in% lab3[c("wald", "ar", "med")], ], aes(x = F, y = value, colour = series), size = 1.8) +
  facet_wrap(~panel, scales = "free_y") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_colour_manual(values = setNames(c(col_purple, col_blue, col_orange, col_true, col_green), unname(lab3))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "solid", "dashed", "dashed"), unname(lab3))) +
  labs(x = "Mean first-stage F across Monte Carlo samples (log scale; dot-dash line at 10)", y = NULL, colour = NULL,
       linetype = NULL) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  base_theme + theme(strip.text = element_text(face = "plain"))
save_fig(p3, "fig_09_weak_iv.svg")
