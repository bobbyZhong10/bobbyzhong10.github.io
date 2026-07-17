# ============================================================================
# 07_figures.R  --  two SVGs for chapter 8
#   (a) fig_08_curve.svg : the demand curve (linear misses, nonparametric
#       recovers) + the bias-variance tradeoff with the CV-picked bandwidth
#   (b) fig_08_semipar.svg : Robinson partially-linear beta + the curse of
#       dimensionality
# ============================================================================

suppressMessages({library(ggplot2); library(patchwork); library(svglite)})
S  <- readRDS("solstice.rds"); truth <- S$truth
d  <- S$demand; grid <- S$grid
rl <- readRDS("res_linear.rds"); rk <- readRDS("res_kernel.rds")
rc <- readRDS("res_cv.rds"); rr <- readRDS("res_robinson.rds"); rcu <- readRDS("res_curse.rds")

fig_dir <- "/Users/bobbyzhong/Library/Mobile Documents/com~apple~CloudDocs/Research/Method/Micro Theory_Microeconometrics/notes_v2/assets/fig"
navy <- "#1a365d"; blue <- "#2b6cb0"; red <- "#c53030"
olive<- "#6b7c3a"; grey <- "#718096"; gold <- "#b7791f"; bg <- "#fafaf9"; ink <- "#2d3748"
base_theme <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = bg, colour = NA),
        panel.background = element_rect(fill = bg, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e3dd", linewidth = 0.3),
        plot.title = element_text(colour = navy, face = "bold", size = 12.5),
        plot.subtitle = element_text(colour = grey, size = 9.5),
        axis.title = element_text(colour = ink), axis.text = element_text(colour = ink),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(colour = ink, size = 9))

# ---- FIGURE A ---------------------------------------------------------------
## A1: the demand curve
lin <- lm(y ~ x, data = d)
curveA <- data.frame(x = grid, truth = truth$g_grid,
                     linear = predict(lin, data.frame(x = grid)),
                     nonpar = rc$sp_grid)
library(tidyr)
cl <- pivot_longer(curveA, -x, names_to = "fit", values_to = "y")
cl$fit <- factor(cl$fit, levels = c("truth","linear","nonpar"),
                 labels = c("true curve","linear fit","nonparametric (GCV spline)"))
pA1 <- ggplot() +
  geom_point(data = d, aes(x, y), colour = "#c9c4b8", size = 0.7, alpha = 0.7) +
  geom_line(data = cl, aes(x, y, colour = fit, linetype = fit), linewidth = 1) +
  scale_colour_manual(values = c("true curve"=red,"linear fit"=gold,
                                 "nonparametric (GCV spline)"=navy)) +
  scale_linetype_manual(values = c("true curve"=1,"linear fit"=2,"nonparametric (GCV spline)"=1)) +
  labs(title = "One line cannot be a demand curve",
       subtitle = "the constant-slope fit misses the flat zones and the steep drop",
       x = "price", y = "log demand") + base_theme

## A2: bias-variance vs bandwidth
bv <- rc$bv; hs <- rc$hs
bvl <- rbind(
  data.frame(h = hs, val = bv["bias2",], part = "bias^2"),
  data.frame(h = hs, val = bv["variance",], part = "variance"),
  data.frame(h = hs, val = bv["mse",], part = "MSE"))
bvl$part <- factor(bvl$part, levels = c("bias^2","variance","MSE"))
pA2 <- ggplot(bvl, aes(h, val, colour = part, linetype = part)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  geom_vline(xintercept = rc$cv_bw, colour = navy, linetype = "dotted", linewidth = 0.7) +
  annotate("text", x = rc$cv_bw, y = max(bvl$val)*0.9,
           label = sprintf("CV picks\nh = %.2f", rc$cv_bw), colour = navy, size = 3, hjust = -0.1, lineheight = 0.9) +
  scale_x_log10() +
  scale_colour_manual(values = c("bias^2"=red,"variance"=blue,"MSE"=navy)) +
  scale_linetype_manual(values = c("bias^2"=2,"variance"=2,"MSE"=1)) +
  labs(title = "The bias-variance tradeoff",
       subtitle = "smoothing trades bias for variance; CV finds the bottom of the U",
       x = "bandwidth h (log scale)", y = "integrated error") + base_theme
figA <- pA1 + pA2
ggsave(file.path(fig_dir, "fig_08_curve.svg"), figA, width = 10.8, height = 4.3, device = svglite::svglite)

# ---- FIGURE B ---------------------------------------------------------------
## B1: Robinson partially-linear beta
dfB1 <- data.frame(
  method = factor(c("y ~ z","y ~ z + x (linear)","y ~ z + poly(x,3)","Robinson"),
                  levels = c("y ~ z","y ~ z + x (linear)","y ~ z + poly(x,3)","Robinson")),
  beta = c(rr$b_naive1, rr$b_naive2, rr$b_naive3, rr$b_rob))
pB1 <- ggplot(dfB1, aes(method, beta, fill = method)) +
  geom_col(width = 0.66) +
  geom_hline(yintercept = 0.80, colour = red, linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", beta)), vjust = -0.4, size = 3.1, colour = ink) +
  annotate("text", x = 0.7, y = 0.83, label = "truth 0.80", colour = red, size = 3, hjust = 0) +
  scale_fill_manual(values = c("y ~ z"=gold,"y ~ z + x (linear)"=olive,
                               "y ~ z + poly(x,3)"=blue,"Robinson"=navy), guide="none") +
  scale_y_continuous(limits = c(0, 0.95), expand = c(0,0)) +
  labs(title = "Partially linear: flexible control recovers beta",
       subtitle = "omitting or linearizing the price effect biases the ad-exposure coefficient",
       x = NULL, y = expression(hat(beta)~"(ad exposure)")) + base_theme +
  theme(axis.text.x = element_text(size = 7.5))

## B2: curse of dimensionality (edge length)
dfB2 <- data.frame(d = rcu$ds, edge = rcu$edge)
pB2 <- ggplot(dfB2, aes(factor(d), edge)) +
  geom_col(width = 0.6, fill = navy) +
  geom_text(aes(label = sprintf("%.2f", edge)), vjust = -0.4, size = 3.2, colour = ink) +
  scale_y_continuous(limits = c(0, 0.92), expand = c(0,0)) +
  labs(title = "The curse of dimensionality",
       subtitle = "edge a neighborhood must span to hold 10% of the data",
       x = "number of regressors d", y = "edge length per axis") + base_theme
figB <- pB1 + pB2 + plot_layout(widths = c(1.25, 1))
ggsave(file.path(fig_dir, "fig_08_semipar.svg"), figB, width = 11.0, height = 4.4, device = svglite::svglite)

if (nzchar(Sys.getenv("PNG_INSPECT"))) {
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "A8.png"), figA, width=10.8, height=4.3)
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "B8.png"), figB, width=11.0, height=4.4)
}
cat("Wrote fig_08_curve.svg and fig_08_semipar.svg\n")
