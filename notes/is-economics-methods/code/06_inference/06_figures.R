# ============================================================================
# 06_figures.R  --  two SVGs for chapter 6
#   (a) fig_06_delta_bootstrap.svg : WTP ratio, Delta-normal vs bootstrap dist,
#       strong denominator (they agree) vs weak denominator (bootstrap skews)
#   (b) fig_06_cluster.svg : the false-positive-rate punchline + CI widths
# ============================================================================

suppressMessages({library(ggplot2); library(patchwork); library(svglite)})
H  <- readRDS("halcyon.rds"); truth <- H$truth
dl <- readRDS("res_delta.rds"); cl <- readRDS("res_cluster.rds")

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

# ---- FIGURE A : Delta vs bootstrap, strong vs weak denominator -------------
panel_ratio <- function(res, title, sub, xlim) {
  ws <- res$wstar; ws <- ws[ws >= xlim[1] & ws <= xlim[2]]
  xs <- seq(xlim[1], xlim[2], length = 400)
  dn <- data.frame(x = xs, d = dnorm(xs, res$wtp_hat, res$se_delta))
  ggplot() +
    geom_histogram(data = data.frame(w = ws), aes(w, after_stat(density)),
                   bins = 50, fill = "#d9d5cc", colour = NA) +
    geom_line(data = dn, aes(x, d), colour = navy, linewidth = 1) +
    geom_vline(xintercept = truth$WTP, colour = red, linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = truth$WTP, y = 0, label = "truth 0.20", colour = red,
             size = 3.1, vjust = -0.5, hjust = -0.05) +
    labs(title = title, subtitle = sub, x = "WTP = b_time / b_fee", y = "density") +
    coord_cartesian(xlim = xlim) + base_theme
}
pA1 <- panel_ratio(dl$main, "Strong denominator: Delta = bootstrap",
                   "navy = Delta normal; bars = bootstrap draws", c(0.16, 0.26))
pA2 <- panel_ratio(dl$weak, "Weak denominator: bootstrap skews right",
                   "the symmetric Delta normal misses the long tail", c(-0.3, 1.4))
figA <- pA1 + pA2
ggsave(file.path(fig_dir, "fig_06_delta_bootstrap.svg"), figA,
       width = 10.5, height = 4.3, device = svglite::svglite)

# ---- FIGURE B : clustering -- false-positive rates + CI widths -------------
rej <- data.frame(
  method = factor(rep(c("naive i.i.d.", "cluster (z)", "cluster t(G-1)"), 2),
                  levels = c("naive i.i.d.", "cluster (z)", "cluster t(G-1)")),
  G = factor(rep(c("G = 40 cities", "G = 10 cities"), each = 3),
             levels = c("G = 40 cities", "G = 10 cities")),
  rate = c(cl$s40["iid"], cl$s40["cr1_z"], cl$s40["cr1_t"],
           cl$s10["iid"], cl$s10["cr1_z"], cl$s10["cr1_t"]))
pB1 <- ggplot(rej, aes(method, rate, fill = method)) +
  geom_col(width = 0.65) + facet_wrap(~G) +
  geom_hline(yintercept = 0.05, colour = red, linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", rate)), vjust = -0.4, size = 3, colour = ink) +
  scale_fill_manual(values = c("naive i.i.d." = gold, "cluster (z)" = blue,
                               "cluster t(G-1)" = navy), guide = "none") +
  scale_y_continuous(limits = c(0, 0.72), expand = c(0,0)) +
  labs(title = "Rejecting a TRUE null at nominal 5%",
       subtitle = "naive i.i.d. inference fires two-thirds of the time; red line is the 5% target",
       x = NULL, y = "false-positive rate") +
  base_theme + theme(axis.text.x = element_text(size = 8))

## CI half-widths for the single G=40 rollout
ciw <- data.frame(
  method = factor(c("naive i.i.d.", "HC1 robust", "cluster CR1"),
                  levels = c("naive i.i.d.", "HC1 robust", "cluster CR1")),
  se = c(cl$se_iid, cl$se_hc1, cl$se_cr1))
pB2 <- ggplot(ciw, aes(x = method)) +
  geom_hline(yintercept = 0, colour = grey, linewidth = 0.3) +
  geom_errorbar(aes(ymin = cl$tau_hat - 1.96*se, ymax = cl$tau_hat + 1.96*se,
                    colour = method), width = 0.18, linewidth = 1) +
  geom_point(aes(y = cl$tau_hat), colour = ink, size = 2) +
  scale_colour_manual(values = c("naive i.i.d." = gold, "HC1 robust" = olive,
                                 "cluster CR1" = navy), guide = "none") +
  labs(title = "One rollout: the same estimate, three CIs",
       subtitle = "heteroskedasticity-robust HC1 barely differs; only clustering widens it",
       x = NULL, y = expression(hat(tau)%+-%1.96~SE)) +
  base_theme + theme(axis.text.x = element_text(size = 8.5))
figB <- pB1 + pB2 + plot_layout(widths = c(1.35, 1))
ggsave(file.path(fig_dir, "fig_06_cluster.svg"), figB,
       width = 11.0, height = 4.4, device = svglite::svglite)

if (nzchar(Sys.getenv("PNG_INSPECT"))) {
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "A6.png"), figA, width=10.5, height=4.3)
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "B6.png"), figB, width=11.0, height=4.4)
}
cat("Wrote fig_06_delta_bootstrap.svg and fig_06_cluster.svg\n")
