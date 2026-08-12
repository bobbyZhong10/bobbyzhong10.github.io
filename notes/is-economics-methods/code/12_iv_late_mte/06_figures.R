# ============================================================================
# 06_figures.R  --  two publication-quality SVGs for the chapter
#   (a) fig_12_compliers.svg : compliance types + the first stage
#   (b) fig_12_mte.svg       : the MTE curve over U_D with ATE/ATT/LATE marked
# ============================================================================

suppressMessages({library(ggplot2); library(patchwork); library(svglite)})

d   <- readRDS("meridian_data.rds")
iv  <- readRDS("res_iv.rds")
mte <- readRDS("res_mte.rds")
par <- d$params; truth <- d$truth

fig_dir <- "/Users/bobbyzhong/Library/Mobile Documents/com~apple~CloudDocs/Research/Method/Micro Theory_Microeconometrics/notes_v2/assets/fig"

## ---- palette (warm navy, matches the notes) --------------------------------
navy  <- "#1a365d"; blue <- "#2b6cb0"; red <- "#c53030"
olive <- "#6b7c3a"; grey <- "#718096"; gold <- "#b7791f"
bg    <- "#fafaf9"; ink  <- "#2d3748"
P0 <- par$P0_main; P1 <- par$P1_main
# empirical adoption rates for the first-stage figure: E[D|Z], whose difference
# is the estimated first stage 0.466 (= what the rest of the chapter reports).
# Figure B (MTE) keeps the DESIGN propensities P0/P1 for the complier window.
pD0a <- mean(d$dat$D[d$dat$Z == 0]); pD1a <- mean(d$dat$D[d$dat$Z == 1])

base_theme <- theme_minimal(base_size = 12) +
  theme(plot.background  = element_rect(fill = bg, colour = NA),
        panel.background = element_rect(fill = bg, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e3dd", linewidth = 0.3),
        plot.title    = element_text(colour = navy, face = "bold", size = 13),
        plot.subtitle = element_text(colour = grey, size = 10),
        axis.title    = element_text(colour = ink),
        axis.text     = element_text(colour = ink),
        legend.position = "top",
        legend.title  = element_blank(),
        legend.text   = element_text(colour = ink, size = 10))

# ============================================================================
# FIGURE A : compliance types and the first stage
# ============================================================================

## Panel A1: adoption rate by encouragement, with the complier gap ------------
barA <- data.frame(
  Z = factor(c("No push (Z=0)", "Push (Z=1)"),
             levels = c("No push (Z=0)", "Push (Z=1)")),
  rate = c(pD0a, pD1a))
pA1 <- ggplot(barA, aes(Z, rate)) +
  geom_col(width = 0.55, fill = c(navy, blue)) +
  geom_text(aes(label = sprintf("%.3f", rate)), vjust = -0.6,
            colour = ink, size = 4) +
  annotate("segment", x = 2, xend = 2, y = pD0a, yend = pD1a,
           colour = red, linewidth = 0.8,
           arrow = arrow(ends = "both", length = unit(0.15, "cm"))) +
  annotate("text", x = 1.62, y = (pD0a + pD1a) / 2,
           label = sprintf("first stage\n= complier share\n= %.3f", pD1a - pD0a),
           colour = red, size = 3.4, hjust = 1, lineheight = 0.95) +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(title = "The first stage",
       subtitle = "Encouragement lifts adoption from E[D|Z=0] to E[D|Z=1]",
       x = NULL, y = "adoption rate  E[D | Z]") +
  base_theme

## Panel A2: the resistance line partitioned into types ----------------------
seg <- data.frame(
  xmin = c(0, pD0a, pD1a), xmax = c(pD0a, pD1a, 1),
  type = factor(c("always-takers", "compliers", "never-takers"),
                levels = c("always-takers", "compliers", "never-takers")),
  share = c(pD0a, pD1a - pD0a, 1 - pD1a))
seg$mid <- (seg$xmin + seg$xmax) / 2
seg$label <- c("always-\ntakers", "compliers", "never-\ntakers")
pA2 <- ggplot(seg) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = type),
            colour = "white", linewidth = 0.6) +
  geom_text(aes(x = mid, y = 0.60, label = label), colour = "white",
            size = 3.2, fontface = "bold", lineheight = 0.9) +
  geom_text(aes(x = mid, y = 0.34, label = sprintf("%.1f%%", 100 * share)),
            colour = "white", size = 3.1) +
  scale_fill_manual(values = c("always-takers" = navy,
                               "compliers"     = red,
                               "never-takers"  = grey), guide = "none") +
  scale_x_continuous(breaks = round(c(0, pD0a, pD1a, 1), 2), limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Who gets moved? The unobserved resistance line",
       subtitle = expression("adopt if resistance " * U[D] * " < propensity;  encouragement raises the cutoff from P(0) to P(1)"),
       x = expression("unobserved resistance  " * U[D]), y = NULL) +
  base_theme +
  theme(axis.text.y = element_blank(), panel.grid = element_blank())

figA <- pA1 + pA2 + plot_layout(widths = c(1, 1.6))
ggsave(file.path(fig_dir, "fig_12_compliers.svg"), figA,
       width = 10.5, height = 4.2, device = svglite::svglite)

# ============================================================================
# FIGURE B : the MTE curve over U_D, with ATE / ATT / LATE
# ============================================================================
g <- mte$grid
curves <- rbind(
  data.frame(u = g, mte = mte$mte_true,    est = "true MTE = 6 - 8u"),
  data.frame(u = g, mte = mte$mte_cf,       est = "control-function estimate"),
  data.frame(u = g, mte = mte$mte_A_line,   est = "local-IV estimate (noisy)"))
curves$est <- factor(curves$est,
  levels = c("true MTE = 6 - 8u", "control-function estimate", "local-IV estimate (noisy)"))

ate  <- mte$ate_int
late <- mte$late_int
att  <- mte$att_int
u_att <- (mte$cf_t0 - att) / mte$cf_t1     # where MTE = ATT (weight centroid)

pB <- ggplot() +
  # complier window shading
  annotate("rect", xmin = P0, xmax = P1, ymin = -Inf, ymax = Inf,
           fill = red, alpha = 0.06) +
  geom_hline(yintercept = 0, colour = grey, linewidth = 0.3) +
  # MTE curves
  geom_line(data = curves, aes(u, mte, colour = est, linetype = est,
                               linewidth = est)) +
  # ATE : uniform-weighted average, full-width dashed line
  geom_segment(aes(x = 0, xend = 1, y = ate, yend = ate),
               colour = gold, linewidth = 0.7, linetype = "dashed") +
  annotate("text", x = 0.02, y = ate + 0.45, hjust = 0, colour = gold,
           size = 3.5, label = sprintf("ATE = %.2f  (uniform average)", ate)) +
  # LATE : average over the complier window
  geom_segment(aes(x = P0, xend = P1, y = late, yend = late),
               colour = red, linewidth = 1.3) +
  annotate("text", x = (P0 + P1) / 2, y = late + 0.48, hjust = 0.5, colour = red,
           size = 3.5, label = sprintf("LATE = %.2f  (complier window)", late)) +
  # ATT : weighted toward low resistance
  annotate("point", x = u_att, y = att, colour = navy, size = 2.6) +
  annotate("segment", x = 0.29, xend = u_att - 0.01, y = 4.35, yend = att + 0.06,
           colour = navy, linewidth = 0.3) +
  annotate("text", x = 0.02, y = 4.35, hjust = 0, colour = navy,
           size = 3.5, label = sprintf("ATT = %.2f  (weights low U_D)", att)) +
  scale_colour_manual(values = c("true MTE = 6 - 8u" = navy,
                                 "control-function estimate" = blue,
                                 "local-IV estimate (noisy)" = grey)) +
  scale_linetype_manual(values = c("true MTE = 6 - 8u" = "solid",
                                   "control-function estimate" = "solid",
                                   "local-IV estimate (noisy)" = "dashed")) +
  scale_linewidth_manual(values = c("true MTE = 6 - 8u" = 1.5,
                                    "control-function estimate" = 0.9,
                                    "local-IV estimate (noisy)" = 0.7)) +
  scale_x_continuous(breaks = round(c(0, P0, 0.5, P1, 1), 2),
                     limits = c(0, 1), expand = expansion(mult = c(0.005, 0.01))) +
  labs(title = "Marginal treatment effect over unobserved resistance",
       subtitle = "sophisticated firms (low resistance) gain most; reluctant firms beyond u ~ 0.75 lose. ATE, ATT, LATE are weighted integrals of one curve.",
       x = expression("unobserved resistance to adoption  " * U[D]),
       y = "effect on engagement  MTE(u)") +
  base_theme +
  theme(legend.position = c(0.80, 0.90),
        legend.background = element_rect(fill = bg, colour = "#e6e3dd"))

ggsave(file.path(fig_dir, "fig_12_mte.svg"), pB,
       width = 9.5, height = 5.6, device = svglite::svglite)

if (nzchar(Sys.getenv("PNG_INSPECT"))) {
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "A.png"), figA,
         width = 10.5, height = 4.2, dpi = 130, device = ragg::agg_png)
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "B.png"), pB,
         width = 9.5, height = 5.6, dpi = 130, device = ragg::agg_png)
}

cat("Wrote:\n  ", file.path(fig_dir, "fig_12_compliers.svg"), "\n  ",
    file.path(fig_dir, "fig_12_mte.svg"), "\n")
cat(sprintf("Figure markers: ATE=%.2f LATE=%.2f ATT=%.2f (u_att=%.2f) window=[%.2f,%.2f]\n",
            ate, late, att, u_att, P0, P1))
