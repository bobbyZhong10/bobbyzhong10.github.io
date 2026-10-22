# ------------------------------------------------------------------
# 06_figures.R -- chapter 4 figures (SVG), same style as 01_consumer_producer_welfare/05_figures.R.
# Ruling P2: written to figures/ in this directory only (the text tasks copy them into assets/fig/).
# Ruling P24: they replace the hand-drawn legacy fig_04_monopoly.svg and fig_04_cournot.svg.
#   fig_04_monopoly.svg       block C0: demand p = 12 - q, MR, MC = 2; consumer surplus, profit and deadweight loss
#   fig_04_cournot.svg        block C0: Cirrus's and Garnet's best responses, the Cournot equilibrium, the
#                             collusive and competitive outputs
#   fig_04_passthrough.svg    pass-through of a common cost shock against the number of firms: Cournot under three
#                             demand curvatures (C0) and logit under two conducts (C1 parameters, 02_logit.R)
#   fig_04_personalized.svg   block C2: profit and consumer surplus under uniform, coarse and fine pricing against
#                             the surplus triangle of Bergemann, Brooks and Morris (2015); coarse and fine lie
#                             inside it, the one price (set on the training half) just below its floor
# Reads truth04.rds, res_logit.rds, res_pd.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

tr <- readRDS("truth04.rds")
c1 <- readRDS("res_logit.rds")
c2 <- readRDS("res_pd.rds")
C0 <- tr$C0

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"; col_purple <- "#CC79A7"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}
num <- function(x) format(round(x, 2), nsmall = 0, trim = TRUE)

## ---- Figure (a): block C0, single-price monopoly ----------------------------------------------------
a <- C0$a; cst <- C0$cost; qm <- C0$qm; pm <- C0$pm; qc <- C0$qc
lab_a <- c(d = sprintf("Demand p = %g - q", a), mr = sprintf("Marginal revenue MR = %g - 2q", a),
           mc = sprintf("Marginal cost MC = %g", cst))
lines_a <- rbind(data.frame(q = c(0, a), p = c(a, 0), series = lab_a[["d"]]),
                 data.frame(q = c(0, a / 2), p = c(a, 0), series = lab_a[["mr"]]),
                 data.frame(q = c(0, a), p = c(cst, cst), series = lab_a[["mc"]]))
lines_a$series <- factor(lines_a$series, levels = unname(lab_a))
areas_a <- rbind(data.frame(q = c(0, 0, qm), p = c(pm, a, pm), area = "cs"),
                 data.frame(q = c(0, qm, qm, 0), p = c(cst, cst, pm, pm), area = "profit"),
                 data.frame(q = c(qm, qm, qc), p = c(pm, cst, cst), area = "dwl"))
fills_a <- c(cs = col_blue, profit = col_green, dwl = col_orange)
p1 <- ggplot() +
  geom_polygon(data = areas_a, aes(x = q, y = p, fill = area, group = area), alpha = 0.22, show.legend = FALSE) +
  geom_segment(aes(x = c(qm, qc), xend = c(qm, qc), y = c(0, 0), yend = c(pm, cst)), linetype = "dotted",
               colour = col_true, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = qm, y = pm, yend = pm), linetype = "dotted", colour = col_true, linewidth = 0.4) +
  geom_line(data = lines_a, aes(x = q, y = p, colour = series, linetype = series), linewidth = 0.8) +
  geom_point(aes(x = c(qm, qc), y = c(pm, cst)), size = 2.2, colour = col_true) +
  annotate("text", x = qm + 0.2, y = pm + 0.35, hjust = 0, size = 3.2, colour = col_true,
           label = sprintf("Monopoly: q = %g, p = %g", qm, pm)) +
  annotate("text", x = qc + 0.2, y = cst + 0.45, hjust = 0, size = 3.2, colour = col_true,
           label = sprintf("Competition: q = %g, p = %g", qc, cst)) +
  annotate("text", x = 3.4, y = 7.8, size = 3.2, colour = col_true, label = sprintf("CS = %s", num(C0$cs_m))) +
  annotate("text", x = qm / 2, y = (pm + cst) / 2, size = 3.2, colour = col_true, label = sprintf("Profit = %s", num(C0$profit_m))) +
  annotate("text", x = qm + (qc - qm) / 3, y = cst + (pm - cst) / 3, size = 3.2, colour = col_true,
           label = sprintf("DWL = %s", num(C0$dwl_m))) +
  scale_colour_manual(values = setNames(c(col_true, col_purple, col_blue), unname(lab_a))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "solid"), unname(lab_a))) +
  scale_fill_manual(values = fills_a) +
  scale_x_continuous(limits = c(0, 13.5), breaks = seq(0, 12, 2), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 12.5), breaks = seq(0, 12, 2), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Operating activity q", y = "Price p (USD per unit of activity)", colour = NULL, linetype = NULL) +
  base_theme
save_fig(p1, "fig_04_monopoly.svg")

## ---- Figure (b): block C0, Cournot best responses -------------------------------------------------
s <- seq(0, 10, length.out = 2)
lab_b <- c(br1 = sprintf("Cirrus's best response q1 = (%g - q2) / 2", a - cst),
           br2 = sprintf("Garnet's best response q2 = (%g - q1) / 2", a - cst),
           mono = sprintf("Joint monopoly output q1 + q2 = %g", qm),
           comp = sprintf("Competitive output q1 + q2 = %g", qc))
lines_b <- rbind(data.frame(q1 = (a - cst - s) / 2, q2 = s, series = lab_b[["br1"]]),
                 data.frame(q1 = s, q2 = (a - cst - s) / 2, series = lab_b[["br2"]]),
                 data.frame(q1 = c(0, qm), q2 = c(qm, 0), series = lab_b[["mono"]]),
                 data.frame(q1 = c(0, qc), q2 = c(qc, 0), series = lab_b[["comp"]]))
lines_b$series <- factor(lines_b$series, levels = unname(lab_b))
pts_b <- data.frame(q1 = c(C0$cournot_q, C0$collusion_q_each, C0$competitive_q_each),
                    q2 = c(C0$cournot_q, C0$collusion_q_each, C0$competitive_q_each),
                    lab = c(sprintf("Cournot-Nash: q1 = q2 = 10/3, p = %s", num(C0$cournot_p)),
                            sprintf("Collusion: q1 = q2 = %g, p = %g", C0$collusion_q_each, C0$pm),
                            sprintf("Competition: q1 = q2 = %g, p = %g", C0$competitive_q_each, C0$cost)),
                    lx = c(4.1, 0.15, 5.6), ly = c(3.9, 0.85, 5.45),          # label anchors (left middle)
                    sx = c(4.1, 1.4, 5.6), sy = c(3.9, 1.0, 5.45))            # where each leader line ends
p2 <- ggplot() +
  geom_line(data = lines_b, aes(x = q1, y = q2, colour = series, linetype = series), linewidth = 0.8) +
  geom_point(data = pts_b, aes(x = q1, y = q2), size = 2.4, colour = col_true) +
  geom_segment(data = pts_b, aes(x = q1, y = q2, xend = sx, yend = sy), colour = "grey55", linewidth = 0.3) +
  geom_label(data = pts_b, aes(x = lx, y = ly, label = lab), hjust = 0, size = 3.2, colour = col_true, fill = "white",
             linewidth = 0, label.padding = grid::unit(0.12, "lines")) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange, "grey55", "grey55"), unname(lab_b))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed", "dotted"), unname(lab_b))) +
  scale_x_continuous(limits = c(0, 10.5), breaks = 0:10, expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 10.5), breaks = 0:10, expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cirrus's quantity q1", y = "Garnet's quantity q2", colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
  base_theme
save_fig(p2, "fig_04_cournot.svg")

## ---- Figure (c): pass-through, curvature and conduct ------------------------------------------------
pc <- c1$pt_curves
cp <- c1$pt_curves_par
panel_q <- "Quantity competition (Cournot): demand curvature"
panel_p <- "Price competition (logit): conduct"
lab_c <- c(cournot_linear = "Cournot, linear demand: n / (n + 1)",
           cournot_exp = "Cournot, exponential demand: 1",
           cournot_ce2 = sprintf("Cournot, constant elasticity %g: %gn / (%gn - 1)", cp$ce_e, cp$ce_e, cp$ce_e),
           logit_bertrand = sprintf("Logit (quality %g, \u03b1 = %g, cost %g), Bertrand (\u03ba = 0)", cp$delta, cp$alpha, cp$c),
           logit_joint = "Logit, one owner of all products (\u03ba = 1): equals s0")
pc$series <- factor(lab_c[pc$model], levels = unname(lab_c))
pc$panel <- factor(ifelse(grepl("^cournot", pc$model), panel_q, panel_p), levels = c(panel_q, panel_p))
p3 <- ggplot(pc, aes(x = n, y = pt, colour = series, linetype = series)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = col_true) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = "dotted", colour = col_true) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  facet_wrap(~panel, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = setNames(c(col_blue, col_true, col_orange, col_green, col_purple), unname(lab_c))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "solid", "solid", "dashed"), unname(lab_c))) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Number of symmetric firms n", y = "Pass-through of a common cost shock dp/dc", colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
  base_theme + theme(strip.text = element_text(size = 9.5))
save_fig(p3, "fig_04_passthrough.svg")

## ---- Figure (d): block C2, coarse and fine inside the BBM surplus triangle, one price below its floor ----
b <- c2$bbm
E <- b$efficient; U <- b$uniform_profit
tri <- data.frame(cs = c(0, E - U, 0), ps = c(U, U, E))
lad <- c2$ladder
reg <- data.frame(cs = unname(c2$cs), ps = unname(c2$profit),
                  lab = c(sprintf("One price (K = 1): profit %s, CS %s", num(c2$profit[["uniform"]]), num(c2$cs[["uniform"]])),
                          sprintf("Coarse, K = %d segments: profit %s, CS %s", c2$K[["coarse"]], num(c2$profit[["coarse"]]), num(c2$cs[["coarse"]])),
                          sprintf("Fine, K = %d segments: profit %s, CS %s", c2$K[["fine"]], num(c2$profit[["fine"]]), num(c2$cs[["fine"]]))),
                  lx = c(14.2, 14.2, 14.2), ly = c(11.2, 17.6, 19.4))
reg$lab <- factor(reg$lab, levels = reg$lab)
p4 <- ggplot() +
  geom_polygon(data = tri, aes(x = cs, y = ps), fill = "grey85", alpha = 0.6) +
  geom_path(data = rbind(tri, tri[1, ]), aes(x = cs, y = ps), colour = col_true, linewidth = 0.4) +
  annotate("text", x = 0.8, y = 21.5, hjust = 0, size = 3.1, colour = col_true,
           label = "Reachable by some\nsegmentation (BBM 2015)") +
  annotate("point", x = c(0, E - U), y = c(E, U), size = 2.4, colour = col_true, shape = 17) +
  annotate("text", x = 0.4, y = E, hjust = 0, vjust = 0.4, size = 3.1, colour = col_true,
           label = sprintf("First-degree discrimination: profit = efficient surplus = %s", num(E))) +
  annotate("text", x = E - U + 0.5, y = U + 1.1, hjust = 0, size = 3.1, colour = col_true,
           label = sprintf("Consumer-optimal segmentation (BBM):\nCS = %s, profit = %s", num(E - U), num(U))) +
  geom_point(data = lad, aes(x = cs, y = profit), size = 1.4, colour = "grey55") +
  geom_segment(data = reg, aes(x = cs, y = ps, xend = lx, yend = ly), colour = "grey55", linewidth = 0.3) +
  geom_point(data = reg, aes(x = cs, y = ps, colour = lab), size = 2.6) +
  geom_label(data = reg, aes(x = lx, y = ly, label = lab, colour = lab), hjust = 0, size = 3, linewidth = 0.2,
             show.legend = FALSE) +
  scale_colour_manual(values = setNames(c(col_true, col_blue, col_orange), levels(reg$lab))) +
  scale_x_continuous(limits = c(0, 31), breaks = seq(0, 30, 5), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(10, 34), breaks = seq(10, 35, 5)) +
  labs(x = "Consumer surplus (USD per merchant per month)", y = "Cirrus's profit (USD per merchant per month)",
       colour = NULL,
       caption = sprintf(paste0("Evaluation half, %s merchants. Grey dots: K = %s segments.\n",
                                "Coarse and fine lie inside the triangle. The one-price point lies %s USD below its floor:\n",
                                "its price (%g USD) was set on the training half; the evaluation half's best single price is %g USD."),
                         format(c2$n_eval, big.mark = ","), paste(lad$K, collapse = ", "), num(b$uniform_gap),
                         c2$prices$uniform, b$uniform_price)) +
  guides(colour = "none") +
  base_theme + theme(plot.caption = element_text(size = 8.5, colour = col_true, hjust = 0))
save_fig(p4, "fig_04_personalized.svg")
