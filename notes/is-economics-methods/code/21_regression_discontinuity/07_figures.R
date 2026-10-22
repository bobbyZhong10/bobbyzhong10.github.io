# =====================================================================
# 21_regression_discontinuity / 07_figures.R
# SVG figures of chapter 21 (ggplot2 + svglite, house style of fig_21_rdd.svg), written to
# ../../assets/fig/:
#   fig_21_rdd.svg       R1: binned means, local-linear fits inside the CCT bandwidth, the
#                        jump (legacy 05_figures.R code; inputs now aster_R1.rds, res_core.rds)
#   fig_21_discrete.svg  R2: means at the 21 support points with the local-linear fits, and
#                        the six 95% intervals for the jump with their Monte Carlo coverage
#   fig_21_kink.svg      R3: the commission schedule and binned log GMV around the USD 1,500
#                        kink, local-quadratic fits, the left fit extended past the kink
#   fig_21_multi.svg     R4: cutoff-specific estimates and the extrapolation to score 75
#                        against the true effect line, and the pooled estimators
#   fig_21_heaping.svg   R1h: the heap at 80.0 in the reported score, naive vs donut estimates
#   fig_21_geo.svg       R4g: the delivery zone with the six boundary points, and the
#                        boundary-point estimates against their truths
# Every label is in English (the figure files are shared with the English edition). Truths,
# cutoffs and design values in labels come from truth21.rds, the block meta or res_*.rds.
# Set PREVIEW_DIR to also write PNG previews there. Run from this directory:
#   Rscript 07_figures.R   (a few seconds)
# =====================================================================

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(patchwork)
}))

in_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(in_dir) == 0 || in_dir == "") in_dir <- "."
fig_dir <- normalizePath(file.path(in_dir, "..", "..", "assets", "fig"),
                         mustWork = FALSE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
preview_dir <- Sys.getenv("PREVIEW_DIR")

theme_set(theme_minimal(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 13),
                legend.position = "top"))

save_fig <- function(p, name, width = 7.2, height = 4.6) {
  ggsave(file.path(fig_dir, name), p, width = width, height = height, device = svglite::svglite)
  if (nzchar(preview_dir))
    ggsave(file.path(preview_dir, sub("\\.svg$", ".png", name)), p, width = width, height = height, dpi = 110)
  cat("wrote", file.path(fig_dir, name), "\n")
}

# =====================================================================
# FIGURE 1 -- RDD (legacy code; only the two input files changed)
# =====================================================================
rdd  <- readRDS(file.path(in_dir, "aster_R1.rds"))
d    <- rdd$dat
c_cut <- rdd$meta$c_cut
res_core <- readRDS(file.path(in_dir, "res_core.rds"))
tr <- readRDS(file.path(in_dir, "truth21.rds"))
bw <- res_core$sharp$h

# binned means (bin width 2 points)
d$bin <- cut(d$score, breaks = seq(45, 100, by = 2), include.lowest = TRUE)
binned <- d %>% group_by(bin) %>%
  summarise(x = mean(score), y = mean(y_sharp), n = n(), .groups = "drop") %>%
  mutate(side = ifelse(x >= c_cut, "above", "below"))

# local-linear fits within the CCT bandwidth on each side
grid_below <- seq(c_cut - bw, c_cut, length.out = 50)
grid_above <- seq(c_cut, c_cut + bw, length.out = 50)
fitL <- lm(y_sharp ~ I(score - c_cut),
           data = subset(d, score >= c_cut - bw & score < c_cut))
fitR <- lm(y_sharp ~ I(score - c_cut),
           data = subset(d, score >= c_cut & score <= c_cut + bw))
lineL <- data.frame(score = grid_below,
                    yhat = predict(fitL, data.frame(score = grid_below)))
lineR <- data.frame(score = grid_above,
                    yhat = predict(fitR, data.frame(score = grid_above)))
jumpL <- predict(fitL, data.frame(score = c_cut))
jumpR <- predict(fitR, data.frame(score = c_cut))

p_rdd <- ggplot() +
  geom_point(data = binned, aes(x, y, size = n, colour = side), alpha = 0.55) +
  geom_line(data = lineL, aes(score, yhat), colour = "#2c7fb8", linewidth = 1.2) +
  geom_line(data = lineR, aes(score, yhat), colour = "#d95f0e", linewidth = 1.2) +
  geom_vline(xintercept = c_cut, linetype = "dashed", colour = "grey40") +
  annotate("segment", x = c_cut, xend = c_cut, y = jumpL, yend = jumpR,
           colour = "black", linewidth = 1.1,
           arrow = arrow(ends = "both", length = unit(0.15, "cm"))) +
  annotate("label", x = c_cut + 6, y = (jumpL + jumpR) / 2,
           label = sprintf("RDD jump = %.3f\n(true %.3f)", res_core$sharp$est, tr$R1$jump),
           size = 3.6, hjust = 0) +
  scale_colour_manual(values = c(below = "#2c7fb8", above = "#d95f0e"),
                      guide = "none") +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(title = sprintf("Preferred badge at the quality-score cutoff (c = %d)", as.integer(c_cut)),
       subtitle = "Binned means; local-linear fits inside the CCT bandwidth on each side",
       x = "Quality score (running variable)", y = "log GMV") +
  coord_cartesian(xlim = c(c_cut - 18, c_cut + 16))

save_fig(p_rdd, "fig_21_rdd.svg")

# =====================================================================
# shared pieces for the new figures
# =====================================================================
col_below <- "#2c7fb8"; col_above <- "#d95f0e"; col_third <- "#1b9e77"
col_cover <- "grey25"; col_miss <- "#7570b3"     # interval plots: covers / misses the truth
pct <- function(v) sprintf("%.1f%%", 100 * v)
# local polynomial fit of order p on one side of c with triangular weights, evaluated on a grid
side_fit <- function(y, x, c, h, side, p = 1, n_grid = 60, extend = 0) {
  keep <- if (side == "left") x < c & x >= c - h else x >= c & x <= c + h
  w <- pmax(1 - abs(x[keep] - c) / h, 0)
  u <- x[keep] - c
  fit <- lm(y[keep] ~ poly(u, p, raw = TRUE), weights = w)
  g <- if (side == "left") seq(c - h, c + extend, length.out = n_grid) else seq(c, c + h, length.out = n_grid)
  data.frame(x = g, yhat = as.numeric(cbind(1, poly(g - c, p, raw = TRUE)) %*% coef(fit)))
}
# horizontal interval plot: one row per estimator, truth as a dashed vertical line, a note
# column on the right (coverage or other)
forest <- function(df, truth, xlab, note_head = NULL, truth_lab = "truth") {
  df$label <- factor(df$label, levels = rev(df$label))
  span <- range(c(df$lo, df$hi, truth))
  x_note <- span[2] + 0.06 * diff(span)
  p <- ggplot(df, aes(y = label)) +
    geom_vline(xintercept = truth, linetype = "dashed", colour = "grey40") +
    geom_segment(aes(x = lo, xend = hi, yend = label, colour = covers), linewidth = 1.1) +
    geom_point(aes(x = est, colour = covers), size = 2.4) +
    scale_colour_manual(values = c(`TRUE` = col_cover, `FALSE` = col_miss),
                        labels = c(`TRUE` = "covers the truth", `FALSE` = "misses the truth"), name = NULL) +
    annotate("text", x = truth, y = nrow(df) + 0.65, label = truth_lab, size = 3.2, colour = "grey30") +
    labs(x = xlab, y = NULL) +
    coord_cartesian(xlim = c(span[1], x_note + 0.22 * diff(span)), ylim = c(0.6, nrow(df) + 0.8), clip = "off")
  if (!is.null(df$note)) {
    p <- p + geom_text(aes(x = x_note, label = note), hjust = 0, size = 3.3, colour = "grey20")
    if (!is.null(note_head))
      p <- p + annotate("text", x = x_note, y = nrow(df) + 0.65, label = note_head, hjust = 0, size = 3.2,
                        fontface = "bold", colour = "grey20")
  }
  p
}
ci_df <- function(label, est, ci, truth, note = NULL) {
  out <- data.frame(label = label, est = est, lo = ci[1], hi = ci[2], covers = ci[1] <= truth & truth <= ci[2],
                    stringsAsFactors = FALSE)
  if (!is.null(note)) out$note <- note
  out
}

# =====================================================================
# FIGURE 2 -- discrete running variable (R2)
# =====================================================================
R2 <- readRDS(file.path(in_dir, "aster_R2.rds"))
rd2 <- readRDS(file.path(in_dir, "res_discrete.rds"))
c2 <- R2$meta$cutoff; tau2 <- R2$meta$jump; h2 <- rd2$conv$h
means2 <- data.frame(x = R2$x, y = R2$y) %>% group_by(x) %>%
  summarise(y = mean(y), n = n(), .groups = "drop") %>%
  mutate(side = ifelse(x >= c2, "above", "below"))
cef2 <- R2$meta$cef %>% mutate(m = m0 + tau, side = ifelse(x >= c2, "above", "below"))
fl2 <- side_fit(R2$y, R2$x, c2, h2, "left"); fr2 <- side_fit(R2$y, R2$x, c2, h2, "right")
xlim2 <- c(3.55, 5.0)
p2a <- ggplot() +
  geom_line(data = subset(cef2, x >= xlim2[1]), aes(x, m, group = side), colour = "grey55", linewidth = 0.6) +
  geom_point(data = subset(cef2, x >= xlim2[1]), aes(x, m), shape = 4, colour = "grey35", size = 2) +
  geom_point(data = subset(means2, x >= xlim2[1]), aes(x, y, size = n, colour = side), alpha = 0.7) +
  geom_line(data = fl2, aes(x, yhat), colour = col_below, linewidth = 1.1) +
  geom_line(data = fr2, aes(x, yhat), colour = col_above, linewidth = 1.1) +
  geom_vline(xintercept = c2, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 3.57, y = max(cef2$m), hjust = 0, vjust = 1, size = 2.7, colour = "grey20",
           label = sprintf("\u00d7 true mean\ndots: sample means\njump %.3f\n(true %.2f)", rd2$conv$est, tau2)) +
  scale_colour_manual(values = c(below = col_below, above = col_above), guide = "none") +
  scale_size_continuous(range = c(1, 4.5), guide = "none") +
  scale_x_continuous(breaks = seq(3.6, 5.0, by = 0.2)) +
  labs(title = NULL, x = sprintf("Average of the last %d review stars", as.integer(R2$meta$par$K)), y = "log GMV")
cov2 <- rd2$mc
rows2 <- rbind(
  ci_df("Conventional (rdrobust)", rd2$conv$est, rd2$conv$ci, tau2, pct(cov2$cover_conv)),
  ci_df("Robust bias-corrected", rd2$robust$est, rd2$robust$ci, tau2, pct(cov2$cover_robust)),
  ci_df("Clustered by support point", rd2$crv$est, rd2$crv$ci, tau2, pct(cov2$cover_crv)),
  ci_df(sprintf("Honest, M = %.2f (true bound)", rd2$honest$M), rd2$honest$est, rd2$honest$ci, tau2, pct(cov2$cover_honest)),
  ci_df(sprintf("Bias-aware, M = %.2f (data lower bound)", rd2$honest_lb$M), rd2$honest_lb$est, rd2$honest_lb$ci, tau2,
        paste0(pct(cov2$cover_honest_lb), "*")),
  ci_df(sprintf("Local randomization [%.1f, %.1f]", rd2$locrand$window[1], rd2$locrand$window[2]), rd2$locrand$est, rd2$locrand$ci, tau2, paste0(pct(cov2$cover_lr), "**")))
p2b <- forest(rows2, tau2, sprintf("95%% interval for the jump at %.2f", c2), note_head = "MC coverage") +
  theme(legend.position = "bottom", plot.margin = margin(5.5, 18, 5.5, 5.5))   # right margin keeps the note column inside the canvas
filter_star <- floor(2 * c2 + 0.5) / 2          # the half-star display level that the cutoff rounds up to
p_disc <- p2a + p2b + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(title = sprintf("Discrete running variable: rating averages at the %.1f-star filter", filter_star),
                  subtitle = sprintf("One sample of %s sellers; coverage from %d simulated samples",
                                     format(length(R2$x), big.mark = ","), cov2$reps),
                  caption = sprintf(paste0("* M re-estimated in each sample.  ** Neyman interval in the fixed window ",
                                           "[%.1f, %.1f]; the plotted interval inverts the randomization test."),
                                    cov2$lr_window[1], cov2$lr_window[2]),
                  theme = theme(plot.title = element_text(face = "bold", size = 13),
                                plot.caption = element_text(size = 8.5, colour = "grey30")))
save_fig(p_disc, "fig_21_discrete.svg", width = 9, height = 4.8)

# =====================================================================
# FIGURE 3 -- regression kink (R3)
# =====================================================================
R3 <- readRDS(file.path(in_dir, "aster_R3.rds"))
rk <- readRDS(file.path(in_dir, "res_kink.rds"))
k0 <- R3$meta$cutoff; hk <- rk$rk$h
win3 <- c(k0 - 500, k0 + 500)                   # both panels show the same window, so the kink lines align
BIN3 <- 25                                      # bin width (USD) of the binned means
usd <- function(v) paste0("USD ", format(v, big.mark = ",", trim = TRUE))
x_scale3 <- scale_x_continuous(breaks = seq(win3[1], win3[2], by = 250), labels = usd)
sched <- data.frame(x = seq(win3[1], win3[2], by = 5)); sched$b <- R3$meta$b_fun(sched$x)
p3a <- ggplot(sched, aes(x, b)) +
  geom_line(colour = "grey20", linewidth = 1) +
  geom_vline(xintercept = k0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = win3[1] + 10, y = min(sched$b) + 0.15 * diff(range(sched$b)), hjust = 0, vjust = 0, size = 3.3,
           label = sprintf("%g%% up to %s, then %.1f pp lower\nper USD 100 (floor %g%% from %s)", R3$meta$rate_base,
                           usd(k0), -100 * R3$meta$slope_right, R3$meta$rate_floor, usd(R3$meta$floor_from))) +
  x_scale3 + coord_cartesian(xlim = win3) +
  labs(x = NULL, y = "Commission rate (%)", title = NULL)
bins3 <- data.frame(x = R3$x, y = R3$y) %>% filter(x >= win3[1], x <= win3[2]) %>%
  mutate(bin = floor((x - k0) / BIN3)) %>% group_by(bin) %>%
  summarise(x = mean(x), y = mean(y), n = n(), .groups = "drop") %>%
  mutate(side = ifelse(x >= k0, "above", "below"))
fl3 <- side_fit(R3$y, R3$x, k0, hk, "left", p = 2, extend = hk)
fr3 <- side_fit(R3$y, R3$x, k0, hk, "right", p = 2)
fl3_in <- subset(fl3, x <= k0); fl3_out <- subset(fl3, x >= k0)
p3b <- ggplot() +
  geom_point(data = bins3, aes(x, y, colour = side), size = 1.6, alpha = 0.75) +
  geom_line(data = fl3_in, aes(x, yhat), colour = col_below, linewidth = 1.1) +
  geom_line(data = fl3_out, aes(x, yhat), colour = col_below, linewidth = 0.8, linetype = "dotted") +
  geom_line(data = fr3, aes(x, yhat), colour = col_above, linewidth = 1.1) +
  geom_vline(xintercept = k0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = k0 + 30, y = min(bins3$y) + 0.02, hjust = 0, vjust = 0, size = 3.3, colour = "grey20",
           label = sprintf("kink estimate %.3f per pp of commission\n(true %.3f); dotted: left fit continued",
                           rk$rk$est, tr$R3$kink_effect)) +
  scale_colour_manual(values = c(below = col_below, above = col_above), guide = "none") +
  x_scale3 + coord_cartesian(xlim = win3) +
  labs(x = "Last quarter's average monthly GMV (running variable)", y = "log GMV this quarter")
p_kink <- p3a / p3b + plot_layout(heights = c(1, 2)) +
  plot_annotation(title = "Regression kink: Aster's volume discount on commissions",
                  subtitle = sprintf("Binned means (USD %d bins); local-quadratic fits inside the bandwidth on each side", BIN3),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig(p_kink, "fig_21_kink.svg", width = 7.2, height = 5.4)

# =====================================================================
# FIGURE 4 -- multiple cutoffs and extrapolation (R4)
# =====================================================================
mm <- readRDS(file.path(in_dir, "res_multi.rds"))
bc <- mm$by_cutoff
# the true effect line through the cutoff effects (linear in gen_R4; the extrapolation truth
# lies on it)
cuts4 <- bc$cutoff; tc4 <- as.numeric(tr$R4$tau_c[as.character(cuts4)])
s4 <- (tc4[2] - tc4[1]) / (cuts4[2] - cuts4[1])
stopifnot(abs(tc4[3] - (tc4[1] + s4 * (cuts4[3] - cuts4[1]))) < 1e-12,
          abs(tr$R4$tau_extrap - (tc4[1] + s4 * (tr$R4$x_extrap - cuts4[1]))) < 1e-12)
tau_line <- data.frame(x = seq(66, 94, by = 0.5)); tau_line$tau <- tc4[1] + s4 * (tau_line$x - cuts4[1])
pts4 <- data.frame(x = bc$cutoff, est = bc$est, lo = vapply(bc$ci, `[`, 1, 1), hi = vapply(bc$ci, `[`, 1, 2),
                   what = "at its cutoff", lab = sprintf("%s\nweight %.2f", bc$category, mm$weights),
                   lx = cuts4 + c(0.8, -0.8, -0.8), ly = c(bc$ci[[1]][2] + 0.01, 0.215, 0.235), hj = c(0, 1, 1),
                   vj = c(0, 0.5, 0))
ex4 <- data.frame(x = tr$R4$x_extrap, est = mm$extrap$est, lo = mm$extrap$ci[1], hi = mm$extrap$ci[2], what = "extrapolated")
p4a <- ggplot() +
  geom_line(data = tau_line, aes(x, tau), colour = "grey55", linewidth = 0.8) +
  geom_segment(data = pts4, aes(x = x, xend = x, y = lo, yend = hi), colour = col_above, linewidth = 1.1) +
  geom_point(data = pts4, aes(x, est), colour = col_above, size = 2.6) +
  geom_segment(data = ex4, aes(x = x, xend = x, y = lo, yend = hi), colour = col_third, linewidth = 1.1) +
  geom_point(data = ex4, aes(x, est), colour = col_third, size = 2.6, shape = 17) +
  geom_text(data = pts4, aes(x = lx, y = ly, label = lab, hjust = hj, vjust = vj), size = 3.1, colour = "grey20") +
  annotate("text", x = 66.5, y = 0.11, hjust = 0, size = 3.1, colour = "grey20",
           label = sprintf("triangle: vintage effect\nextrapolated to score %g\n(parallel Y(0) across categories)",
                           tr$R4$x_extrap)) +
  annotate("text", x = 66.5, y = -0.03, hjust = 0, size = 3.1, colour = "grey35", label = "grey line: true effect at each score") +
  scale_x_continuous(breaks = seq(65, 95, by = 5)) +
  coord_cartesian(xlim = c(66, 95), ylim = c(-0.05, 0.52)) +
  labs(x = "Quality score", y = "Badge effect on log GMV")
rows4 <- rbind(
  ci_df(sprintf("Simple mean of the %d cutoff estimates", nrow(bc)), mm$simple_mean, c(mm$simple_mean, mm$simple_mean), mm$pooled$target),
  ci_df("Lemma 1 weighted estimator", mm$pooled$est, mm$pooled$ci, mm$pooled$target),
  ci_df("Pooled RD on x - cutoff", mm$pooled_norm$est, mm$pooled_norm$ci, mm$pooled$target))
rows4$note <- c(sprintf("%.3f", mm$simple_mean), sprintf("%.3f", mm$pooled$est), sprintf("%.3f", mm$pooled_norm$est))
p4b <- forest(rows4, mm$pooled$target, "Pooled effect", truth_lab = sprintf("weighted target %.3f", mm$pooled$target)) +
  geom_vline(xintercept = mean(tr$R4$tau_c), linetype = "dotted", colour = "grey40") +
  annotate("text", x = mean(tr$R4$tau_c), y = 0.55, label = sprintf("simple mean %.2f", mean(tr$R4$tau_c)),
           size = 3.0, colour = "grey30") +
  theme(legend.position = "bottom")
p_multi <- p4a + p4b + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(title = "Category-specific badge cutoffs: three effects, one pooled number",
                  subtitle = sprintf(paste0("Robust bias-corrected 95%% intervals, except the extrapolation ",
                                            "(conventional, fixed h = %g)\n",
                                            "Weights = share of sellers near their own cutoff"), mm$extrap$h),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig(p_multi, "fig_21_multi.svg", width = 9, height = 4.6)

# =====================================================================
# FIGURE 5 -- heaping and the donut (R1h)
# =====================================================================
H <- readRDS(file.path(in_dir, "aster_R1h.rds"))
xe <- readRDS(file.path(in_dir, "res_extra.rds"))
c5 <- H$meta$cutoff
hist5 <- data.frame(x = H$x) %>% filter(x >= c5 - 6, x <= c5 + 6) %>% count(x) %>%
  mutate(where = ifelse(abs(x - c5) <= xe$donut$radius, "dropped by the donut", "kept"))
p5a <- ggplot(hist5, aes(x, n, fill = where)) +
  geom_col(width = 0.08) +
  geom_vline(xintercept = c5, linetype = "dashed", colour = "grey40") +
  annotate("text", x = c5 + 0.3, y = xe$heap$n_at_80, hjust = 0, vjust = 1, size = 3.2, colour = "grey20",
           label = sprintf("%d sellers at %.1f\ndensity test p = %.2e", xe$heap$n_at_80, c5, xe$heap$density_p)) +
  scale_fill_manual(values = c(`dropped by the donut` = col_above, kept = "grey60"), name = NULL) +
  labs(x = "Reported quality score (one decimal)", y = "Sellers per 0.1 point")
tr_j <- tr$R1h$jump
rows5 <- rbind(
  ci_df("All sellers (naive)", xe$heap_naive$est, xe$heap_naive$ci, tr_j),
  ci_df(sprintf("Donut: drop |x - %g| <= %.1f", c5, xe$donut$radius), xe$donut$est, xe$donut$ci, tr_j),
  ci_df(sprintf("Drop x = %.1f only", c5), xe$donut_heap$est, xe$donut_heap$ci, tr_j))
rows5$note <- sprintf("%.3f", rows5$est)
p5b <- forest(rows5, tr_j, "Badge effect (robust 95% interval)", note_head = "estimate") +
  theme(legend.position = "bottom")
p_heap <- p5a + p5b + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(title = sprintf("Heaping at the cutoff: sellers who push their score to exactly %.1f", c5),
                  subtitle = sprintf("The heapers are weaker sellers sitting at the cutoff; true jump %.2f", tr_j),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig(p_heap, "fig_21_heaping.svg", width = 9, height = 4.6)

# =====================================================================
# FIGURE 6 -- geographic RD (R4g)
# =====================================================================
G <- readRDS(file.path(in_dir, "aster_R4g.rds"))
geo <- mm$geo
set.seed(2106)                                   # which sellers are drawn on the map (display only)
idx <- sample.int(length(G$y), 600)
map6 <- data.frame(x1 = G$x1[idx], x2 = G$x2[idx], zone = ifelse(G$d[idx] == 1, "inside the zone", "outside"))
R6 <- G$meta$radius
circ <- data.frame(t = seq(0, 2 * pi, length.out = 361)); circ$x1 <- R6 * cos(circ$t); circ$x2 <- R6 * sin(circ$t)
p6a <- ggplot() +
  geom_point(data = map6, aes(x1, x2, colour = zone), size = 0.7, alpha = 0.5) +
  geom_path(data = circ, aes(x1, x2), colour = "grey20", linewidth = 0.8) +
  geom_point(data = geo, aes(b1, b2), shape = 21, fill = "white", colour = "black", size = 5, stroke = 0.8) +
  geom_text(data = geo, aes(b1, b2, label = point), size = 3.2, fontface = "bold") +
  annotate("point", x = 0, y = 0, shape = 3, size = 3) +
  annotate("text", x = 0, y = -5, label = "hub", size = 3.0, fontface = "bold") +
  scale_colour_manual(values = c(`inside the zone` = col_above, outside = col_below), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  coord_equal(xlim = c(-40, 40), ylim = c(-40, 40)) +
  labs(x = "km east of the hub", y = "km north of the hub")
stopifnot(isTRUE(all.equal(geo$truth, as.numeric(tr$R4g$tau_b))))
rows6 <- do.call(rbind, lapply(seq_len(nrow(geo)), function(k)
  ci_df(sprintf("Point %d (%d deg)", geo$point[k], as.integer(geo$angle_deg[k])), geo$est[k], geo$ci[[k]], geo$truth[k])))
rows6$truth <- geo$truth
p6b <- ggplot(transform(rows6, label = factor(label, levels = rev(label))), aes(y = label)) +
  geom_segment(aes(x = lo, xend = hi, yend = label, colour = covers), linewidth = 1.1) +
  geom_point(aes(x = est, colour = covers), size = 2.4) +
  geom_point(aes(x = truth), shape = 4, size = 2.6, colour = "black") +
  scale_colour_manual(values = c(`TRUE` = col_cover, `FALSE` = col_miss),
                      labels = c(`TRUE` = "covers the truth", `FALSE` = "misses the truth"), name = NULL) +
  labs(x = "Same-day effect (robust 95% CI; \u00d7 = truth)", y = NULL) +
  theme(legend.position = "bottom")
p_geo <- p6a + p6b + plot_layout(widths = c(1, 1.1)) +
  plot_annotation(title = "Geographic RD at the edge of the same-day delivery zone",
                  subtitle = sprintf("%g km zone; effects rise toward the east; %d of %d intervals cover",
                                     R6, sum(geo$covers), nrow(geo)),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig(p_geo, "fig_21_geo.svg", width = 9, height = 4.8)
