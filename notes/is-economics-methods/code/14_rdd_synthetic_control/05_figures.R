# =====================================================================
# 14_rdd_synthetic_control / 05_figures.R
# Publication-quality SVGs for Module 12:
#   fig_14_rdd.svg   -- binned means + local-linear fits + the jump.
#   fig_14_synth.svg -- treated vs synthetic path + placebo-gaps panel.
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

theme_set(theme_minimal(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 13),
                legend.position = "top"))

# =====================================================================
# FIGURE 1 -- RDD
# =====================================================================
rdd  <- readRDS(file.path(in_dir, "data_rdd.rds"))
d    <- rdd$dat
c_cut <- rdd$meta$c_cut
res_rdd <- readRDS(file.path(in_dir, "res_rdd.rds"))
bw <- res_rdd$bw_cct

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
           label = sprintf("RDD jump = %.3f\n(true 0.200)", res_rdd$rd_conv_est),
           size = 3.6, hjust = 0) +
  scale_colour_manual(values = c(below = "#2c7fb8", above = "#d95f0e"),
                      guide = "none") +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(title = "Preferred badge at the quality-score cutoff (c = 80)",
       subtitle = "Binned means; local-linear fits inside the CCT bandwidth on each side",
       x = "Quality score (running variable)", y = "log GMV") +
  coord_cartesian(xlim = c(c_cut - 18, c_cut + 16))

ggsave(file.path(fig_dir, "fig_14_rdd.svg"), p_rdd,
       width = 7.2, height = 4.6, device = svglite::svglite)
cat("wrote", file.path(fig_dir, "fig_14_rdd.svg"), "\n")

# =====================================================================
# FIGURE 2 -- Synthetic control
# =====================================================================
res_sc <- readRDS(file.path(in_dir, "res_sc.rds"))
traj   <- res_sc$traj
T0     <- res_sc$T0
tall   <- res_sc$traj_all

# ---- panel A: treated vs synthetic path ----------------------------
pathA <- data.frame(
  month = rep(traj$time_unit, 2),
  y     = c(traj$real_y, traj$synth_y),
  series = rep(c("Aster-Metro (observed)", "Synthetic Aster-Metro"),
               each = nrow(traj)))
pA <- ggplot(pathA, aes(month, y, colour = series, linetype = series)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = T0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = T0 - 0.5, y = min(pathA$y), label = "fee reform",
           hjust = 1, size = 3.2, colour = "grey30") +
  scale_colour_manual(values = c("Aster-Metro (observed)" = "#d95f0e",
                                 "Synthetic Aster-Metro" = "#2c7fb8")) +
  scale_linetype_manual(values = c("Aster-Metro (observed)" = "solid",
                                   "Synthetic Aster-Metro" = "22")) +
  labs(title = "Treated metro vs its synthetic control",
       x = NULL, y = "log GMV", colour = NULL, linetype = NULL)

# ---- panel B: gap for treated + in-space placebos ------------------
# Standard Abadie et al. practice: drop placebo donors whose pre-period
# fit is far worse than the treated unit's (their large gaps are just
# bad fit, not evidence). Keep those with pre-RMSPE <= 5x treated's.
tall <- tall %>% mutate(gap = real_y - synth_y,
                        is_treated = .id == "Aster-Metro")
tr_pre_rmspe <- tall %>%
  filter(is_treated, time_unit < T0) %>%
  summarise(r = sqrt(mean(gap^2))) %>% pull(r)
pre_fit <- tall %>% filter(time_unit < T0) %>%
  group_by(.id) %>% summarise(r = sqrt(mean(gap^2)), .groups = "drop")
keep_ids <- pre_fit %>% filter(r <= 5 * tr_pre_rmspe) %>% pull(.id)
tall <- tall %>% filter(.id %in% keep_ids)

pB <- ggplot() +
  geom_line(data = filter(tall, !is_treated),
            aes(time_unit, gap, group = .id),
            colour = "grey75", linewidth = 0.4, alpha = 0.8) +
  geom_line(data = filter(tall, is_treated),
            aes(time_unit, gap), colour = "#d95f0e", linewidth = 1.2) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = T0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = T0 + 0.5, y = max(tall$gap),
           label = "treated gap = effect", hjust = 0, size = 3.2,
           colour = "#d95f0e") +
  labs(title = "Gap: treated (orange) vs donor placebos (grey)",
       x = "Month", y = "log GMV gap")

p_sc <- pA / pB + plot_layout(heights = c(1, 1))
ggsave(file.path(fig_dir, "fig_14_synth.svg"), p_sc,
       width = 7.2, height = 6.4, device = svglite::svglite)
cat("wrote", file.path(fig_dir, "fig_14_synth.svg"), "\n")
