# ------------------------------------------------------------------
# 05_figures.R -- Chapter figures (SVG) for Module 11.
#   fig_13_eventstudy.svg : naive TWFE ES vs Sun-Abraham vs truth
#   fig_13_raw_trends.svg : raw mean log GMV by cohort, calendar time
# ------------------------------------------------------------------
library(data.table)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
fig_dir  <- normalizePath(file.path(code_dir, "..", "..", "assets", "fig"))

dt     <- readRDS(file.path(code_dir, "kestrel_panel.rds"))
truth  <- readRDS(file.path(code_dir, "truth.rds"))
r_twfe <- readRDS(file.path(code_dir, "res_twfe.rds"))
r_mod  <- readRDS(file.path(code_dir, "res_modern.rds"))

## ---- Figure (a): event-study comparison --------------------------
e_lo <- -12; e_hi <- 18

truth_es <- rbind(
  data.table(e = e_lo:-1, estimate = 0, se = NA_real_),
  truth$true_theta_e[event_time <= e_hi,
                     .(e = event_time, estimate = theta, se = NA_real_)]
)
truth_es[, series := "True dynamic effect"]

naive_es <- copy(r_twfe$naive_es)[e >= e_lo & e <= e_hi]
naive_es[, series := "Naive TWFE event study"]
naive_es <- rbind(naive_es[, .(e, estimate, se, series)],
                  data.table(e = -1L, estimate = 0, se = NA_real_,
                             series = "Naive TWFE event study"))

sa_es <- copy(r_mod$sa_es)[e >= e_lo & e <= e_hi]
sa_es[, series := "Sun-Abraham (sunab)"]
sa_es <- rbind(sa_es[, .(e, estimate, se, series)],
               data.table(e = -1L, estimate = 0, se = NA_real_,
                          series = "Sun-Abraham (sunab)"))

es <- rbind(truth_es[, .(e, estimate, se, series)], naive_es, sa_es)
es[, series := factor(series, levels = c("True dynamic effect",
                                         "Naive TWFE event study",
                                         "Sun-Abraham (sunab)"))]
es[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]

pd <- position_dodge(width = 0.55)
p1 <- ggplot(es[order(e)], aes(x = e, y = estimate, colour = series, shape = series)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = -0.5, linewidth = 0.3, colour = "grey40", linetype = "dashed") +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = pd, width = 0,
                linewidth = 0.35, na.rm = TRUE) +
  geom_line(data = es[series == "True dynamic effect"],
            aes(group = series), linewidth = 0.5, alpha = 0.7) +
  geom_point(position = pd, size = 1.6, na.rm = TRUE) +
  scale_colour_manual(values = c("True dynamic effect"    = "grey25",
                                 "Naive TWFE event study" = "#D55E00",
                                 "Sun-Abraham (sunab)"    = "#0072B2")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_x_continuous(breaks = seq(e_lo, e_hi, 3)) +
  labs(x = "Event time (months since commission cut)",
       y = "Effect on log GMV",
       colour = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_13_eventstudy.svg"), p1,
       width = 8, height = 4.8, device = svglite::svglite)

## ---- Figure (b): raw cohort trends -------------------------------
trends <- dt[, .(mean_lgmv = mean(lgmv)), by = .(cohort, month)]
trends[, cohort := factor(cohort, levels = c("g13", "g19", "g25", "never"),
                          labels = c("Cohort g = 13", "Cohort g = 19",
                                     "Cohort g = 25", "Never treated"))]
vlines <- data.table(g = c(13, 19, 25),
                     cohort = factor(c("Cohort g = 13", "Cohort g = 19", "Cohort g = 25"),
                                     levels = levels(trends$cohort)))

p2 <- ggplot(trends, aes(x = month, y = mean_lgmv, colour = cohort)) +
  geom_vline(data = vlines, aes(xintercept = g - 0.5, colour = cohort),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  scale_colour_manual(values = c("Cohort g = 13" = "#D55E00",
                                 "Cohort g = 19" = "#0072B2",
                                 "Cohort g = 25" = "#009E73",
                                 "Never treated" = "grey40")) +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  labs(x = "Calendar month",
       y = "Mean log GMV (city average)",
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_13_raw_trends.svg"), p2,
       width = 8, height = 4.8, device = svglite::svglite)

cat("Saved figures:\n")
cat(file.path(fig_dir, "fig_13_eventstudy.svg"), "\n")
cat(file.path(fig_dir, "fig_13_raw_trends.svg"), "\n")
