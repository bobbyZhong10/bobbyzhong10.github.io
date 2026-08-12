# ------------------------------------------------------------------
# 03_figures.R -- Chapter 21 figures (SVG).
#   fig_25_prominence.svg : P(buy) by (randomized) position -- the
#                           causal value of prominence / advertising.
#   fig_25_depth.svg      : distribution of search depth (listings
#                           inspected) -- limited consideration.
#   fig_25_estimates.svg  : three price-coefficient estimates vs truth.
# ------------------------------------------------------------------
library(data.table)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
fig_dir  <- normalizePath(file.path(code_dir, "..", "..", "assets", "fig"))
d      <- readRDS(file.path(code_dir, "lumen_log.rds"))
truth  <- readRDS(file.path(code_dir, "truth25.rds"))
res    <- readRDS(file.path(code_dir, "res25_estimate.rds"))
log_dt <- d$log; setDT(log_dt)

okabe <- c(orange = "#D55E00", blue = "#0072B2", green = "#009E73", grey = "grey40")

## ---- Figure 1: prominence (advertising) effect by position -------
prom <- log_dt[, .(buy = mean(bought), insp = mean(inspected)), by = position]
p1 <- ggplot(prom, aes(position)) +
  geom_col(aes(y = insp), fill = "grey80", width = 0.7) +
  geom_line(aes(y = buy), colour = okabe["orange"], linewidth = 0.8) +
  geom_point(aes(y = buy), colour = okabe["orange"], size = 2) +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(name = "Probability",
                     sec.axis = sec_axis(~ ., name = "(bars: P(inspect); line: P(buy))")) +
  labs(x = "Randomized position in the ranking",
       title = NULL) +
  annotate("text", x = 6.2, y = 0.22, label = "P(buy) by position",
           colour = okabe["orange"], size = 3.4, hjust = 0) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_25_prominence.svg"), p1,
       width = 7.4, height = 4.4, device = svglite::svglite)

## ---- Figure 2: search-depth distribution -------------------------
depth <- data.table(n = d$n_insp)
dd <- depth[, .N, by = n][order(n)]
dd[, share := N / sum(N)]
p2 <- ggplot(dd, aes(factor(n), share)) +
  geom_col(fill = okabe["blue"], width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", share)), vjust = -0.5, size = 3.2) +
  geom_vline(xintercept = truth$mean_depth, linetype = "dashed",
             colour = okabe["grey"], linewidth = 0.4) +
  labs(x = "Number of listings inspected (of 8)",
       y = "Share of consumers") +
  annotate("text", x = 5.5, y = 0.55,
           label = sprintf("mean depth = %.2f of 8\nfull-info assumes all 8", truth$mean_depth),
           hjust = 0, size = 3.3, colour = okabe["grey"]) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(fig_dir, "fig_25_depth.svg"), p2,
       width = 7.4, height = 4.4, device = svglite::svglite)

## ---- Figure 3: three price-coefficient estimates -----------------
est <- data.table(
  method = factor(c("Naive full-info logit", "Inspected-set logit", "Structural search (SMM)"),
                  levels = rev(c("Naive full-info logit", "Inspected-set logit", "Structural search (SMM)"))),
  alpha  = c(res$a_naive, res$a_insp, res$a_hat),
  se     = c(res$se_naive, res$se_insp, NA))
est[, `:=`(lo = alpha - 1.96 * se, hi = alpha + 1.96 * se)]

p3 <- ggplot(est, aes(alpha, method)) +
  geom_vline(xintercept = truth$alpha, colour = okabe["grey"],
             linetype = "dashed", linewidth = 0.5) +
  geom_text(data = data.table(alpha = truth$alpha,
              method = factor("Structural search (SMM)", levels = levels(est$method))),
            label = "true alpha = 0.030", hjust = -0.05, vjust = -1.4,
            colour = okabe["grey"], size = 3.2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.15,
                colour = okabe["blue"], linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 2.6, colour = okabe["blue"]) +
  geom_text(aes(label = sprintf("%.4f", alpha)), vjust = -1.1, size = 3.2) +
  scale_x_continuous(limits = c(0.022, 0.038)) +
  labs(x = "Estimated price coefficient  alpha", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_25_estimates.svg"), p3,
       width = 7.4, height = 3.4, device = svglite::svglite)

cat("Saved figures: fig_25_prominence.svg, fig_25_depth.svg, fig_25_estimates.svg\n")
cat(sprintf("P(buy) by position: %s\n", paste(round(prom$buy, 3), collapse = ", ")))
