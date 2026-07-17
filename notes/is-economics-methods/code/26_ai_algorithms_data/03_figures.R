# ------------------------------------------------------------------
# 03_figures.R -- Chapter 22 figures (SVG).
#   fig_26_genreg.svg    : three estimates of the LLM-measured effect
#                          (naive attenuated / gold-only / corrected).
#   fig_26_collusion.svg : Q-learning price path converging above the
#                          Bertrand-Nash price toward monopoly.
# ------------------------------------------------------------------
library(data.table)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
fig_dir  <- normalizePath(file.path(code_dir, "..", "..", "assets", "fig"))
gr  <- readRDS(file.path(code_dir, "res26_genreg.rds"))
co  <- readRDS(file.path(code_dir, "res26_collusion.rds"))

okabe <- c(orange = "#D55E00", blue = "#0072B2", green = "#009E73", grey = "grey40")

## ---- Figure 1: generated-regressor estimates ---------------------
est <- data.table(
  method = factor(c("Naive OLS on LLM score", "Gold-standard only (n=400)",
                    "Reliability-corrected"),
                  levels = rev(c("Naive OLS on LLM score", "Gold-standard only (n=400)",
                                 "Reliability-corrected"))),
  beta = c(gr$b_naive, gr$b_gold, gr$b_eiv),
  se   = c(gr$se_naive, gr$se_gold, gr$se_eiv))
est[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]

p1 <- ggplot(est, aes(beta, method)) +
  geom_vline(xintercept = gr$b1, colour = okabe["grey"],
             linetype = "dashed", linewidth = 0.5) +
  geom_text(data = data.table(beta = gr$b1,
              method = factor("Reliability-corrected", levels = levels(est$method))),
            label = "true beta = 2.0", hjust = -0.05, vjust = 2.0,
            colour = okabe["grey"], size = 3.2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.16,
                colour = okabe["blue"], linewidth = 0.5) +
  geom_point(size = 2.6, colour = okabe["blue"]) +
  geom_text(aes(label = sprintf("%.2f", beta)), vjust = -1.1, size = 3.3) +
  scale_x_continuous(limits = c(0.85, 2.45)) +
  labs(x = "Estimated effect of sentiment on sales  beta", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_26_genreg.svg"), p1,
       width = 7.4, height = 3.3, device = svglite::svglite)

## ---- Figure 2: Q-learning price path -----------------------------
## downsample + rolling mean of the average price path
path <- co$path1
idx <- seq(1, length(path), by = 500)
roll <- function(x, k = 40) as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))
dt <- data.table(t = idx, p = roll(path[idx]))
dt <- dt[!is.na(p)]

p2 <- ggplot(dt, aes(t / 1000, p)) +
  geom_hline(yintercept = co$pMon, colour = okabe["orange"],
             linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = co$pNash, colour = okabe["green"],
             linetype = "dashed", linewidth = 0.5) +
  geom_line(colour = okabe["blue"], linewidth = 0.5) +
  annotate("text", x = max(dt$t) / 1000 * 0.62, y = co$pMon + 0.012,
           label = "monopoly price", colour = okabe["orange"], size = 3.2, hjust = 0) +
  annotate("text", x = max(dt$t) / 1000 * 0.62, y = co$pNash - 0.012,
           label = "Bertrand-Nash price", colour = okabe["green"], size = 3.2, hjust = 0) +
  labs(x = "Learning periods (thousands)",
       y = "Average price of the two Q-learners") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_26_collusion.svg"), p2,
       width = 7.4, height = 4.2, device = svglite::svglite)

cat("Saved fig_26_genreg.svg, fig_26_collusion.svg\n")
cat(sprintf("mean collusion Delta = %.3f (range %.3f-%.3f)\n",
            co$Delta_mean, min(co$Delta_vec), max(co$Delta_vec)))
