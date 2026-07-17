# ------------------------------------------------------------------
# 06_figures.R -- Chapter figures (SVG) for Module 13.
#   fig_16_mc.svg     : sampling distribution, naive ML vs DML (MC)
#   fig_16_cate.svg   : GATES calibration, estimated vs true by quintile
#   fig_16_policy.svg : value vs fraction treated (targeting curve)
# ------------------------------------------------------------------
suppressMessages({library(data.table); library(ggplot2)})

code_dir <- dirname(normalizePath(sub("--file=", "",
              grep("--file=", commandArgs(FALSE), value = TRUE))))
fig_dir  <- normalizePath(file.path(code_dir, "..", "..", "assets", "fig"))

mc    <- as.data.table(readRDS(file.path(code_dir, "mc_frame.rds")))
cf    <- readRDS(file.path(code_dir, "res_cf.rds"))
cfr   <- as.data.table(readRDS(file.path(code_dir, "cate_frame.rds")))
sc    <- readRDS(file.path(code_dir, "dr_scores.rds"))
rp    <- readRDS(file.path(code_dir, "res_policy.rds"))
truth <- readRDS(file.path(code_dir, "truth.rds"))

orange <- "#D55E00"; blue <- "#0072B2"; green <- "#009E73"; grey <- "grey30"

## ---- Figure 1: MC sampling distribution --------------------------
sate0 <- mean(mc$sate)
mlong <- rbind(
  data.table(est = mc$naive, method = "Naive ML plug-in"),
  data.table(est = mc$dml,   method = "DML (cross-fit)")
)
mlong[, method := factor(method, levels = c("Naive ML plug-in",
                                            "DML (cross-fit)"))]
p1 <- ggplot(mlong, aes(x = est, fill = method, colour = method)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = sate0, linewidth = 0.5, colour = grey,
             linetype = "dashed") +
  annotate("text", x = sate0, y = Inf, label = " true SATE",
           hjust = 0, vjust = 1.4, size = 3.2, colour = grey) +
  scale_fill_manual(values = c("Naive ML plug-in" = orange,
                               "DML (cross-fit)" = blue)) +
  scale_colour_manual(values = c("Naive ML plug-in" = orange,
                                 "DML (cross-fit)" = blue)) +
  labs(x = "ATE estimate across Monte Carlo draws", y = "Density",
       fill = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_16_mc.svg"), p1,
       width = 8, height = 4.6, device = svglite::svglite)

## ---- Figure 2: GATES calibration ---------------------------------
g <- as.data.table(cf$gates)
setnames(g, c("est", "se", "tau_true"))
g[, quint := factor(1:5, labels = c("Q1\n(lowest CATE)", "Q2", "Q3",
                                    "Q4", "Q5\n(highest CATE)"))]
p2 <- ggplot(g, aes(x = quint)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_errorbar(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                width = 0.16, linewidth = 0.4, colour = blue) +
  geom_point(aes(y = est, colour = "AIPW group ATE"), size = 2.4) +
  geom_point(aes(y = tau_true, colour = "True group mean"),
             shape = 17, size = 2.4) +
  scale_colour_manual(values = c("AIPW group ATE" = blue,
                                 "True group mean" = orange)) +
  labs(x = "Merchants sorted into quintiles by estimated CATE",
       y = "Treatment effect on log GMV growth", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_16_cate.svg"), p2,
       width = 8, height = 4.6, device = svglite::svglite)

## ---- Figure 3: value vs fraction treated (targeting curve) -------
cost <- rp$cost
ord  <- order(sc$tau_hat, decreasing = TRUE)
tt   <- sc$tau_true[ord]
n    <- length(tt)
frac <- (1:n) / n
# value gain over treat-none when treating the top-k by estimated CATE
cum_gain <- cumsum(tt - cost) / n
curve <- data.table(frac = frac, gain = cum_gain)
idx   <- unique(c(seq(1, n, length.out = 250), n))   # thin for SVG size
curve <- curve[round(idx)]
# key operating points
tree_share <- rp$share["tree"]
tree_gain  <- rp$true_gain["tree"]
all_gain   <- rp$true_gain["all"]
oracle_gain <- rp$true_gain["oracle"]
oracle_share <- rp$share["oracle"]
p3 <- ggplot(curve, aes(x = frac, y = gain)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.7, colour = blue) +
  geom_point(aes(x = 1, y = all_gain), colour = orange, size = 2.6) +
  annotate("text", x = 1, y = all_gain, label = "treat all ",
           hjust = 1.05, vjust = -0.6, size = 3.1, colour = orange) +
  geom_point(aes(x = tree_share, y = tree_gain), colour = green,
             size = 2.6, shape = 15) +
  annotate("text", x = tree_share, y = tree_gain, label = " policy tree",
           hjust = 0, vjust = 1.8, size = 3.1, colour = green) +
  labs(x = "Fraction of merchants treated (ranked by estimated CATE)",
       y = "Net value gain over treating no one") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_16_policy.svg"), p3,
       width = 8, height = 4.6, device = svglite::svglite)

cat("Saved: fig_16_mc.svg, fig_16_cate.svg, fig_16_policy.svg\n")
