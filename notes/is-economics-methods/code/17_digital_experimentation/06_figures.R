# ------------------------------------------------------------------
# 06_figures.R -- Chapter figures (SVG) for Module 14.
#   fig_17_power.svg   : power vs sample size, plain vs CUPED
#   fig_17_peeking.svg : type-I error vs number of looks, naive vs mSPRT
#   fig_17_interf.svg  : design estimates vs true global ATE
# ------------------------------------------------------------------
suppressMessages({library(data.table); library(ggplot2)})

code_dir <- dirname(normalizePath(sub("--file=", "",
              grep("--file=", commandArgs(FALSE), value = TRUE))))
fig_dir  <- normalizePath(file.path(code_dir, "..", "..", "assets", "fig"))

p   <- readRDS(file.path(code_dir, "lumen_params.rds"))
pw  <- readRDS(file.path(code_dir, "res_power.rds"))
cu  <- readRDS(file.path(code_dir, "res_cuped.rds"))
pk  <- readRDS(file.path(code_dir, "res_peek.rds"))
it  <- readRDS(file.path(code_dir, "res_interf.rds"))

orange <- "#D55E00"; blue <- "#0072B2"; green <- "#009E73"; grey <- "grey30"

## ---- Figure 1: power vs n, plain vs CUPED ------------------------
za <- qnorm(0.975); zb <- qnorm(0.80)
sig <- p$sd_Y0; sig_cv <- sig * sqrt(cu$one_minus_rho2); del <- p$delta
ns  <- seq(3000, 120000, by = 1000)
powf <- function(n, s) { se <- s * sqrt(2 / n); nc <- del / se
                         pnorm(nc - za) + pnorm(-nc - za) }
pcv <- data.table(
  n = rep(ns, 2),
  power = c(sapply(ns, powf, s = sig), sapply(ns, powf, s = sig_cv)),
  design = rep(c("Plain difference in means", "CUPED-adjusted"), each = length(ns))
)
pcv[, design := factor(design, levels = c("Plain difference in means",
                                          "CUPED-adjusted"))]
p1 <- ggplot(pcv, aes(x = n / 1000, y = power, colour = design)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.3,
             colour = "grey55") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(x = 15, y = pw$power_run), colour = orange, size = 2.4) +
  annotate("text", x = 15, y = pw$power_run, label = " opening run",
           hjust = 0, vjust = -0.8, size = 3.1, colour = orange) +
  scale_colour_manual(values = c("Plain difference in means" = orange,
                                 "CUPED-adjusted" = blue)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Sample size per arm (thousands)",
       y = "Power to detect the true 2% lift", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_17_power.svg"), p1,
       width = 8, height = 4.6, device = svglite::svglite)

## ---- Figure 2: peeking type-I error ------------------------------
nc <- as.data.table(pk$naive_curve)
p2 <- ggplot(nc, aes(x = looks, y = fpr)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 0.3,
             colour = "grey55") +
  annotate("text", x = 1, y = 0.05, label = "nominal 0.05", hjust = 0,
           vjust = -0.6, size = 3, colour = "grey40") +
  geom_line(aes(colour = "Naive fixed-0.05 peeking"), linewidth = 0.8) +
  geom_point(aes(colour = "Naive fixed-0.05 peeking"), size = 2) +
  geom_hline(aes(yintercept = pk$msprt_typeI, colour = "mSPRT (always-valid)"),
             linewidth = 0.8) +
  scale_colour_manual(values = c("Naive fixed-0.05 peeking" = orange,
                                 "mSPRT (always-valid)" = blue)) +
  labs(x = "Number of looks at the data", y = "Type-I error rate (A/A test)",
       colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_17_peeking.svg"), p2,
       width = 8, height = 4.6, device = svglite::svglite)

## ---- Figure 3: interference / design comparison ------------------
di <- data.table(
  design = c("User-level A/B", "Switchback", "Cluster"),
  est = c(it$userAB["mean"], it$switchback["mean"], it$cluster["mean"]),
  sd  = c(it$userAB["sd"], it$switchback["sd"], it$cluster["sd"])
)
di[, design := factor(design, levels = c("User-level A/B", "Switchback",
                                         "Cluster"))]
p3 <- ggplot(di, aes(x = design, y = est)) +
  geom_hline(aes(yintercept = it$true_ate, linetype = "True global ATE"),
             colour = green, linewidth = 0.6) +
  geom_errorbar(aes(ymin = est - 1.96 * sd, ymax = est + 1.96 * sd),
                width = 0.14, linewidth = 0.4, colour = blue) +
  geom_point(colour = blue, size = 2.8) +
  scale_linetype_manual(values = c("True global ATE" = "dashed")) +
  labs(x = NULL, y = "Estimated treatment effect", linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_17_interf.svg"), p3,
       width = 8, height = 4.4, device = svglite::svglite)

cat("Saved: fig_17_power.svg, fig_17_peeking.svg, fig_17_interf.svg\n")
