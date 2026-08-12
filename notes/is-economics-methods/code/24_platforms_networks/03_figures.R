# ------------------------------------------------------------------
# 03_figures.R -- Chapter 20 figures (SVG).
#   fig_24_tipping.svg   : fulfilled-expectations demand, multiple
#                          equilibria and critical mass (network eff.)
#   fig_24_inertia.svg   : persistence decomposition (counterfactual
#                          stay rates: neither / hetero / state-dep / both)
#   fig_24_estimates.svg : delta estimates (naive / RE / CRE) vs truth
# ------------------------------------------------------------------
library(data.table)
library(ggplot2)

code_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
fig_dir  <- normalizePath(file.path(code_dir, "..", "..", "assets", "fig"))
truth <- readRDS(file.path(code_dir, "truth24.rds"))
res   <- readRDS(file.path(code_dir, "res24_estimate.rds"))

okabe <- c(orange = "#D55E00", blue = "#0072B2", green = "#009E73", grey = "grey40")

## ---- Figure 1: tipping / critical mass ---------------------------
## Fulfilled-expectations inverse demand p(n) = (1 - n) * (v0 + v1 * n):
##  marginal consumer's WTP = (untaken fraction proxy) times network
##  benefit v0 + v1 n. Hump-shaped in n -> a horizontal price can
##  cross it twice (low unstable, high stable equilibrium).
v0 <- 0.15; v1 <- 1.35
pfun <- function(n) (1 - n) * (v0 + v1 * n)
ngrid <- seq(0, 1, length.out = 400)
demand <- data.table(n = ngrid, p = pfun(ngrid))
p_star <- 0.22                                   # an illustrative price
## equilibria solve (1-n)(v0+v1 n) = p_star
roots <- polyroot(c(v0 - p_star, v1 - v0, -v1))  # -v1 n^2 + (v1-v0) n + (v0 - p*)
roots <- sort(Re(roots[abs(Im(roots)) < 1e-8]))
eq <- data.table(n = roots, p = p_star,
                 kind = c("unstable (critical mass)", "stable"))

p1 <- ggplot(demand, aes(n, p)) +
  geom_line(colour = okabe["blue"], linewidth = 0.8) +
  geom_hline(yintercept = p_star, colour = okabe["orange"],
             linewidth = 0.5, linetype = "dashed") +
  geom_point(data = eq, aes(n, p, shape = kind), size = 2.6,
             colour = okabe["orange"]) +
  annotate("text", x = 0.02, y = p_star + 0.015, hjust = 0,
           label = "price p", colour = okabe["orange"], size = 3.4) +
  scale_shape_manual(values = c("stable" = 16, "unstable (critical mass)" = 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x = "Installed base / expected network size  n",
       y = "Fulfilled-expectations willingness to pay  p(n)",
       shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_24_tipping.svg"), p1,
       width = 7.2, height = 4.6, device = svglite::svglite)

## ---- Figure 2: persistence decomposition -------------------------
dec <- data.table(
  world = factor(c("Fees + noise\nonly", "Heterogeneity\nonly (delta=0)",
                   "State dependence\nonly (sigma=0)", "Both channels"),
                 levels = c("Fees + noise\nonly", "Heterogeneity\nonly (delta=0)",
                            "State dependence\nonly (sigma=0)", "Both channels")),
  stay  = c(truth$stay_neither, truth$stay_no_sd,
            truth$stay_no_hetero, truth$stay_full))

p2 <- ggplot(dec, aes(world, stay, fill = world)) +
  geom_col(width = 0.66, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", stay)), vjust = -0.5, size = 3.5) +
  geom_hline(yintercept = truth$stay_neither, colour = okabe["grey"],
             linetype = "dotted", linewidth = 0.4) +
  scale_fill_manual(values = unname(c(okabe["grey"], okabe["blue"],
                               okabe["green"], okabe["orange"]))) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.2)) +
  labs(x = NULL, y = "Share repeating last month's platform (stay rate)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(fig_dir, "fig_24_inertia.svg"), p2,
       width = 7.6, height = 4.6, device = svglite::svglite)

## ---- Figure 3: delta estimates vs truth --------------------------
tab <- res$tab[model != "True (DGP)"]
tab[, model := factor(model, levels = rev(c("Naive pooled logit",
                                            "RE logit (no IC fix)",
                                            "Wooldridge CRE logit")))]
tab[, `:=`(lo = delta - 1.96 * delta_se, hi = delta + 1.96 * delta_se)]

p3 <- ggplot(tab, aes(delta, model)) +
  geom_vline(xintercept = truth$delta, colour = okabe["grey"],
             linetype = "dashed", linewidth = 0.5) +
  geom_text(data = data.table(delta = truth$delta, model = factor("Wooldridge CRE logit",
              levels = levels(tab$model))),
            label = "true delta = 0.90", hjust = -0.05, vjust = 2.0,
            colour = okabe["grey"], size = 3.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.16,
                colour = okabe["blue"], linewidth = 0.5) +
  geom_point(size = 2.6, colour = okabe["blue"]) +
  geom_text(aes(label = sprintf("%.3f", delta)), vjust = -1.0, size = 3.4) +
  scale_x_continuous(limits = c(0.8, 1.45)) +
  labs(x = "Estimated state dependence  delta  (switching cost)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_24_estimates.svg"), p3,
       width = 7.4, height = 3.4, device = svglite::svglite)

cat("Saved figures:\n",
    file.path(fig_dir, "fig_24_tipping.svg"), "\n",
    file.path(fig_dir, "fig_24_inertia.svg"), "\n",
    file.path(fig_dir, "fig_24_estimates.svg"), "\n")
cat(sprintf("Tipping equilibria at p=%.2f: n = %s\n", p_star,
            paste(round(roots, 3), collapse = ", ")))
