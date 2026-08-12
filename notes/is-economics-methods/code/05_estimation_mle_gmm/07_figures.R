# ============================================================================
# 07_figures.R  --  two publication-quality SVGs for chapter 5
#   (a) fig_05_identification.svg : the criterion peak (identification) +
#       strong-vs-weak objective slice (identification strength)
#   (b) fig_05_sandwich.svg : the QMLE sampling distribution of b1_hat vs the
#       too-narrow Fisher-information normal and the correct sandwich normal
# ============================================================================

suppressMessages({library(ggplot2); library(patchwork); library(svglite)})

oracle <- readRDS("oracle.rds")
qm     <- readRDS("res_qmle.rds")
wk     <- readRDS("res_weak.rds")
dat    <- readRDS("helix_data.rds")

fig_dir <- "/Users/bobbyzhong/Library/Mobile Documents/com~apple~CloudDocs/Research/Method/Micro Theory_Microeconometrics/notes_v2/assets/fig"

navy <- "#1a365d"; blue <- "#2b6cb0"; red <- "#c53030"
olive<- "#6b7c3a"; grey <- "#718096"; gold <- "#b7791f"
bg   <- "#fafaf9"; ink  <- "#2d3748"

base_theme <- theme_minimal(base_size = 12) +
  theme(plot.background  = element_rect(fill = bg, colour = NA),
        panel.background = element_rect(fill = bg, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e3dd", linewidth = 0.3),
        plot.title    = element_text(colour = navy, face = "bold", size = 13),
        plot.subtitle = element_text(colour = grey, size = 10),
        axis.title    = element_text(colour = ink),
        axis.text     = element_text(colour = ink),
        legend.position = "top", legend.title = element_blank(),
        legend.text   = element_text(colour = ink, size = 9.5))

# ============================================================================
# FIGURE A : identification as a well-separated peak
# ============================================================================
## Panel A1: the logit average log-likelihood as a function of b1 (b0 at MLE),
## a smooth concave peak at the truth -- the "extremum" of an extremum estimator
x <- dat$x; D <- dat$D
b0hat <- coef(glm(D ~ x, family = binomial()))[1]
grid1 <- seq(0.2, 2.0, length = 121)
ll <- sapply(grid1, function(bb) mean(D*plogis(b0hat+bb*x, log.p=TRUE) +
                                      (1-D)*plogis(-(b0hat+bb*x), log.p=TRUE)))
dA1 <- data.frame(b1 = grid1, ll = ll)
b1max <- grid1[which.max(ll)]
pA1 <- ggplot(dA1, aes(b1, ll)) +
  geom_line(colour = navy, linewidth = 1) +
  geom_vline(xintercept = oracle$a1, colour = red, linetype = "dashed", linewidth=0.6) +
  annotate("point", x = b1max, y = max(ll), colour = blue, size = 2.6) +
  annotate("text", x = oracle$a1, y = min(ll)+0.02, label = "truth b1 = 1.10",
           colour = red, size = 3.3, hjust = -0.05) +
  annotate("text", x = b1max, y = max(ll), label = "  argmax = MLE",
           colour = blue, size = 3.3, hjust = 0, vjust = 1.6) +
  labs(title = "An extremum estimator maximizes a criterion",
       subtitle = "logit average log-likelihood peaks at the truth (identification)",
       x = expression(b[1]), y = "average log-likelihood") +
  base_theme

## Panel A2: strong vs weak GMM objective slice (normalized) -----------------
sl <- rbind(transform(wk$sl_strong, id = "strong (regressor varies)"),
            transform(wk$sl_weak,   id = "weak (regressor barely varies)"))
pA2 <- ggplot(sl, aes(b1, Q, colour = id)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = oracle$b1, colour = red, linetype = "dashed", linewidth=0.6) +
  scale_colour_manual(values = c("strong (regressor varies)" = navy,
                                 "weak (regressor barely varies)" = gold)) +
  annotate("text", x = oracle$b1, y = 0.02, label = "truth", colour = red,
           size = 3.2, hjust = -0.1) +
  labs(title = "Identification strength = curvature of the criterion",
       subtitle = "a flat objective cannot tell one theta from another",
       x = expression(b[1]), y = "GMM objective (normalized)") +
  base_theme
figA <- pA1 + pA2 + plot_layout(widths = c(1,1))
ggsave(file.path(fig_dir, "fig_05_identification.svg"), figA,
       width = 10.5, height = 4.3, device = svglite::svglite)

# ============================================================================
# FIGURE B : the QMLE sandwich -- which SE tells the truth
# ============================================================================
## Reconstruct the MC sampling distribution of b1_hat (Poisson QMLE, NB truth)
set.seed(313)
b0 <- oracle$b0; b1 <- oracle$b1; th <- oracle$theta_nb; n <- nrow(dat)
R <- 2000
b1_mc <- replicate(R, {
  xx <- rnorm(n); mu <- exp(b0 + b1*xx); yy <- rnbinom(n, size=th, mu=mu)
  coef(glm(yy ~ xx, family = poisson()))[2]
})
dens <- data.frame(b1 = b1_mc)
xs <- seq(b1 - 4*qm$sd_true, b1 + 4*qm$sd_true, length = 400)
curv <- rbind(
  data.frame(b1 = xs, d = dnorm(xs, b1, qm$mean_se_mod),  which = "naive Fisher-info SE (too narrow)"),
  data.frame(b1 = xs, d = dnorm(xs, b1, qm$mean_se_sand), which = "sandwich SE (correct)"))
pB <- ggplot() +
  geom_histogram(data = dens, aes(b1, after_stat(density)), bins = 45,
                 fill = "#d9d5cc", colour = NA) +
  geom_line(data = curv, aes(b1, d, colour = which), linewidth = 1) +
  geom_vline(xintercept = b1, colour = red, linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = c("naive Fisher-info SE (too narrow)" = gold,
                                 "sandwich SE (correct)" = navy)) +
  annotate("text", x = b1, y = 0, label = "truth b1 = 0.50", colour = red,
           size = 3.3, vjust = -0.5, hjust = -0.05) +
  labs(title = "Under overdispersion, the naive SE lies",
       subtitle = "Monte-Carlo sampling distribution of the Poisson-QMLE slope vs the two SE stories",
       x = expression(hat(b)[1]), y = "density") +
  base_theme
ggsave(file.path(fig_dir, "fig_05_sandwich.svg"), pB,
       width = 9.5, height = 5.0, device = svglite::svglite)

if (nzchar(Sys.getenv("PNG_INSPECT"))) {
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "A.png"), figA, width=10.5, height=4.3)
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "B.png"), pB, width=9.5, height=5.0)
}
cat("Wrote fig_05_identification.svg and fig_05_sandwich.svg\n")
