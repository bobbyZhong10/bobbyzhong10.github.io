# ============================================================================
# 05_figures.R  —  two publication-quality SVGs into assets/fig/
#   (a) fig_09_overlap.svg    : propensity-score distributions by group
#   (b) fig_09_estimators.svg : each estimator's point estimate + 95% CI
#                               against the true ATT and ATE
# ============================================================================
suppressMessages({ library(ggplot2); library(svglite) })

d   <- readRDS("northwind.rds")
obs <- d$obs
obs$D <- obs$D_cia
ps_fit <- glm(D ~ size + dig + eff + ind, data = obs, family = binomial)
obs$e  <- predict(ps_fit, type = "response")

figdir <- normalizePath(file.path("..", "..", "assets", "fig"), mustWork = TRUE)

# warm-navy palette consistent with the series
col_treat <- "#2b6cb0"   # link blue  -> adopters
col_ctrl  <- "#b7791f"   # warm ochre -> non-adopters
navy      <- "#1a365d"
slate     <- "#2d3748"
red       <- "#c53030"

theme_nw <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e3df"),
        plot.background  = element_rect(fill = "#fafaf9", colour = NA),
        panel.background = element_rect(fill = "#fafaf9", colour = NA),
        axis.title = element_text(colour = slate),
        axis.text  = element_text(colour = slate),
        plot.title = element_text(colour = navy, face = "bold"),
        legend.position = "top",
        legend.title = element_blank())

# ---- (a) overlap -----------------------------------------------------------
pd <- data.frame(e = obs$e,
                 grp = factor(ifelse(obs$D == 1, "Adopters (D=1)",
                                     "Non-adopters (D=0)"),
                              levels = c("Adopters (D=1)","Non-adopters (D=0)")))
p1 <- ggplot(pd, aes(e, fill = grp, colour = grp)) +
  geom_density(alpha = 0.35, adjust = 1.1) +
  scale_fill_manual(values = c(col_treat, col_ctrl)) +
  scale_colour_manual(values = c(col_treat, col_ctrl)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Propensity-score overlap by adoption status",
       x = "Estimated propensity score  e(X) = P(adopt | X)",
       y = "Density") +
  theme_nw
ggsave(file.path(figdir, "fig_09_overlap.svg"), p1,
       width = 7.2, height = 4.4, device = svglite)

# ---- (b) estimators vs truth ----------------------------------------------
re  <- readRDS("res_experiment.rds")
rn  <- readRDS("res_naive_reg.rds")
rm4 <- readRDS("res_matching_ipw_aipw.rds")
tr  <- d$truths
z <- qnorm(0.975)

est <- data.frame(
  label = c("RCT diff-in-means", "Naive diff-in-means", "OLS adjustment",
            "NN matching (AI SE)", "IPW (Hajek)", "AIPW"),
  est   = c(re$dim_est, rn$naive, rn$ols_att,
            rm4$match_ai_pt, rm4$ipw_hajek, rm4$aipw),
  se    = c(re$se, rn$se_naive, rn$se_ols_att,
            rm4$match_ai_se, rm4$se_hajek, rm4$se_aipw),
  study = c("Study A (RCT)", rep("Study B (observational)", 5)))
est$lo <- est$est - z*est$se
est$hi <- est$est + z*est$se
est$label <- factor(est$label, levels = rev(est$label))

lab_att <- data.frame(x = tr$true_ATT_cia,
                      label = factor("RCT diff-in-means", levels = levels(est$label)),
                      txt = sprintf("true ATT = %.3f", tr$true_ATT_cia))
lab_ate <- data.frame(x = tr$true_ATE_pop,
                      label = factor("Naive diff-in-means", levels = levels(est$label)),
                      txt = sprintf("true ATE = %.3f", tr$true_ATE_pop))

p2 <- ggplot(est, aes(est, label, colour = study)) +
  geom_vline(xintercept = tr$true_ATT_cia, linetype = "dashed", colour = navy) +
  geom_vline(xintercept = tr$true_ATE_pop, linetype = "dotted", colour = red) +
  geom_text(data = lab_att, aes(x = x, y = label, label = txt),
            colour = navy, hjust = -0.04, vjust = -0.9, size = 3.4,
            inherit.aes = FALSE) +
  geom_text(data = lab_ate, aes(x = x, y = label, label = txt),
            colour = red, hjust = 1.04, vjust = -0.9, size = 3.4,
            inherit.aes = FALSE) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.22, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("Study A (RCT)" = col_ctrl,
                                 "Study B (observational)" = col_treat)) +
  scale_x_continuous(expand = expansion(mult = c(0.13, 0.05))) +
  labs(title = "Estimates vs. the truth: naive bias and its repair",
       x = "Estimated effect on log tickets resolved per agent",
       y = NULL) +
  theme_nw
ggsave(file.path(figdir, "fig_09_estimators.svg"), p2,
       width = 7.6, height = 4.6, device = svglite)

# ---- (c) raw contrast (Section 1 hook): naive (Study B) vs RCT (Study A) ---
raw <- data.frame(
  label = c("Study B  naive difference", "Study A  randomized pilot"),
  est   = c(rn$naive, re$dim_est),
  grp   = c("Study B (observational)", "Study A (RCT)"))
raw$label <- factor(raw$label, levels = rev(raw$label))   # Study B on top
p3 <- ggplot(raw, aes(est, label, fill = grp)) +
  geom_col(width = 0.5) +
  geom_vline(xintercept = tr$true_ATT_cia, linetype = "dashed",
             colour = navy, linewidth = 0.7) +
  annotate("text", x = tr$true_ATT_cia, y = 2.62,
           label = sprintf("true ATT = %.3f", tr$true_ATT_cia),
           colour = navy, hjust = -0.05, vjust = 0.5, size = 3.6) +
  geom_text(aes(label = sprintf("%.3f", est)), hjust = -0.25,
            colour = slate, size = 4.2) +
  scale_fill_manual(values = c("Study A (RCT)" = col_ctrl,
                               "Study B (observational)" = col_treat),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 0.66), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
  labs(title = "One product, two answers: naive vs. randomized",
       x = "Estimated effect on log tickets resolved per agent",
       y = NULL) +
  theme_nw
ggsave(file.path(figdir, "fig_09_raw_contrast.svg"), p3,
       width = 7.2, height = 3.4, device = svglite)

cat("wrote:\n  ", file.path(figdir, "fig_09_overlap.svg"),
    "\n  ", file.path(figdir, "fig_09_estimators.svg"),
    "\n  ", file.path(figdir, "fig_09_raw_contrast.svg"), "\n")
