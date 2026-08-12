# ============================================================================
# 04_figures.R  --  two SVGs for chapter 7
#   (a) fig_07_estimators.svg : static (pool/RE/FE) and dynamic (pool/FE/AB/BB)
#       estimator comparisons -- biased naive vs consistent, truth bracketed
#   (b) fig_07_nickell.svg : Nickell bias shrinking as 1/T; and rho->1 weak IV
# ============================================================================

suppressMessages({library(ggplot2); library(patchwork); library(svglite)})
st <- readRDS("res_static.rds"); dy <- readRDS("res_dynamic.rds")

fig_dir <- "/Users/bobbyzhong/Library/Mobile Documents/com~apple~CloudDocs/Research/Method/Micro Theory_Microeconometrics/notes_v2/assets/fig"
navy <- "#1a365d"; blue <- "#2b6cb0"; red <- "#c53030"
olive<- "#6b7c3a"; grey <- "#718096"; gold <- "#b7791f"; bg <- "#fafaf9"; ink <- "#2d3748"
base_theme <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = bg, colour = NA),
        panel.background = element_rect(fill = bg, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e3dd", linewidth = 0.3),
        plot.title = element_text(colour = navy, face = "bold", size = 12.5),
        plot.subtitle = element_text(colour = grey, size = 9.5),
        axis.title = element_text(colour = ink), axis.text = element_text(colour = ink),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(colour = ink, size = 9))

# ---- FIGURE A : estimator comparisons --------------------------------------
dfA1 <- data.frame(
  est = factor(c("pooled OLS","random effects","fixed effects"),
               levels = c("pooled OLS","random effects","fixed effects")),
  val = c(unname(st$pool), unname(st$re), unname(st$fe)))
pA1 <- ggplot(dfA1, aes(est, val, fill = est)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0.50, colour = red, linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", val)), vjust = -0.4, size = 3.2, colour = ink) +
  annotate("text", x = 0.7, y = 0.53, label = "truth 0.50", colour = red, size = 3.1, hjust = 0) +
  scale_fill_manual(values = c("pooled OLS"=gold,"random effects"=olive,"fixed effects"=navy), guide="none") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0,0)) +
  labs(title = "Static panel: effect of promotions",
       subtitle = "OLS/RE inflated by merchant quality; only FE recovers the truth",
       x = NULL, y = expression(hat(beta))) + base_theme +
  theme(axis.text.x = element_text(size = 8.5))

dfA2 <- data.frame(
  est = factor(c("pooled OLS","fixed effects","Arellano-Bond","Blundell-Bond"),
               levels = c("pooled OLS","fixed effects","Arellano-Bond","Blundell-Bond")),
  val = c(unname(dy$pool), unname(dy$fe), unname(dy$ab), unname(dy$bb)))
pA2 <- ggplot(dfA2, aes(est, val, fill = est)) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = 0.70, colour = red, linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", val)), vjust = -0.4, size = 3.2, colour = ink) +
  annotate("text", x = 0.7, y = 0.73, label = "truth 0.70", colour = red, size = 3.1, hjust = 0) +
  scale_fill_manual(values = c("pooled OLS"=gold,"fixed effects"=blue,
                               "Arellano-Bond"=navy,"Blundell-Bond"=olive), guide="none") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0,0)) +
  labs(title = "Dynamic panel: persistence of GMV",
       subtitle = "OLS biased up, FE biased down (Nickell); GMM recovers rho",
       x = NULL, y = expression(hat(rho))) + base_theme +
  theme(axis.text.x = element_text(size = 8))
figA <- pA1 + pA2
ggsave(file.path(fig_dir, "fig_07_estimators.svg"), figA,
       width = 10.8, height = 4.3, device = svglite::svglite)

# ---- FIGURE B : Nickell vs T, and rho->1 weak IV ---------------------------
nk <- dy$nick; Ts <- as.numeric(colnames(nk))
dfB1 <- rbind(
  data.frame(T = Ts, rho = nk["fe",], series = "within / FE"),
  data.frame(T = Ts, rho = nk["ab",], series = "Arellano-Bond"),
  data.frame(T = Ts, rho = nk["nickell_approx",], series = "Nickell -(1+rho)/(T-1)"))
pB1 <- ggplot(dfB1, aes(T, rho, colour = series, linetype = series)) +
  geom_hline(yintercept = 0.70, colour = red, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  scale_x_continuous(breaks = Ts) +
  scale_colour_manual(values = c("within / FE"=blue,"Arellano-Bond"=navy,
                                 "Nickell -(1+rho)/(T-1)"=grey)) +
  scale_linetype_manual(values = c("within / FE"=1,"Arellano-Bond"=1,
                                   "Nickell -(1+rho)/(T-1)"=2)) +
  annotate("text", x = 26, y = 0.72, label = "truth 0.70", colour = red, size = 3) +
  labs(title = "Nickell bias vanishes like 1/T",
       subtitle = "FE attenuates rho at short T and creeps up; GMM is flat",
       x = "panel length T", y = expression(hat(rho))) + base_theme

dfB2 <- data.frame(
  method = factor(c("difference GMM","system GMM"), levels=c("difference GMM","system GMM")),
  rho = c(unname(dy$ab9), unname(dy$bb9)),
  se  = c(unname(dy$se_ab9), unname(dy$se_bb9)))
pB2 <- ggplot(dfB2, aes(method, rho, colour = method)) +
  geom_hline(yintercept = 0.90, colour = red, linetype = "dashed", linewidth = 0.6) +
  geom_errorbar(aes(ymin = rho-1.96*se, ymax = rho+1.96*se), width = 0.16, linewidth = 1) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.2f", rho)), hjust = -0.35, size = 3.2, colour = ink) +
  annotate("text", x = 0.6, y = 0.915, label = "truth 0.90", colour = red, size = 3, hjust = 0) +
  scale_colour_manual(values = c("difference GMM"=gold,"system GMM"=navy), guide="none") +
  labs(title = "As rho approaches 1, difference GMM goes weak",
       subtitle = "system GMM's level moments restore precision",
       x = NULL, y = expression(hat(rho)%+-%1.96~SE)) + base_theme +
  theme(axis.text.x = element_text(size = 8.5))
figB <- pB1 + pB2 + plot_layout(widths = c(1.25, 1))
ggsave(file.path(fig_dir, "fig_07_nickell.svg"), figB,
       width = 11.0, height = 4.4, device = svglite::svglite)

if (nzchar(Sys.getenv("PNG_INSPECT"))) {
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "A7.png"), figA, width=10.8, height=4.3)
  ggsave(file.path(Sys.getenv("PNG_INSPECT"), "B7.png"), figB, width=11.0, height=4.4)
}
cat("Wrote fig_07_estimators.svg and fig_07_nickell.svg\n")
