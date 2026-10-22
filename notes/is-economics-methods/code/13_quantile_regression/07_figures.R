# ------------------------------------------------------------------
# 07_figures.R -- chapter 13 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R
# and 20_did_beyond_binary/05_figures.R.
#   fig_13_qte.svg     block Q1 (and Q1r): estimated QTE with 95% intervals against the true QTE (the same
#                      in Q1 and Q1r), the OLS mean effect, and the true quantiles of the individual effects
#                      in Q1 (rank invariance) and in Q1r (re-paired): the QTE is not the distribution of
#                      individual effects
#   fig_13_uqpe.svg    block Q2: the RIF-regression coefficient on FBM (unconditional quantile partial effect)
#                      and the conditional QR coefficient, with 95% intervals and their true values
#   fig_13_cic.svg     block Q5: CIC and quantile DiD estimates of the QTT with 95% intervals against the
#                      true QTT
# Reads truth13.rds, res_qte.rds, res_extra.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- normalizePath(file.path("..", "..", "assets", "fig"))

tr <- readRDS("truth13.rds"); q <- readRDS("res_qte.rds"); x <- readRDS("res_extra.rds")
P <- tr$params

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"; col_purple <- "#CC79A7"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}
tau_breaks <- c(0.1, 0.25, 0.5, 0.75, 0.9)
dodge <- 0.012

## ---- Figure (a): block Q1, QTE against the distribution of individual effects ------------------
P1 <- P$Q1
tt <- seq(0.02, 0.98, by = 0.005)
q0 <- exp(P1$m + sqrt(P1$v) * qnorm(tt))
lab1 <- c(est = "QTE estimate, Q1 data (95% CI)", qte = "True QTE (Q1 and Q1r)",
          ind1 = "Quantiles of individual effects, Q1",
          ind1r = "Quantiles of individual effects, Q1r",
          ols = "OLS mean effect (95% CI)")
cv <- tr$Q1$curve                                                       # truth on the grid 0.02, 0.03, ..., 0.98
lines1 <- rbind(data.frame(tau = tt, value = exp(P1$logc) * q0^P1$k - q0, series = lab1[["qte"]]),
                data.frame(tau = cv$tau, value = cv$ind_Q1, series = lab1[["ind1"]]),
                data.frame(tau = cv$tau, value = cv$ind_Q1r, series = lab1[["ind1r"]]),
                data.frame(tau = c(0.02, 0.98), value = q$ols$est, series = lab1[["ols"]]))
lines1$series <- factor(lines1$series, levels = unname(lab1[c("qte", "ind1", "ind1r", "ols")]))
est1 <- q$qte_fine
est1$series <- factor(lab1[["est"]], levels = lab1[["est"]])
p1 <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = q$ols$ci[1], ymax = q$ols$ci[2], fill = col_green, alpha = 0.12) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(data = lines1, aes(x = tau, y = value, colour = series, linetype = series), linewidth = 0.75) +
  geom_errorbar(data = est1, aes(x = tau, ymin = lo, ymax = hi), colour = col_blue, width = 0.012, linewidth = 0.45) +
  geom_point(data = est1, aes(x = tau, y = est, shape = series), colour = col_blue, size = 2) +
  scale_colour_manual(values = setNames(c(col_true, col_orange, col_purple, col_green), levels(lines1$series))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "dotdash", "dotted"), levels(lines1$series))) +
  scale_shape_manual(values = setNames(16, lab1[["est"]])) +
  scale_x_continuous(breaks = tau_breaks, limits = c(0.02, 0.98)) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(x = expression("Quantile level " * tau * " (for individual effects: level u)"),
       y = "Effect of AutoPrice on monthly revenue\n(USD thousands)", colour = NULL, linetype = NULL, shape = NULL) +
  guides(colour = guide_legend(ncol = 2, order = 2), linetype = guide_legend(ncol = 2, order = 2),
         shape = guide_legend(order = 1)) +
  base_theme + theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
save_fig(p1, "fig_13_qte.svg")

## ---- Figure (b): block Q2, conditional against unconditional quantile effects -----------------
g2 <- x$Q2$grid
P2 <- P$Q2
lab2 <- c(rif = "RIF regression coefficient on FBM (95% CI)", cqr = "Conditional QR coefficient on FBM (95% CI)",
          uqpe = "True unconditional effect (UQPE)", cond = "True conditional effect (0.5 at every level)")
pts2 <- rbind(data.frame(tau = g2$tau - dodge / 2, est = g2$rif_est, lo = g2$rif_lo, hi = g2$rif_hi, series = lab2[["rif"]]),
              data.frame(tau = g2$tau + dodge / 2, est = g2$cqr_est, lo = g2$cqr_lo, hi = g2$cqr_hi, series = lab2[["cqr"]]))
pts2$series <- factor(pts2$series, levels = unname(lab2[c("rif", "cqr")]))
mixq <- function(t, p) uniroot(function(qq) (p * pnorm(qq, P2$mu + P2$g, sqrt(P2$bx^2 + P2$s1^2)) +
                                               (1 - p) * pnorm(qq, P2$mu, sqrt(P2$bx^2 + P2$s0^2))) - t, P2$mu + c(-10, 10), tol = 1e-12)$root
tt2 <- seq(0.03, 0.97, by = 0.005)
uq <- vapply(tt2, function(t) {
  qq <- mixq(t, P2$p)
  f <- P2$p * dnorm(qq, P2$mu + P2$g, sqrt(P2$bx^2 + P2$s1^2)) + (1 - P2$p) * dnorm(qq, P2$mu, sqrt(P2$bx^2 + P2$s0^2))
  (pnorm(qq, P2$mu, sqrt(P2$bx^2 + P2$s0^2)) - pnorm(qq, P2$mu + P2$g, sqrt(P2$bx^2 + P2$s1^2))) / f
}, numeric(1))
lines2 <- rbind(data.frame(tau = tt2, value = uq, series = lab2[["uqpe"]]),
                data.frame(tau = tt2, value = P2$g + (P2$s1 - P2$s0) * qnorm(tt2), series = lab2[["cond"]]))
lines2$series <- factor(lines2$series, levels = unname(lab2[c("uqpe", "cond")]))
p2 <- ggplot() +
  geom_line(data = lines2, aes(x = tau, y = value, linetype = series), colour = col_true, linewidth = 0.7) +
  geom_errorbar(data = pts2, aes(x = tau, ymin = lo, ymax = hi, colour = series), width = 0.01, linewidth = 0.45) +
  geom_point(data = pts2, aes(x = tau, y = est, colour = series, shape = series), size = 2) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), levels(pts2$series))) +
  scale_shape_manual(values = setNames(c(16, 17), levels(pts2$series))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), levels(lines2$series))) +
  scale_x_continuous(breaks = tau_breaks) +
  labs(x = expression("Quantile level " * tau), y = "Effect of FBM on log monthly revenue\n(log points)",
       colour = NULL, shape = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 2, order = 1), shape = guide_legend(nrow = 2, order = 1),
         linetype = guide_legend(nrow = 2, order = 2)) +
  base_theme
save_fig(p2, "fig_13_uqpe.svg")

## ---- Figure (c): block Q5, CIC against quantile DiD ------------------------------------------
g5 <- x$Q5$grid
P5 <- P$Q5
lab3 <- c(cic = "Changes-in-changes (95% CI)", qdid = "Quantile DiD (95% CI)", qtt = "True QTT")
pts3 <- rbind(data.frame(tau = g5$tau - dodge / 2, est = g5$cic_est, lo = g5$cic_lo, hi = g5$cic_hi, series = lab3[["cic"]]),
              data.frame(tau = g5$tau + dodge / 2, est = g5$qdid_est, lo = g5$qdid_lo, hi = g5$qdid_hi, series = lab3[["qdid"]]))
pts3$series <- factor(pts3$series, levels = unname(lab3[c("cic", "qdid")]))
tt3 <- seq(0.03, 0.97, by = 0.005)
line3 <- data.frame(tau = tt3, value = P5$delta * P5$A * exp(P5$g + P5$s * (1 + P5$kap) * (P5$mu1 + P5$sd1 * qnorm(tt3))),
                    series = factor(lab3[["qtt"]], levels = lab3[["qtt"]]))
p3 <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(data = line3, aes(x = tau, y = value, linetype = series), colour = col_true, linewidth = 0.75) +
  geom_errorbar(data = pts3, aes(x = tau, ymin = lo, ymax = hi, colour = series), width = 0.01, linewidth = 0.45) +
  geom_point(data = pts3, aes(x = tau, y = est, colour = series, shape = series), size = 2) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), levels(pts3$series))) +
  scale_shape_manual(values = setNames(c(16, 17), levels(pts3$series))) +
  scale_linetype_manual(values = setNames("solid", lab3[["qtt"]])) +
  scale_x_continuous(breaks = tau_breaks) +
  labs(x = expression("Quantile level " * tau), y = "Effect of One-Click Checkout on monthly revenue\n(USD thousands)",
       colour = NULL, shape = NULL, linetype = NULL) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  base_theme
save_fig(p3, "fig_13_cic.svg")
