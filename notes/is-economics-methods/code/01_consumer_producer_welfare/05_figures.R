# ------------------------------------------------------------------
# 05_figures.R -- chapter 1 figures (SVG), same style as 05_m_estimation_nonlinear_regression/06_figures.R.
# Ruling P2: written to figures/ in this directory only (the text tasks copy them into assets/fig/).
#   fig_01_budget.svg          block W0: budget lines at p1 = 2 and 4 USD, indifference curves u0 and u1,
#                              the compensated budget at p1 = 4 (income e(p1 = 4, u0)) and the Hicks decomposition
#   fig_01_welfare_areas.svg   block W0: Marshallian demand and the two Hicksian demands; the areas to their
#                              left between 2 and 4 USD are dCS, CV and EV
#   fig_01_extrapolation.svg   block W1: the three fitted demand forms (averaged over subscribers' budgets)
#                              against the true mean demand, inside and far outside the 2 to 4 USD support
#   fig_01_ccei.svg            block W2: distribution of the Afriat efficiency index (CCEI) by type, with the
#                              Bronars random choosers
# Reads truth01.rds, res_welfare.rds, res_rp.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

tr <- readRDS("truth01.rds")
w1 <- readRDS("res_welfare.rds")
rp <- readRDS("res_rp.rds")
W0 <- tr$W0

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}
usd <- function(x) sprintf("%.2f", x)

## ---- Figure (a): block W0, budget lines, indifference curves, Hicks decomposition -----------
a <- W0$a
icurve <- function(u, x1) (u / x1^a)^(1 / (1 - a))                    # x2 on the indifference curve u
xs <- seq(4, 52, length.out = 400)
lab_a <- c(b0 = "Budget at p1 = 2 USD and indifference curve u0",
           b1 = "Budget at p1 = 4 USD and indifference curve u1",
           bc = sprintf("Compensated budget at p1 = 4 USD (income e(p1 = 4, u0) = %s USD)", usd(W0$e_p1_u0)))
lines_a <- rbind(
  data.frame(x1 = c(0, W0$w / W0$p1_0), x2 = c(W0$w, 0), series = lab_a[["b0"]], what = "budget"),
  data.frame(x1 = c(0, W0$w / W0$p1_1), x2 = c(W0$w, 0), series = lab_a[["b1"]], what = "budget"),
  data.frame(x1 = c(0, W0$e_p1_u0 / W0$p1_1), x2 = c(W0$e_p1_u0, 0), series = lab_a[["bc"]], what = "budget"))
curves_a <- rbind(data.frame(x1 = xs, x2 = icurve(W0$u0, xs), series = lab_a[["b0"]]),
                  data.frame(x1 = xs, x2 = icurve(W0$u1, xs), series = lab_a[["b1"]]))
curves_a <- curves_a[curves_a$x2 <= 140, ]
lines_a$series <- factor(lines_a$series, levels = unname(lab_a))
curves_a$series <- factor(curves_a$series, levels = unname(lab_a))
pts_a <- data.frame(x1 = c(W0$x1_0, W0$h1_p1_u0, W0$x1_1), x2 = c(W0$x2_0, W0$h2_p1_u0, W0$x2_1),
                    lab = c(sprintf("x(p1 = 2) = (%g, %g)", W0$x1_0, W0$x2_0),          # p1 is good 1's price throughout
                            sprintf("h(p1 = 4, u0) = (%s, %s)", usd(W0$h1_p1_u0), usd(W0$h2_p1_u0)),
                            sprintf("x(p1 = 4) = (%g, %g)", W0$x1_1, W0$x2_1)),
                    hj = c(-0.12, -0.1, 1.1), vj = c(-0.6, -0.5, 1.6))
arr <- data.frame(x = c(W0$x1_0, W0$h1_p1_u0), y = c(W0$x2_0, W0$h2_p1_u0),
                  xend = c(W0$h1_p1_u0, W0$x1_1), yend = c(W0$h2_p1_u0, W0$x2_1))
p1 <- ggplot() +
  geom_line(data = lines_a, aes(x = x1, y = x2, colour = series, linetype = series), linewidth = 0.7) +
  geom_line(data = curves_a, aes(x = x1, y = x2, colour = series), linewidth = 0.7, show.legend = FALSE) +
  geom_segment(data = arr, aes(x = x, y = y, xend = xend, yend = yend), colour = col_true, linewidth = 0.4,
               arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")) +
  geom_point(data = pts_a, aes(x = x1, y = x2), size = 2.2, colour = col_true) +
  geom_text(data = pts_a, aes(x = x1, y = x2, label = lab, hjust = hj, vjust = vj), size = 3.2, colour = col_true) +
  annotate("text", x = 17.2, y = 73.5, label = "substitution", size = 3.1, colour = col_true, hjust = 0) +
  annotate("text", x = 11.2, y = 72, label = "income", size = 3.1, colour = col_true, hjust = 1) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange, col_green), unname(lab_a))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed"), unname(lab_a))) +
  scale_x_continuous(limits = c(0, 52), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 140), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Premium content x1 (hours per month)", y = "Everything else x2 (USD per month)",
       colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  base_theme
save_fig(p1, "fig_01_budget.svg")

## ---- Figure (b): block W0, the three demand curves and the CV, dCS and EV areas -------------
ps <- seq(1.6, 5, length.out = 300)
h_u0 <- function(p) W0$h1_0 * (p / W0$p1_0)^(a - 1)                   # Hicksian demand at u0 (Shephard)
h_u1 <- function(p) W0$x1_1 * (p / W0$p1_1)^(a - 1)                   # Hicksian demand at u1
x_m <- function(p) a * W0$w / p                                         # Marshallian demand at w = 100
lab_b <- c(cv = sprintf("Hicksian demand at u0: area to its left from 2 to 4 USD is CV = %s USD", usd(W0$cv)),
           dcs = sprintf("Marshallian demand: area to its left from 2 to 4 USD is \u0394CS = %s USD", usd(W0$dcs)),
           ev = sprintf("Hicksian demand at u1: area to its left from 2 to 4 USD is EV = %s USD", usd(W0$ev)))
fns <- list(cv = h_u0, dcs = x_m, ev = h_u1)
curves_b <- do.call(rbind, lapply(names(fns), function(k) data.frame(p = ps, x = fns[[k]](ps), series = lab_b[[k]])))
pa <- seq(W0$p1_0, W0$p1_1, length.out = 200)
areas_b <- do.call(rbind, lapply(names(fns), function(k)
  data.frame(p = c(pa, rev(pa)), x = c(fns[[k]](pa), rep(0, length(pa))), series = lab_b[[k]])))
curves_b$series <- factor(curves_b$series, levels = unname(lab_b))
areas_b$series <- factor(areas_b$series, levels = unname(lab_b))
cols_b <- setNames(c(col_blue, col_true, col_orange), unname(lab_b))
p2 <- ggplot() +
  geom_polygon(data = areas_b, aes(x = x, y = p, fill = series), alpha = 0.16) +
  geom_hline(yintercept = c(W0$p1_0, W0$p1_1), linewidth = 0.3, linetype = "dotted", colour = col_true) +
  geom_path(data = curves_b, aes(x = x, y = p, colour = series), linewidth = 0.8) +
  annotate("text", x = 0.4, y = c(W0$p1_0, W0$p1_1), label = c("p1 = 2 USD", "p1 = 4 USD"), vjust = -0.5, hjust = 0,
           size = 3.2, colour = col_true) +
  scale_colour_manual(values = cols_b) +
  scale_fill_manual(values = cols_b) +
  scale_x_continuous(limits = c(0, 26), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Premium content x1 (hours per month)", y = "Premium price p1 (USD per hour)", colour = NULL, fill = NULL) +
  guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
  base_theme
save_fig(p2, "fig_01_welfare_areas.svg")

## ---- Figure (c): block W1, the fitted forms inside and far outside the support ----------------
fm <- w1$W1$forms
nm_c <- c(truth = "True mean demand", linear = "Linear form", loglinear = "Log-linear form",
          income = "Income-term form")
cv_lab <- function(cv_in, cv_out) sprintf("CV 2 to 4 USD = %s, CV 2 to 12 USD = %s", usd(cv_in), usd(cv_out))
lab_c <- c(truth = paste0(nm_c[["truth"]], ": ", cv_lab(tr$W1$cv_mean, tr$W1$cv12_mean)),
           linear = paste0(nm_c[["linear"]], ": ", cv_lab(fm$linear$cv, fm$linear$cv_out)),
           loglinear = paste0(nm_c[["loglinear"]], ": ", cv_lab(fm$loglinear$cv, fm$loglinear$cv_out)),
           income = paste0(nm_c[["income"]], ": ", cv_lab(fm$income$cv, fm$income$cv_out)))
dc <- tr$W1$demand_curve
curves_c <- rbind(data.frame(p = dc$p1, x = dc$x1_mean, series = lab_c[["truth"]]),
                  data.frame(p = w1$pgrid, x = fm$linear$demand_curve, series = lab_c[["linear"]]),
                  data.frame(p = w1$pgrid, x = fm$loglinear$demand_curve, series = lab_c[["loglinear"]]),
                  data.frame(p = w1$pgrid, x = fm$income$demand_curve, series = lab_c[["income"]]))
curves_c$series <- factor(curves_c$series, levels = unname(lab_c))
arms <- w1$W1$arms
p3 <- ggplot() +
  annotate("rect", xmin = 2, xmax = 4, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  annotate("text", x = 3, y = 0.5, label = "experiment's prices", size = 3.1, colour = col_true) +
  geom_line(data = curves_c, aes(x = p, y = x, colour = series, linetype = series), linewidth = 0.75) +
  geom_point(data = arms, aes(x = p1, y = x1_mean), size = 2, colour = col_true) +
  scale_colour_manual(values = setNames(c(col_true, col_orange, col_blue, col_green), unname(lab_c))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed", "dotdash"), unname(lab_c))) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
  labs(x = "Premium price p1 (USD per hour)", y = "Mean premium hours per subscriber", colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  base_theme
save_fig(p3, "fig_01_extrapolation.svg")

## ---- Figure (d): block W2, CCEI by type ----------------------------------------------------------
u <- rp$users; b <- rp$bronars
pct <- function(x) sprintf("%.1f%%", 100 * x)
lab_d <- c(rational = sprintf("CES optimisers (%d users): GARP pass rate %s", rp$n_rational, pct(rp$pass_rational)),
           noisy = sprintf("Noisy users (%d): pass rate %s", rp$n_noisy, pct(rp$pass_noisy)),
           bronars = sprintf("Uniform random choosers (%s, Bronars): pass rate %s",
                             format(rp$n_bronars, big.mark = ","), pct(rp$pass_bronars)))
cc <- rbind(data.frame(ccei = u$ccei[u$type == "rational"], series = lab_d[["rational"]]),
            data.frame(ccei = u$ccei[u$type == "noisy"], series = lab_d[["noisy"]]),
            data.frame(ccei = b$ccei, series = lab_d[["bronars"]]))
cc$series <- factor(cc$series, levels = unname(lab_d))
p4 <- ggplot(cc, aes(x = ccei, colour = series, linetype = series)) +
  stat_ecdf(linewidth = 0.75, pad = TRUE) +
  scale_colour_manual(values = setNames(c(col_green, col_blue, col_orange), unname(lab_d))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "solid"), unname(lab_d))) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(0.6, 1), ylim = c(0, 1)) +             # a few Bronars choosers fall below 0.6
  labs(x = "Afriat efficiency index (CCEI; 1 = passes GARP)", y = "Share with CCEI at or below x",
       colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  base_theme
save_fig(p4, "fig_01_ccei.svg")
