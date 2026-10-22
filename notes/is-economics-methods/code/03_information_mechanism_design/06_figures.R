# ------------------------------------------------------------------
# 06_figures.R -- chapter 3 figures (SVG), same style as chapter 1's 05_figures.R.
# Ruling P2: written to figures/ in this directory only (the text tasks copy them into assets/fig/).
# Ruling P24: R versions replace the legacy hand-drawn fig_03_position.svg and fig_03_screening.svg.
#   fig_03_position.svg    block V0: values and per-click prices in Vireo's two-slot auction (VCG, GSP with
#                          truthful bids), and the revenue of each rule (1240, 760, 680)
#   fig_03_screening.svg   block V2, two types: participation lines t = theta q, the high type's indifference
#                          line through both bundles, the low type's quantity pushed down, the information rent
#   fig_03_menu.svg        block V2, continuous types: the optimal menu q(theta) with and without ironing
#                          against the efficient q = theta / cost2, for the regular and the bimodal density
#   fig_03_revenue.svg     block V1: revenue per auction by bidder count for the four formats against the
#                          truth, and the gain from the reserve r*, with Monte Carlo 2.5% to 97.5% ranges
#   fig_03_gpv.svg         block V1: first-price bid quantiles, the GPV-recovered values and the true value
#                          quantiles against tau, for n = 2, 4, 6, with the Monte Carlo band
# Reads truth03.rds, res_auction.rds, res_menu.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

tr <- readRDS("truth03.rds")
a1 <- readRDS("res_auction.rds")
s2 <- readRDS("res_menu.rds")
V0 <- tr$V0

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

## ---- Figure (a): block V0, the position auction ---------------------------------------------------
adv <- c(sprintf("Advertiser 1\ntop slot, %g clicks", V0$ctr[1]), sprintf("Advertiser 2\nsecond slot, %g clicks", V0$ctr[2]),
         "Advertiser 3\nno slot")
panel_a <- "Value and price per click (USD)"
panel_b <- "Revenue per month (USD)"
lab_val <- "Value per click"; lab_vcg <- "VCG price per click"; lab_gsp <- "GSP price per click, truthful bids"
what_lv <- c(lab_val, "Not an equilibrium", lab_gsp, lab_vcg)  # explicit level order, not locale-dependent sort()
what_f <- function(x) factor(x, levels = what_lv)
bars_a <- data.frame(panel = panel_a, x = factor(adv, levels = adv), y = V0$v, what = what_f(lab_val))
pts_a <- rbind(data.frame(panel = panel_a, x = factor(adv[1:2], levels = adv), y = c(V0$vcg_per_click_1, V0$vcg_per_click_2),
                          what = what_f(lab_vcg), lab = c(sprintf("%g, pays %g", V0$vcg_per_click_1, V0$vcg_pay1),
                                                  sprintf("%g, pays %g", V0$vcg_per_click_2, V0$vcg_pay2)), vj = 1.9),
               data.frame(panel = panel_a, x = factor(adv[1:2], levels = adv), y = c(V0$v[2], V0$v[3]),
                          what = what_f(lab_gsp), lab = c(sprintf("%g, pays %g", V0$v[2], V0$v[2] * V0$ctr[1]),
                                                  sprintf("%g, pays %g", V0$v[3], V0$v[3] * V0$ctr[2])), vj = -0.9))
rules <- c("Pay your bid,\ntruthful reports", "GSP,\ntruthful bids", "GSP lowest envy-free\n= VCG")
bars_b <- data.frame(panel = panel_b, x = factor(rules, levels = rules), y = c(V0$pm_naive, V0$gsp_truth, V0$vcg_rev),
                     what = what_f(c("Not an equilibrium", lab_gsp, lab_vcg)))
fill_map <- setNames(c("grey80", "grey60", col_orange, col_blue), c(lab_val, "Not an equilibrium", lab_gsp, lab_vcg))
bars_a$panel <- factor(bars_a$panel, levels = c(panel_a, panel_b)); bars_b$panel <- factor(bars_b$panel, levels = c(panel_a, panel_b))
pts_a$panel <- factor(pts_a$panel, levels = c(panel_a, panel_b))
p_a <- ggplot() +
  geom_col(data = bars_a, aes(x = x, y = y, fill = what), width = 0.55) +
  geom_col(data = bars_b, aes(x = x, y = y, fill = what), width = 0.55) +
  geom_text(data = bars_b, aes(x = x, y = y, label = sprintf("%g", y)), vjust = -0.4, size = 3.3, colour = col_true) +
  geom_point(data = pts_a, aes(x = x, y = y, colour = what), size = 3,
             position = position_dodge(width = 0.35)) +
  geom_text(data = pts_a, aes(x = x, y = y, label = lab, colour = what, vjust = vj), size = 3,
            position = position_dodge(width = 0.35), show.legend = FALSE) +
  facet_wrap(~panel, scales = "free") +
  scale_fill_manual(values = fill_map, breaks = c(lab_val, "Not an equilibrium")) +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), c(lab_vcg, lab_gsp))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = NULL, fill = NULL, colour = NULL) +
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  base_theme
save_fig(p_a, "fig_03_position.svg")

## ---- Figure (b): block V2, two-type screening ------------------------------------------------------
tt <- s2$two_type
qs <- seq(0, 9, length.out = 200)
lab_L <- sprintf("Low type: participation line t = %g q", tt$theta_L)
lab_H <- sprintf("High type: participation line t = %g q", tt$theta_H)
lab_I <- sprintf("High type: indifference line t = %g q - %g (through both bundles)", tt$theta_H, tt$rent_H)
lines_b <- rbind(data.frame(q = qs, t = tt$theta_L * qs, series = lab_L),
                 data.frame(q = qs, t = tt$theta_H * qs, series = lab_H),
                 data.frame(q = qs[qs >= 0.8], t = tt$theta_H * qs[qs >= 0.8] - tt$rent_H, series = lab_I))
lines_b$series <- factor(lines_b$series, levels = c(lab_L, lab_H, lab_I))
pts_b <- data.frame(q = c(tt$q_L, tt$q_H, tt$q_L_fb, tt$q_H_fb), t = c(tt$t_L, tt$t_H, tt$t_L_fb, tt$t_H_fb),
                    kind = factor(c("Optimal menu", "Optimal menu", "First best (types observed)", "First best (types observed)"),
                                  levels = c("First best (types observed)", "Optimal menu")),
                    lab = c(sprintf("B: low type (q = %g, t = %g)", tt$q_L, tt$t_L), sprintf("A: high type (q = %g, t = %g)", tt$q_H, tt$t_H),
                            sprintf("low type, first best (q = %g)", tt$q_L_fb), sprintf("high type, first best (t = %g)", tt$t_H_fb)),
                    hj = c(-0.08, -0.08, -0.08, 1.06), vj = c(1.6, 1.6, 1.6, -0.5))
p_b <- ggplot() +
  geom_line(data = lines_b, aes(x = q, y = t, colour = series, linetype = series), linewidth = 0.7) +
  geom_segment(aes(x = tt$q_H, xend = tt$q_H, y = tt$t_H, yend = tt$t_H_fb), colour = col_true, linewidth = 0.5,
               arrow = grid::arrow(length = grid::unit(0.15, "cm"), ends = "both", type = "closed")) +
  annotate("text", x = tt$q_H + 0.12, y = (tt$t_H + tt$t_H_fb) / 2, hjust = 0, size = 3.2, colour = col_true,
           label = sprintf("information rent %g", tt$rent_H)) +
  geom_segment(aes(x = tt$q_L_fb, xend = tt$q_L + 0.05, y = -3, yend = -3), colour = col_true,
               linewidth = 0.4, arrow = grid::arrow(length = grid::unit(0.15, "cm"), type = "closed")) +
  annotate("text", x = (tt$q_L + tt$q_L_fb) / 2, y = -5.6, size = 3.1, colour = col_true,
           label = sprintf("low type's quantity pushed down from %g to %g", tt$q_L_fb, tt$q_L)) +
  geom_point(data = pts_b, aes(x = q, y = t, shape = kind), size = 2.6, colour = col_true) +
  geom_text(data = pts_b, aes(x = q, y = t, label = lab, hjust = hj, vjust = vj), size = 3, colour = col_true) +
  scale_colour_manual(values = setNames(c(col_orange, col_blue, col_blue), c(lab_L, lab_H, lab_I))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed"), c(lab_L, lab_H, lab_I))) +
  scale_shape_manual(values = c("Optimal menu" = 16, "First best (types observed)" = 1)) +
  scale_x_continuous(limits = c(0, 10.2), breaks = 0:5 * 2, expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(-8, 64), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Usage tier q (thousand calls a month)", y = "Monthly fee t (USD)", colour = NULL, linetype = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 3, order = 1), linetype = guide_legend(nrow = 3, order = 1),
         shape = guide_legend(nrow = 2, order = 2)) +
  base_theme
save_fig(p_b, "fig_03_screening.svg")

## ---- Figure (c): block V2, continuous types, with and without ironing ----------------------------------
fac <- c(regular = "Regular type density", irregular = "Bimodal type density: ironing")
lab_eff <- "Efficient q = theta"; lab_opt <- "Optimal menu (ironed virtual value)"
lab_raw <- "Without ironing (not incentive compatible)"
mk <- function(nm) {
  m <- s2[[nm]]$menu; m0 <- s2[[nm]]$menu_naive
  rbind(data.frame(theta = m$theta, q = m$q_eff, series = lab_eff, panel = fac[[nm]]),
        data.frame(theta = m$theta, q = m$q, series = lab_opt, panel = fac[[nm]]),
        data.frame(theta = m0$theta, q = m0$q, series = lab_raw, panel = fac[[nm]]))
}
dc <- rbind(mk("regular"), mk("irregular"))
dc$series <- factor(dc$series, levels = c(lab_eff, lab_opt, lab_raw)); dc$panel <- factor(dc$panel, levels = unname(fac))
ir <- s2$irregular
rect_c <- data.frame(panel = factor(fac[["irregular"]], levels = unname(fac)), xmin = ir$bunch_lo, xmax = ir$bunch_hi)
txt_c <- data.frame(panel = factor(fac[["irregular"]], levels = unname(fac)), x = (ir$bunch_lo + ir$bunch_hi) / 2, y = 9,
                    lab = sprintf("bunch: types %s to %s\nall get q = %s", fmt(ir$bunch_lo), fmt(ir$bunch_hi), fmt(ir$bunch_q)))
p_c <- ggplot(dc, aes(x = theta, y = q, colour = series, linetype = series)) +
  geom_rect(data = rect_c, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
            fill = "grey90", colour = NA) +
  geom_line(linewidth = 0.75) +
  geom_text(data = txt_c, aes(x = x, y = y, label = lab), inherit.aes = FALSE, size = 3, colour = col_true) +
  facet_wrap(~panel) +
  scale_colour_manual(values = setNames(c("grey45", col_blue, col_orange), c(lab_eff, lab_opt, lab_raw))) +
  scale_linetype_manual(values = setNames(c("dashed", "solid", "dotted"), c(lab_eff, lab_opt, lab_raw))) +
  labs(x = "Type theta (USD per thousand calls)", y = "Usage tier q (thousand calls a month)", colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  base_theme
save_fig(p_c, "fig_03_menu.svg")

## ---- Figure (d): block V1, revenue by format and bidder count; the reserve gain ---------------------------
bn <- a1$by_n; mcb <- a1$mc$by_n; gbn <- a1$mc$gain_by_n
pan1 <- "Revenue per auction (USD per click)"; pan2 <- "Gain from the reserve r* (USD per click)"
lab_t0 <- "Expected revenue without a reserve (both formats)"
lab_t1 <- "Expected revenue with r* (both formats); right: expected gain"
f_lab <- c(rev_fp = "First price", rev_sp = "Second price", rev_fp_res = "First price, reserve r*",
           rev_sp_res = "Second price, reserve r*")
pts_d <- rbind(
  data.frame(panel = pan1, n = rep(bn$n, 4), format = rep(names(f_lab), each = nrow(bn)),
             y = c(bn$rev_fp, bn$rev_sp, bn$rev_fp_res, bn$rev_sp_res)),
  data.frame(panel = pan2, n = bn$n, format = "rev_sp_res", y = bn$gain_reserve))
rng_d <- rbind(
  data.frame(panel = pan1, n = mcb$n, format = mcb$format, lo = mcb$lo, hi = mcb$hi),
  data.frame(panel = pan2, n = gbn$n, format = "rev_sp_res", lo = gbn$gain_lo, hi = gbn$gain_hi))
truth_d <- rbind(
  data.frame(panel = pan1, n = tr$V1$rev_true$n, y = tr$V1$rev_true$none, series = lab_t0),
  data.frame(panel = pan1, n = tr$V1$rev_true$n, y = tr$V1$rev_true$reserve, series = lab_t1),
  data.frame(panel = pan2, n = tr$V1$rev_true$n, y = tr$V1$rev_true$reserve - tr$V1$rev_true$none, series = lab_t1))
for (d in c("pts_d", "rng_d")) { x <- get(d); x$format <- factor(f_lab[x$format], levels = unname(f_lab)); assign(d, x) }
pts_d$panel <- factor(pts_d$panel, levels = c(pan1, pan2)); rng_d$panel <- factor(rng_d$panel, levels = c(pan1, pan2))
truth_d$panel <- factor(truth_d$panel, levels = c(pan1, pan2))
dodge <- position_dodge(width = 0.5)
p_d <- ggplot() +
  geom_line(data = truth_d, aes(x = n, y = y, linetype = series), colour = "grey45", linewidth = 0.6) +
  geom_errorbar(data = rng_d, aes(x = n, ymin = lo, ymax = hi, colour = format), width = 0.25, position = dodge) +
  geom_point(data = pts_d, aes(x = n, y = y, colour = format, shape = format), size = 2, position = dodge) +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = setNames(c(col_blue, col_orange, col_blue, col_orange), unname(f_lab))) +
  scale_shape_manual(values = setNames(c(16, 16, 17, 17), unname(f_lab))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), c(lab_t0, lab_t1)), breaks = c(lab_t0, lab_t1)) +
  labs(x = "Bidders per auction n", y = NULL, colour = "Sample (2,000 auctions per n); bars: Monte Carlo 2.5% to 97.5%",
       shape = "Sample (2,000 auctions per n); bars: Monte Carlo 2.5% to 97.5%", linetype = NULL) +
  guides(colour = guide_legend(nrow = 2, title.position = "top", order = 1),
         shape = guide_legend(nrow = 2, title.position = "top", order = 1),
         linetype = guide_legend(nrow = 2, order = 2)) +
  base_theme + theme(legend.title = element_text(size = 9))
save_fig(p_d, "fig_03_revenue.svg")

## ---- Figure (e): block V1, GPV in quantile form ----------------------------------------------------------
gp <- a1$gpv; gb <- a1$mc$gpv_band
keep <- c(2L, 4L, 6L)
lab_v <- "True value quantile"; lab_bid <- "First-price bid quantile (observed)"
lab_g <- "GPV-recovered value (this sample)"; lab_band <- "GPV, Monte Carlo 2.5% to 97.5%"
de <- rbind(data.frame(n = gp$n, tau = gp$tau, y = gp$v_true, series = lab_v),
            data.frame(n = gp$n, tau = gp$tau, y = gp$b, series = lab_bid),
            data.frame(n = gp$n, tau = gp$tau, y = gp$v_hat, series = lab_g))
de <- de[de$n %in% keep, ]; de$series <- factor(de$series, levels = c(lab_v, lab_bid, lab_g))
de$panel <- factor(paste0("n = ", de$n), levels = paste0("n = ", keep))
be <- gb[gb$n %in% keep, ]; be$panel <- factor(paste0("n = ", be$n), levels = paste0("n = ", keep)); be$fill <- lab_band
p_e <- ggplot() +
  geom_ribbon(data = be, aes(x = tau, ymin = lo, ymax = hi, fill = fill), alpha = 0.25) +
  geom_line(data = de, aes(x = tau, y = y, colour = series, linetype = series), linewidth = 0.7) +
  facet_wrap(~panel) +
  scale_colour_manual(values = setNames(c(col_true, col_orange, col_blue), c(lab_v, lab_bid, lab_g))) +
  scale_linetype_manual(values = setNames(c("dashed", "solid", "solid"), c(lab_v, lab_bid, lab_g))) +
  scale_fill_manual(values = setNames(col_blue, lab_band)) +
  labs(x = "Quantile level tau", y = "USD per click", colour = NULL, linetype = NULL, fill = NULL) +
  guides(colour = guide_legend(nrow = 2, order = 1), linetype = guide_legend(nrow = 2, order = 1),
         fill = guide_legend(order = 2)) +
  base_theme
save_fig(p_e, "fig_03_gpv.svg")
