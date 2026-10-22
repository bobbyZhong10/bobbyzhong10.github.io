# ------------------------------------------------------------------
# 06_figures.R -- chapter 2 figures (SVG), same style as 01_consumer_producer_welfare/05_figures.R.
# Ruling P2: written to figures/ in this directory only (the text tasks copy them into assets/fig/).
#   fig_02_qre_branches.svg    block G1: logit QRE of the compatibility game (x = 0) against the precision chi,
#                              Osprey's and Boreal's probability of S1; the Nash limits; the two precisions
#                              of block G1 and the point where the principal branch splits
#   fig_02_price_wars.svg      block G2: share of market-months at war against demand noise sigma (stationary
#                              share of the equilibrium and the simulated panels), and the trigger, collusive
#                              and Cournot prices
#   fig_02_regret_matching.svg block G3: CE violation of the empirical joint distribution against the number of
#                              periods, mean and largest over 100 regret-matching runs per game
#   fig_02_global_game.svg     block G4: share of merchants adopting against the fundamental theta for the
#                              three noise levels, simulated populations and the equilibrium curves
#   fig_02_gametree.svg        block G0 (Ruling P24, replaces the legacy hand-drawn figure): the sequential
#                              compatibility game, Osprey first, Boreal after observing; payoffs at the leaves;
#                              the backward-induction choices bold, the equilibrium path coloured
#   fig_02_folk.svg            block G0 (Ruling P24, replaces the legacy hand-drawn figure): the pricing game's
#                              feasible payoff set (convex hull of the four outcomes), both minmax lines, and the
#                              individually rational part shaded
# Reads truth02.rds, res_qre.rds, res_gp.rds, res_learn.rds, res_gg.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
suppressPackageStartupMessages(library(ggplot2))
source("games_fns.R")
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

tr <- readRDS("truth02.rds")
gp <- readRDS("res_gp.rds"); ln <- readRDS("res_learn.rds"); gg <- readRDS("res_gg.rds")
P1 <- tr$params$G1

col_true <- "grey25"; col_blue <- "#0072B2"; col_orange <- "#D55E00"; col_green <- "#009E73"
base_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name) {
  ggsave(file.path(fig_dir, name), p, width = 8, height = 4.8, device = svglite::svglite)
  cat(file.path(fig_dir, name), "\n")
}
chi <- "\u03c7"; sig <- "\u03c3"; th <- "\u03b8"

## ---- Figure (a): QRE branches against the precision ----------------------------------------
grid_chi <- seq(0.02, 5, by = 0.02)
bif <- tr$G1$bifurcation
lab_a <- c(principal = "Principal branch (from 1/2 at precision 0)", S1 = "Branch favouring S1",
           S2 = "Branch favouring S2")
rows <- lapply(grid_chi, function(l) {
  r <- qre_all(P1$A0, P1$B0, l)
  br <- if (nrow(r) == 1L) "principal" else c("S2", "principal", "S1")
  data.frame(chi = l, branch = br, p = r[, "p"], q = r[, "q"])
})
qd <- do.call(rbind, rows)
start <- data.frame(chi = c(bif[["lambda"]], bif[["lambda"]]), branch = c("S1", "S2"), p = bif[["p"]], q = bif[["q"]])
qd <- rbind(qd, start)
qd <- qd[order(qd$branch, qd$chi, method = "radix"), ]   # radix: C-locale string order, the same in every shell
long_a <- rbind(data.frame(chi = qd$chi, branch = qd$branch, who = "Osprey: Pr(S1) = p", prob = qd$p),
                data.frame(chi = qd$chi, branch = qd$branch, who = "Boreal: Pr(S1) = q", prob = qd$q))
long_a$who <- factor(long_a$who, levels = c("Osprey: Pr(S1) = p", "Boreal: Pr(S1) = q"))
long_a$branch <- factor(lab_a[long_a$branch], levels = unname(lab_a))
nash_a <- data.frame(who = factor(rep(levels(long_a$who), each = 3), levels = levels(long_a$who)),
                     y = c(1, 3 / 4, 0, 1, 1 / 4, 0))
vl <- data.frame(x = c(P1$lambda[["unique"]], P1$lambda[["high"]]),
                 lab = c(sprintf("%s = %g: unique-QRE block", chi, P1$lambda[["unique"]]),
                         sprintf("%s = %g: high-precision block", chi, P1$lambda[["high"]])))
p1 <- ggplot(long_a) +
  geom_hline(data = nash_a, aes(yintercept = y), linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_vline(data = vl, aes(xintercept = x), linetype = "dotted", linewidth = 0.4, colour = col_true) +
  geom_line(aes(x = chi, y = prob, colour = branch), linewidth = 0.8) +
  geom_text(data = vl, aes(x = x, y = 0.5, label = lab), angle = 90, vjust = -0.5, size = 2.9, colour = col_true) +
  facet_wrap(~who, nrow = 1) +
  scale_colour_manual(values = setNames(c(col_true, col_blue, col_orange), unname(lab_a))) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  labs(x = sprintf("Precision %s (per thousand USD)", chi), y = "Probability of choosing S1", colour = NULL,
       caption = sprintf("Dashed lines: Nash equilibria (the limits as %s grows). The principal branch splits at %s = %.3f.",
                         chi, chi, bif[["lambda"]])) +
  base_theme
save_fig(p1, "fig_02_qre_branches.svg")

## ---- Figure (b): price wars against demand noise ------------------------------------------------
lab_b <- c(theory = "Stationary share of the equilibrium",
           sim = sprintf("Simulated: %d markets x %d months (95%% interval)", gp$K, gp$T))
wd <- data.frame(sigma = gp$sigma, theory = gp$war_share_theory, sim = gp$war_freq, se = gp$war_freq_se)
lab_c <- c(trig = "Trigger price", coll = "Expected collusive price", cour = "Cournot price (price war)")
pd <- rbind(data.frame(sigma = gp$sigma, price = gp$trigger, s = lab_c[["trig"]]),
            data.frame(sigma = gp$sigma, price = gp$price_collusive, s = lab_c[["coll"]]),
            data.frame(sigma = gp$sigma, price = gp$cournot_price, s = lab_c[["cour"]]))
# two panels side by side (facets of one long table)
long_b <- rbind(data.frame(sigma = wd$sigma, y = wd$theory, lo = NA_real_, hi = NA_real_, series = lab_b[["theory"]],
                           panel = "Share of market-months at war"),
                data.frame(sigma = wd$sigma, y = wd$sim, lo = wd$sim - 1.96 * wd$se, hi = wd$sim + 1.96 * wd$se,
                           series = lab_b[["sim"]], panel = "Share of market-months at war"),
                data.frame(sigma = pd$sigma, y = pd$price, lo = NA_real_, hi = NA_real_, series = as.character(pd$s),
                           panel = "Prices (USD per unit)"))
long_b$panel <- factor(long_b$panel, levels = c("Share of market-months at war", "Prices (USD per unit)"))
ser_b <- c(unname(lab_b), unname(lab_c))
long_b$series <- factor(long_b$series, levels = ser_b)
cols_b <- setNames(c(col_true, col_orange, col_green, col_blue, "grey55"), ser_b)
ltys_b <- setNames(c("solid", "blank", "solid", "solid", "dashed"), ser_b)
p2 <- ggplot(long_b, aes(x = sigma, y = y, colour = series)) +
  geom_line(data = long_b[long_b$series != lab_b[["sim"]], ], aes(linetype = series), linewidth = 0.7) +
  geom_errorbar(data = long_b[!is.na(long_b$lo), ], aes(ymin = lo, ymax = hi), width = 0.01, show.legend = FALSE) +
  geom_point(size = 1.9) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = cols_b) +
  scale_linetype_manual(values = ltys_b, guide = "none") +
  scale_x_continuous(breaks = gp$sigma) +
  labs(x = sprintf("Demand noise %s (USD)", sig), y = NULL, colour = NULL) +
  guides(colour = guide_legend(ncol = 2)) +
  base_theme
save_fig(p2, "fig_02_price_wars.svg")

## ---- Figure (c): regret matching, CE violation against periods -------------------------------------
pa <- ln$paths
lab_g <- c(compat = "Compatibility game", chicken = "Chicken")
lab_s <- c(mean = sprintf("Mean over %d runs", ln$R), max = sprintf("Largest of %d runs", ln$R))
long_c <- rbind(data.frame(t = pa$t, v = pa$viol_reps_mean, game = lab_g[pa$game], stat = lab_s[["mean"]]),
                data.frame(t = pa$t, v = pa$viol_reps_max, game = lab_g[pa$game], stat = lab_s[["max"]]))
long_c$game <- factor(long_c$game, levels = unname(lab_g)); long_c$stat <- factor(long_c$stat, levels = unname(lab_s))
p3 <- ggplot(long_c, aes(x = t, y = v, colour = game, linetype = stat)) +
  geom_line(linewidth = 0.7) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 50000), labels = c("10", "100", "1,000", "10,000", "50,000")) +
  scale_y_log10() +
  scale_colour_manual(values = setNames(c(col_blue, col_orange), unname(lab_g))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), unname(lab_s))) +
  labs(x = "Periods of play t (log scale)", y = "Largest gain from disobeying (million USD, log scale)",
       colour = NULL, linetype = NULL) +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +   # fixed legend order (byte-stable SVG)
  base_theme
save_fig(p3, "fig_02_regret_matching.svg")

## ---- Figure (d): global game, adoption against the fundamental ----------------------------------------
lab_d <- setNames(sprintf("%s = %g: threshold %s* = %.4f", sig, gg$sigma, th, gg$theta_star), as.character(gg$sigma))
tg <- seq(0.25, 0.85, by = 0.001)
curves_d <- do.call(rbind, lapply(seq_along(gg$sigma), function(k)
  data.frame(theta = tg, share = stats::pnorm((tg - gg$theta_star[k]) / gg$sigma[k]), s = lab_d[[k]])))
pts_d <- data.frame(theta = gg$pop$theta, share = gg$pop$share, s = lab_d[as.character(gg$pop$sigma)])
curves_d$s <- factor(curves_d$s, levels = unname(lab_d)); pts_d$s <- factor(pts_d$s, levels = unname(lab_d))
p4 <- ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.4, colour = col_true) +
  annotate("text", x = 0.497, y = 0.92, label = sprintf("limit %s* = 1/2", th), hjust = 1, size = 3.1, colour = col_true) +
  geom_line(data = curves_d, aes(x = theta, y = share, colour = s), linewidth = 0.7) +
  geom_point(data = pts_d, aes(x = theta, y = share, colour = s), size = 1.1, alpha = 0.8) +
  scale_colour_manual(values = setNames(c(col_orange, col_green, col_blue), unname(lab_d))) +
  scale_x_continuous(breaks = seq(0.3, 0.8, by = 0.1)) +
  labs(x = sprintf("Fundamental %s", th), y = sprintf("Share of %s merchants adopting", format(gg$n_pop, big.mark = ",")),
       colour = NULL, caption = sprintf("Public signal %s ~ N(%g, %g^2). Lines: equilibrium share; points: simulated populations.",
                                        th, gg$prior_mean, gg$prior_sd)) +
  guides(colour = guide_legend(ncol = 1)) +
  base_theme
save_fig(p4, "fig_02_global_game.svg")

## ---- Figure (e): the sequential compatibility game, backward induction ----------------------------------
G0 <- tr$G0
Acp <- G0$compat_A; Bcp <- G0$compat_B; follow <- G0$seq_boreal_follow; first <- G0$seq_osprey_choice
std <- c("S1", "S2")
nodes <- data.frame(x = c(0, -2, 2), y = c(3, 1.8, 1.8), lab = c("Osprey", "Boreal", "Boreal"))
leaf_x <- c(-3, -1, 1, 3)
edges <- rbind(
  data.frame(x = 0, y = 3, xend = c(-2, 2), yend = 1.8, lab = std, who = "Osprey", i = 1:2, j = NA_integer_),
  data.frame(x = rep(c(-2, 2), each = 2), y = 1.8, xend = leaf_x, yend = 0.6, lab = rep(std, 2), who = "Boreal",
             i = rep(1:2, each = 2), j = rep(1:2, 2)))
edges$bi <- ifelse(edges$who == "Osprey", edges$i == first, edges$j == follow[edges$i])   # backward-induction choice
edges$path <- ifelse(edges$who == "Osprey", edges$i == first, edges$i == first & edges$j == follow[edges$i])
lab_e <- c(path = sprintf("Equilibrium path: Osprey %s, Boreal follows (payoffs %g, %g)", std[first],
                          G0$seq_payoff_osprey, G0$seq_payoff_boreal),
           bi = "Backward-induction choice off the path", other = "Choice not taken")
edges$kind <- factor(ifelse(edges$path, lab_e[["path"]], ifelse(edges$bi, lab_e[["bi"]], lab_e[["other"]])), levels = unname(lab_e))
edges$mx <- (edges$x + edges$xend) / 2 + ifelse(edges$xend < edges$x, -0.22, 0.22)
edges$my <- (edges$y + edges$yend) / 2 + 0.08
leaves <- data.frame(x = leaf_x, y = 0.6, i = rep(1:2, each = 2), j = rep(1:2, 2))
leaves$lab <- sprintf("(%g, %g)", Acp[cbind(leaves$i, leaves$j)], Bcp[cbind(leaves$i, leaves$j)])
p5 <- ggplot() +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend, colour = kind, linewidth = kind)) +
  geom_text(data = edges, aes(x = mx, y = my, label = lab), size = 3.4, colour = col_true) +
  geom_point(data = nodes, aes(x = x, y = y), size = 3, colour = col_true) +
  geom_label(data = nodes, aes(x = x, y = y + 0.3, label = lab), size = 3.4, colour = col_true, linewidth = 0) +
  geom_point(data = leaves, aes(x = x, y = y), size = 1.6, colour = col_true) +
  geom_text(data = leaves, aes(x = x, y = y - 0.28, label = lab), size = 3.6, colour = col_true) +
  scale_colour_manual(values = setNames(c(col_blue, col_true, "grey65"), unname(lab_e))) +
  scale_linewidth_manual(values = setNames(c(1.7, 1.7, 0.6), unname(lab_e))) +
  scale_x_continuous(limits = c(-3.6, 3.6)) + scale_y_continuous(limits = c(0.1, 3.5)) +
  labs(colour = NULL, linewidth = NULL,
       caption = "Payoffs (Osprey, Boreal), million USD a month. Boreal sees Osprey's standard before choosing.") +
  guides(colour = guide_legend(ncol = 1), linewidth = guide_legend(ncol = 1)) +
  theme_void(base_size = 11) + theme(legend.position = "bottom")
save_fig(p5, "fig_02_gametree.svg")

## ---- Figure (f): the folk theorem in the pricing game -------------------------------------------------------
hull <- as.data.frame(G0$folk_hull); ir <- as.data.frame(G0$folk_ir); mm <- G0$folk_minmax
names(hull) <- names(ir) <- c("u1", "u2")
lab_f <- c(hull = "Feasible payoffs: convex hull of the four outcomes",
           ir = sprintf("Feasible and above both minmax levels (%g, %g): sustainable for patient firms", mm[1], mm[2]))
outc <- data.frame(u1 = as.vector(G0$pd_A), u2 = as.vector(G0$pd_B),
                   lab = c("(High, High)", "(Low, High)", "(High, Low)", "(Low, Low)"))
outc$lab <- sprintf("%s: (%g, %g)", outc$lab, outc$u1, outc$u2)
outc$hj <- c(-0.08, 1.08, 0.4, 1.08); outc$vj <- c(-0.6, 1.6, -0.9, 1.4)   # (High, Low) above its point, clear of the minmax line
hull$s <- factor(lab_f[["hull"]], levels = unname(lab_f)); ir$s <- factor(lab_f[["ir"]], levels = unname(lab_f))
p6 <- ggplot() +
  geom_polygon(data = hull, aes(x = u1, y = u2, fill = s), alpha = 0.35, colour = "grey55", linewidth = 0.4) +
  geom_polygon(data = ir, aes(x = u1, y = u2, fill = s), alpha = 0.45) +
  geom_vline(xintercept = mm[1], linetype = "dashed", colour = col_orange, linewidth = 0.5) +
  geom_hline(yintercept = mm[2], linetype = "dashed", colour = col_orange, linewidth = 0.5) +
  annotate("text", x = mm[1] + 0.15, y = 15.3, label = sprintf("Osprey's minmax %g", mm[1]), hjust = 0, size = 3.1, colour = col_orange) +
  annotate("text", x = 0.15, y = mm[2] + 0.35, label = sprintf("Boreal's minmax %g", mm[2]), hjust = 0, size = 3.1, colour = col_orange) +
  geom_point(data = outc, aes(x = u1, y = u2), size = 2.2, colour = col_true) +
  geom_text(data = outc, aes(x = u1, y = u2, label = lab, hjust = hj, vjust = vj), size = 3.1, colour = col_true) +
  scale_fill_manual(values = setNames(c("grey80", col_blue), unname(lab_f))) +
  scale_x_continuous(limits = c(0, 16), breaks = seq(0, 16, 2), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 2), expand = expansion(mult = 0)) +
  coord_equal() +
  labs(x = "Osprey's payoff (million USD a month)", y = "Boreal's payoff (million USD a month)", fill = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  base_theme
save_fig(p6, "fig_02_folk.svg")
