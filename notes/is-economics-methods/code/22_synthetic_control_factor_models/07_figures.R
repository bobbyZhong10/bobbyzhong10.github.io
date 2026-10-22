# =====================================================================
# 22_synthetic_control_factor_models / 07_figures.R
# The chapter 22 figures (SVG, ggplot2 + svglite), written to
# ../../assets/fig/:
#   fig_22_synth.svg      S1: treated vs synthetic path + placebo-gaps panel
#                         (the inherited figure; the plotting code is the
#                         legacy 05_figures.R code, reading res_sc_core.rds,
#                         so the file stays byte-identical)
#   fig_22_scpi.svg       S1: scpi prediction intervals for each post-month
#                         effect against the true effect path
#   fig_22_hull.svg       S2: Aster-Capital outside the donor hull; SC, ASC
#                         and SDID counterfactual paths, and their gaps
#                         against the true effect
#   fig_22_staggered.svg  S3: estimates and 95% intervals vs the true 0.06,
#                         and coverage across the 24 screened draws
#   fig_22_design.svg     S4: synthetic treated vs synthetic control in the
#                         history, and estimates over the 500 replicate pilots
# Run after 02 to 06 (and 05b), from this directory:  Rscript 07_figures.R
# =====================================================================

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(patchwork)
}))

in_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(in_dir) == 0 || in_dir == "") in_dir <- "."
fig_dir <- normalizePath(file.path(in_dir, "..", "..", "assets", "fig"),
                         mustWork = FALSE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 13),
                legend.position = "top"))

# =====================================================================
# FIGURE 2 -- Synthetic control   (legacy 05_figures.R, unchanged except
# that it reads res_sc_core.rds, written by 02_core.R)
# =====================================================================
res_sc <- readRDS(file.path(in_dir, "res_sc_core.rds"))
traj   <- res_sc$traj
T0     <- res_sc$T0
tall   <- res_sc$traj_all

# ---- panel A: treated vs synthetic path ----------------------------
pathA <- data.frame(
  month = rep(traj$time_unit, 2),
  y     = c(traj$real_y, traj$synth_y),
  series = rep(c("Aster-Metro (observed)", "Synthetic Aster-Metro"),
               each = nrow(traj)))
pA <- ggplot(pathA, aes(month, y, colour = series, linetype = series)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = T0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = T0 - 0.5, y = min(pathA$y), label = "fee reform",
           hjust = 1, size = 3.2, colour = "grey30") +
  scale_colour_manual(values = c("Aster-Metro (observed)" = "#d95f0e",
                                 "Synthetic Aster-Metro" = "#2c7fb8")) +
  scale_linetype_manual(values = c("Aster-Metro (observed)" = "solid",
                                   "Synthetic Aster-Metro" = "22")) +
  labs(title = "Treated metro vs its synthetic control",
       x = NULL, y = "log GMV", colour = NULL, linetype = NULL)

# ---- panel B: gap for treated + in-space placebos ------------------
# Standard Abadie et al. practice: drop placebo donors whose pre-period
# fit is far worse than the treated unit's (their large gaps are just
# bad fit, not evidence). Keep those with pre-RMSPE <= 5x treated's.
tall <- tall %>% mutate(gap = real_y - synth_y,
                        is_treated = .id == "Aster-Metro")
tr_pre_rmspe <- tall %>%
  filter(is_treated, time_unit < T0) %>%
  summarise(r = sqrt(mean(gap^2))) %>% pull(r)
pre_fit <- tall %>% filter(time_unit < T0) %>%
  group_by(.id) %>% summarise(r = sqrt(mean(gap^2)), .groups = "drop")
keep_ids <- pre_fit %>% filter(r <= 5 * tr_pre_rmspe) %>% pull(.id)
tall <- tall %>% filter(.id %in% keep_ids)

pB <- ggplot() +
  geom_line(data = filter(tall, !is_treated),
            aes(time_unit, gap, group = .id),
            colour = "grey75", linewidth = 0.4, alpha = 0.8) +
  geom_line(data = filter(tall, is_treated),
            aes(time_unit, gap), colour = "#d95f0e", linewidth = 1.2) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = T0, linetype = "dashed", colour = "grey40") +
  annotate("text", x = T0 + 0.5, y = max(tall$gap),
           label = "treated gap = effect", hjust = 0, size = 3.2,
           colour = "#d95f0e") +
  labs(title = "Gap: treated (orange) vs donor placebos (grey)",
       x = "Month", y = "log GMV gap")

p_sc <- pA / pB + plot_layout(heights = c(1, 1))
ggsave(file.path(fig_dir, "fig_22_synth.svg"), p_sc,
       width = 7.2, height = 6.4, device = svglite::svglite)
cat("wrote", file.path(fig_dir, "fig_22_synth.svg"), "\n")

# =====================================================================
# New figures (same theme, the same orange for the treated metro / the
# truth side and blue for synthetic controls)
# =====================================================================
theme_update(plot.title.position = "plot", plot.subtitle = element_text(size = 10, colour = "grey30"))
col_tr  <- "#d95f0e"   # treated metro, DiD family
col_sc  <- "#2c7fb8"   # synthetic control
col_asc <- "#31a354"   # augmented SC
col_sd  <- "#756bb1"   # synthetic DiD
tr22 <- readRDS(file.path(in_dir, "truth22.rds"))

# ---- fig_22_scpi.svg: S1 prediction intervals per post month --------
s1 <- readRDS(file.path(in_dir, "res_s1_ext.rds"))
dS <- data.frame(month = 25:36, est = s1$scpi$est, lo = s1$scpi$lo, hi = s1$scpi$hi,
                 truth = tr22$S1$tau_path)
# labels built from the stored interval level, so a change of u_alpha or e_alpha cannot leave
# the figure contradicting numbers.json (coverage at least 1 - u_alpha - e_alpha, CFT 2021)
lab_pi <- sprintf("scpi prediction interval (at least %g%%)", 100 * (1 - s1$scpi$u_alpha - s1$scpi$e_alpha))
sub_pi <- sprintf("SC estimate (dot) and scpi interval: simplex weights, Gaussian bounds, %g%% + %g%% error budget",
                  100 * s1$scpi$u_alpha, 100 * s1$scpi$e_alpha)
p_scpi <- ggplot(dS, aes(month)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi, colour = lab_pi), width = 0.35, linewidth = 0.8) +
  geom_point(aes(y = est, colour = lab_pi), size = 2) +
  geom_line(aes(y = truth, linetype = "true effect"), colour = col_tr, linewidth = 1.1) +
  scale_colour_manual(values = setNames(col_sc, lab_pi)) +
  scale_linetype_manual(values = c("true effect" = "solid")) +
  scale_x_continuous(breaks = 25:36) +
  guides(linetype = guide_legend(order = 1), colour = guide_legend(order = 2)) +   # fixed legend order
  labs(title = "Aster-Metro: effect of the fee reform, month by month",
       subtitle = sub_pi,
       x = "Month", y = "effect on log GMV", colour = NULL, linetype = NULL)
ggsave(file.path(fig_dir, "fig_22_scpi.svg"), p_scpi, width = 7.2, height = 4.2, device = svglite::svglite)
cat("wrote", file.path(fig_dir, "fig_22_scpi.svg"), "\n")

# ---- fig_22_hull.svg: S2 paths and gaps ------------------------------
s2 <- readRDS(file.path(in_dir, "res_hull.rds"))
S2 <- readRDS(file.path(in_dir, "aster_S2.rds"))
mo <- 1:36
lab_tr <- "Aster-Capital (observed)"
lab_sc <- "SC (convex weights)"; lab_asc <- "Augmented SC (ridge)"; lab_sd <- "Synthetic DiD"
dH <- data.frame(month = rep(mo, 4),
                 y = c(s2$y1, s2$sc$synth, s2$asc$synth, s2$sdid$synth),
                 series = factor(rep(c(lab_tr, lab_sc, lab_asc, lab_sd), each = 36),
                                 levels = c(lab_tr, lab_sc, lab_asc, lab_sd)))
donors_long <- data.frame(month = rep(mo, each = 20), y = as.vector(S2$Y[-1, ]),
                          id = rep(S2$unit_names[-1], times = 36))
cols_h <- setNames(c(col_tr, col_sc, col_asc, col_sd), c(lab_tr, lab_sc, lab_asc, lab_sd))
ltys_h <- setNames(c("solid", "22", "42", "12"), c(lab_tr, lab_sc, lab_asc, lab_sd))
pH1 <- ggplot() +
  geom_line(data = donors_long, aes(month, y, group = id), colour = "grey85", linewidth = 0.35) +
  geom_line(data = dH, aes(month, y, colour = series, linetype = series), linewidth = 1) +
  geom_vline(xintercept = 25, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 24.5, y = min(donors_long$y), label = "fee reform", hjust = 1, size = 3.2, colour = "grey30") +
  scale_colour_manual(values = cols_h) + scale_linetype_manual(values = ltys_h) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  labs(title = sprintf("Aster-Capital (orange) among the %d donors (grey)", nrow(S2$Y) - 1L),
       x = NULL, y = "log GMV", colour = NULL, linetype = NULL)
dG <- data.frame(month = rep(mo, 3),
                 gap = c(s2$sc$gap, s2$asc$gap, s2$y1 - s2$sdid$synth),
                 series = factor(rep(c(lab_sc, lab_asc, lab_sd), each = 36), levels = c(lab_sc, lab_asc, lab_sd)))
dT <- data.frame(month = mo, truth = c(rep(0, 24), tr22$S2$tau_path))
pH2 <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(data = dG, aes(month, gap, colour = series, linetype = series), linewidth = 0.9) +
  geom_line(data = dT, aes(month, truth), colour = col_tr, linewidth = 1.2) +
  geom_vline(xintercept = 25, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 36, y = max(dT$truth) + 0.02, label = "true effect", hjust = 1, size = 3.2, colour = col_tr) +
  scale_colour_manual(values = cols_h[-1]) + scale_linetype_manual(values = ltys_h[-1]) +
  guides(colour = "none", linetype = "none") +
  labs(title = "Gap (observed minus counterfactual) vs the true effect (orange)",
       x = "Month", y = "log GMV gap")
p_hull <- pH1 / pH2 + plot_layout(heights = c(1.1, 1))
ggsave(file.path(fig_dir, "fig_22_hull.svg"), p_hull, width = 7.2, height = 6.8, device = svglite::svglite)
cat("wrote", file.path(fig_dir, "fig_22_hull.svg"), "\n")

# ---- fig_22_staggered.svg: S3 estimates and the seed screen ----------
s3 <- readRDS(file.path(in_dir, "res_stag.rds"))
keys <- c("twfe", "cs", "cs_nyt", "ife", "mc_default", "mc", "multisynth", "sdid_cohort")
labs3 <- c(twfe = "TWFE", cs = "Callaway-Sant'Anna (never-treated)", cs_nyt = "Callaway-Sant'Anna (not-yet-treated)",
           ife = "Interactive FE (fect, CV r)", mc_default = "Matrix completion (fect default)",
           mc = "Matrix completion (refit, tight tol)", multisynth = "Staggered SC (multisynth)",
           sdid_cohort = "Synthetic DiD by cohort")
fam <- c(twfe = "DiD family", cs = "DiD family", cs_nyt = "DiD family", ife = "factor / SC family",
         mc_default = "factor / SC family", mc = "factor / SC family", multisynth = "factor / SC family",
         sdid_cohort = "factor / SC family")
d3 <- data.frame(key = keys, label = factor(labs3[keys], levels = rev(labs3[keys])), family = fam[keys],
                 est = vapply(keys, function(k) s3[[k]]$est, 0),
                 lo = vapply(keys, function(k) s3[[k]]$ci[1], 0), hi = vapply(keys, function(k) s3[[k]]$ci[2], 0))
cols_f <- c("DiD family" = col_tr, "factor / SC family" = col_sc)
# The level of every S3 interval, in one place for the labels. stag_fns.R builds them with
# qnorm(0.975) (CS, SDID by cohort) or level 0.95 (TWFE; fect and multisynth defaults); the
# normal-based ones are checked here, so a change of level there cannot leave a stale label.
level_s3 <- 0.95
for (k in c("cs", "sdid_cohort"))
  stopifnot(abs(diff(s3[[k]]$ci) - 2 * qnorm(1 - (1 - level_s3) / 2) * s3[[k]]$se) < 1e-9)
lvl_s3 <- sprintf("%g%%", 100 * level_s3)
# multisynth's default horizon (augsynth n_leads = NULL): the last cohort's post window, so its
# estimate averages each treated metro's first ms_leads post months, not every treated metro-month
ms_leads <- s3$design$n_months - max(s3$design$waves) + 1
p31 <- ggplot(d3, aes(y = label, colour = family)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(xintercept = tr22$S3$att, linetype = "dashed", colour = "grey20") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.3, linewidth = 0.8) +
  geom_point(aes(x = est), size = 2.4) +
  annotate("text", x = tr22$S3$att, y = length(keys) + 0.6, label = sprintf("true effect %s", format(tr22$S3$att)), hjust = -0.05,
           size = 3.2, colour = "grey20") +
  scale_colour_manual(values = cols_f) +
  coord_cartesian(clip = "off") +
  labs(title = sprintf("Staggered seller support program: estimate and %s interval", lvl_s3),
       x = sprintf("average effect on log GMV (treated metro-months)\nmultisynth: each treated metro's first %d post months", ms_leads),
       y = NULL, colour = NULL)
plots3 <- p31
if (!is.null(s3$screen)) {
  sk <- c("twfe", "cs", "ife", "mc_default", "mc", "multisynth")
  dS3 <- data.frame(label = factor(labs3[sk], levels = rev(labs3[keys])), family = fam[sk],
                    covered = vapply(sk, function(k) s3$screen$covered[[k]], 0))
  p32 <- ggplot(dS3, aes(y = label, x = covered, fill = family)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = paste0(covered, "/", s3$screen$n_seeds)), hjust = -0.15, size = 3.2, colour = "grey20") +
    scale_fill_manual(values = cols_f, guide = "none") +
    scale_x_continuous(limits = c(0, s3$screen$n_seeds + 3), breaks = c(0, 6, 12, 18, 24)) +
    labs(title = sprintf("Intervals covering the truth across %d simulated draws", s3$screen$n_seeds),
         x = sprintf("draws (%d resamples) whose %s interval covers %s", s3$screen$n_seeds,
                     lvl_s3, format(s3$screen$truth)), y = NULL)
  plots3 <- p31 / p32 + plot_layout(heights = c(1.5, 1))
}
ggsave(file.path(fig_dir, "fig_22_staggered.svg"), plots3, width = 7.6, height = 7.0, device = svglite::svglite)
cat("wrote", file.path(fig_dir, "fig_22_staggered.svg"), "\n")

# ---- fig_22_design.svg: S4 synthetic design vs random assignment -----
s4 <- readRDS(file.path(in_dir, "res_design.rds"))
S4 <- readRDS(file.path(in_dir, "aster_S4.rds"))
sy <- s4$synthetic
w_full <- setNames(numeric(s4$n_metros), S4$unit_names); w_full[names(sy$w)] <- sy$w
v_full <- setNames(numeric(s4$n_metros), S4$unit_names); v_full[names(sy$v)] <- sy$v
Yobs <- S4$Y0
Yobs[names(sy$w), 25:30] <- Yobs[names(sy$w), 25:30] + s4$tau            # replicate 1, pilot applied
lab_pop <- sprintf("All %d metros (average)", s4$n_metros)
lab_st <- sprintf("Synthetic treated (%d pilot metros)", sy$n_treated); lab_sco <- "Synthetic control"
dD <- data.frame(month = rep(1:30, 3),
                 y = c(colMeans(S4$Y0), drop(w_full %*% Yobs), drop(v_full %*% Yobs)),
                 series = factor(rep(c(lab_pop, lab_st, lab_sco), each = 30), levels = c(lab_pop, lab_st, lab_sco)))
pD1 <- ggplot(dD, aes(month, y, colour = series, linetype = series)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 24.5, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 24, y = min(dD$y), label = "pilot starts", hjust = 1, size = 3.2, colour = "grey30") +
  scale_colour_manual(values = setNames(c("grey45", col_tr, col_sc), c(lab_pop, lab_st, lab_sco))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "22"), c(lab_pop, lab_st, lab_sco))) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  labs(title = sprintf("Pilot metros chosen to match the %d-metro average", s4$n_metros),
       x = "Month", y = "log GMV", colour = NULL, linetype = NULL)
lab_d <- c("Synthetic design", sprintf("Random %d, diff. in means", s4$m_pilot), sprintf("Random %d, DiD", s4$m_pilot))
rmse_d <- c(s4$rmse_synthetic_design, s4$rmse_random_design, s4$rmse_random_did)
dE <- data.frame(design = factor(rep(lab_d, each = s4$reps), levels = rev(lab_d)),
                 est = s4$tau + c(s4$err_syn, s4$err_random, s4$err_random_did))
pD2 <- ggplot(dE, aes(x = est, y = design)) +
  geom_vline(xintercept = s4$tau, linetype = "dashed", colour = "grey20") +
  geom_boxplot(fill = "grey92", colour = "grey30", outlier.size = 0.6, outlier.colour = "grey55", width = 0.55) +
  annotate("text", x = s4$tau, y = 3.55, label = sprintf("true effect %s", format(s4$tau)), hjust = -0.05, size = 3.2, colour = "grey20") +
  annotate("text", x = max(dE$est), y = 3:1 + 0.38, label = sprintf("RMSE %.4f", rmse_d), hjust = 1, size = 3.2,
           colour = "grey20") +
  coord_cartesian(clip = "off") +
  labs(title = sprintf("Estimated pilot effect over %d replicate pilot periods", s4$reps),
       x = sprintf("estimated effect on log GMV (%d-month average)", s4$n_pilot_months), y = NULL)
p_design <- pD1 / pD2 + plot_layout(heights = c(1, 1))
ggsave(file.path(fig_dir, "fig_22_design.svg"), p_design, width = 7.2, height = 7.0, device = svglite::svglite)
cat("wrote", file.path(fig_dir, "fig_22_design.svg"), "\n")
