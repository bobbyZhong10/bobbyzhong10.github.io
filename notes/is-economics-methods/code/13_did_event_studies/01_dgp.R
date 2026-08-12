# ------------------------------------------------------------------
# 01_dgp.R -- Kestrel Marketplace simulation DGP
#
# Staggered commission cut (12% -> 8%) across cities; outcome is
# city-month mean log GMV of sellers. Three treatment cohorts
# (months 13, 19, 25) plus a never-treated group. Treatment effects
# are heterogeneous and GROW in event time, which is exactly the
# setting where static TWFE fails.
#
# Design note: the never-treated group is deliberately small (8 of
# 60 cities). With a large never-treated group (say 20 of 60) almost
# all TWFE weight falls on clean treated-vs-untreated comparisons and
# the static TWFE estimate cannot fall below roughly half of the true
# ATT no matter how steep the dynamics are. Shrinking the clean
# control group shifts weight onto forbidden timing comparisons and
# produces the stark attenuation the chapter opens with.
# ------------------------------------------------------------------
library(data.table)

set.seed(11)

## ---- Panel dimensions -------------------------------------------
T_months     <- 36
cohort_g     <- c(13, 19, 25)   # first treated month per cohort
cohort_size  <- c(17, 17, 18)   # cities per cohort
n_never      <- 8               # never-treated cities
n_cities     <- sum(cohort_size) + n_never   # 60

## ---- Treatment-effect profile -----------------------------------
## Piecewise-linear ramp in event time e = month - g:
##   e in [0,5]  : slow start, tau0 -> tau5
##   e in (5,12] : rapid growth, tau5 -> plateau
##   e > 12      : plateau + slope_late * (e - 12), mild further drift
## plus a cohort boost (later cohorts slightly larger) and a
## city-specific multiplier (effect heterogeneity across cities).
tau0         <- 0.015
tau5         <- 0.020
plateau      <- 0.140
slope_late   <- 0.004
cohort_boost <- c(0.000, 0.004, 0.008)
mult_sd      <- 0.10

tau_profile <- function(e) {
  ifelse(e <= 5, tau0 + (tau5 - tau0) * e / 5,
  ifelse(e <= 12, tau5 + (plateau - tau5) * (e - 5) / 7,
         plateau + slope_late * (e - 12)))
}

## ---- Baseline structure -----------------------------------------
## Levels selection: earlier-treated cities are bigger markets
## (higher city FE). Harmless under parallel trends.
level_shift  <- c(0.50, 0.25, 0.10)
sd_cityfe    <- 0.35
trend_common <- 0.010            # common monthly growth in month FE
season_amp   <- 0.05
ar_rho       <- 0.35
ar_sd        <- 0.030

## ---- Build city table -------------------------------------------
cities <- data.table(
  city_id = 1:n_cities,
  g = c(rep(cohort_g, times = cohort_size), rep(0L, n_never))  # 0 = never treated
)
cities[, cohort := ifelse(g == 0, "never", paste0("g", g))]
cities[, city_fe := rnorm(n_cities, 0, sd_cityfe) +
         ifelse(g == 0, 0, level_shift[match(g, cohort_g)])]
cities[, tau_mult := pmax(rnorm(n_cities, 1, mult_sd), 0.5)]

## ---- Month fixed effects ----------------------------------------
months <- data.table(month = 1:T_months)
months[, month_fe := trend_common * month + season_amp * sin(2 * pi * month / 12)]

## ---- Assemble panel ---------------------------------------------
dt <- CJ(city_id = 1:n_cities, month = 1:T_months)
dt <- merge(dt, cities, by = "city_id")
dt <- merge(dt, months, by = "month")
setkey(dt, city_id, month)

## AR(1) noise within city
dt[, eps := {
  e <- numeric(.N)
  e[1] <- rnorm(1, 0, ar_sd / sqrt(1 - ar_rho^2))
  for (s in 2:.N) e[s] <- ar_rho * e[s - 1] + rnorm(1, 0, ar_sd)
  e
}, by = city_id]

## Event time and realized treatment effect
dt[, event_time := ifelse(g > 0, month - g, NA_integer_)]
dt[, D := as.integer(g > 0 & month >= g)]
dt[, tau_true := 0]
dt[D == 1, tau_true := (tau_profile(event_time) +
                          cohort_boost[match(g, cohort_g)]) * tau_mult]

## Outcome: log GMV
dt[, lgmv := 8 + city_fe + month_fe + tau_true + eps]

## ---- True estimands (from realized effects) ---------------------
true_att_simple <- dt[D == 1, mean(tau_true)]
true_att_by_g   <- dt[D == 1, .(att = mean(tau_true), n_cells = .N), keyby = g]
true_theta_e    <- dt[D == 1, .(theta = mean(tau_true), n_cells = .N), keyby = event_time]

cat("== Kestrel Marketplace DGP ==\n")
cat(sprintf("Cities: %d (cohorts g=13/19/25: %d/%d/%d, never-treated: %d); months: %d\n",
            n_cities, cohort_size[1], cohort_size[2], cohort_size[3], n_never, T_months))
cat(sprintf("True simple ATT (mean realized effect over treated post cells): %.4f\n",
            true_att_simple))
cat("True ATT by cohort:\n"); print(true_att_by_g[, .(g, att = round(att, 4), n_cells)])
cat("True dynamic effect at selected horizons:\n")
print(true_theta_e[event_time %in% c(0, 3, 6, 9, 12, 18, 23),
                   .(event_time, theta = round(theta, 4))])

## ---- Save --------------------------------------------------------
out_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
saveRDS(dt, file.path(out_dir, "kestrel_panel.rds"))
fwrite(dt, file.path(out_dir, "kestrel_panel.csv"))
saveRDS(list(true_att_simple = true_att_simple,
             true_att_by_g = true_att_by_g,
             true_theta_e = true_theta_e),
        file.path(out_dir, "truth.rds"))
cat("Saved kestrel_panel.rds / kestrel_panel.csv / truth.rds\n")
