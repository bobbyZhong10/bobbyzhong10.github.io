# ------------------------------------------------------------------
# 07_mc_hausman.R -- block P7: size of the classical and the robust (Mundlak) Hausman tests when RE.1
# holds but RE.3 fails (heteroskedastic alpha_i, heteroskedastic AR(1) eps_it; Wooldridge 2010 §10.7.3,
# p. 333). One sample (nimbus_P7.rds, also checked against plm::phtest) and a Monte Carlo of R_mc = 2000
# samples (seed P7_mc). Classical: hausman_classic() on the FE and Swamy-Arora RE fits with their
# non-robust variances (= plm::phtest). Robust: hausman_robust(mundlak_fit()), cluster-robust Wald of xi = 0.
# Reads nimbus_P7.rds, truth10.rds; writes res_p7.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("panel_fns.R")
suppressPackageStartupMessages(library(plm))
tr <- readRDS("truth10.rds"); P7 <- tr$params$P7

both_tests <- function(d) {
  fe <- fe_fit(d$y, d$x, d$id); re <- re_fit(d$y, d$x, d$id)
  hc <- hausman_classic(fe$b, fe$V_classic, re$b["x"], re$V_classic["x", "x"])
  hr <- hausman_robust(mundlak_fit(d$y, d$x, d$id))
  c(fe = unname(fe$b), re = unname(re$b["x"]), H_classic = hc$stat, p_classic = hc$p, H_robust = hr$stat,
    p_robust = hr$p, vd_classic = unname(fe$V_classic[1, 1] - re$V_classic["x", "x"]), varphi = re$varphi)
}

## ---- one sample, checked against plm -------------------------------------------------------
d7 <- readRDS("nimbus_P7.rds")
s1 <- both_tests(d7)
pd <- pdata.frame(d7, index = c("id", "t"))
ph <- phtest(plm(y ~ x, data = pd, model = "within"), plm(y ~ x, data = pd, model = "random"))
stopifnot(abs(unname(ph$statistic) - s1[["H_classic"]]) < 1e-6)
sample7 <- as.list(s1)

## ---- Monte Carlo ------------------------------------------------------------------------------
set.seed(tr$params$seed$P7_mc)
mc <- t(vapply(seq_len(P7$R_mc), function(r) both_tests(gen_P7(P7)), numeric(8)))
levels <- c(0.01, 0.02, 0.05, 0.10, 0.15, 0.20)
curve <- data.frame(nominal = levels,
                    classic = vapply(levels, function(a) mean(mc[, "p_classic"] < a), numeric(1)),
                    robust = vapply(levels, function(a) mean(mc[, "p_robust"] < a), numeric(1)))
res7 <- list(R = nrow(mc),
             size_classic = mean(mc[, "p_classic"] < 0.05), size_robust = mean(mc[, "p_robust"] < 0.05),
             size_classic_01 = mean(mc[, "p_classic"] < 0.01), size_robust_01 = mean(mc[, "p_robust"] < 0.01),
             size_classic_10 = mean(mc[, "p_classic"] < 0.10), size_robust_10 = mean(mc[, "p_robust"] < 0.10),
             mcse_size = sqrt(0.05 * 0.95 / nrow(mc)),
             mean_H_classic = mean(mc[, "H_classic"]), mean_H_robust = mean(mc[, "H_robust"]),
             sd_diff = sd(mc[, "fe"] - mc[, "re"]), mean_sd_diff_classic = mean(sqrt(mc[, "vd_classic"])),
             mean_fe = mean(mc[, "fe"]), mean_re = mean(mc[, "re"]), sd_fe = sd(mc[, "fe"]), sd_re = sd(mc[, "re"]),
             mean_varphi = mean(mc[, "varphi"]),
             curve = curve, sample = sample7, phtest_sample = unname(ph$statistic))
res7$ratio_sd_diff <- res7$sd_diff / res7$mean_sd_diff_classic
saveRDS(list(P7 = res7), "res_p7.rds")
cat(sprintf("P7 Monte Carlo (R = %d): size at 5%%: classical %.4f, robust %.4f (MC se %.4f); at 1%%: %.4f %.4f; at 10%%: %.4f %.4f\n",
            res7$R, res7$size_classic, res7$size_robust, res7$mcse_size, res7$size_classic_01, res7$size_robust_01,
            res7$size_classic_10, res7$size_robust_10))
cat(sprintf("   SD of FE - RE across samples %.4f vs the classical formula %.4f (ratio %.2f); mean FE %.4f, RE %.4f\n",
            res7$sd_diff, res7$mean_sd_diff_classic, res7$ratio_sd_diff, res7$mean_fe, res7$mean_re))
cat(sprintf("   one sample: H classical %.3f (p %.3f), robust %.3f (p %.3f)\n",
            s1[["H_classic"]], s1[["p_classic"]], s1[["H_robust"]], s1[["p_robust"]]))
