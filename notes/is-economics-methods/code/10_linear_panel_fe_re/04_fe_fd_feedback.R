# ------------------------------------------------------------------
# 04_fe_fd_feedback.R -- blocks P2 and P3.
#   P2  FE vs FD with iid errors (FE.3) and random-walk errors (FD.3): one sample (estimates,
#       cluster-robust SEs, the AR(1) coefficient of the FD residuals, plm's pwartest and pwfdtest)
#       and a Monte Carlo of R_mc = 1000 samples for the sampling SDs (Wooldridge 2010 §10.7.1).
#   P3  feedback from last month's shock: FE and FD against their probability limits, the lead test
#       of strict exogeneity (§10.7.1, p. 325) and Chamberlain's overidentification test (§11.1.2).
# Reads nimbus_P2.rds, nimbus_P3.rds, truth10.rds; writes res_p23.rds.
# "se" is the cluster-robust (by merchant, CR0) standard error.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("panel_fns.R")
suppressPackageStartupMessages(library(plm))
tr <- readRDS("truth10.rds"); P2 <- tr$params$P2; P3 <- tr$params$P3

row_fe <- function(f, truth) {
  e <- unname(f$b[1]); se <- sqrt(f$V[1, 1])
  list(est = e, se = se, ci = ci95(e, se), covers_truth = as.numeric(ci95(e, se)[1] <= truth && truth <= ci95(e, se)[2]))
}
row_fd <- function(f, truth) {
  e <- unname(f$b["x"]); se <- sqrt(f$V["x", "x"])
  list(est = e, se = se, ci = ci95(e, se), covers_truth = as.numeric(ci95(e, se)[1] <= truth && truth <= ci95(e, se)[2]))
}

## =================================================================
## P2: one sample
## =================================================================
d2 <- readRDS("nimbus_P2.rds")
S2 <- list()
for (k in c("iid", "rw")) {
  y <- d2[[paste0("y_", k)]]
  ff <- fe_fit(y, d2$x, d2$id); fdf <- fd_fit(y, d2$x, d2$id, d2$t)
  pd <- pdata.frame(data.frame(id = d2$id, t = d2$t, y = y, x = d2$x), index = c("id", "t"))
  S2[[k]] <- list(fe = row_fe(ff, P2$beta), fd = row_fd(fdf, P2$beta), fd_ar1 = fd_ar1(fdf),
                  pwartest_p = pwartest(y ~ x, data = pd)$p.value,
                  pwfdtest_fd_p = pwfdtest(y ~ x, data = pd, h0 = "fd")$p.value,
                  pwfdtest_fe_p = pwfdtest(y ~ x, data = pd, h0 = "fe")$p.value)
}

## =================================================================
## P2: Monte Carlo, R_mc samples (seed P2_mc)
## =================================================================
set.seed(tr$params$seed$P2_mc)
mc <- t(vapply(seq_len(P2$R_mc), function(r) {
  d <- gen_P2(P2)
  out <- numeric(0)
  for (k in c("iid", "rw")) {
    y <- d[[paste0("y_", k)]]
    ff <- fe_fit(y, d$x, d$id); fdf <- fd_fit(y, d$x, d$id, d$t)
    out <- c(out, ff$b[1], sqrt(ff$V[1, 1]), fdf$b["x"], sqrt(fdf$V["x", "x"]))
  }
  out
}, numeric(8)))
colnames(mc) <- c("fe_iid", "se_fe_iid", "fd_iid", "se_fd_iid", "fe_rw", "se_fe_rw", "fd_rw", "se_fd_rw")
cover <- function(est, se) mean(abs(est - P2$beta) <= qnorm(0.975) * se)
res2 <- list(R = nrow(mc),
             sd_fe_iid = sd(mc[, "fe_iid"]), sd_fd_iid = sd(mc[, "fd_iid"]),
             sd_fe_rw = sd(mc[, "fe_rw"]), sd_fd_rw = sd(mc[, "fd_rw"]),
             mean_fe_iid = mean(mc[, "fe_iid"]), mean_fd_iid = mean(mc[, "fd_iid"]),
             mean_fe_rw = mean(mc[, "fe_rw"]), mean_fd_rw = mean(mc[, "fd_rw"]),
             mean_se_fe_iid = mean(mc[, "se_fe_iid"]), mean_se_fd_iid = mean(mc[, "se_fd_iid"]),
             mean_se_fe_rw = mean(mc[, "se_fe_rw"]), mean_se_fd_rw = mean(mc[, "se_fd_rw"]),
             cover_fe_iid = cover(mc[, "fe_iid"], mc[, "se_fe_iid"]), cover_fd_iid = cover(mc[, "fd_iid"], mc[, "se_fd_iid"]),
             cover_fe_rw = cover(mc[, "fe_rw"], mc[, "se_fe_rw"]), cover_fd_rw = cover(mc[, "fd_rw"], mc[, "se_fd_rw"]),
             sample = S2)
res2$ratio_fd_fe_iid <- res2$sd_fd_iid / res2$sd_fe_iid
res2$ratio_fe_fd_rw <- res2$sd_fe_rw / res2$sd_fd_rw

## =================================================================
## P3: feedback
## =================================================================
d3 <- readRDS("nimbus_P3.rds")
f3 <- fe_fit(d3$y, d3$x, d3$id); fd3 <- fd_fit(d3$y, d3$x, d3$id, d3$t)
lt <- lead_test(d3$y, d3$x, d3$id, d3$t)
ch <- chamberlain_md(d3$y, d3$x, d3$id, d3$t)
res3 <- list(fe = row_fe(f3, P3$beta), fd = row_fd(fd3, P3$beta),
             lead_p = as.numeric(lt), lead_stat = attr(lt, "stat"), lead_delta = unname(attr(lt, "delta")),
             lead_se = unname(attr(lt, "se")), lead_n = attr(lt, "n"),
             chamberlain = c(ch[c("beta", "se", "stat", "df", "p")], list(ci = ci95(ch$beta, ch$se))))
res3$fe$covers_plim <- as.numeric(res3$fe$ci[1] <= tr$P3$fe_plim && tr$P3$fe_plim <= res3$fe$ci[2])
res3$fd$covers_plim <- as.numeric(res3$fd$ci[1] <= tr$P3$fd_plim && tr$P3$fd_plim <= res3$fd$ci[2])

saveRDS(list(P2 = res2, P3 = res3), "res_p23.rds")
cat(sprintf("P2 Monte Carlo (R = %d): SD iid FE %.4f FD %.4f; SD random walk FE %.4f FD %.4f\n",
            res2$R, res2$sd_fe_iid, res2$sd_fd_iid, res2$sd_fe_rw, res2$sd_fd_rw))
cat(sprintf("P2 sample: iid FE %.4f (%.4f) FD %.4f (%.4f), FD residual rho %.3f; rw FE %.4f (%.4f) FD %.4f (%.4f), rho %.3f\n",
            S2$iid$fe$est, S2$iid$fe$se, S2$iid$fd$est, S2$iid$fd$se, S2$iid$fd_ar1$rho,
            S2$rw$fe$est, S2$rw$fe$se, S2$rw$fd$est, S2$rw$fd$se, S2$rw$fd_ar1$rho))
cat(sprintf("P3: FE %.4f (plim %.4f) FD %.4f (plim %.4f); lead test p %.3g (delta %.4f); Chamberlain overid %.1f on %d df, p %.3g\n",
            res3$fe$est, tr$P3$fe_plim, res3$fd$est, tr$P3$fd_plim, res3$lead_p, res3$lead_delta,
            ch$stat, ch$df, ch$p))
