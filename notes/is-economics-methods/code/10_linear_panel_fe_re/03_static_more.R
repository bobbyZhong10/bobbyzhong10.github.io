# ------------------------------------------------------------------
# 03_static_more.R -- block P1 (the inherited static Nimbus panel) beyond the legacy 02:
#   between and FD estimates, the RE variance components and varphi, RE as a matrix-weighted
#   average of between and within, dummy-variable regression = FE, the Mundlak regression
#   (one-way and with time effects) and the robust Hausman test, the classical Hausman statistic
#   from hausman_classic (= the legacy phtest), the lead test of strict exogeneity, the AR(1)
#   coefficient of the FD residuals, and Chamberlain's minimum distance estimator.
# Reads nimbus.rds and res_static.rds (02_static_fe_re.R); writes res_p1.rds.
# Standard errors: "se" is the cluster-robust (by merchant, CR0) standard error unless named
# se_classic (the textbook non-robust formula).
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("panel_fns.R")
suppressPackageStartupMessages(library(plm))

d <- readRDS("nimbus.rds")$static; st <- readRDS("res_static.rds"); tr <- readRDS("truth10.rds")
pd <- pdata.frame(d, index = c("id", "t"))
cl <- function(m) vcovHC(m, method = "arellano", type = "HC0", cluster = "group")
est_row <- function(m, nm = "x") {
  e <- unname(coef(m)[nm]); se <- sqrt(cl(m)[nm, nm]); sc <- sqrt(vcov(m)[nm, nm])
  list(est = e, se = se, se_classic = sc, ci = ci95(e, se), covers_truth = as.numeric(covers_ci(ci95(e, se), tr$P1$beta)))
}
covers_ci <- function(ci, v) ci[1] <= v && v <= ci[2]

pool <- plm(y ~ x, data = pd, model = "pooling")
re <- plm(y ~ x, data = pd, model = "random")
fe <- plm(y ~ x, data = pd, model = "within")
be <- plm(y ~ x, data = pd, model = "between")
fd <- plm(y ~ x, data = pd, model = "fd")
tw <- plm(y ~ x, data = pd, model = "within", effect = "twoways")
stopifnot(abs(coef(fe) - st$fe) < 1e-12, abs(coef(re)["x"] - st$re) < 1e-12, abs(coef(pool)["x"] - st$pool) < 1e-12)

out <- list(pool = est_row(pool), re = est_row(re), fe = est_row(fe), fd = est_row(fd))
out$be <- list(est = unname(coef(be)["x"]), se_classic = sqrt(vcov(be)["x", "x"]))   # one row per merchant
out$be$ci <- ci95(out$be$est, out$be$se_classic)
out$twfe <- unname(coef(tw))

## ---- RE variance components and quasi-demeaning (Wooldridge 2010 §10.7.2) -------------------
ec <- ercomp(re)
out$varphi <- unname(ec$theta)
out$sigma2_eps <- unname(ec$sigma2[["idios"]]); out$sigma2_alpha <- unname(ec$sigma2[["id"]])
out$share_alpha <- out$sigma2_alpha / (out$sigma2_alpha + out$sigma2_eps)
# RE = w * between + (1 - w) * within, w = (1 - varphi)^2 S_B / (S_W + (1 - varphi)^2 S_B)
TT <- 6
xw <- within_tf(d$y, d$x, d$id)$X[, 1]
xb <- tapply(d$x, d$id, mean); SB <- TT * sum((xb - mean(xb))^2); SW <- sum(xw^2)
w <- (1 - out$varphi)^2 * SB / (SW + (1 - out$varphi)^2 * SB)
out$re_decomp <- list(w = w, S_W = SW, S_B = SB, check = w * out$be$est + (1 - w) * st$fe - st$re)

## ---- least-squares dummy-variable regression (Wooldridge 2010 §10.5.3) -----------------------
lsdv <- lm(y ~ x + factor(id), data = d)
out$lsdv <- unname(coef(lsdv)["x"])
out$lsdv_se_classic <- sqrt(vcov(lsdv)["x", "x"])     # = the FE classical SE (same residuals, same df)

## ---- Mundlak regression and Hausman tests (Wooldridge 2010 §10.7.3) --------------------------
mk <- mundlak_fit(d$y, d$x, d$id)
mk_tw <- mundlak_fit(d$y, d$x, d$id, time_fe = TRUE, t = d$t)
out$mundlak <- mk; out$mundlak_tw <- mk_tw
out$mundlak_se_x <- sqrt(mk$V["x", "x"]); out$mundlak_se_xbar <- sqrt(mk$V["xbar_x", "xbar_x"])
out$haus_classic <- hausman_classic(coef(fe), vcov(fe), coef(re)["x"], vcov(re)["x", "x"])
out$haus_robust <- hausman_robust(mk)
out$haus_robust_tw <- hausman_robust(mk_tw)
aux <- phtest(y ~ x, data = pd, method = "aux", vcov = function(m) cl(m))      # RE + Xbar, robust Wald (plm)
out$haus_aux <- list(stat = unname(aux$statistic), df = unname(aux$parameter), p = aux$p.value)

## ---- strict exogeneity: lead test (Wooldridge 2010 §10.7.1, p. 325) --------------------------
lt <- lead_test(d$y, d$x, d$id, d$t)
out$lead_p <- as.numeric(lt); out$lead_stat <- attr(lt, "stat"); out$lead_delta <- unname(attr(lt, "delta"))
out$lead_se <- unname(attr(lt, "se")); out$lead_n <- attr(lt, "n")

## ---- FD residuals: AR(1) coefficient (Wooldridge 2010 §10.6.3) -------------------------------
out$fd_ar1 <- fd_ar1(fd_fit(d$y, d$x, d$id, d$t))

## ---- Chamberlain's method by minimum distance (Wooldridge 2010 §11.1.2) ----------------------
out$chamberlain <- chamberlain_md(d$y, d$x, d$id, d$t)
out$chamberlain$ci <- ci95(out$chamberlain$beta, out$chamberlain$se)

saveRDS(out, "res_p1.rds")
cat(sprintf("P1: pooled %.4f  RE %.4f  FE %.4f  BE %.4f  FD %.4f  TWFE %.4f  (truth %.2f)\n",
            out$pool$est, out$re$est, out$fe$est, out$be$est, out$fd$est, out$twfe, tr$P1$beta))
cat(sprintf("    varphi %.4f, w(between) %.4f; Mundlak b_x %.4f, xi %.4f (BE - FE %.4f)\n",
            out$varphi, w, mk$b_x, mk$b_xbar, out$be$est - st$fe))
cat(sprintf("    Hausman classic %.2f (p %.2g); robust %.2f (p %.2g); plm aux %.2f\n",
            out$haus_classic$stat, out$haus_classic$p, out$haus_robust$stat, out$haus_robust$p, out$haus_aux$stat))
cat(sprintf("    lead test p %.4f (delta %.4f, se %.4f); FD residual rho %.4f; Chamberlain beta %.4f (se %.4f), overid %.2f on %d df, p %.3f\n",
            out$lead_p, out$lead_delta, out$lead_se, out$fd_ar1$rho, out$chamberlain$beta, out$chamberlain$se,
            out$chamberlain$stat, out$chamberlain$df, out$chamberlain$p))
