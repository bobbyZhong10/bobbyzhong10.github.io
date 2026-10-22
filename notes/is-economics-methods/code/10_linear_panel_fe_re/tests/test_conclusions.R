# ------------------------------------------------------------------
# tests/test_conclusions.R -- the teaching conclusions of chapter 10.
# Run from the chapter code directory after 01 to 08 (tests/run_tests.R does this).
# The first block is the chapter card's assertions, verbatim; the second block adds the
# algebraic identities the text states (two-way Mundlak, the classical Hausman of the legacy
# code, RE as a matrix-weighted average, dummy-variable regression = FE).
# ------------------------------------------------------------------
source("../../tools/R/chapter_kit.R")

## ---- chapter card (verbatim) ---------------------------------------------------
tr <- readRDS("truth10.rds"); x <- readRDS("res_extra.rds"); r <- readRDS("res_static.rds")
stopifnot(abs(x$mundlak$b_x - r$fe) < 1e-8, abs(x$mundlak_unbal$b_x - x$fe_unbal) < 1e-8)    # Mundlak 恒等式
stopifnot(x$P7$size_classic > 0.10, abs(x$P7$size_robust - 0.05) < 0.02)
stopifnot(x$P2$sd_fd_rw < x$P2$sd_fe_rw, x$P2$sd_fe_iid < x$P2$sd_fd_iid)                  # 谁更有效看误差序列相关
stopifnot(x$P3$lead_p < 0.01, x$P1$lead_p > 0.05)
stopifnot(covers(x$P4$ht$ci, tr$P4$gamma), !covers(x$P4$re$ci, tr$P4$gamma))
stopifnot(covers(x$P5$fe$ci, tr$P5$fe_plim), abs(tr$P5$fe_plim - tr$P5$avg_slope) > 0.05, covers(x$P5b$fe$ci, tr$P5$avg_slope))
stopifnot(covers(x$P6$attr_alpha$ci, tr$P6$beta), !covers(x$P6$attr_eps$ci, tr$P6$beta))

## ---- identities the text states (beyond the card) ------------------------------
# two-way Mundlak (time dummies, X-bar and the unit averages of the time dummies) = TWFE,
# balanced (P1) and unbalanced (P6, attrition on alpha); without the averaged time dummies the
# identity fails in the unbalanced panel (Wooldridge 2019)
stopifnot(abs(x$mundlak_tw$b_x - x$twfe) < 1e-8, abs(x$mundlak_tw_unbal$b_x - x$twfe_unbal) < 1e-8,
          abs(x$mundlak_tw_naive_unbal - x$twfe_unbal) > 1e-6)
# hausman_classic on the plm FE and RE fits = the legacy phtest statistic; the robust (Mundlak) test also rejects
stopifnot(abs(x$P1$haus_classic$stat - r$haus_stat) < 1e-6, x$P1$haus_classic$df == 1, x$P1$haus_robust$p < 1e-6)
# Mundlak: xi_hat = between minus within (balanced panel)
stopifnot(abs(x$mundlak$b_xbar - (x$P1$be$est - r$fe)) < 1e-8)
# RE = w * between + (1 - w) * within, w from varphi_hat (balanced panel, one regressor)
stopifnot(abs(x$P1$re_decomp$w * x$P1$be$est + (1 - x$P1$re_decomp$w) * r$fe - r$re) < 1e-8,
          x$P1$re_decomp$w > 0, x$P1$re_decomp$w < 1)
# least-squares dummy-variable regression = FE (Wooldridge 2010 §10.5.3)
stopifnot(abs(x$P1$lsdv - r$fe) < 1e-8)
