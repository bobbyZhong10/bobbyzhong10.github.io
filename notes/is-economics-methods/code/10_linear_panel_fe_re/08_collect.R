# ------------------------------------------------------------------
# 08_collect.R -- res_extra.rds: the results of 03 to 07 under the field names of the chapter card
# (tests/test_conclusions.R reads it together with res_static.rds and truth10.rds).
#   mundlak, mundlak_tw, twfe                      block P1 (balanced): Mundlak = FE, two-way Mundlak = TWFE
#   mundlak_unbal, fe_unbal, mundlak_tw_unbal,     block P6 with attrition on alpha_i (unbalanced): the same
#   twfe_unbal, mundlak_tw_naive_unbal             identities; the naive two-way version without the averaged
#                                                  month dummies
#   P1 ... P7                                      per-block results (03 to 07)
# Reads res_p1.rds, res_p23.rds, res_p45.rds, res_p6.rds, res_p7.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
p1 <- readRDS("res_p1.rds"); p23 <- readRDS("res_p23.rds"); p45 <- readRDS("res_p45.rds")
p6 <- readRDS("res_p6.rds"); p7 <- readRDS("res_p7.rds")
x <- list(mundlak = p1$mundlak, mundlak_tw = p1$mundlak_tw, twfe = p1$twfe,
          mundlak_unbal = p6$mundlak_unbal, fe_unbal = p6$fe_unbal,
          mundlak_tw_unbal = p6$mundlak_tw_unbal, twfe_unbal = p6$twfe_unbal,
          mundlak_tw_naive_unbal = p6$mundlak_tw_naive_unbal,
          P1 = p1[setdiff(names(p1), c("mundlak", "mundlak_tw", "twfe"))],
          P2 = p23$P2, P3 = p23$P3, P4 = p45$P4, P5 = p45$P5, P5b = p45$P5b, P6 = p6$P6, P7 = p7$P7)
saveRDS(x, "res_extra.rds")
cat("Saved res_extra.rds\n")
