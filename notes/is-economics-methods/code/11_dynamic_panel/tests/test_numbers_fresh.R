# ------------------------------------------------------------------
# tests/test_numbers_fresh.R -- numbers.json must not be older than the res_*.rds files
# (Phase 2 Task 7): its keys and values equal collect_numbers11() to 1e-12.
# Run from the chapter code directory (tests/run_tests.R does this).
# ------------------------------------------------------------------
source("numbers_def.R"); j <- jsonlite::read_json("numbers.json"); n <- collect_numbers11()
stopifnot(setequal(names(j), names(n)), all(abs(unlist(j[names(n)]) - unlist(n)) < 1e-12))
# One rule for tiny p-values (numbers_def.R header): a p-value key (ends in _p, _p_full or _p_collapsed, not a
# mean or a share) is stored as is when it is at least 1e-4; a smaller one is stored as the flag *_p_lt_0001 = 1
# plus log10 where it is positive (a p-value that underflows to exactly 0, e.g. m_1 with Z = -40 in some of
# the D3 and D4 fits, has the flag only); none is dropped or stored as 0; the plm one-step Sargan p is not stored.
k <- names(n)
isp <- grepl("(^|_)p$|_p_(full|collapsed)$", k) & !grepl("mean|share", k)
stopifnot(all(unlist(n[isp]) >= 1e-4), all(unlist(n[isp]) <= 1))
flag <- k[grepl("_p_lt_0001$", k)]
lg <- k[grepl("_p_log10$", k)]
stopifnot(all(unlist(n[flag]) == 1), length(flag) > 0,                      # m1_p is below 1e-4 in every GMM fit here
          all(sub("_p_log10$", "_p_lt_0001", lg) %in% flag),                 # a log10 comes with its flag ...
          all(unlist(n[lg]) < -4), length(lg) > 0,                           # ... and is below log10(1e-4)
          !any(sub("_p_lt_0001$", "_p", flag) %in% k),                       # a flagged p-value is not stored as a number as well
          !any(grepl("sargan_onestep", k)), !any(grepl("(^|_)u_sd|_u_it", k)))
# every GMM fit has its m1 p-value (flag or number): none was silently dropped
m1 <- unique(sub("_m1$", "", k[grepl("_m1$", k)]))
stopifnot(length(m1) > 0, all(paste0(m1, "_m1_p") %in% k | paste0(m1, "_m1_p_lt_0001") %in% k))
# fix round 2 (11.3a E1, 11.3b C1 and I2): the one-step large-sample D4 limits, the seed and the summaries of the
# Anderson-Hsiao (differenced instrument) Monte Carlo at rho = 0.90
stopifnot(all(c("dgp_D4_large_bb_onestep", "dgp_D4_large_ab_onestep", "dgp_D4_large_bb", "dgp_seed_D1p_ah_mc",
                "d1p_ah_diff_mc_R", "d1p_ah_diff_mc_cover", "d1p_ah_diff_mc_cover_mcse", "d1p_ah_diff_mc_median",
                "d1p_ah_diff_mc_iqr", "d1p_ah_diff_mc_share_ci_below_truth") %in% k),
          !any(c("d1p_ah_diff_mc_sd", "d1p_ah_diff_mc_mean") %in% k),                                            # no SD: a just-identified IV has no finite moments
          n$d1p_ah_diff_mc_R == 1000, n$dgp_seed_D1p_ah_mc == 1161)
# numbers_cited.json is a JSON object (empty for this chapter until Task 11.2 adds cited numbers)
stopifnot(is.list(jsonlite::read_json("numbers_cited.json")))
