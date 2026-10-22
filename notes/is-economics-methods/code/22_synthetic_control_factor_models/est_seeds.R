# ------------------------------------------------------------------
# est_seeds.R -- the estimation seeds of chapter 22, in one place (Ruling S9, review
# Minor 6; the data-generating seeds are seeds22 in 01_dgp.R). Sourced by 03, 04, 05,
# 05b, 05c and 06. Each seed fixes one random step (a bootstrap, a permutation set, a search
# start, a random assignment), so every result in numbers.json is reproducible; changing a
# value here changes the numbers that step produces and nothing else.
# ------------------------------------------------------------------
est_seeds22 <- c(
  s1_sdid_placebo   = 22101,   # 03: synthdid placebo variance of SDID (500 replications)
  s1_sc_placebo     = 22102,   # 03: placebo variance of synthdid's SC
  s1_did_placebo    = 22103,   # 03: placebo variance of synthdid's DiD
  s1_asc_summary    = 22104,   # 03: augsynth summary (conformal inference)
  s1_conformal_iid  = 22105,   # 03: the 9999 iid permutations of conformal_p
  s1_scpi           = 22106,   # 03: scpi simulations (sims = 200)
  s2_pen_orderings  = 22201,   # 04: the 20 random donor orderings (uniqueness check)
  s2_asc_summary    = 22202,   # 04: augsynth summary
  s2_sdid_placebo   = 22203,   # 04: synthdid placebo variance of SDID
  s3_fits           = 1,       # 05: fect, CS, multisynth and SDID-by-cohort bootstraps on the stored draw
  s3_screen         = 1,       # 05b: the same recipes on gen_S3(1..24) (Task 9's recorded seed = 1)
  s3_scpi           = 22301,   # 05c: scpi simulations on the stored draw (sims = 200; Ruling S12)
  s4_search_starts  = 22401,   # 06: the 100 random starts of the design search
  s4_random_draws   = 22402    # 06: the random 5-metro assignments, one per replicate
)
