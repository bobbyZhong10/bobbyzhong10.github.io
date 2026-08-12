# ============================================================================
# 03_naive_and_reg.R  —  Study B (observational)
#   (1) naive difference-in-means  -> shows selection bias (upward)
#   (2) OLS regression adjustment  -> recovers ATT under CIA (heterogeneous
#       effect fully interacted, ATT read off at treated-firm covariates)
# CIA-holds arm uses D_cia / Y_cia ; robust (HC1) SEs via sandwich.
# ============================================================================
suppressMessages(library(sandwich))
suppressMessages(library(lmtest))

d   <- readRDS("northwind.rds")
obs <- d$obs
obs$D <- obs$D_cia
obs$Y <- obs$Y_cia

# ---- (1) naive difference in means ----------------------------------------
y1 <- obs$Y[obs$D==1]; y0 <- obs$Y[obs$D==0]
naive <- mean(y1) - mean(y0)
se_naive <- sqrt(var(y1)/length(y1) + var(y0)/length(y0))

# ---- (2a) plain additive OLS: Y ~ D + X  (constant-effect spec) ------------
m_add <- lm(Y ~ D + size + dig + eff + ind, data = obs)
V_add <- vcovHC(m_add, type = "HC1")
b_add <- coef(m_add)["D"]; se_add <- sqrt(V_add["D","D"])

# ---- (2b) fully interacted OLS -> ATT ------------------------------------
# center covariates at TREATED means so the D coefficient = ATT.
Xv <- c("size","dig","eff")
tm <- colMeans(obs[obs$D==1, Xv])
oc <- obs
for (v in Xv) oc[[paste0(v,"_c")]] <- oc[[v]] - tm[v]
# industry: use treated-group proportions via interaction; center dummies too
ind_dum <- model.matrix(~ ind, obs)[, -1, drop=FALSE]     # indB, indC, indD
ind_tm  <- colMeans(ind_dum[obs$D==1, , drop=FALSE])
ind_c   <- sweep(ind_dum, 2, ind_tm, "-")
colnames(ind_c) <- paste0(colnames(ind_dum), "_c")
oc <- cbind(oc, ind_c)

form <- as.formula(paste0(
  "Y ~ D * (size_c + dig_c + eff_c + ",
  paste(colnames(ind_c), collapse = " + "), ")"))
m_int <- lm(form, data = oc)
V_int <- vcovHC(m_int, type = "HC1")
b_att <- coef(m_int)["D"]; se_att <- sqrt(V_int["D","D"])

cat("==== Study B : naive + regression adjustment (CIA arm) ====\n")
cat(sprintf("naive diff-in-means : %.4f  (SE %.4f)\n", naive, se_naive))
cat(sprintf("OLS additive  (D)   : %.4f  (SE %.4f)\n", b_add, se_add))
cat(sprintf("OLS interacted ATT  : %.4f  (SE %.4f)\n", b_att, se_att))
cat(sprintf("true ATT (CIA)      : %.4f\n", d$truths$true_ATT_cia))
cat(sprintf("true ATE            : %.4f\n", d$truths$true_ATE_pop))

saveRDS(list(naive=naive, se_naive=se_naive,
             ols_add=b_add, se_ols_add=se_add,
             ols_att=b_att, se_ols_att=se_att),
        "res_naive_reg.rds")
