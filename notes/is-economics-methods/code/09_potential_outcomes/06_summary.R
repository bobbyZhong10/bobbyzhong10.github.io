# ============================================================================
# 06_summary.R  —  collect every quotable number for Module 09.
#   Truths, RCT, naive, OLS, matching, IPW (HT + Hajek), AIPW, and the
#   CIA-fails variant's residual bias. Printed to console, rounded to 3 dp.
# ============================================================================
suppressMessages({ library(sandwich) })

d   <- readRDS("northwind.rds")
tr  <- d$truths
re  <- readRDS("res_experiment.rds")
rn  <- readRDS("res_naive_reg.rds")
rm4 <- readRDS("res_matching_ipw_aipw.rds")

r3 <- function(x) formatC(round(x, 3), format = "f", digits = 3)

cat("################  MODULE 09  Northwind AI assistant  ################\n\n")

cat("== TRUTHS (from individual effects) ==\n")
cat(sprintf("  true ATE (population)     : %s\n", r3(tr$true_ATE_pop)))
cat(sprintf("  true ATE (RCT sample)     : %s\n", r3(tr$true_ATE_rct)))
cat(sprintf("  true ATT (CIA arm)        : %s\n", r3(tr$true_ATT_cia)))
cat(sprintf("  true ATT (CIA-fails arm)  : %s\n", r3(tr$true_ATT_fail)))
cat(sprintf("  adoption rate CIA / fail  : %s / %s\n",
            r3(tr$adopt_rate_cia), r3(tr$adopt_rate_fail)))

cat("\n== STUDY A : randomized pilot ==\n")
cat(sprintf("  diff-in-means (ATE)       : %s   SE %s\n", r3(re$dim_est), r3(re$se)))
cat(sprintf("  95%% CI                    : [%s, %s]\n", r3(re$ci[1]), r3(re$ci[2])))
cat(sprintf("  permutation p-value       : %s\n", r3(re$p_perm)))

cat("\n== STUDY B : observational (CIA arm) ==\n")
cat(sprintf("  naive diff-in-means       : %s   SE %s   (bias vs ATT = %s)\n",
            r3(rn$naive), r3(rn$se_naive), r3(rn$naive - tr$true_ATT_cia)))
cat(sprintf("  OLS regression adjustment : %s   SE %s\n", r3(rn$ols_att), r3(rn$se_ols_att)))
cat(sprintf("  NN matching (MatchIt pt)  : %s\n", r3(rm4$match_pt)))
cat(sprintf("  NN matching (AI SE)       : %s   AI-SE %s  (naive OLS-SE %s)\n",
            r3(rm4$match_ai_pt), r3(rm4$match_ai_se), r3(rm4$match_naive_se)))
cat(sprintf("  IPW Horvitz-Thompson      : %s   SE %s\n", r3(rm4$ipw_ht), r3(rm4$se_ht)))
cat(sprintf("  IPW Hajek                 : %s   SE %s\n", r3(rm4$ipw_hajek), r3(rm4$se_hajek)))
cat(sprintf("  AIPW (doubly robust)      : %s   SE %s\n", r3(rm4$aipw), r3(rm4$se_aipw)))
cat(sprintf("  propensity range          : [%s, %s]\n",
            r3(rm4$e_range[1]), r3(rm4$e_range[2])))

# ---- CIA-fails variant : residual bias after full adjustment ---------------
obs <- d$obs
D <- obs$D_fail; Y <- obs$Y_fail
ATT_f <- tr$true_ATT_fail
naive_f <- mean(Y[D==1]) - mean(Y[D==0])
ps <- glm(D ~ size+dig+eff+ind, data = cbind(obs, D=D), family = binomial)
e  <- predict(ps, type = "response"); w <- e/(1-e)
dff <- cbind(obs, D=D, Y=Y)
m0 <- predict(lm(Y ~ size+dig+eff+ind, data = dff[D==0,]), dff)
pp <- mean(D)
aipw_f <- mean((D*(Y-m0) - (1-D)*w*(Y-m0))/pp)
# OLS interacted ATT on the fail arm
Xv <- c("size","dig","eff"); tm <- colMeans(dff[D==1, Xv]); oc <- dff
for (v in Xv) oc[[paste0(v,"_c")]] <- oc[[v]] - tm[v]
idm <- model.matrix(~ind, dff)[,-1,drop=FALSE]; itm <- colMeans(idm[D==1,,drop=FALSE])
idc <- sweep(idm,2,itm,"-"); colnames(idc) <- paste0(colnames(idc),"_c"); oc <- cbind(oc, idc)
fo <- as.formula(paste0("Y~D*(size_c+dig_c+eff_c+",paste(colnames(idc),collapse="+"),")"))
ols_f <- coef(lm(fo, oc))["D"]

cat("\n== STUDY B : CIA-FAILS variant (selection on unobserved U) ==\n")
cat(sprintf("  true ATT                  : %s\n", r3(ATT_f)))
cat(sprintf("  naive diff-in-means       : %s   (bias = %s)\n",
            r3(naive_f), r3(naive_f - ATT_f)))
cat(sprintf("  OLS adjustment            : %s   (residual bias = %s)\n",
            r3(ols_f), r3(ols_f - ATT_f)))
cat(sprintf("  AIPW (doubly robust)      : %s   (residual bias = %s)\n",
            r3(aipw_f), r3(aipw_f - ATT_f)))
cat("  -> adjustment on X removes only part of the bias; the unobserved\n")
cat("     confounder leaves a stubborn residual gap (motivates IV, ch.10).\n")

cat("\n####################################################################\n")
